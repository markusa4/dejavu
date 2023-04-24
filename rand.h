// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_RAND_H
#define DEJAVU_RAND_H

#include "ir.h"
#include "groups.h"

namespace dejavu::search_strategy {

    /**
     * \brief IR search using random walks.
     *
     * Performs random walks of the IR shared_tree, sifting resulting automorphisms into the given Schreier structure. If the
     * Schreier structure is complete with respect to the base, or the probabilistic abort criterion satisfied, the
     * process terminates. The algorithm guarantees to find all automorphisms up to the specified error bound.
     *
     * Alternatively, a limit for the amount of discovered differing leaves can be set.
     *
     * Objects contain a local workspace which contain datastructures only used for random_ir, as well as settings and
     * for statistics of the search.
     */
    class random_ir {
        //shared_leaves leaf_storage; /**< stores all the leaves */
        std::default_random_engine generator; /**< random number generator */

        groups::automorphism_workspace automorphism;
        groups::schreier_workspace     w;
        std::vector<int>               heuristic_reroll;

        /**
         * Loads the leaf into the \p local_state. Only works if base of leaf was stored, i.e., whenever
         * `full_save = false` was used in the constructor.
         *
         * @param leaf The leaf.
         * @param g The graph.
         * @param local_state Local state in which the leaf will be stored.
         * @param start_from State from which the walk to the leaf is performed.
         */
        static void load_state_from_leaf(sgraph *g, ir::controller &local_state, ir::reduced_save &start_from,
                                  ir::stored_leaf *leaf) {
            assert(leaf->get_store_type() == ir::stored_leaf::STORE_BASE);
            std::vector<int> base;
            base.reserve(leaf->get_lab_or_base_size());
            for(int i = 0; i < leaf->get_lab_or_base_size(); ++i) base.push_back(leaf->get_lab_or_base()[i]);
            local_state.walk(g, start_from, base);
        }

        /**
         * Co-routine which adds a leaf to leaf_storage, and sifts resulting automorphism into a given group.
         *
         * @param g
         * @param hook
         * @param group
         * @param leaf_storage
         * @param local_state
         * @param root_save
         * @param uniform
         * @return
         */
        bool add_leaf_to_storage_and_group(sgraph *g, dejavu_hook *hook, groups::shared_schreier &group,
                                           ir::shared_leaves &leaf_storage, ir::controller &local_state,
                                           ir::reduced_save &root_save, bool uniform) {
            bool used_load = false;

            assert(g->v_size == local_state.c->cells);

            // Outer loop to handle hash collisions -- at most h_hash_col_limit will be checked and stored
            for(int hash_offset = 0; hash_offset < h_hash_col_limit; ++hash_offset) {
                // After first hash collision, we write a stronger invariant
                if(hash_offset == 1) local_state.write_strong_invariant(g); // TODO might lead to discarding leaves... should actually save whether cert failed, and then just do this before the loop if that flag is set

                // First, test whether leaf with same hash has already been stored
                const long hash_c = local_state.T->get_hash() + hash_offset; // '+hash_offset' is for hash collisions
                auto other_leaf = leaf_storage.lookup_leaf(hash_c);

                // If not, add leaf to leaf_storage
                if (other_leaf == nullptr) {
                    s_rolling_success = (9.0 * s_rolling_success + 0.0) / 10.0;
                    leaf_storage.add_leaf(hash_c, *local_state.c, local_state.base_vertex);
                    return false;
                }

                // If there is a leaf with the same hash, load the leaf and test automorphism
                automorphism.reset();

                // Sometimes, a lab array is already stored for the leaf
                if(other_leaf->get_store_type() == ir::stored_leaf::STORE_LAB) {
                    const int* lab = other_leaf->get_lab_or_base();
                    automorphism.write_color_diff(local_state.c->vertex_to_col, lab);
                } else {
                    // If not, we need to use a more involved loading procedure which changes the local_state
                    memcpy(w.scratch_apply1.get_array(), local_state.c->vertex_to_col, g->v_size * sizeof(int));
                    load_state_from_leaf(g, local_state, root_save, other_leaf);
                    automorphism.write_color_diff(w.scratch_apply1.get_array(), local_state.c->lab);
                    used_load = true;
                }


                const bool cert = local_state.certify_automorphism(g, automorphism);

                if (cert) {
                    // We found an automorphism!
                    s_rolling_success = (9.0 * s_rolling_success + 1.0) / 10.0;

                    // Output automorphism
                    if(hook) (*hook)(0, automorphism.perm(), automorphism.nsupport(), automorphism.support());

                    // Sift into Schreier structure
                    bool sift = group.sift(w, automorphism, uniform);
                    automorphism.reset();

                    if(sift && h_sift_random && s_paths > h_sift_random_lim) {
                        int fail = 3;
                        while(fail >= 0) fail -= !group.sift_random(w, automorphism, generator);
                    }

                    automorphism.reset();
                    return sift;
                } else {
                    if(hash_offset > 0) std::cout << "cert fail " << other_leaf->get_store_type() << "/" << ir::stored_leaf::STORE_LAB << std::endl;
                    // If we used the more involved loading procedure, break for now
                    if(used_load) break; // TODO maybe there is better fix here...
                    continue;
                }
            }
            automorphism.reset();
            return false;
        }

    public:
        // stats
        double    s_rolling_success = 0;                  /**< rolling probability how many random paths succeed     */
        double    s_rolling_first_level_success  = 1.0;   /**< rolling probability how many random paths succeed on the
                                                            *  first level*/

        int       s_paths      = 0;                       /**< how many total paths have been computed */
        int       s_paths_fail = 0;                       /**< how many total paths failed             */

        // settings for heuristics
        int       h_leaf_limit = 0;                       /**< limit to how many leaves can be stored         */
        bool      h_look_close = false;                   /**< whether to use trace early out on first level  */
        const int h_hash_col_limit = 32;                  /**< limit for how many hash collisions are allowed */
        bool      h_sift_random     = true;
        int       h_sift_random_lim = 128;

        random_ir(const int domain_size) : automorphism(domain_size), w(domain_size) {}

        void setup(int error, int leaf_store_limit, bool look_close = false) {
            h_leaf_limit = leaf_store_limit;
            h_look_close = look_close;
        }

        void reset() {
            s_paths      = 0;
            s_paths_fail =
            h_leaf_limit = 0;
            s_rolling_success = 0;
            s_rolling_first_level_success  = 1.0;
        }

        bool h_almost_done(groups::shared_schreier &group) {
            return group.s_consecutive_success >= 1;
        }

        void
        specific_walk(sgraph *g, ir::shared_tree &ir_tree, ir::controller &local_state, std::vector<int> &base_vertex) {
            local_state.walk(g, *ir_tree.pick_node_from_level(0,0)->get_save(), base_vertex);
            auto other_leaf = ir_tree.stored_leaves.lookup_leaf(local_state.T->get_hash());
            if(other_leaf == nullptr) {
                ir_tree.stored_leaves.add_leaf(local_state.T->get_hash(), *local_state.c, local_state.base_vertex);
            }
        }

        // TODO I should have one random_walk method that is used by the methods below, i think
        // TODO random_walks_from_tree should be the only outward facing method? base_aligned should just happen automatically, or flag

        void random_walks(sgraph *g, dejavu_hook *hook, std::function<ir::type_selector_hook> *selector,
                          ir::shared_tree &ir_tree, groups::shared_schreier &group, ir::controller &local_state) {
            local_state.use_reversible(false);
            local_state.use_trace_early_out(false);
            ir::reduced_save my_own_save;

            ir::reduced_save* root_save  = ir_tree.pick_node_from_level(0,0)->get_save();
            ir::reduced_save* start_from = root_save;

            while(!group.probabilistic_abort_criterion() && !group.deterministic_abort_criterion()
                  && (ir_tree.stored_leaves.s_leaves <= h_leaf_limit || (h_almost_done(group) && ir_tree.stored_leaves.s_leaves <= 2*h_leaf_limit))) { // * s_rolling_first_level_success
                local_state.load_reduced_state(*start_from); // TODO can load more efficiently, this uses copy_force

                int could_start_from = group.finished_up_to_level();

                // can start from below the root if we finished Schreier table at the current root
                if(local_state.s_base_pos < could_start_from) {
                    while (local_state.s_base_pos <= could_start_from) {
                        local_state.move_to_child(g, group.base_point(local_state.s_base_pos));
                    }
                    assert(local_state.T->trace_equal());
                    local_state.save_reduced_state(my_own_save);
                    start_from = &my_own_save;
                }

                const int start_from_base_pos = local_state.s_base_pos;
                int base_pos                  = local_state.s_base_pos;

                // track whether current walk is base-aligned and/or uniform
                bool base_aligned = true;
                bool uniform      = true;

                //walk down the tree as long as we are not in a leaf
                while (g->v_size != local_state.c->cells) {
                    const int col = (*selector)(local_state.c, base_pos);
                    const int col_sz = local_state.c->ptn[col] + 1;
                    const int rand = ((int) generator()) % col_sz;
                    int v = local_state.c->lab[col + rand];

                    if(group.finished_up_to_level() + 1 == base_pos && group.is_in_base_orbit(base_pos, v)  && ir_tree.stored_leaves.s_leaves <= 1) {
                        heuristic_reroll.clear();
                        for(int i = 0; i < col_sz; ++i) {
                            heuristic_reroll.push_back(local_state.c->lab[col + i]);
                        }
                        group.reduce_to_unfinished(w, heuristic_reroll, base_pos);
                        if(!heuristic_reroll.empty()) {
                            const int rand = ((int) generator()) % heuristic_reroll.size();
                            v = heuristic_reroll[rand];
                        }
                    }

                    if(group.is_in_base_orbit(base_pos, v) && base_aligned && ir_tree.stored_leaves.s_leaves <= 1) {
                        v = group.base_point(local_state.s_base_pos);
                        assert(local_state.c->vertex_to_col[v] == col);
                        uniform = false;
                    } else {
                        base_aligned = false;
                    }
                    assert(local_state.c->vertex_to_col[v] == col);
                    local_state.move_to_child(g, v);

                    if(base_pos == start_from_base_pos) {
                        s_rolling_first_level_success =
                                (9.0 * s_rolling_first_level_success + (local_state.T->trace_equal())) / 10.0;
                    }

                    ++base_pos;
                }

                ++s_paths;

                add_leaf_to_storage_and_group(g, hook, group, ir_tree.stored_leaves, local_state, *root_save, uniform);
            }
        }

        // TODO implement dejavu strategy, more simple
        // TODO depends on ir_tree, selector, and given base (no need to sift beyond base!)
        // TODO: swap out ir_reduced to weighted IR shared_tree later? or just don't use automorphism pruning on BFS...?
        void random_walks_from_tree(sgraph *g, dejavu_hook *hook, std::function<ir::type_selector_hook> *selector,
                                    ir::shared_tree &ir_tree, groups::shared_schreier &group,
                                    ir::controller &local_state) {
            local_state.use_reversible(false);
            local_state.use_trace_early_out(false);
            s_rolling_first_level_success = 1;
            const int pick_from_level = ir_tree.get_finished_up_to();

            while(!group.probabilistic_abort_criterion() && !group.deterministic_abort_criterion()
                  && (ir_tree.stored_leaves.s_leaves <= h_leaf_limit || (h_almost_done(group) && ir_tree.stored_leaves.s_leaves <= 2*h_leaf_limit))
                  && (ir_tree.stored_leaves.s_leaves <= h_leaf_limit/4 || s_rolling_success > 0.001 || s_rolling_first_level_success > 0.1) // re-consider...
                  && s_rolling_first_level_success > 0.001) { //  * s_rolling_first_level_success

                auto node = ir_tree.pick_node_from_level(pick_from_level, (int) generator());
                local_state.load_reduced_state(*node->get_save());

                int base_pos                  = local_state.s_base_pos;
                const int start_from_base_pos = base_pos;

                while (g->v_size != local_state.c->cells) {
                    const int col    = (*selector)(local_state.c, base_pos);
                    const int col_sz = local_state.c->ptn[col] + 1;
                    const int rand = ((int) generator()) % col_sz;
                    int v = local_state.c->lab[col + rand];
                    local_state.use_trace_early_out((base_pos == start_from_base_pos) && !h_look_close);
                    local_state.move_to_child(g, v);

                    if(base_pos == start_from_base_pos) {
                        s_rolling_first_level_success =
                                (9.0 * s_rolling_first_level_success + (local_state.T->trace_equal())) / 10.0;
                        if(!h_look_close && !local_state.T->trace_equal()) break;
                    }
                    ++base_pos;
                }

                ++s_paths;
                if(base_pos == start_from_base_pos) {
                    ++s_paths_fail;
                    continue;
                }

                add_leaf_to_storage_and_group(g, hook, group, ir_tree.stored_leaves, local_state,
                                              *ir_tree.pick_node_from_level(0, 0)->get_save(), true);
            }
        }
    };
}

#endif //DEJAVU_RAND_H
