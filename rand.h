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
    private:
        std::default_random_engine generator; /**< random number generator */
        std::vector<int> heuristic_reroll;

        groups::schreier_workspace*     gws_schreierw;
        groups::automorphism_workspace* gws_automorphism;

        /**
         * Loads the leaf into the \p local_state. Only works if base of leaf was stored, i.e., whenever
         * `full_save = false` was used in the constructor.
         *
         * @param leaf The leaf.
         * @param g The graph.
         * @param local_state Local state in which the leaf will be stored.
         * @param start_from State from which the walk to the leaf is performed.
         */
        static void load_state_from_leaf(sgraph *g, ir::controller &local_state, ir::limited_save &start_from,
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
                                           ir::limited_save &root_save, bool uniform) {
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
                    ++s_leaves;
                    ++s_paths_failany;
                    s_rolling_success = (9.0 * s_rolling_success + 0.0) / 10.0;
                    leaf_storage.add_leaf(hash_c, *local_state.c, local_state.base_vertex);
                    //return false;
                    break;
                }

                // If there is a leaf with the same hash, load the leaf and test automorphism
                gws_automorphism->reset();

                // Sometimes, a lab array is already stored for the leaf
                if(other_leaf->get_store_type() == ir::stored_leaf::STORE_LAB) {
                    const int* lab = other_leaf->get_lab_or_base();
                    gws_automorphism->write_color_diff(local_state.c->vertex_to_col, lab);
                } else {
                    // If not, we need to use a more involved loading procedure which changes the local_state
                    memcpy(gws_schreierw->scratch_apply1.get_array(), local_state.c->vertex_to_col, g->v_size * sizeof(int));
                    load_state_from_leaf(g, local_state, root_save, other_leaf);
                    gws_automorphism->write_color_diff(gws_schreierw->scratch_apply1.get_array(), local_state.c->lab);
                    used_load = true;
                }


                const bool cert = local_state.certify(g, *gws_automorphism);

                if (cert) {
                    // We found an automorphism!
                    s_rolling_success = (9.0 * s_rolling_success + 1.0) / 10.0;
                    ++s_succeed;

                    // Output automorphism
                    if(hook) (*hook)(0, gws_automorphism->perm(), gws_automorphism->nsupport(), gws_automorphism->support());

                    // Sift into Schreier structure
                    bool sift = group.sift(*gws_schreierw, *gws_automorphism, uniform);
                    gws_automorphism->reset();

                    if(sift && h_sift_random && s_paths > h_sift_random_lim) {
                        int fail = 3;
                        while(fail >= 0) fail -= !group.sift_random(*gws_schreierw, *gws_automorphism, generator);
                    }

                    gws_automorphism->reset();
                    return sift;
                } else {
                    if(hash_offset > 0) std::cout << "cert fail " << other_leaf->get_store_type() << "/" << ir::stored_leaf::STORE_LAB << std::endl;
                    // If we used the more involved loading procedure, break for now
                    if(used_load) {
                        ++s_paths_failany;
                        break; // TODO maybe there is better fix here...
                    }
                    continue;
                }
            }
            gws_automorphism->reset();
            return false;
        }

    public:
        // stats
        double    s_rolling_success = 0;                  /**< rolling probability how many random paths succeed     */
        double    s_rolling_first_level_success  = 1.0;   /**< rolling probability how many random paths succeed on the
                                                            *  first level*/

        long      s_trace_cost1   = 0;                    /**< total cost incurred on first level         */
        int       s_paths         = 0;                    /**< how many total paths have been computed    */
        int       s_paths_fail1   = 0;                    /**< how many total paths failed on first level */
        int       s_paths_failany = 0;                    /**< how many total paths have failed           */
        int       s_succeed       = 0;                    /**< how many total paths have succeeded        */
        int       s_leaves        = 0;                    /**< how many leaves were added                 */

        // settings for heuristics
        bool      h_look_close      = false;              /**< whether to use trace early out on first level   */
        const int h_hash_col_limit  = 32;                 /**< limit for how many hash collisions are allowed  */
        bool      h_sift_random     = true;               /**< sift random elements into Schreier structure    */
        int       h_sift_random_lim = 8;                  /**< after how many paths random elements are sifted */

        void setup(bool look_close = false) {
            h_look_close = look_close;
        }

        void link_to_workspace(groups::schreier_workspace* schreier, groups::automorphism_workspace* automorphism) {
            gws_automorphism = automorphism;
            gws_schreierw    = schreier;
        }

        void reset_statistics() {
            s_paths            = 0;
            s_paths_fail1      = 0;
            s_trace_cost1      = 0;
            s_paths_failany    = 0;
            s_leaves           = 0;
            s_succeed          = 0;
            s_rolling_success = 0;
            s_rolling_first_level_success  = 1.0;
        }

        static bool h_almost_done(groups::shared_schreier &group) {
            return group.s_consecutive_success >= 1;
        }

        static void
        specific_walk(sgraph *g, ir::shared_tree &ir_tree, ir::controller &local_state, std::vector<int> &base_vertex) {
            local_state.walk(g, *ir_tree.pick_node_from_level(0,0)->get_save(), base_vertex);
            auto other_leaf = ir_tree.stored_leaves.lookup_leaf(local_state.T->get_hash());
            if(other_leaf == nullptr) {
                ir_tree.stored_leaves.add_leaf(local_state.T->get_hash(), *local_state.c, local_state.base_vertex);
            }
        }

        // TODO setting to escape base aligned search for more preprocessing? could even return skiplevel-coloring
        // TODO preprocess using that, reset schreier structure, continue... ("in-in-processing" lol)
        void random_walks(sgraph *g, dejavu_hook *hook, std::function<ir::type_selector_hook> *selector,
                          ir::shared_tree &ir_tree, groups::shared_schreier &group, ir::controller &local_state,
                          int fail_limit) {
            local_state.use_reversible(false);
            local_state.use_trace_early_out(false);
            ir::limited_save my_own_save;

            ir::limited_save* root_save  = ir_tree.pick_node_from_level(0, 0)->get_save();
            ir::limited_save* start_from = root_save;

            while(!group.probabilistic_abort_criterion() && !group.deterministic_abort_criterion() &&
                    s_paths_failany < fail_limit) {
                local_state.load_reduced_state(*start_from); // TODO can load more efficiently, this uses copy_force

                int could_start_from = group.finished_up_to_level();

                // can start from below the root if we finished Schreier table at the current root
                if(local_state.s_base_pos < could_start_from) {
                    while (local_state.s_base_pos <= could_start_from) {
                        local_state.move_to_child(g, group.base_point(local_state.s_base_pos));
                    }
                    assert(local_state.T->trace_equal());
                    local_state.save_reduced_state(my_own_save); // from now on, we start from this save!
                    start_from = &my_own_save;
                    //std::cout << "start from: " << local_state.s_base_pos << ", cells: " << local_state.c->cells << std::endl;
                }

                const int start_from_base_pos = local_state.s_base_pos;
                int base_pos                  = local_state.s_base_pos;

                // track whether current walk is base-aware and/or uniform
                bool base_aligned = true;
                bool uniform      = true;

                //walk down the tree as long as we are not in a leaf
                while (g->v_size != local_state.c->cells) {
                    const int col = (*selector)(local_state.c, base_pos);
                    const int col_sz = local_state.c->ptn[col] + 1;
                    const int rand = ((int) generator()) % col_sz;
                    int v = local_state.c->lab[col + rand];

                    // base-aware search: if we are still walking along the base, and the vertex we picked is in the
                    // same orbit as the base -- we might as well keep walking on the base, or choose a different vertex
                    if(group.finished_up_to_level() + 1 == base_pos && group.is_in_base_orbit(base_pos, v)
                         && ir_tree.stored_leaves.s_leaves <= 1) {
                        heuristic_reroll.clear();
                        for(int i = 0; i < col_sz; ++i) {
                            heuristic_reroll.push_back(local_state.c->lab[col + i]);
                        }
                        group.reduce_to_unfinished(*gws_schreierw, heuristic_reroll, base_pos);
                        if(!heuristic_reroll.empty()) {
                            const int rand = ((int) generator()) % heuristic_reroll.size();
                            v = heuristic_reroll[rand];
                        }
                    }

                    if(group.is_in_base_orbit(base_pos, v) && base_aligned && ir_tree.stored_leaves.s_leaves <= 1) {
                        v = group.base_point(local_state.s_base_pos);
                        assert(local_state.c->vertex_to_col[v] == col);
                        uniform = false; // not uniform anymore!
                    } else {
                        base_aligned = false; // we stopped walking along the base
                    }
                    assert(local_state.c->vertex_to_col[v] == col);
                    const int trace_pos_pre = local_state.T->get_position();
                    local_state.use_trace_early_out((base_pos == start_from_base_pos) && !h_look_close);
                    local_state.move_to_child(g, v);

                    // keep track of some statistics for the first individualization (these statistics are used for
                    // decisions concerning breadth-first search)
                    if(base_pos == start_from_base_pos) {
                        s_trace_cost1 += local_state.T->get_position() - trace_pos_pre;
                        s_paths_fail1 += !local_state.T->trace_equal();
                        s_rolling_first_level_success =
                                (9.0 * s_rolling_first_level_success + (local_state.T->trace_equal())) / 10.0;
                        if(!h_look_close && !local_state.T->trace_equal()) break;
                    }

                    ++base_pos;
                }

                ++s_paths;
                if(base_pos == start_from_base_pos) { // did not arrive in a leaf
                    ++s_paths_failany;
                    continue;
                }

                // we arrived at a leaf... let's check whether we already have an equivalent leaf to form an
                // automorphism -- or add this leaf to the storage
                add_leaf_to_storage_and_group(g, hook, group, ir_tree.stored_leaves, local_state,
                                              *root_save, uniform);
            }
        }

        // TODO implement dejavu strategy, more simple
        // TODO depends on ir_tree, selector, and given base (no need to sift beyond base!)
        // TODO: swap out ir_reduced to weighted IR shared_tree later? or just don't use automorphism pruning on BFS...?
        void random_walks_from_tree(sgraph *g, dejavu_hook *hook, std::function<ir::type_selector_hook> *selector,
                                    ir::shared_tree &ir_tree, groups::shared_schreier &group,
                                    ir::controller &local_state, int fail_limit) {
            local_state.use_reversible(false);
            local_state.use_trace_early_out(false);
            s_rolling_first_level_success = 1;
            const int pick_from_level = ir_tree.get_finished_up_to();

            while(!group.probabilistic_abort_criterion() && !group.deterministic_abort_criterion() &&
                    s_paths_failany < fail_limit) {
                auto node = ir_tree.pick_node_from_level(pick_from_level, (int) generator());
                local_state.load_reduced_state(*node->get_save());

                int base_pos                  = local_state.s_base_pos;
                const int start_from_base_pos = base_pos;

                while (g->v_size != local_state.c->cells) {
                    const int col    = (*selector)(local_state.c, base_pos);
                    const int col_sz = local_state.c->ptn[col] + 1;
                    const int rand = ((int) generator()) % col_sz;
                    int v = local_state.c->lab[col + rand];
                    const int trace_pos_pre = local_state.T->get_position();
                    local_state.use_trace_early_out((base_pos == start_from_base_pos) && !h_look_close);
                    local_state.move_to_child(g, v);

                    if(base_pos == start_from_base_pos) {
                        s_trace_cost1 += local_state.T->get_position() - trace_pos_pre;
                        s_paths_fail1 += !local_state.T->trace_equal();
                        s_rolling_first_level_success =
                                (9.0 * s_rolling_first_level_success + (local_state.T->trace_equal())) / 10.0;
                        if(!h_look_close && !local_state.T->trace_equal()) break;
                    }
                    ++base_pos;
                }

                ++s_paths;
                if(base_pos == start_from_base_pos) {
                    ++s_paths_failany;
                    continue;
                }

                add_leaf_to_storage_and_group(g, hook, group, ir_tree.stored_leaves, local_state,
                                              *ir_tree.pick_node_from_level(0, 0)->get_save(), true);
            }
        }
    };
}

#endif //DEJAVU_RAND_H
