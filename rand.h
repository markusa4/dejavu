#ifndef DEJAVU_RAND_H
#define DEJAVU_RAND_H

#include "ir.h"
#include "groups.h"

namespace dejavu {
    namespace search_strategy {
        class stored_leaf {
            coloring store_c;
        public:
            stored_leaf(coloring& c, std::vector<int>& base) {
                // TODO make the stored leaf more properly
                store_c.copy(&c);
            }

            coloring* get_coloring() {
                return &store_c;
            }
        };

        class stored_leafs {
            std::mutex lock;
            std::unordered_multimap<long, stored_leaf*> leaf_store;
            std::vector<stored_leaf*> garbage_collector;

        public:
            int s_leaves = 0;

            ~stored_leafs() {
                for(int i = 0; i < garbage_collector.size(); ++i) {
                    delete garbage_collector[i];
                }
            }

            /**
             * Lookup whether a leaf with the given hash already exists.
             *
             * @param hash
             * @return
             */
            stored_leaf* lookup_leaf(long hash) {
                lock.lock();
                auto find = leaf_store.find(hash);
                if(find != leaf_store.end()) {
                    auto result = find->second;
                    lock.unlock();
                    return result;
                } else {
                    lock.unlock();
                    return nullptr;
                }
            }

            /**
             * Add leaf with the given hash. Does not add the leaf, if a leaf with the given hash already exists.
             *
             * @param hash
             * @param ptr
             */
            void add_leaf(long hash, coloring& c, std::vector<int>& base) {
                lock.lock();
                if(!leaf_store.contains(hash)) {
                    auto new_leaf = new stored_leaf(c, base);
                    leaf_store.insert(std::pair<long, stored_leaf*>(hash, new_leaf));
                    garbage_collector.push_back(new_leaf);
                    ++s_leaves;
                    lock.unlock();
                } else {
                    lock.unlock();
                }
            }

            void clear() {
                s_leaves = 0;
                leaf_store.clear();
            }
        };

        /**
         * \brief IR search using random walks.
         *
         * Performs random walks of the IR shared_tree, sifting resulting automorphisms into the given Schreier structure. If the
         * Schreier structure is complete with respect to the base, or the probabilistic abort criterion satisfied, the
         * process terminates. The algorithm guarantees to find all automorphisms up to the specified error bound.
         *
         * Alternatively, a limit for the amount of discovered differing leafs can be set.
         */
        class random_ir {
            stored_leafs leaf_storage; /**< stores all the leaves */

            std::default_random_engine generator; /**< random number generator */
            int consecutive_success = 0; /**< track repeated sifts for probabilistic abort criterion */
            int error_bound = 5;         /**< determines error probability */

            int    h_leaf_limit = 0;
            double h_rolling_success = 0;
            double h_rolling_first_level_success  = 1.0;
            double h_required_first_level_success = 0;

        public:
            void setup(int error, int leaf_store_limit, double required_first_level_success) {
                h_leaf_limit = leaf_store_limit;
                error_bound  = error;
                h_required_first_level_success = required_first_level_success;
            }

            void reset() {
                h_leaf_limit = 0;
                h_rolling_success = 0;
                h_rolling_first_level_success  = 1.0;
                h_required_first_level_success = 0;
            }

            void clear_leaves() {
                leaf_storage.clear();
            }

            double get_rolling_sucess_rate() {
                return h_rolling_success;
            }

            double get_rolling_first_level_success_rate() {
                return h_rolling_first_level_success;
            }

            /**
             * Records a sift result for the probabilistic abort criterion.
             *
             * @param changed Whether the sift was successful or not.
             */
            void record_sift_result(const bool changed) {
                if(!changed) {
                    ++consecutive_success;
                } else {
                    if(consecutive_success > 0) {
                        ++error_bound;
                        consecutive_success = 0;
                    }
                }
            }

            /**
             * @return Whether the probabilistic abort criterion allows termination or not.
             */
            bool probabilistic_abort_criterion() {
                return (consecutive_success > error_bound);
            }

            /**
             * @param group The Schreier structure of the group we are considering.
             * @return Whether the deterministic abort criterion allows termination or not.
             */
            bool deterministic_abort_criterion(groups::schreier &group) {
                return (group.finished_up_to_level() + 1 == group.base_size());
            }

            /**
             * @return How many leaves were stored during random search.
             */
            int stat_leaves() {
                return leaf_storage.s_leaves;
            }

            void specific_walk(refinement &R, std::vector<int>& base_vertex, sgraph *g,
                               groups::schreier &group, ir::controller &local_state, ir::reduced_save &start_from) {
                local_state.load_reduced_state(start_from);
                int base_pos = local_state.base_pos;

                while (g->v_size != local_state.c->cells) {
                    int v = base_vertex[local_state.base_pos];
                    local_state.move_to_child(&R, g, v);
                    ++base_pos;
                }

                auto other_leaf = leaf_storage.lookup_leaf(local_state.T->get_hash());
                if(other_leaf == nullptr) {
                    leaf_storage.add_leaf(local_state.T->get_hash(), *local_state.c, local_state.base_vertex);
                }
            }

            void random_walk(refinement &R, std::function<ir::type_selector_hook> *selector, sgraph *g,
                             groups::schreier &group, ir::controller &local_state, ir::reduced_save &start_from) {
                groups::automorphism_workspace automorphism(g->v_size);
                groups::schreier_workspace w(g->v_size, &R, g);
                std::vector<int> heuristic_reroll;

                local_state.load_reduced_state(start_from);
                int could_start_from = group.finished_up_to_level();
                if(local_state.base_pos < could_start_from) {
                    while (local_state.base_pos <= could_start_from) {
                        //std::cout << local_state.base_pos << ", " << group.base_point(local_state.base_pos) << std::endl;
                        local_state.move_to_child(&R, g, group.base_point(local_state.base_pos));
                    }
                    local_state.save_reduced_state(start_from);
                }

                int base_pos = local_state.base_pos;

                bool base_aligned = true;
                bool uniform = true;

                while (g->v_size != local_state.c->cells) {
                    const int col = (*selector)(local_state.c, base_pos);
                    const int col_sz = local_state.c->ptn[col] + 1;
                    const int rand = ((int) generator()) % col_sz;
                    int v = local_state.c->lab[col + rand];

                    if(group.finished_up_to_level() + 1 == base_pos && group.is_in_base_orbit(base_pos, v)  && leaf_storage.s_leaves <= 1) {
                        heuristic_reroll.clear();
                        for(int i = 0; i < col_sz; ++i) {
                            heuristic_reroll.push_back(local_state.c->lab[col + i]);
                        }
                        group.reduce_to_unfinished(w, heuristic_reroll, base_pos);
                        if(heuristic_reroll.size() > 0) {
                            const int rand = ((int) generator()) % heuristic_reroll.size();
                            v = heuristic_reroll[rand];
                            //std::cout << group.is_finished(base_pos) << "re-roll!"
                            //          << group.is_in_base_orbit(base_pos, v) << std::endl;
                        }
                    }

                    if(group.is_in_base_orbit(base_pos, v) && base_aligned && leaf_storage.s_leaves <= 1) {
                        v = group.base_point(local_state.base_pos);
                        assert(local_state.c->vertex_to_col[v] == col);
                        uniform = false;
                    } else {
                        base_aligned = false;
                    }

                    local_state.move_to_child(&R, g, v);
                    ++base_pos;
                }

                auto other_leaf = leaf_storage.lookup_leaf(local_state.T->get_hash());
                if(other_leaf == nullptr) {
                    //std::cout << "adding leaf " << local_state.T->get_hash() << std::endl;
                    leaf_storage.add_leaf(local_state.T->get_hash(), *local_state.c, local_state.base_vertex);
                } else {
                    //std::cout << "reading leaf " << local_state.T->get_hash() << std::endl;
                    automorphism.write_color_diff(local_state.c->vertex_to_col, other_leaf->get_coloring()->lab);
                    const bool cert = R.certify_automorphism_sparse(g, automorphism.perm(), automorphism.nsupport(), automorphism.support());
                    if(cert) {
                        //std::cout << "found automorphism, hash " << local_state.T->get_hash() << " support " << automorphism.nsupport() << std::endl;
                        const bool sift = group.sift(w, g, &R, automorphism);
                        if(uniform) record_sift_result(sift);
                    }

                    automorphism.reset();
                }
            }

            // TODO implement dejavu strategy, more simple
            // TODO depends on ir_tree, selector, and given base (no need to sift beyond base!)
            // TODO: swap out ir_reduced to weighted IR shared_tree later? or just don't use automorphism pruning on BFS...?
            void random_walks(refinement &R, std::function<ir::type_selector_hook> *selector, sgraph *g,
                              groups::schreier &group, ir::controller &local_state, ir::reduced_save &start_from) {
                groups::automorphism_workspace automorphism(g->v_size);
                groups::schreier_workspace w(g->v_size, &R, g);
                std::vector<int> heuristic_reroll;

                local_state.use_trace_early_out(false);

                while(!probabilistic_abort_criterion() && !deterministic_abort_criterion(group)
                      && leaf_storage.s_leaves <= h_leaf_limit) { // * h_rolling_first_level_success
                    local_state.load_reduced_state(start_from);
                    int could_start_from = group.finished_up_to_level();
                    if(local_state.base_pos < could_start_from) {
                        while (local_state.base_pos <= could_start_from) {
                            //std::cout << local_state.base_pos << ", " << group.base_point(local_state.base_pos) << std::endl;
                            local_state.move_to_child(&R, g, group.base_point(local_state.base_pos));
                        }
                        local_state.save_reduced_state(start_from);
                    }

                    const int start_from_base_pos = local_state.base_pos;
                    int base_pos                  = local_state.base_pos;

                    bool base_aligned = true;
                    bool uniform = true;

                    while (g->v_size != local_state.c->cells) {
                        const int col = (*selector)(local_state.c, base_pos);
                        const int col_sz = local_state.c->ptn[col] + 1;
                        const int rand = ((int) generator()) % col_sz;
                        int v = local_state.c->lab[col + rand];

                        if(group.finished_up_to_level() + 1 == base_pos && group.is_in_base_orbit(base_pos, v)  && leaf_storage.s_leaves <= 1) {
                            heuristic_reroll.clear();
                            for(int i = 0; i < col_sz; ++i) {
                                heuristic_reroll.push_back(local_state.c->lab[col + i]);
                            }
                            group.reduce_to_unfinished(w, heuristic_reroll, base_pos);
                            if(heuristic_reroll.size() > 0) {
                                const int rand = ((int) generator()) % heuristic_reroll.size();
                                v = heuristic_reroll[rand];
                                //std::cout << group.is_finished(base_pos) << "re-roll!"
                                //          << group.is_in_base_orbit(base_pos, v) << std::endl;
                            }
                        }

                        if(group.is_in_base_orbit(base_pos, v) && base_aligned && leaf_storage.s_leaves <= 1) {
                            v = group.base_point(local_state.base_pos);
                            assert(local_state.c->vertex_to_col[v] == col);
                            uniform = false;
                        } else {
                            base_aligned = false;
                        }

                        local_state.move_to_child(&R, g, v);

                        if(base_pos == start_from_base_pos) {
                            h_rolling_first_level_success =
                                    (9.0 * h_rolling_first_level_success + (local_state.T->trace_equal())) / 10.0;
                        }

                        ++base_pos;
                    }

                    // local_state.write_strong_invariant(g);

                    auto other_leaf = leaf_storage.lookup_leaf(local_state.T->get_hash());
                    if(other_leaf == nullptr) {
                        //std::cout << "adding leaf " << local_state.T->get_hash() << std::endl;
                        leaf_storage.add_leaf(local_state.T->get_hash(), *local_state.c, local_state.base_vertex);
                        h_rolling_success = (9.0*h_rolling_success + 0.0) / 10.0;
                    } else {
                        //std::cout << "reading leaf " << local_state.T->get_hash() << std::endl;
                        automorphism.write_color_diff(local_state.c->vertex_to_col, other_leaf->get_coloring()->lab);
                        const bool cert = R.certify_automorphism_sparse(g, automorphism.perm(), automorphism.nsupport(), automorphism.support());
                        if(cert) {
                            h_rolling_success = (9.0*h_rolling_success + 1.0) / 10.0;
                            //std::cout << "found automorphism, hash " << local_state.T->get_hash() << " support " << automorphism.nsupport() << std::endl;
                            const bool sift = group.sift(w, g, &R, automorphism);
                            if(uniform) record_sift_result(sift);
                        } else {
                            std::cout << "cert fail " << std::endl;
                        }
                        automorphism.reset();
                    }
                }
            }

            // TODO implement dejavu strategy, more simple
            // TODO depends on ir_tree, selector, and given base (no need to sift beyond base!)
            // TODO: swap out ir_reduced to weighted IR shared_tree later? or just don't use automorphism pruning on BFS...?
            void __attribute__ ((noinline)) random_walks_from_tree(refinement &R, std::function<ir::type_selector_hook> *selector, sgraph *g,
                                                                   groups::schreier &group, ir::controller &local_state, ir::shared_tree &ir_tree) {
                groups::automorphism_workspace automorphism(g->v_size);
                groups::schreier_workspace w(g->v_size, &R, g);
                std::vector<int> heuristic_reroll;
                local_state.use_trace_early_out(false);

                h_rolling_first_level_success = 1;
                const int pick_from_level = ir_tree.get_finished_up_to();

                while(!probabilistic_abort_criterion() && !deterministic_abort_criterion(group)
                      && leaf_storage.s_leaves <= h_leaf_limit
                      && (leaf_storage.s_leaves <= h_leaf_limit/4 || h_rolling_success > 0.001 || h_rolling_first_level_success > 0.1) // re-consider...
                        ) { //  * h_rolling_first_level_success
                    auto node = ir_tree.pick_node_from_level(pick_from_level, (int) generator());
                    local_state.load_reduced_state(*node->get_save());

                    long started_from_hash = local_state.T->get_hash();

                    int base_pos                  = local_state.base_pos;
                    const int start_from_base_pos = base_pos;

                    bool uniform = true;

                    while (g->v_size != local_state.c->cells) {
                        const int col    = (*selector)(local_state.c, base_pos);
                        const int col_sz = local_state.c->ptn[col] + 1;
                        const int rand = ((int) generator()) % col_sz;
                        int v = local_state.c->lab[col + rand];
                        local_state.move_to_child(&R, g, v);

                        if(base_pos == start_from_base_pos) {
                            h_rolling_first_level_success =
                                    (9.0 * h_rolling_first_level_success + (local_state.T->trace_equal())) / 10.0;
                        }
                        ++base_pos;
                    }

                    //local_state.write_strong_invariant(g);

                    for(int hashf = 0; hashf < 32; ++hashf) {
                        long hash_c = local_state.T->get_hash()+hashf;
                        auto other_leaf = leaf_storage.lookup_leaf(hash_c);
                        if (other_leaf == nullptr) {
                            h_rolling_success = (9.0 * h_rolling_success + 0.0) / 10.0;
                            leaf_storage.add_leaf(hash_c, *local_state.c, local_state.base_vertex);
                            break;
                        } else {
                            automorphism.reset();
                            for (int i = 0; i < g->v_size; ++i) {
                                assert(i == automorphism.perm()[i]);
                            }
                            automorphism.write_color_diff(local_state.c->vertex_to_col,
                                                          other_leaf->get_coloring()->lab);
                            const bool cert = R.certify_automorphism_sparse(g, automorphism.perm(),
                                                                            automorphism.nsupport(),
                                                                            automorphism.support());
                            if (cert) {
                                //std::cout << "cert success" << hash_c << std::endl;
                                h_rolling_success = (9.0 * h_rolling_success + 1.0) / 10.0;
                                //std::cout << "found automorphism, hash " << local_state.T->get_hash() << " support " << automorphism.nsupport() << std::endl;
                                const bool sift = group.sift(w, g, &R, automorphism);
                                if (uniform) record_sift_result(sift);
                                automorphism.reset();
                                break;
                            } else {
                                std::cout << "cert fail " << std::endl;
                                //std::cout << "cert fail" << hash_c << " / " << local_state.T->trace_equal()  << "/" << started_from_hash << std::endl;
                                automorphism.reset();
                                continue;
                            }
                        }
                    }
                }
            }
        };
    }
}

#endif //DEJAVU_RAND_H
