#ifndef DEJAVU_RAND_H
#define DEJAVU_RAND_H

#include "ir.h"
#include "groups.h"

namespace dejavu {
    namespace search_strategy {
        // TODO: this should be a method of local_state?
        static void specific_walk(ir::refinement &R, std::vector<int>& base_vertex, sgraph *g,
                           ir::controller &local_state, ir::reduced_save &start_from) {
            local_state.load_reduced_state(start_from);

            while (g->v_size != local_state.c->cells) {
                int v = base_vertex[local_state.base_pos];
                local_state.move_to_child(&R, g, v);
            }
        }


        class stored_leaf {
            work_list lab;
            std::vector<int> base;
            bool full_save;
        public:
            stored_leaf(coloring& c, std::vector<int>& base, bool full_save = true) {
                // TODO make the stored leaf more properly
                this->full_save = full_save;
                if(full_save) {
                    //store_c.copy(&c);
                    lab.initialize(c.lab_sz);
                    memcpy(lab.get_array(), c.lab, c.lab_sz * sizeof(int));
                } else {
                    this->base = base;
                }
            }

            void load_state(ir::refinement &R, sgraph *g, ir::controller &local_state, ir::reduced_save &start_from) {
                specific_walk(R, base, g, local_state, start_from);
            }

            /*coloring* get_coloring_if_possible() {
                if(full_save) return &store_c;
                return nullptr;
            }*/

            int* get_lab_if_possible() {
                if(full_save) return lab.get_array();
                return nullptr;
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
                    auto new_leaf = new stored_leaf(c, base, s_leaves < 100);
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
            int    h_paths = 0;
            bool   h_look_close = false;

        public:
            void setup(int error, int leaf_store_limit, double required_first_level_success, bool look_close = false) {
                h_leaf_limit = leaf_store_limit;
                error_bound  = error;
                h_required_first_level_success = required_first_level_success;
                h_look_close = look_close;
            }

            void reset() {
                h_paths      = 0;
                h_leaf_limit = 0;
                h_rolling_success = 0;
                h_rolling_first_level_success  = 1.0;
                h_required_first_level_success = 0;
            }

            void reset_prob() {
                consecutive_success = 0;
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

            bool h_almost_done() {
                return consecutive_success >= 1;
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
            bool deterministic_abort_criterion(groups::shared_schreier &group) {
                return (group.finished_up_to_level() + 1 == group.base_size());
            }

            /**
             * @return How many leaves were stored during random search.
             */
            int stat_leaves() {
                return leaf_storage.s_leaves;
            }

            void specific_walk(ir::refinement &R, std::vector<int>& base_vertex, sgraph *g,
                               groups::shared_schreier &group, ir::controller &local_state, ir::reduced_save &start_from) {
                local_state.load_reduced_state(start_from);

                while (g->v_size != local_state.c->cells) {
                    //std::cout << local_state.T->get_position() <<  ", " << local_state.c->cells << ", " << local_state.T->trace_equal() << std::endl;
                    int v = base_vertex[local_state.base_pos];
                    local_state.move_to_child(&R, g, v);
                    assert(local_state.T->trace_equal());
                }

                local_state.write_strong_invariant(g);

                auto other_leaf = leaf_storage.lookup_leaf(local_state.T->get_hash());
                if(other_leaf == nullptr) {
                    //std::cout << "added" << local_state.T->get_hash() << std::endl;
                    leaf_storage.add_leaf(local_state.T->get_hash(), *local_state.c, local_state.base_vertex);
                }
            }

            // TODO implement dejavu strategy, more simple
            // TODO depends on ir_tree, selector, and given base (no need to sift beyond base!)
            // TODO: swap out ir_reduced to weighted IR shared_tree later? or just don't use automorphism pruning on BFS...?
            void random_walks(ir::refinement &R, std::function<ir::type_selector_hook> *selector, sgraph *g,
                              groups::shared_schreier &group, ir::controller &local_state, ir::reduced_save* start_from) {
                groups::automorphism_workspace automorphism(g->v_size);
                groups::schreier_workspace w(g->v_size);
                std::vector<int> heuristic_reroll;

                local_state.use_trace_early_out(false);
                ir::reduced_save my_own_save;

                while(!probabilistic_abort_criterion() && !deterministic_abort_criterion(group)
                      && (leaf_storage.s_leaves <= h_leaf_limit || (h_almost_done() && leaf_storage.s_leaves <= 2*h_leaf_limit))) { // * h_rolling_first_level_success

                    local_state.load_reduced_state(*start_from);

                    int could_start_from = group.finished_up_to_level();
                    //int could_start_from = 0;
                    if(local_state.base_pos < could_start_from) {
                        while (local_state.base_pos <= could_start_from) {
                            //std::cout << local_state.base_pos << ", " << group.base_point(local_state.base_pos) << std::endl;
                            local_state.move_to_child(&R, g, group.base_point(local_state.base_pos));
                        }
                        assert(local_state.T->trace_equal());
                        local_state.save_reduced_state(my_own_save);
                        start_from = &my_own_save;
                    }

                    const int start_from_base_pos = local_state.base_pos;
                    int base_pos                  = local_state.base_pos;

                    bool base_aligned = true;
                    bool uniform = true;

                    while (g->v_size != local_state.c->cells) {
                        //std::cout << local_state.T->get_position() <<  ", " << local_state.c->cells << ", " << local_state.T->trace_equal() << std::endl;
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
                            if(!heuristic_reroll.empty()) {
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
                        assert(local_state.c->vertex_to_col[v] == col);
                        local_state.move_to_child(&R, g, v);

                        if(base_pos == start_from_base_pos) {
                            h_rolling_first_level_success =
                                    (9.0 * h_rolling_first_level_success + (local_state.T->trace_equal())) / 10.0;
                        }

                        ++base_pos;
                    }

                    ++h_paths;

                    local_state.write_strong_invariant(g);

                    auto other_leaf = leaf_storage.lookup_leaf(local_state.T->get_hash());
                    if(other_leaf == nullptr) {
                        //std::cout << "added" << local_state.T->get_hash() << std::endl;
                        //std::cout << "adding leaf " << local_state.T->get_hash() << std::endl;
                        leaf_storage.add_leaf(local_state.T->get_hash(), *local_state.c, local_state.base_vertex);
                        h_rolling_success = (9.0*h_rolling_success + 0.0) / 10.0;
                    } else {
                        //std::cout << "reading leaf " << local_state.T->get_hash() << std::endl;
                        int* lab = other_leaf->get_lab_if_possible();
                        if(lab != nullptr) {
                            automorphism.write_color_diff(local_state.c->vertex_to_col, lab);
                        } else {
                            memcpy(w.scratch_apply1.get_array(), local_state.c->vertex_to_col, g->v_size * sizeof(int));
                            other_leaf->load_state(R, g, local_state, *start_from);
                            automorphism.write_color_diff(w.scratch_apply1.get_array(), local_state.c->lab);
                        }
                        const bool cert = R.certify_automorphism_sparse(g, automorphism.perm(), automorphism.nsupport(), automorphism.support());
                        if(cert) {
                            //std::cout << "certified" << local_state.T->get_hash() << std::endl;
                            h_rolling_success = (9.0*h_rolling_success + 1.0) / 10.0;
                            //std::cout << "found automorphism, hash " << local_state.T->get_hash() << " support " << automorphism.nsupport() << std::endl;
                            const bool sift = group.sift(w, automorphism);
                            if(uniform) record_sift_result(sift);
                        } else {
                            //std::cout << "testing" << local_state.T->get_hash() << std::endl;
                            std::cout << "cert fail " << std::endl;
                            assert(false);
                        }
                        automorphism.reset();
                    }
                }
            }

            // TODO implement dejavu strategy, more simple
            // TODO depends on ir_tree, selector, and given base (no need to sift beyond base!)
            // TODO: swap out ir_reduced to weighted IR shared_tree later? or just don't use automorphism pruning on BFS...?
            void __attribute__ ((noinline)) random_walks_from_tree(ir::refinement &R, std::function<ir::type_selector_hook> *selector, sgraph *g,
                                                                   groups::shared_schreier &group, ir::controller &local_state, ir::shared_tree &ir_tree) {
                groups::automorphism_workspace automorphism(g->v_size);
                groups::schreier_workspace w(g->v_size);
                std::vector<int> heuristic_reroll;
                local_state.use_trace_early_out(false);

                h_rolling_first_level_success = 1;
                const int pick_from_level = ir_tree.get_finished_up_to();

                while(!probabilistic_abort_criterion() && !deterministic_abort_criterion(group)
                      && (leaf_storage.s_leaves <= h_leaf_limit || (h_almost_done() && leaf_storage.s_leaves <= 2*h_leaf_limit))
                      && (leaf_storage.s_leaves <= h_leaf_limit/4 || h_rolling_success > 0.001 || h_rolling_first_level_success > 0.1) // re-consider...
                      && h_rolling_first_level_success > 0.001) { //  * h_rolling_first_level_success
                    auto node = ir_tree.pick_node_from_level(pick_from_level, (int) generator());
                    local_state.load_reduced_state(*node->get_save());

                    int base_pos                  = local_state.base_pos;
                    const int start_from_base_pos = base_pos;

                    while (g->v_size != local_state.c->cells) {
                        const int col    = (*selector)(local_state.c, base_pos);
                        const int col_sz = local_state.c->ptn[col] + 1;
                        const int rand = ((int) generator()) % col_sz;
                        int v = local_state.c->lab[col + rand];
                        local_state.use_trace_early_out((base_pos == start_from_base_pos) && !h_look_close);
                        local_state.move_to_child(&R, g, v);

                        if(base_pos == start_from_base_pos) {
                            h_rolling_first_level_success =
                                    (9.0 * h_rolling_first_level_success + (local_state.T->trace_equal())) / 10.0;
                            if(!h_look_close && !local_state.T->trace_equal()) break;
                        }
                        ++base_pos;
                    }

                    ++h_paths;

                    if(base_pos == start_from_base_pos) {continue;}

                    local_state.write_strong_invariant(g);

                    bool used_load = false;

                    for(int hashf = 0; hashf < 32; ++hashf) {
                        long hash_c = local_state.T->get_hash()+hashf;
                        auto other_leaf = leaf_storage.lookup_leaf(hash_c);
                        if (other_leaf == nullptr) {
                            h_rolling_success = (9.0 * h_rolling_success + 0.0) / 10.0;
                            leaf_storage.add_leaf(hash_c, *local_state.c, local_state.base_vertex);
                            //std::cout << "add leaf" << std::endl;
                            break;
                        } else {
                            automorphism.reset();
                            for (int i = 0; i < g->v_size; ++i) {
                                assert(i == automorphism.perm()[i]);
                            }
                            int* lab = other_leaf->get_lab_if_possible();
                            if(lab != nullptr) {
                                automorphism.write_color_diff(local_state.c->vertex_to_col, lab);
                            } else {
                                memcpy(w.scratch_apply1.get_array(), local_state.c->vertex_to_col, g->v_size * sizeof(int));
                                other_leaf->load_state(R, g, local_state, *ir_tree.pick_node_from_level(0,0)->get_save());
                                automorphism.write_color_diff(w.scratch_apply1.get_array(), local_state.c->lab);
                                used_load = true;
                            }

                            const bool cert = R.certify_automorphism_sparse(g, automorphism.perm(),
                                                                            automorphism.nsupport(),
                                                                            automorphism.support());
                            if (cert) {
                                //std::cout << "cert success" << hash_c << std::endl;
                                h_rolling_success = (9.0 * h_rolling_success + 1.0) / 10.0;
                                //std::cout << "found automorphism, hash " << local_state.T->get_hash() << " support " << automorphism.nsupport() << std::endl;
                                const bool sift = group.sift(w, automorphism);
                                record_sift_result(sift);
                                automorphism.reset();
                                break;
                            } else {
                                std::cout << "cert fail " << std::endl;
                                //std::cout << "cert fail" << hash_c << " / " << local_state.T->trace_equal()  << "/" << started_from_hash << std::endl;
                                automorphism.reset();
                                if(used_load) break; // TODO maybe there is better fix here...
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
