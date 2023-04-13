#ifndef DEJAVU_DFS_H
#define DEJAVU_DFS_H

#include <random>
#include <chrono>
#include "refinement.h"
#include "bijection.h"
#include "coloring.h"
#include "sgraph.h"
#include "trace.h"
#include "groups.h"
#include "ir.h"

namespace dejavu {

    static void progress_print_header() {

        PRINT("________________________________________________________________");
        PRINT(std::setw(16) << std::left <<"T (ms)"                                  << std::setw(16) << "proc"  << std::setw(16) << "P1"        << std::setw(16)        << "P2");
        PRINT("________________________________________________________________");
        PRINT(std::setw(16) << std::left << 0 << std::setw(16) << "start" << std::setw(16) << "_" << std::setw(16) << "_" );
    }

    static void progress_print(const std::string
    proc, const std::string p1, const std::string p2) {
        static std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
        PRINT(std::setw(16) << std::left << (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count()) / 1000000.0  << std::setw(16) << proc << std::setw(16) << p1 << std::setw(16) << p2);
    }

    /**
     * \brief High-level IR search strategies.
     *
     * Contains different high-level individualization-refinement search strategies. In particular, depth-first search,
     * breadth-first search as well as different methods of random search.
     */
    namespace search_strategy {

        typedef std::pair<ir::tree_node*, int> bfs_todo;

        class bfs_queue {
        private:
            std::mutex lock;
            std::vector<bfs_todo> queue;
        public:
            void add(ir::tree_node* base, int v) {
                lock.lock();
                queue.emplace_back(base, v);
                lock.unlock();
            }

            void reserve(int n) {
                lock.lock();
                queue.reserve(n);
                lock.unlock();
            }

            bool empty() {
                return queue.empty();
            }

            bfs_todo pop() {
                lock.lock();
                auto element = queue.back();
                queue.pop_back();
                lock.unlock();
                return element;
            }

        };

        /**
         * \brief Breadth-first search.
         */
        class bfs_ir {
            search_strategy::bfs_queue bfs_todo;
        public:
            void do_a_level(refinement* R, sgraph* g, ir::tree& ir_tree, ir::controller& local_state, std::function<ir::type_selector_hook> *selector) {
                int current_level = ir_tree.get_finished_up_to();
                queue_up_level(bfs_todo, selector, ir_tree, current_level);
                work_on_todo(bfs_todo, R, g, &ir_tree, local_state);
                ir_tree.set_finished_up_to(current_level + 1);
            }

            int next_level_estimate(ir::tree& ir_tree, std::function<ir::type_selector_hook> *selector) {
                const int base_pos = ir_tree.get_finished_up_to();
                const auto start_node = ir_tree.get_level(base_pos);
                const auto level_size = ir_tree.get_level_size(base_pos);
                auto next_node_save = start_node->get_save();
                auto c = next_node_save->get_coloring();
                auto base_pos_   = next_node_save->get_base_position();
                int col = (*selector)(c, base_pos_);
                return level_size * (c->ptn[col] + 1);
            }

            void queue_up_level(bfs_queue& queue, std::function<ir::type_selector_hook> *selector, ir::tree& ir_tree, int base_pos) {
                const auto start_node = ir_tree.get_level(base_pos);
                const auto level_size = ir_tree.get_level_size(base_pos);
                auto next_node = start_node;
                bool reserve = false;

                do {
                    auto next_node_save = next_node->get_save();
                    auto c = next_node_save->get_coloring();
                    auto base_pos    = next_node_save->get_base_position();
                    int col = (*selector)(c, base_pos);
                    if(!reserve && col >= 0) {
                        queue.reserve((c->ptn[col] + 1) * level_size);
                        reserve = true;
                    }

                    if(col >= 0) {
                        for (int i = 0; i < c->ptn[col] + 1; ++i) queue.add(next_node, c->lab[col + i]);
                    }
                    next_node = next_node->get_next();
                } while(next_node != start_node);
            }

            void compute_node(refinement* R, sgraph* g, ir::tree* ir_tree,  ir::controller& local_state, ir::tree_node* node, const int v, ir::reduced_save* last_load) {
                auto next_node_save = node->get_save();

                if(next_node_save != last_load) {
                    local_state.load_reduced_state(*next_node_save);
                } else {
                    local_state.move_to_parent();
                    local_state.load_reduced_state_without_coloring(*next_node_save);
                }
                local_state.reset_trace_equal();
                local_state.use_limited_reversible_for_next();
                local_state.use_trace_early_out(true);
                local_state.move_to_child(R, g, v);

                if(local_state.T->trace_equal()) { // TODO: what if leaf?
                    auto new_save = new ir::reduced_save();
                    local_state.save_reduced_state(*new_save);
                    ir_tree->add_node(local_state.base_pos, new_save);
                } else {
                }
            }

            void work_on_todo(bfs_queue& queue, refinement* R, sgraph* g, ir::tree* ir_tree, ir::controller& local_state) {
                ir::reduced_save* last_load = nullptr;
                while(!queue.empty()) {
                    const auto todo = queue.pop();
                    compute_node(R, g, ir_tree, local_state, todo.first, todo.second, last_load);
                    last_load = todo.first->get_save();
                }
            }
        };

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
                    ++s_leaves;
                    lock.unlock();
                } else {
                    lock.unlock();
                }
            }
        };

        /**
         * \brief IR search using random walks.
         *
         * Performs random walks of the IR tree, sifting resulting automorphisms into the given Schreier structure. If the
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

            double get_rolling_sucess_rate() {
                return h_rolling_success;
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
            // TODO: swap out ir_reduced to weighted IR tree later? or just don't use automorphism pruning on BFS...?
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
                        }
                        automorphism.reset();
                    }
                }
            }

            // TODO implement dejavu strategy, more simple
            // TODO depends on ir_tree, selector, and given base (no need to sift beyond base!)
            // TODO: swap out ir_reduced to weighted IR tree later? or just don't use automorphism pruning on BFS...?
            void random_walks_from_tree(refinement &R, std::function<ir::type_selector_hook> *selector, sgraph *g,
                              groups::schreier &group, ir::controller &local_state, ir::tree &ir_tree) {
                groups::automorphism_workspace automorphism(g->v_size);
                groups::schreier_workspace w(g->v_size, &R, g);
                std::vector<int> heuristic_reroll;

                local_state.use_trace_early_out(false);

                h_rolling_first_level_success = 1;
                const int pick_from_level = ir_tree.get_finished_up_to();

                while(!probabilistic_abort_criterion() && !deterministic_abort_criterion(group)
                      && leaf_storage.s_leaves <= h_leaf_limit) { //  * h_rolling_first_level_success
                    auto node = ir_tree.pick_node_from_level(pick_from_level, (int) generator());
                    local_state.load_reduced_state(*node->get_save());

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

                    auto other_leaf = leaf_storage.lookup_leaf(local_state.T->get_hash());
                    if(other_leaf == nullptr) {
                        h_rolling_success = (9.0*h_rolling_success + 0.0) / 10.0;
                        //std::cout << "adding leaf " << local_state.T->get_hash() << std::endl;
                        leaf_storage.add_leaf(local_state.T->get_hash(), *local_state.c, local_state.base_vertex);
                    } else {
                        //std::cout << "reading leaf " << local_state.T->get_hash() << std::endl;
                        automorphism.write_color_diff(local_state.c->vertex_to_col, other_leaf->get_coloring()->lab);
                        const bool cert = R.certify_automorphism_sparse(g, automorphism.perm(), automorphism.nsupport(), automorphism.support());
                        if(cert) {
                            h_rolling_success = (9.0*h_rolling_success + 1.0) / 10.0;
                            //std::cout << "found automorphism, hash " << local_state.T->get_hash() << " support " << automorphism.nsupport() << std::endl;
                            const bool sift = group.sift(w, g, &R, automorphism);
                            if(uniform) record_sift_result(sift);
                        }

                        automorphism.reset();
                    }
                }
            }
        };




        /**
         * \brief Depth-first search without backtracking.
         *
         * Depth-first IR search which does not backtrack. Can parallelize along the base.
         * Due to the search not back-tracking, this module can not deal with difficult parts of combinatorial graphs.
         */
        class dfs_ir {
            int fail_cnt = 0;
            int threads = 1;
            refinement *R = nullptr;
            ir::splitmap *SM = nullptr;
            int cost_snapshot = 0; /**< used to track cost-based abort criterion */
            ir::trace compare_T; // TODO should not be inside dfs_ir

        public:
            long double grp_sz_man = 1.0; /**< group size mantissa */
            int         grp_sz_exp = 0;   /**< group size exponent */

            /**
             * Setup the DFS module.
             *
             * @param threads number of threads we are allowed to dispatch
             * @param R refinement workspace
             */
            void setup(int threads, refinement *R) {
                this->R = R;
                this->threads = threads;
            }

            std::pair<bool, bool> recurse_to_equal_leaf(sgraph *g, int *initial_colors, ir::controller *state,
                                                        groups::automorphism_workspace &automorphism) {
                bool prev_fail = false;
                int prev_fail_pos = -1;
                int cert_pos = 0;

                while ((size_t) state->base_pos < state->compare_base_color.size()) {
                    const int col = state->compare_base_color[state->base_pos];
                    const int col_sz = state->c->ptn[col] + 1;
                    if (col_sz < 2)
                        return {false, false};
                    const int ind_v = state->c->lab[col];

                    state->move_to_child(R, g, ind_v);
                    automorphism.write_singleton(&state->compare_singletons, &state->singletons,
                                                 state->base_singleton_pt[state->base_singleton_pt.size() - 1],
                                                 state->singletons.size());
                    //bool found_auto = R->certify_automorphism_sparse(g, initial_colors->get_array(), automorphism->get_array(),
                    //                                                 automorphism_supp->cur_pos, automorphism_supp->get_array());

                    bool prev_cert = true;

                    assert(state->h_hint_color_is_singleton_now ? state->h_last_refinement_singleton_only : true);

                    if (prev_fail && state->h_last_refinement_singleton_only && state->h_hint_color_is_singleton_now) {
                        prev_cert = R->check_single_failure(g, initial_colors, automorphism.perm(), prev_fail_pos);
                    }

                    //if(state->c->cells == g->v_size) {
                    if (prev_cert && state->h_last_refinement_singleton_only && state->h_hint_color_is_singleton_now) {
                        // TODO: add better heuristic to not always do this check, too expensive!
                        auto cert_res = R->certify_automorphism_sparse_report_fail_resume(g, initial_colors,
                                                                                          automorphism.perm(),
                                                                                          automorphism.nsupport(),
                                                                                          automorphism.support(),
                                                                                          cert_pos);
                        cert_pos = std::get<2>(cert_res);
                        if (std::get<0>(cert_res)) {
                            cert_res = R->certify_automorphism_sparse_report_fail_resume(g, initial_colors,
                                                                                         automorphism.perm(),
                                                                                         automorphism.nsupport(),
                                                                                         automorphism.support(),
                                                                                         0);
                        }

                        if (std::get<0>(cert_res)) {
                            return {true, true};
                        } else {
                            prev_fail = true;
                            prev_fail_pos = std::get<1>(cert_res);
                        }
                    }
                }
                if (state->c->cells == g->v_size) {
                    return {true, false};
                } else {
                    return {false, false};
                }
            }

            // TODO: maybe this should go into the controller?
            void make_leaf_snapshot(ir::controller *state) {
                cost_snapshot = state->T->get_position();

                compare_T.set_compare(true);
                compare_T.set_record(false);
                compare_T.set_compare_trace(state->T);
                compare_T.set_position(state->T->get_position());
                state->T = &compare_T;

                state->compare_base_color.resize(state->base_color.size());
                std::copy(state->base_color.begin(), state->base_color.end(), state->compare_base_color.begin());
                state->compare_base_cells.resize(state->base_cells.size());
                std::copy(state->base_cells.begin(), state->base_cells.end(), state->compare_base_cells.begin());
                state->compare_singletons.resize(state->singletons.size());
                std::copy(state->singletons.begin(), state->singletons.end(), state->compare_singletons.begin());

                state->mode = ir::IR_MODE_COMPARE_TRACE;
            }

            /**
             * Performs DFS from a given leaf node. Does not backtrack and returns the level up to which DFS succeeded.
             *
             * @param g The graph.
             * @param initial_colors The initial coloring of the graph \p g.
             * @param local_state The state from which DFS will be performed. Must be a leaf node of the IR tree.
             * @return The level up to which DFS succeeded.
             */
            int do_dfs(sgraph *g, int *initial_colors, ir::controller &local_state) {
                // orbit algorithm structure
                groups::orbit orbs;
                orbs.initialize(g->v_size);

                // automorphism workspace
                groups::automorphism_workspace pautomorphism(g->v_size);

                // make a snapshot of the leaf to compare to!
                make_leaf_snapshot(&local_state);
                coloring leaf_color;
                leaf_color.copy(local_state.get_coloring());

                // abort criteria
                double recent_cost_snapshot = 0;
                bool   fail = false;

                // loop that serves to optimize Tinhofer graphs
                while (recent_cost_snapshot < 0.25 && local_state.base_pos > 0 && !fail) {
                    local_state.move_to_parent();
                    const int col = local_state.base_color[local_state.base_pos];
                    const int col_sz = local_state.c->ptn[col] + 1;
                    const int vert = local_state.base_vertex[local_state.base_pos];

                    // iterate over current color class
                    for (int i = col_sz - 1; i >= 0; --i) {
                        const int ind_v = leaf_color.lab[col + i];
                        if (ind_v == vert || !orbs.represents_orbit(ind_v)) continue;
                        if (orbs.are_in_same_orbit(ind_v, vert))        continue;

                        // track cost of this refinement for whatever is to come
                        const int cost_start = local_state.T->get_position();

                        // perform individualization-refinement
                        const int prev_base_pos = local_state.base_pos;
                        local_state.T->reset_trace_equal(); // probably need to re-adjust position?
                        local_state.move_to_child(R, g, ind_v);
                        const int wr_pos_st = local_state.base_singleton_pt[local_state.base_singleton_pt.size() - 1];
                        const int wr_pos_end= (int) local_state.singletons.size();

                        // check for sparse automorphism
                        pautomorphism.write_singleton(&local_state.compare_singletons, &local_state.singletons,
                                                      wr_pos_st,wr_pos_end);
                        bool found_auto = R->certify_automorphism_sparse(g, initial_colors, pautomorphism.perm(),
                                                                    pautomorphism.nsupport(),
                                                                    pautomorphism.support());
                        assert(pautomorphism.perm()[vert] == ind_v);

                        // if no luck with sparse automorphism, try more proper walk to leaf node
                        if (!found_auto) {
                            auto rec_succeeded = recurse_to_equal_leaf(g, initial_colors, &local_state, pautomorphism);
                            found_auto = (rec_succeeded.first && rec_succeeded.second);
                            if (rec_succeeded.first && !rec_succeeded.second) {
                                pautomorphism.reset();
                                pautomorphism.write_color_diff(local_state.c->vertex_to_col, leaf_color.lab);
                                found_auto = R->certify_automorphism_sparse(g, initial_colors,
                                                                            pautomorphism.perm(),
                                                                            pautomorphism.nsupport(),
                                                                            pautomorphism.support());
                            }
                        }

                        // track cost-based abort criterion
                        const int cost_end = local_state.T->get_position();
                        double cost_partial  = (cost_end - cost_start) / (cost_snapshot*1.0);
                        recent_cost_snapshot = (cost_partial + recent_cost_snapshot * 3) / 4;

                        // if we found automorphism, add to abort, (and TODO: call hook)
                        if (found_auto) {
                            assert(pautomorphism.perm()[vert] == ind_v);
                            orbs.add_automorphism_to_orbit(pautomorphism);
                        }
                        pautomorphism.reset();

                        // move state back up where we started in this iteration
                        while (prev_base_pos < local_state.base_pos) {
                            local_state.move_to_parent();
                        }

                        // if no automorphism could be determined we would have to backtrack -- so stop!
                        if (!found_auto) {
                            fail = true; // TODO: need to track global failure
                            break;
                        }

                        // if orbit size equals color size now, we are done on this DFS level
                        if (orbs.orbit_size(vert) == col_sz) {
                            break;
                        }
                    }

                    // if we did not fail, accumulate size of current level to group size
                    if (!fail) {
                        grp_sz_man *= col_sz;
                        while (grp_sz_man > 10) {
                            grp_sz_man /= 10;
                            grp_sz_exp += 1;
                        }
                    }
                }

                // if DFS failed on current level, we did not finish the current level -- has to be accounted for
                return local_state.base_pos + (fail);
            }
        };
    }

    /**
     * \brief The dejavu solver.
     *
     * Contains the high-level strategy of the dejavu solver, controlling the interactions between different modules
     * of the solver.
     */
    class dejavu2 {
    private:
        // high-level modules of the algorithm
        sassy::preprocessor m_prep;        /**< preprocessor */
        search_strategy::dfs_ir    m_dfs;  /**< depth-first search */
        search_strategy::bfs_ir    m_bfs;  /**< breadth-first search */
        search_strategy::random_ir m_rand; /**< randomized search */

        // utility tools used by other modules
        refinement m_refinement;          /**< workspace for color refinement and other utilities */
        ir::selector_factory m_selectors; /**< cell selector creation */
        groups::schreier m_schreier; /**< Schreier-Sims algorithm */
        ir::tree m_tree;             /**< IR-tree */

        // TODO: should not be necessary in the end!
        void transfer_sgraph_to_sassy_sgraph(sgraph* g, sassy::sgraph* gg) {
            gg->v = g->v;
            gg->d = g->d;
            gg->e = g->e;
            gg->v_size = g->v_size;
            gg->d_size = g->d_size;
            gg->e_size = g->e_size;
        }
        // TODO: should not be necessary in the end!
        void transfer_sassy_sgraph_to_sgraph(sgraph* g, sassy::sgraph* gg) {
            g->v = gg->v;
            g->d = gg->d;
            g->e = gg->e;
            g->v_size = gg->v_size;
            g->d_size = gg->d_size;
            g->e_size = gg->e_size;
        }

        void man_to_exp_group_size() {
            while(grp_sz_man >= 10.0) {
                grp_sz_exp += 1;
                grp_sz_man = grp_sz_man / 10;
            }
        }

        void add_to_group_size(std::pair<long double, int> add_grp_sz) {
            add_to_group_size(add_grp_sz.first, add_grp_sz.second);
        }

        void add_to_group_size(long double add_grp_sz_man, int add_grp_sz_exp) {
            while (add_grp_sz_man >= 10.0) {
                grp_sz_exp += 1;
                add_grp_sz_man = add_grp_sz_man / 10;
            }
            grp_sz_exp += add_grp_sz_exp;
            grp_sz_man *= add_grp_sz_man;
            man_to_exp_group_size();
        }

    public:
        long double grp_sz_man = 1.0; /**< group size mantissa, see also \a grp_sz_exp */
        int         grp_sz_exp = 0;   /**< group size exponent, see also \a grp_sz_man  */

        /**
         * Compute the automorphisms of the graph \p g colored with vertex colors \p colmap. Automorphisms are returned
         * using the function pointer \p hook.
         *
         * @param g The graph.
         * @param colmap The vertex coloring of \p g. A null pointer is admissible as the trivial coloring.
         * @param hook The hook used for returning automorphisms. A null pointer is admissible if this is not needed.
         *
         */
        void automorphisms(sgraph* g, int* colmap = nullptr, dejavu_hook* hook = nullptr) {
            // TODO high-level strategy
            //  - full restarts, change selector, but make use of previously computed automorphisms (maybe inprocess) -- goal:
            //    stable, single-threaded performance
            //  - exploit strengths of DFS, random walks, BFS, and their synergies
            //  - try to use no-pruning DFS to fill up available leafs more efficiently

            // TODO miscellaneous SAT-inspired stuff
            //      - consider calloc when appropriate
            //      - weigh variables for earlier individualization in restarts (maybe just non-sense, have to see how well
            //        restarts work out first anyway) -- could weigh variables more that turn out conflicting in fails

            // control values
            int h_limit_leaf = 0;
            int h_limit_fail = 0;
            int h_error_prob = 10;
            int h_restarts   = 0;

            // TODO facilities for restarts

            // preprocess the graph using sassy
            std::cout << "preprocessing..." << std::endl;
            sassy::sgraph gg;
            transfer_sgraph_to_sassy_sgraph(g, &gg);
            m_prep.reduce(&gg, colmap, hook);
            add_to_group_size(m_prep.base, m_prep.exp);
            transfer_sassy_sgraph_to_sgraph(g, &gg);

            std::cout << std::endl << "solving..." << std::endl;
            progress_print_header();

            // root state for random and BFS search, as well as restarts
            ir::reduced_save root_save;

            // initialize a coloring using colors of preprocessed graph
            coloring local_coloring;
            g->initialize_coloring(&local_coloring, colmap);

            // set up a local state for IR computations
            ir::trace local_trace;
            ir::controller local_state(&local_coloring, &local_trace);

            // save root state for random and BFS search
            local_state.save_reduced_state(root_save);

            // loop to enable restarts
            while(g->v_size > 0) {
                ++h_restarts; // Dry land is not a myth, I've seen it!

                grp_sz_man = 1.0;
                grp_sz_exp = 0;
                add_to_group_size(m_prep.base, m_prep.exp);

                // find a selector, moves local_state to a leaf of IR tree
                m_selectors.find_sparse_optimized_base(&m_refinement, g, &local_state);
                std::function<ir::type_selector_hook> *current_selector = m_selectors.get_selector_hook();
                progress_print("selector", std::to_string(local_state.base_pos),
                               std::to_string(local_trace.get_position()));
                int base_size = local_state.base_pos;

                std::vector<int> base = local_state.base_vertex;
                std::vector<int> base_sizes = local_state.base_color_sz;

                // TODO fix this, should work...
                //local_state.T->update_blueprint_hash();
                //std::cout << local_state.T->get_hash() << std::endl;

                // depth-first search starting from the computed leaf in local_state
                m_dfs.setup(0, &m_refinement);
                //m_dfs.make_leaf_snapshot(&local_state);

                const int dfs_reached_level = m_dfs.do_dfs(g, colmap, local_state);
                progress_print("dfs", std::to_string(base_size) + "-" + std::to_string(dfs_reached_level),
                               std::to_string(m_dfs.grp_sz_man) + "*10^" + std::to_string(m_dfs.grp_sz_exp));
                add_to_group_size(m_dfs.grp_sz_man, m_dfs.grp_sz_exp);

                if(dfs_reached_level == 0) break;

                // set up schreier structure
                m_schreier.initialize(g->v_size, base, base_sizes, dfs_reached_level);

                // set up IR tree
                ir::tree ir_tree;
                ir_tree.initialize(base_size, &root_save);

                int bfs_cost_estimate;
                int leaf_store_limit;
                bfs_cost_estimate = m_bfs.next_level_estimate(ir_tree, current_selector);
                leaf_store_limit  = 1 + bfs_cost_estimate / 20;

                // special code to tune heuristic for regular graph
                if(root_save.get_coloring()->cells <= 2) {
                    leaf_store_limit = 1;
                    bfs_cost_estimate = bfs_cost_estimate / 5;
                }

                m_rand.specific_walk(m_refinement, base, g, m_schreier, local_state, root_save);

                while(true) {
                    // random automorphisms, single origin
                    m_rand.setup(h_error_prob, leaf_store_limit, 0.1);

                    if(ir_tree.get_finished_up_to() == 0) {
                        m_rand.random_walks(m_refinement, current_selector, g, m_schreier, local_state, root_save);
                    } else {
                        m_rand.random_walks_from_tree(m_refinement, current_selector, g, m_schreier, local_state,
                                                      ir_tree);
                    }
                    // TODO: random walks should record how many "would fail" on first next level to get better bfs_cost_estimate

                    progress_print("urandom", std::to_string(m_rand.stat_leaves()), std::to_string(m_rand.get_rolling_sucess_rate()));
                    progress_print("schreier", "s" + std::to_string(m_schreier.stat_sparsegen()) + "/d" +
                                               std::to_string(m_schreier.stat_densegen()), "_");

                    if(m_rand.deterministic_abort_criterion(m_schreier) || m_rand.probabilistic_abort_criterion())
                        break;

                    // TODO: if rolling sucess very high, and base large, should start skipping levels to fill Schreier faster! (CFI)

                    if(m_rand.get_rolling_sucess_rate() > 0.25) {
                        leaf_store_limit *= 2;
                        continue;
                    }

                    m_bfs.do_a_level(&m_refinement, g, ir_tree, local_state, current_selector);
                    progress_print("bfs", "0-" + std::to_string(ir_tree.get_finished_up_to()) + "(" + std::to_string(bfs_cost_estimate) + ")",
                                   std::to_string(ir_tree.get_level_size(ir_tree.get_finished_up_to())));

                    bfs_cost_estimate = m_bfs.next_level_estimate(ir_tree, current_selector);
                    leaf_store_limit  += (1 + bfs_cost_estimate / 5);

                    if(ir_tree.get_finished_up_to() == base_size)
                        break;
                }
                // breadth-first search & random automorphisms

                break;
            }

            add_to_group_size(m_schreier.compute_group_size());

            std::cout << "#symmetries: " << grp_sz_man << "*10^" << grp_sz_exp << std::endl;
        }
    };
}

#endif //DEJAVU_DFS_H
