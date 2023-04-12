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

    static void progress_print(const std::string proc, const std::string p1, const std::string p2) {
        static std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
        PRINT(std::setw(16) << std::left << (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count()) / 1000000.0  << std::setw(16) << proc << std::setw(16) << p1 << std::setw(16) << p2);
    }

    namespace search_strategy {
        /**
         * \brief Breadth-first search.
         */
        class bfs_ir {
            // TODO implement new, simpler bfs
            // TODO depends on ir_tree, and selector
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
            std::default_random_engine generator;
            stored_leafs leaf_storage;

            int consecutive_success = 0;
            int error_bound = 5;

        public:
            void setup(int error, int leaf_store_limit) {
            }

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

            bool probabilistic_abort_criterion() {
                return (consecutive_success > error_bound);
            }

            bool deterministic_abort_criterion(groups::schreier &group) {
                return (group.finished_up_to_level() + 1 == group.base_size());
            }

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

            // TODO implement dejavu strategy, more simple
            // TODO depends on ir_tree, selector, and given base (no need to sift beyond base!)
            // TODO: swap out ir_reduced to weighted IR tree later? or just don't use automorphism pruning on BFS...?
            void random_walks(refinement &R, std::function<ir::type_selector_hook> *selector, sgraph *g,
                              groups::schreier &group, ir::controller &local_state, ir::reduced_save &start_from) {
                groups::automorphism_workspace automorphism(g->v_size);
                groups::schreier_workspace w(g->v_size, &R, g);
                std::vector<int> heuristic_reroll;

                while(!probabilistic_abort_criterion() && !deterministic_abort_criterion(group)) {
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
            ir::trace compare_T;

        public:
            long double grp_sz_man = 1.0; /**< group size mantissa */
            int         grp_sz_exp = 0;   /**< group size exponent */

            int         cost_snapshot = 0;

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
                            //std::cout << "deep_automorphism" << "@" << state->base_pos << "/" << state->compare_base_color.size() << "(" << state->c->cells << "/" << g->v_size << ")" << std::endl;
                            return {true, true};
                        } else {
                            prev_fail = true;
                            prev_fail_pos = std::get<1>(cert_res);
                        }
                    }
                }
                if (state->c->cells == g->v_size) {
                    //std::cout << "fail" << std::endl;
                    return {true, false};
                } else {
                    return {false, false};
                }
            }

            // adds a given automorphism to the tiny_orbit structure
            void add_automorphism_to_orbit(tiny_orbit *orbit, int *automorphism, int nsupp, int *supp) {
                for (int i = 0; i < nsupp; ++i) {
                    orbit->combine_orbits(automorphism[supp[i]], supp[i]);
                }
            }

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

            // returns base-level reached (from leaf)
            int do_dfs(sgraph *g, int *initial_colors, ir::controller &local_state) {
                // TODO should depend on given selector
                int gens = 0;

                tiny_orbit orbs;
                orbs.initialize(g->v_size);

                //work_list automorphism(g->v_size);
                //work_list automorphism_supp(g->v_size);
                groups::automorphism_workspace pautomorphism(g->v_size);

                //controller local_state(c, &T);

                // start DFS from a leaf!
                //move_to_leaf(g, &local_state);
                // TODO move this outside DFS, only give state in leaf node (we dont even need a selector here!)

                // make a snapshot of the leaf to compare to!
                make_leaf_snapshot(&local_state);

                coloring leaf_color;
                leaf_color.copy(local_state.get_coloring());

                std::vector<int> individualize;

                double recent_cost_snapshot = 0;

                bool fail = false;

                // loop that serves to optimize Tinhofer graphs
                while (recent_cost_snapshot < 0.25 && local_state.base_pos > 0 && !fail) {
                    local_state.move_to_parent();
                    const int col = local_state.base_color[local_state.base_pos]; // TODO: detect stack of "same color"?
                    const int col_sz = local_state.c->ptn[col] + 1;
                    const int vert = local_state.base_vertex[local_state.base_pos];

                    int count_leaf = 0;
                    int count_orb = 0;

                    for (int i = col_sz - 1; i >= 0; --i) {
                        int cost_start = local_state.T->get_position();

                        const int ind_v = leaf_color.lab[col + i];
                        if (ind_v == vert || !orbs.represents_orbit(ind_v))
                            continue;
                        if (orbs.are_in_same_orbit(ind_v, vert)) { // TODO somehow skip to ones not in same orbit?
                            ++count_orb;
                            continue;
                        }

                        ++count_leaf;

                        // call actual DFS with limited fails
                        // Tinhofer graphs will never fail, and have recursion width = 1
                        const int prev_base_pos = local_state.base_pos;
                        local_state.T->reset_trace_equal(); // probably need to re-adjust position?
                        local_state.move_to_child(R, g, ind_v);
                        bool found_auto = false;
                        //if(local_state.h_last_refinement_singleton_only) {
                        pautomorphism.write_singleton(&local_state.compare_singletons, &local_state.singletons,
                                                      local_state.base_singleton_pt[
                                                              local_state.base_singleton_pt.size() - 1],
                                                      local_state.singletons.size());
                        found_auto = R->certify_automorphism_sparse(g, initial_colors, pautomorphism.perm(),
                                                                    pautomorphism.nsupport(),
                                                                    pautomorphism.support());
                        assert(pautomorphism.perm()[vert] == ind_v);
                        //}
                        //std::cout << found_auto << ", " << automorphism_supp.cur_pos << std::endl;
                        //reset_automorphism(automorphism.get_array(), &automorphism_supp);

                        // try proper recursion
                        if (!found_auto) {
                            // TODO: heuristic that once this becomes common, we start to copy away the state and recover it, and don't track
                            // TODO: touched stuff
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

                        int cost_end = local_state.T->get_position();
                        double cost_partial = (cost_end - cost_start) / (cost_snapshot*1.0);
                        recent_cost_snapshot = (cost_partial + recent_cost_snapshot * 3) / 4;

                        if (found_auto) {
                            assert(pautomorphism.perm()[vert] == ind_v);
                            ++gens;
                            add_automorphism_to_orbit(&orbs, pautomorphism.perm(),
                                                      pautomorphism.nsupport(), pautomorphism.support());
                        }
                        pautomorphism.reset();

                        while (prev_base_pos < local_state.base_pos) {
                            local_state.move_to_parent();
                        }

                        if (!found_auto) {
                            std::cout << ind_v << "(F) " << std::endl;
                            fail = true;
                            break;
                        }

                        if (orbs.orbit_size(vert) == col_sz) {
                            break;
                        }
                    }

                    if (!fail) {
                        grp_sz_man *= col_sz;
                        while (grp_sz_man > 10) {
                            grp_sz_man /= 10;
                            grp_sz_exp += 1;
                        }
                    }
                }

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
        refinement       m_refinement; /**< workspace for color refinement and other utilities */
        ir::selector_factory m_selectors;  /**< cell selector creation */
        groups::schreier m_schreier; /**< Schreier-Sims algorithm */
        ir::ir_tree   m_tree;              /**< IR-tree */

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
            //        restarts work out first anyway)

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
            while(true) {
                ++h_restarts; // Dry land is not a myth, I've seen it!

                grp_sz_man = 1.0; /**< group size mantissa, see also \a grp_sz_exp */
                grp_sz_exp = 0;   /**< group size exponent, see also \a grp_sz_man  */
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

                // set up schreier structure
                m_schreier.setup(g->v_size, base, base_sizes, dfs_reached_level);

                // intertwined random automorphisms and breadth-first search

                // random automorphisms, single origin
                m_rand.setup(h_error_prob, h_limit_leaf);
                m_rand.specific_walk(m_refinement, base, g, m_schreier, local_state, root_save);
                m_rand.random_walks(m_refinement, current_selector, g, m_schreier, local_state, root_save);

                progress_print("urandom", std::to_string(m_rand.stat_leaves()), "_");
                progress_print("schreier", "s" + std::to_string(m_schreier.stat_sparsegen()) + "/d" +
                                           std::to_string(m_schreier.stat_densegen()), "_");

                add_to_group_size(m_schreier.compute_group_size());

                // breadth-first search & random automorphisms

                break;
            }

            std::cout << "#symmetries: " << grp_sz_man << "*10^" << grp_sz_exp << std::endl;
        }
    };
}

#endif //DEJAVU_DFS_H
