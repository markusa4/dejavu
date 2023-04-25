// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_DEJAVU_H
#define DEJAVU_DEJAVU_H

#include "dfs.h"
#include "bfs.h"
#include "rand.h"
#include "sassy/preprocessor.h"

extern dejavu::ir::refinement test_r;
extern sgraph _test_graph;
extern int*   _test_col;

namespace dejavu {
    static bool set_first = false;
    static std::chrono::high_resolution_clock::time_point first    = std::chrono::high_resolution_clock::now();
    static std::chrono::high_resolution_clock::time_point previous = std::chrono::high_resolution_clock::now();

    void test_hook(int n, const int *p, int nsupp, const int *supp) {
        std::cout << "certifying..." << std::endl;
        assert(test_r.certify_automorphism_sparse(&_test_graph, _test_col, p, nsupp, supp));
    }

    static void progress_print_split() {
        PRINT("______________________________________________________________");
    }

    static void progress_print_header() {
        progress_print_split();
        PRINT(std::setw(9) << std::left <<"T (ms)" << std::setw(9) << "δ (ms)" << std::setw(16) << "proc"  << std::setw(16) << "P1"        << std::setw(16)        << "P2");
        progress_print_split();
        PRINT(std::setw(9) << std::left << 0 << std::setw(9) << 0 << std::setw(16) << "start" << std::setw(16) << "_" << std::setw(16) << "_" );
    }

    static void progress_print(const std::string
                               proc, const std::string p1, const std::string p2) {
        if(!set_first) {
            first     = std::chrono::high_resolution_clock::now();
            previous  = first;
            set_first = true;
        }
        auto now = std::chrono::high_resolution_clock::now();
        PRINT(std::fixed << std::setprecision(2) << std::setw(9) << std::left << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - first).count()) / 1000000.0  << std::setw(9) << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - previous).count()) / 1000000.0  << std::setw(16) << proc << std::setw(16) << p1 << std::setw(16) << p2);
        previous = now;
    }

    static void progress_print(const std::string
                                 proc, const double p1, const double p2) {
        if(!set_first) {
            first     = std::chrono::high_resolution_clock::now();
            previous  = first;
            set_first = true;
        }

        auto now = std::chrono::high_resolution_clock::now();
        PRINT(std::fixed << std::setprecision(2) << std::setw(9) << std::left << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - first).count()) / 1000000.0  << std::setw(9) << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - previous).count()) / 1000000.0  << std::setw(16) << proc << std::setw(16) << p1 << std::setw(16) << p2);
        previous = now;
    }

    /**
     * \brief The dejavu solver.
     *
     * Contains the high-level strategy of the dejavu solver, controlling the interactions between different modules
     * of the solver.
     */
    class dejavu2 {
    private:
        // shared modules
        groups::shared_schreier* m_schreier = nullptr;  /**< Schreier-Sims structure */
        ir::shared_tree*         m_tree     = nullptr;  /**< IR tree                 */

        // TODO: should not be necessary in the end!
        void transfer_sgraph_to_sassy_sgraph(sgraph* g, sassy::sgraph* gg) {
            gg->v = g->v;
            gg->d = g->d;
            gg->e = g->e;
            gg->v_size = g->v_size;
            gg->d_size = g->v_size;
            gg->e_size = g->e_size;
        }
        // TODO: should not be necessary in the end!
        void transfer_sassy_sgraph_to_sgraph(sgraph* g, sassy::sgraph* gg) {
            g->v = gg->v;
            g->d = gg->d;
            g->e = gg->e;
            g->v_size = gg->v_size;
            g->e_size = gg->e_size;
        }

        void man_to_exp_group_size() {
            while(s_grp_sz_man >= 10.0) {
                s_grp_sz_exp += 1;
                s_grp_sz_man = s_grp_sz_man / 10;
            }
        }

        void add_to_group_size(std::pair<long double, int> add_grp_sz) {
            add_to_group_size(add_grp_sz.first, add_grp_sz.second);
        }

        void add_to_group_size(long double add_grp_sz_man, int add_grp_sz_exp) {
            while (add_grp_sz_man >= 10.0) {
                s_grp_sz_exp += 1;
                add_grp_sz_man = add_grp_sz_man / 10;
            }
            s_grp_sz_exp += add_grp_sz_exp;
            s_grp_sz_man *= add_grp_sz_man;
            man_to_exp_group_size();
        }

    public:
        // TODO outward-facing, global settings and statistics of solver go here, feed to-and-from sub-modules
        // settings
        int h_error_prob        = 10; /**< error probability is below `(1/2)^h_error_prob`, default value of 10 thus
                                       * misses an automorphism with probabiltiy at most (1/2)^10 < 0.098% */
        int h_limit_leaf        = 0;  /**< limit for the amount of IR leaves to be stored*/
        int h_limit_fail        = 0;  /**< limit for the amount of backtracking allowed*/

        // statistics
        long double s_grp_sz_man = 1.0; /**< group size mantissa, group size is `s_grp_sz_man^s_grp_sz_exp`
                                          * \sa s_grp_sz_exp */
        int         s_grp_sz_exp = 0;   /**< group size exponent, group size is `s_grp_sz_man^s_grp_sz_exp`
                                          * \sa s_grp_sz_man  */

        /**
         * Compute the automorphisms of the graph \p g colored with vertex colors \p colmap. Automorphisms are returned
         * using the function pointer \p hook.
         *
         * @param g The graph.
         * @param colmap The vertex coloring of \p g. A null pointer is admissible as the trivial coloring.
         * @param hook The hook used for returning automorphisms. A null pointer is admissible if this is not needed.
         *
         * \sa A description of the graph format can be found in sgraph.
         */
        void automorphisms(sgraph* g, int* colmap = nullptr, dejavu_hook* hook = nullptr) {
            // TODO high-level strategy
            //  - try to use no-pruning DFS to fill up available leafs more efficiently

            // TODO miscellaneous SAT-inspired stuff
            //      - consider calloc when appropriate
            //      - weigh variables for earlier individualization in restarts (maybe just non-sense, have to see how well
            //        restarts work out first anyway) -- could weigh variables more that turn out conflicting in fails

            sassy::preprocessor m_prep;        /*< preprocessor */

            // preprocess the graph using sassy
            std::cout << "preprocessing..." << std::endl;
            sassy::sgraph gg;
            transfer_sgraph_to_sassy_sgraph(g, &gg);
            m_prep.reduce(&gg, colmap, hook);
            add_to_group_size(m_prep.base, m_prep.exp);

            if(gg.v_size <= 1) return;

            // TODO fix up the hook

            // print that we are solving now...
            std::cout << std::endl << "solving..." << std::endl;
            progress_print_header();

            // we first need to set up quite a few datastructures and modules
            // we start by transferring the graph
            transfer_sassy_sgraph_to_sgraph(g, &gg);

            // heuristic settings
            int  h_restarts         = -1;
            int  h_budget           = 1;
            int  h_budget_inc_fact  = 5;

            // statistics
            int  s_cost             = 0;
            bool s_large_base       = false;
            bool s_short_base       = false;
            int  s_leaves_added_this_restart = 0;
            bool s_inprocessed = false; /*< flag which says that we've inprocessed the graph in the last restart */

            // some data we want to keep track of at high-level
            std::vector<int>                 inproc_fixed_points;      /*< vertices fixed by inprocessing      */
            std::vector<std::pair<int, int>> inproc_can_individualize; /*< vertices that can be individualized */
            std::vector<int> base;       /*< current base of vertices       */
            std::vector<int> base_sizes; /*< current size of colors on base */

            // local modules and workspace, to be used by other modules
            ir::refinement       m_refinement; /*< workspace for color refinement and other utilities */
            ir::selector_factory m_selectors;  /*< cell selector creation */
            groups::automorphism_workspace automorphism(g->v_size); /*< workspace to keep an automorphism */
            groups::schreier_workspace     schreierw(g->v_size);    /*< workspace for Schreier-Sims       */

            // high-level search strategies
            search_strategy::dfs_ir    m_dfs;  /*< depth-first search */
            search_strategy::bfs_ir    m_bfs;  /*< breadth-first search */
            search_strategy::random_ir m_rand; /*< randomized search */

            // link high-level strategies to local workspace
            m_dfs.link_to_workspace(&automorphism);
            m_bfs.link_to_workspace(&automorphism);
            m_rand.link_to_workspace(&schreierw, &automorphism);

            // initialize a coloring using colors of preprocessed graph
            coloring local_coloring;
            g->initialize_coloring(&local_coloring, colmap);
            const bool s_regular = local_coloring.cells == 1; // is the graph regular?

            // set up a local state for IR computations
            ir::controller local_state(&m_refinement, &local_coloring);

            // save root state for random and BFS search, as well as restarts
            ir::reduced_save root_save;
            local_state.save_reduced_state(root_save);
            int s_last_base_size = g->v_size + 1;

            // now that we are set up, let's start solving the graph
            // loop for restarts
            while (true) {
                ++h_restarts; // Dry land is not a myth, I've seen it!
                if (h_restarts >= 1) {
                    local_state.load_reduced_state(root_save);
                    progress_print_split();
                    if(h_restarts >= 3) {
                        h_budget *= h_budget_inc_fact;
                    } else {
                        h_budget *= 2;
                    }
                    s_cost = 0;

                    if(s_inprocessed && m_tree) m_tree->clear_leaves(); // TODO if I link up the corresponding root to each leaf, I can keep using them... also "from BFS" would be more efficient in leaf recovery?
                    m_rand.reset();

                    s_leaves_added_this_restart = 0;
                }

                // keep vertices which are save to individualize for in-processing
                inproc_can_individualize.clear();

                // find a selector, moves local_state to a leaf of IR shared_tree
                m_selectors.find_base(h_restarts, &m_refinement, g, &local_state);
                auto selector = m_selectors.get_selector_hook();
                progress_print("selector" + std::to_string(h_restarts), local_state.s_base_pos,
                               local_state.T->get_position());
                int base_size = local_state.s_base_pos;
                // immediately discard this base if deemed too large
                if (base_size > 5 * s_last_base_size) {
                    progress_print("skip" + std::to_string(h_restarts), local_state.s_base_pos,
                                   s_last_base_size);
                    continue; // TODO: re-consider this
                }
                if(base_size > sqrt(g->v_size)) {
                    s_large_base = true;
                } else if(base_size <= 2) {
                    s_short_base = true;
                }
                // TODO: here needs to go an algorithm to determine whether we want to continue with this base
                // TODO: ... and which parts we want to use of previous restart
                s_last_base_size = base_size;


                // make snapshot of trace and leaf for all following search
                local_state.compare_to_this();
                base       = local_state.base_vertex;   // we want to keep this for later
                base_sizes = local_state.base_color_sz; // used for upper bounds in Schreier structure

                // depth-first search starting from the computed leaf in local_state
                m_dfs.h_recent_cost_snapshot_limit = s_large_base ? 0.33 : 0.25; // set up DFS heuristic
                const int dfs_level =
                        m_dfs.do_dfs(hook, g, root_save.get_coloring(), local_state,&inproc_can_individualize);
                progress_print("dfs", std::to_string(base_size) + "-" + std::to_string(dfs_level),
                               "~"+std::to_string((int)m_dfs.grp_sz_man) + "*10^" + std::to_string(m_dfs.grp_sz_exp));
                if (dfs_level == 0) break; // DFS finished the graph -- we are done!

                // next, we go into the random path + BFS algorithm, so we set up a Schreier structure and IR tree
                // first, the schreier structure
                if(m_schreier)  {
                    // reset the Schreier structure
                    const bool can_keep_previous = (h_restarts >= 3) && !s_inprocessed;
                    const bool reset_prob =
                            m_schreier->reset(schreierw, base, base_sizes, dfs_level, can_keep_previous, inproc_fixed_points);
                    if(reset_prob) {m_schreier->reset_probabilistic_criterion();}
                } else {
                    // create the Schreier structure
                    m_schreier = new groups::shared_schreier();
                    m_schreier->initialize(g->v_size, base, base_sizes, dfs_level);
                }

                // then, the IR tree
                if(m_tree) {
                    // reset the IR tree
                    const bool can_keep_previous = (h_restarts >= 3) && !s_inprocessed;
                    m_tree->reset(base, &root_save, can_keep_previous);
                } else {
                    // create the IR tree
                    m_tree = new ir::shared_tree(g->v_size);
                    m_tree->initialize(base, &root_save);
                }
                s_inprocessed = false;

                int bfs_cost_estimate;
                int leaf_store_limit;
                bfs_cost_estimate = m_bfs.next_level_estimate(*m_tree, selector);
                leaf_store_limit = std::min(1 + bfs_cost_estimate / 20, h_budget);
                //std::cout << bfs_cost_estimate << ", " << leaf_store_limit << std::endl;

                // special code to tune heuristic for regular graph
                bool h_skip_random_paths = false;
                if ((root_save.get_coloring()->cells <= 2)) { // && base_size <= 2) || (root_save.get_coloring()->cells == 1
                    leaf_store_limit = 1;
                    h_skip_random_paths = true;
                }

                bool fail = false;

                // TODO find a better way to add canonical leaf to leaf store
                m_rand.specific_walk(g, *m_tree, local_state, base);

                //if(inproc_can_individualize.size() > 10) fail = true;

                const int flat_leaf_store_inc = m_tree->stat_leaves();
                leaf_store_limit = std::max(std::min(leaf_store_limit, h_budget/2), 2);

                while (!fail) {
                    m_rand.setup(h_error_prob, flat_leaf_store_inc+leaf_store_limit,
                                 (m_rand.s_rolling_first_level_success > 0.5) && !s_short_base); // TODO: look_close does not seem worth it

                    // TODO: if no node was pruned, use m_rand.random_walks instead of "from BFS"! should fix CFI

                    //std::cout << "budget: " << h_budget << " leafs: " << flat_leaf_store_inc << " add " << leaf_store_limit << std::endl;

                    const int leaves_pre = m_tree->stat_leaves();
                    if (m_tree->get_finished_up_to() == 0) {
                        // random automorphisms, single origin
                        m_rand.random_walks(g, hook, selector, *m_tree, *m_schreier,local_state);
                    } else {
                        // random automorphisms, sampled from shared_tree
                        m_rand.random_walks_from_tree(g, hook, selector, *m_tree, *m_schreier,local_state);
                    }
                    s_leaves_added_this_restart = m_tree->stat_leaves() - leaves_pre;

                    progress_print("urandom", m_tree->stat_leaves(),m_rand.s_rolling_success);
                    progress_print("schreier", "s" + std::to_string(m_schreier->s_sparsegen()) + "/d" +
                                               std::to_string(m_schreier->s_densegen()), "_");

                    if (m_schreier->deterministic_abort_criterion() || m_schreier->probabilistic_abort_criterion()) {
                        //std::cout << m_rand.deterministic_abort_criterion(*m_schreier) << ", " << m_rand.probabilistic_abort_criterion() << std::endl;
                        fail = false;
                        break;
                    }

                    s_cost += s_leaves_added_this_restart; // TODO: should only add diff
                    //std::cout << s_cost << "/" << h_budget << std::endl;
                    if (s_cost > h_budget) {
                        fail = true;
                        progress_print("restart", std::to_string(s_cost), std::to_string(h_budget));
                        break;
                    }
                    // TODO: if rolling sucess very high, and base large, should start skipping levels to fill Schreier faster! (CFI)

                    if(m_rand.s_rolling_first_level_success > 0.99) h_skip_random_paths = false;

                    if(s_short_base && s_regular && m_rand.s_rolling_success < 0.01) h_skip_random_paths = true;
                    else if (!(((bfs_cost_estimate <= g->v_size ||
                          bfs_cost_estimate * m_rand.s_rolling_first_level_success < 12)) && m_rand.s_rolling_first_level_success < 0.99) && leaf_store_limit <= (bfs_cost_estimate*m_rand.s_rolling_first_level_success)/2) {
                        if (m_rand.s_rolling_success > 0.1 && !h_skip_random_paths) {
                            leaf_store_limit *= 2;
                            continue;
                        }
                        if (m_rand.s_rolling_first_level_success * bfs_cost_estimate > leaf_store_limit &&
                            !h_skip_random_paths) {
                            //std::cout << "increased here " << m_rand.get_rolling_first_level_success_rate() << ", "
                            //          << bfs_cost_estimate << std::endl;
                            leaf_store_limit = std::min(std::max(
                                    (int) m_rand.s_rolling_first_level_success * bfs_cost_estimate,
                                    ((int) (1.5 * (leaf_store_limit + 1.0)))), h_budget+1);
                            continue;
                        }
                    }

                    if(h_restarts < 2 && s_cost + bfs_cost_estimate > h_budget && !h_skip_random_paths) {
                        fail = true;
                        progress_print("restart", std::to_string(s_cost), std::to_string(h_budget));
                        break;
                    }

                    h_skip_random_paths = false;

                    bfs_cost_estimate = m_tree->get_finished_up_to() < base_size ? m_bfs.next_level_estimate(*m_tree, selector) : -1;
                    m_bfs.do_a_level(g, *m_tree, local_state, selector);
                    progress_print("bfs", "0-" + std::to_string(m_tree->get_finished_up_to()) + "(" +
                                          std::to_string(bfs_cost_estimate) + ")",
                                   std::to_string(m_tree->get_level_size(m_tree->get_finished_up_to())));

                    if (m_tree->get_finished_up_to() == base_size) {
                       //add_to_group_size(m_bfs.top_level_orbit.orbit_size(base[0]), 0);
                       if(base_size == 2) {
                           add_to_group_size(m_tree->h_bfs_top_level_orbit.orbit_size(base[0])*m_tree->h_bfs_automorphism_pw, 0);
                       } else {
                           add_to_group_size(m_tree->get_level_size(m_tree->get_finished_up_to()), 0);
                       }
                        break;
                    }

                    // want to inprocess if BFS was successfull in pruning
                    if(m_tree->get_finished_up_to() == 1 && m_tree->get_level_size(m_tree->get_finished_up_to()) < bfs_cost_estimate / 2) {
                        fail = true;
                        break;
                    }

                    bfs_cost_estimate = m_bfs.next_level_estimate(*m_tree, selector);
                    leaf_store_limit += m_tree->get_level_size(m_tree->get_finished_up_to());
                }
                // breadth-first search & random automorphisms

                if (!fail) {
                    break;
                } else {
                    local_state.load_reduced_state(root_save);
                    const int cell_prev = root_save.get_coloring()->cells;

                    if(m_tree->get_finished_up_to() >= 1) {
                        // TODO: improve saving hash for pruned nodes, then propagate this to first level
                        // TODO: could fix "usr" maybe?
                        // TODO: use "pruned" flag, etc.
                        // TODO: hashing actually makes this worse!
                        work_list hash(g->v_size);
                        mark_set  is_pruned(g->v_size);
                        m_tree->mark_first_level(is_pruned);

                        for (int i = 0; i < g->v_size; ++i) {
                            if(is_pruned.get(i)) hash[i] = 1;
                            else                      hash[i] = 0;
                            hash[i] += (int) (*m_tree->get_node_invariant())[i] % 16;
                        }
                        for (int i = 0; i < g->v_size; ++i) {
                            hash[i] = 17 * local_state.c->vertex_to_col[i] + hash[i];
                        }
                        g->initialize_coloring(local_state.c, hash.get_array());
                        const int cell_after = local_state.c->cells;
                        progress_print("inproc_bfs", std::to_string(cell_prev),
                                       std::to_string(cell_after));
                        if (cell_after != cell_prev) {
                            local_state.refine(g);
                            progress_print("inproc_ref", std::to_string(cell_after),
                                           std::to_string(local_state.c->cells));
                        }
                    }

                    m_schreier->determine_save_to_individualize(&inproc_can_individualize, local_state.get_coloring());

                    if(!inproc_can_individualize.empty()) {
                        for (auto & i : inproc_can_individualize) {
                            const int ind_v   = i.first;
                            const int test_col_sz = i.second;
                            const int ind_col = local_state.c->vertex_to_col[ind_v];
                            assert(test_col_sz == local_state.c->ptn[ind_col]+1);
                            m_prep.multiply_to_group_size(local_state.c->ptn[ind_col]+1);
                            local_state.move_to_child_no_trace(g, ind_v);
                            inproc_fixed_points.push_back(ind_v);
                        }
                        progress_print("inproc_ind", "_",
                                       std::to_string(local_state.c->cells));
                        inproc_can_individualize.clear();
                    }

                    inproc_can_individualize.clear();

                    local_state.save_reduced_state(root_save);
                    // TODO use orbit partition for 1 individualization
                    // TODO could run sassy if individualization succeeds...

                    s_inprocessed = cell_prev != local_state.c->cells;
                }
            }

            // We are done! Let's add up the total group size from all the different modules.
            s_grp_sz_man = 1.0;
            s_grp_sz_exp = 0;
            add_to_group_size(m_prep.base, m_prep.exp);
            add_to_group_size(m_dfs.grp_sz_man, m_dfs.grp_sz_exp);
            if(m_schreier) {
                add_to_group_size(m_schreier->compute_group_size());
                delete m_schreier;
            }
            if(m_tree) delete m_tree;
        }
    };
}

#endif //DEJAVU_DEJAVU_H
