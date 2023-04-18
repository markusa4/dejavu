#ifndef DEJAVU_DEJAVU_H
#define DEJAVU_DEJAVU_H

#include "dfs.h"
#include "bfs.h"
#include "rand.h"
#include "sassy/preprocessor.h"

namespace dejavu {

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
        static std::chrono::high_resolution_clock::time_point first    = std::chrono::high_resolution_clock::now();
        static std::chrono::high_resolution_clock::time_point previous = std::chrono::high_resolution_clock::now();

        auto now = std::chrono::high_resolution_clock::now();
        PRINT(std::fixed << std::setprecision(2) << std::setw(9) << std::left << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - first).count()) / 1000000.0  << std::setw(9) << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - previous).count()) / 1000000.0  << std::setw(16) << proc << std::setw(16) << p1 << std::setw(16) << p2);
        previous = now;
    }

    static void progress_print(const std::string
                                 proc, const double p1, const double p2) {
        static std::chrono::high_resolution_clock::time_point first    = std::chrono::high_resolution_clock::now();
        static std::chrono::high_resolution_clock::time_point previous = std::chrono::high_resolution_clock::now();

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
        // high-level modules of the algorithm
        sassy::preprocessor m_prep;        /**< preprocessor */
        search_strategy::dfs_ir    m_dfs;  /**< depth-first search */
        search_strategy::bfs_ir    m_bfs;  /**< breadth-first search */
        search_strategy::random_ir m_rand; /**< randomized search */

        // utility tools used by other modules
        refinement           m_refinement;/**< workspace for color refinement and other utilities */
        ir::selector_factory m_selectors; /**< cell selector creation */
        groups::schreier*    m_schreier = nullptr;  /**< Schreier-Sims algorithm */
        ir::shared_tree m_tree;           /**< IR-shared_tree */

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
            int h_restarts   = -1;
            int h_budget     = 1;
            int h_cost       = 0;
            int h_budget_inc_fact = 5;
            bool h_large_base = false;
            bool h_short_base = false;

            int h_leaves_added_this_restart = 0;

            bool inprocessed = false;

            // TODO facilities for restarts

            // preprocess the graph using sassy
            std::cout << "preprocessing..." << std::endl;
            sassy::sgraph gg;
            transfer_sgraph_to_sassy_sgraph(g, &gg);
            m_prep.reduce(&gg, colmap, hook);
            add_to_group_size(m_prep.base, m_prep.exp);
            if(gg.v_size > 0) {
                transfer_sassy_sgraph_to_sgraph(g, &gg);

                std::cout << std::endl << "solving..." << std::endl;
                progress_print_header();

                // initialize a coloring using colors of preprocessed graph
                coloring local_coloring;
                g->initialize_coloring(&local_coloring, colmap);

                bool h_regular = local_coloring.cells == 1;

                // set up a local state for IR computations
                ir::controller local_state(&local_coloring);

                // save root state for random and BFS search, as fwell as restarts
                ir::reduced_save root_save;
                local_state.save_reduced_state(root_save);

                int last_base_size = g->v_size + 1;

                // loop to enable restarts
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
                        h_cost = 0;

                        if(inprocessed) m_rand.clear_leaves();
                        m_rand.reset();

                        h_leaves_added_this_restart = 0;

                        // TODO: notice similarities to previous base, don't throw away bfs etc. if same up to certain point!
                        // TODO: inprocess, sift previous automorphisms if graph hard and we got some, make use of restarts more
                        // TODO: record what happens during a run, make next choices based on that -- prevent "obvious"
                        //       bad choices (bases that are much larger, etc. ...)
                        // TODO: prevent restarts if "almost done", or something
                        // TODO: after we tested all 3 cell selectors once, strike out some of them which are stupid
                        // TODO: ... maybe don't increment budget for first 3 too much, but then go more aggressive and add randomness
                        // TODO: also, if it seems like "graph is hard" and sifting is not going to be of any issue, just start
                        // TODO: keeping the schreier structure, doesn't matter that it's a different base then... also sift
                        // TODO dfs elements into the structure!

                        // TODO use trace distances to measure average cost of "~1 experimental path" or  "1 element in bfs"
                    }

                    std::vector<std::pair<int, int>> save_to_individualize;

                    grp_sz_man = 1.0;
                    grp_sz_exp = 0;
                    add_to_group_size(m_prep.base, m_prep.exp);

                    // find a selector, moves local_state to a leaf of IR shared_tree
                    m_selectors.find_base(h_restarts, &m_refinement, g, &local_state);
                    //m_selectors.find_sparse_optimized_base(&m_refinement, g, &local_state);
                    std::function<ir::type_selector_hook> *current_selector = m_selectors.get_selector_hook();
                    progress_print("selector" + std::to_string(h_restarts), std::to_string(local_state.base_pos),
                                   std::to_string(local_state.T->get_position()));
                    int base_size = local_state.base_pos;
                    if (base_size > 5 * last_base_size) {
                        progress_print("skip" + std::to_string(h_restarts), std::to_string(local_state.base_pos),
                                       std::to_string(last_base_size));
                        continue; // TODO: re-consider this
                    }
                    if(base_size > sqrt(g->v_size)) {
                        h_large_base = true;
                    } else if(base_size <= 2) {
                        h_short_base = true;
                    }
                    // TODO: here needs to go an algorithm to determine whether we want to continue with this base
                    // TODO: ... and which parts we want to use of previous restart
                    last_base_size = base_size;


                    // make snapshot of trace and leaf for all following search
                    local_state.compare_to_this();

                    std::vector<int> base = local_state.base_vertex;
                    std::vector<int> base_sizes = local_state.base_color_sz;

                    // depth-first search starting from the computed leaf in local_state
                    // TODO I think there should be no setup functions at all
                    m_dfs.setup(0, &m_refinement, h_large_base?0.33:0.25);

                    const int dfs_reached_level = m_dfs.do_dfs(g, root_save.get_coloring(), local_state, &save_to_individualize);
                    progress_print("dfs", std::to_string(base_size) + "-" + std::to_string(dfs_reached_level),
                                   "~"+std::to_string((int)m_dfs.grp_sz_man) + "*10^" + std::to_string(m_dfs.grp_sz_exp));
                    add_to_group_size(m_dfs.grp_sz_man, m_dfs.grp_sz_exp);
                    //std::cout << "dfs:#symmetries: " << std::to_string(grp_sz_man) << "*10^" << grp_sz_exp << std::endl;

                    if (dfs_reached_level == 0) break;

                    // set up schreier structure
                    if(m_schreier)  {
                        const bool reset_prob = m_schreier->reset(base, base_sizes, dfs_reached_level, (h_restarts >= 3) && !inprocessed);
                        if(reset_prob) {m_rand.reset_prob();}
                    } else {
                        m_schreier = new groups::schreier(); // TODO bad memory leak, make a reset function
                        m_schreier->initialize(g->v_size, base, base_sizes, dfs_reached_level);
                    }
                    // set up IR shared_tree
                    // TODO: make a reset function -- maybe make the reset function smart and have it consider the base
                    // TODO: also consider saving older trees...
                    ir::shared_tree ir_tree;
                    ir_tree.initialize(base_size, &root_save);

                    inprocessed = false;

                    int bfs_cost_estimate;
                    int leaf_store_limit;
                    bfs_cost_estimate = m_bfs.next_level_estimate(ir_tree, current_selector);
                    leaf_store_limit = std::min(1 + bfs_cost_estimate / 20, h_budget);
                    //std::cout << bfs_cost_estimate << ", " << leaf_store_limit << std::endl;

                    // special code to tune heuristic for regular graph
                    bool h_skip_random_paths = false;
                    if ((root_save.get_coloring()->cells <= 2)) { // && base_size <= 2) || (root_save.get_coloring()->cells == 1
                        leaf_store_limit = 1;
                        //std::cout << "1" << std::endl;
                        h_skip_random_paths = true;
                    }

                    bool fail = false;

                    // TODO find a better way to add canonical leaf to leaf store
                    m_rand.specific_walk(m_refinement, base, g, *m_schreier, local_state, root_save);

                    //if(save_to_individualize.size() > 10) fail = true;

                    const int flat_leaf_store_inc = m_rand.stat_leaves();

                    leaf_store_limit = std::max(std::min(leaf_store_limit, h_budget/2), 2);
                    while (!fail) {
                        m_rand.setup(h_error_prob, flat_leaf_store_inc+leaf_store_limit, 0.1, (m_rand.get_rolling_first_level_success_rate() > 0.5) && !h_short_base); // TODO: look_close does not seem worth it

                        // TODO: if no node was pruned, use m_rand.random_walks instead of "from BFS"! should fix CFI

                        //std::cout << "budget: " << h_budget << " leafs: " << flat_leaf_store_inc << " add " << leaf_store_limit << std::endl;

                        const int leaves_pre = m_rand.stat_leaves();
                        if (ir_tree.get_finished_up_to() == 0) {
                            // random automorphisms, single origin
                            m_rand.random_walks(m_refinement, current_selector, g, *m_schreier, local_state, &root_save);
                        } else {
                            // random automorphisms, sampled from shared_tree
                            m_rand.random_walks_from_tree(m_refinement, current_selector, g, *m_schreier,
                                                          local_state,ir_tree);
                        }
                        h_leaves_added_this_restart = m_rand.stat_leaves() - leaves_pre;

                        // TODO: random walks should record how many "would fail" on first next level to get better bfs_cost_estimate
                        // TODO: leaf store should not be in m_rand, but m_rand can contain the schreier workspace! (m_rand is a local workspace, such as the other modules)
                        // TODO: we should pass m_refinement in the constructor of the other modules
                        // TODO: if everything has reset functions, we could pass the other structures as well...
                        progress_print("urandom", m_rand.stat_leaves(),
                                       m_rand.get_rolling_sucess_rate());
                        progress_print("schreier", "s" + std::to_string(m_schreier->stat_sparsegen()) + "/d" +
                                                   std::to_string(m_schreier->stat_densegen()), "_");

                        if (m_rand.deterministic_abort_criterion(*m_schreier) || m_rand.probabilistic_abort_criterion()) {
                            //std::cout << m_rand.deterministic_abort_criterion(*m_schreier) << ", " << m_rand.probabilistic_abort_criterion() << std::endl;
                            fail = false;
                            break;
                        }

                        h_cost += h_leaves_added_this_restart; // TODO: should only add diff
                        //std::cout << h_cost << "/" << h_budget << std::endl;
                        if (h_cost > h_budget) {
                            fail = true;
                            progress_print("restart", std::to_string(h_cost), std::to_string(h_budget));
                            break;
                        }
                        // TODO: if rolling sucess very high, and base large, should start skipping levels to fill Schreier faster! (CFI)

                        if(m_rand.get_rolling_first_level_success_rate() > 0.99) h_skip_random_paths = false;

                        if(h_short_base && h_regular && m_rand.get_rolling_sucess_rate() < 0.01) h_skip_random_paths = true;
                        else if (!(((bfs_cost_estimate <= g->v_size ||
                              bfs_cost_estimate * m_rand.get_rolling_first_level_success_rate() < 12)) && m_rand.get_rolling_first_level_success_rate() < 0.99) && leaf_store_limit <= (bfs_cost_estimate*m_rand.get_rolling_first_level_success_rate())/2) {
                            if (m_rand.get_rolling_sucess_rate() > 0.1 && !h_skip_random_paths) {
                                leaf_store_limit *= 2;
                                continue;
                            }
                            if (m_rand.get_rolling_first_level_success_rate() * bfs_cost_estimate > leaf_store_limit &&
                                !h_skip_random_paths) {
                                //std::cout << "increased here " << m_rand.get_rolling_first_level_success_rate() << ", "
                                //          << bfs_cost_estimate << std::endl;
                                leaf_store_limit = std::min(std::max(
                                        (int) m_rand.get_rolling_first_level_success_rate() * bfs_cost_estimate,
                                        ((int) (1.5 * (leaf_store_limit + 1.0)))), h_budget+1);
                                continue;
                            }
                        }

                        if(h_restarts < 2 && h_cost + bfs_cost_estimate > h_budget && !h_skip_random_paths) {
                            fail = true;
                            progress_print("restart", std::to_string(h_cost), std::to_string(h_budget));
                            break;
                        }

                        h_skip_random_paths = false;

                        m_bfs.do_a_level(&m_refinement, g, ir_tree, local_state, current_selector);
                        progress_print("bfs", "0-" + std::to_string(ir_tree.get_finished_up_to()) + "(" +
                                              std::to_string(bfs_cost_estimate) + ")",
                                       std::to_string(ir_tree.get_level_size(ir_tree.get_finished_up_to())));

                        if (ir_tree.get_finished_up_to() == base_size)
                            break;

                        bfs_cost_estimate = m_bfs.next_level_estimate(ir_tree, current_selector);
                        leaf_store_limit += ir_tree.get_level_size(ir_tree.get_finished_up_to());
                    }
                    // breadth-first search & random automorphisms

                    if (!fail) {
                        break;
                    } else {
                        local_state.load_reduced_state(root_save);
                        const int cell_prev = root_save.get_coloring()->cells;

                        if(ir_tree.get_finished_up_to() >= 1) {
                            // TODO: improve saving hash for pruned nodes, then propagate this to first level
                            // TODO: use "pruned" flag, etc.
                            work_list hash(g->v_size);
                            //mark_set  is_pruned(g->v_size);
                            //ir_tree.mark_first_level(is_pruned);

                            for (int i = 0; i < g->v_size; ++i) {
                                //if(is_pruned.get(i)) hash[i] = 1;
                                //else                      hash[i] = 0;
                                hash[i] = (int) (*ir_tree.get_node_invariant())[i] % 16;
                            }
                            for (int i = 0; i < g->v_size; ++i) {
                                hash[i] = 16 * local_state.c->vertex_to_col[i] + hash[i];
                            }
                            g->initialize_coloring(local_state.c, hash.get_array());
                            const int cell_after = local_state.c->cells;
                            progress_print("inproc_bfs", std::to_string(cell_prev),
                                           std::to_string(cell_after));
                            if (cell_after != cell_prev) {
                                local_state.refine(&m_refinement, g);
                                progress_print("inproc_ref", std::to_string(cell_after),
                                               std::to_string(local_state.c->cells));
                            }
                        }

                        m_schreier->determine_save_to_individualize(&save_to_individualize, local_state.get_coloring());

                        if(!save_to_individualize.empty()) {
                            for (int i = 0; i < save_to_individualize.size(); ++i) {
                                const int ind_v   = save_to_individualize[i].first;
                                const int test_col_sz = save_to_individualize[i].second;
                                //std::cout << test_col_sz << std::endl;
                                const int ind_col = local_state.c->vertex_to_col[ind_v];
                                assert(test_col_sz == local_state.c->ptn[ind_col]+1);
                                m_prep.multiply_to_group_size(local_state.c->ptn[ind_col]+1);
                                local_state.move_to_child_no_trace(&m_refinement, g, ind_v);
                            }
                            progress_print("inproc_ind", "_",
                                           std::to_string(local_state.c->cells));
                            save_to_individualize.clear();
                        }

                        local_state.save_reduced_state(root_save);
                        // TODO use orbit partition for 1 individualization
                        // TODO could run sassy if individualization succeeds...

                        inprocessed = cell_prev != local_state.c->cells;
                    }
                }
                if(m_schreier) {
                    add_to_group_size(m_schreier->compute_group_size());
                    //std::cout << "schreier:#symmetries: " << std::to_string(grp_sz_man) << "*10^" << grp_sz_exp << std::endl;
                }
            }

            std::cout << "#symmetries: " << std::to_string(grp_sz_man) << "*10^" << grp_sz_exp << std::endl;
        }
    };
}

#endif //DEJAVU_DEJAVU_H
