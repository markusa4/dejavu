// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_DEJAVU_H
#define DEJAVU_DEJAVU_H

#include "dfs.h"
#include "bfs.h"
#include "rand.h"
#include "sassy/preprocessor.h"
#include "inprocess.h"

// structures for testing
extern        dejavu::ir::refinement test_r;
extern sgraph _test_graph;
extern int*   _test_col;


namespace dejavu {
    void test_hook(int n, const int *p, int nsupp, const int *supp) {
        std::cout << "certifying..." << std::endl;
        assert(test_r.certify_automorphism_sparse(&_test_graph, _test_col, p, nsupp, supp));
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
        groups::shared_schreier* sh_schreier = nullptr;  /**< Schreier-Sims structure */
        ir::shared_tree*         sh_tree     = nullptr;  /**< IR tree                 */

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

        void multiply_to_group_size(std::pair<long double, int> add_grp_sz) {
            multiply_to_group_size(add_grp_sz.first, add_grp_sz.second);
        }

        void multiply_to_group_size(long double add_grp_sz_man, int add_grp_sz_exp) {
            while (add_grp_sz_man >= 10.0) {
                s_grp_sz_exp += 1;
                add_grp_sz_man = add_grp_sz_man / 10;
            }
            s_grp_sz_exp += add_grp_sz_exp;
            s_grp_sz_man *= add_grp_sz_man;
            man_to_exp_group_size();
        }

    public:
        // settings
        int h_error_bound        = 10; /**< error probability is below `(1/2)^h_error_bound`, default value of 10 thus
                                       * misses an automorphism with probabiltiy at most (1/2)^10 < 0.098% */
        int h_limit_leaf        = 0;  /**< limit for the amount of IR leaves to be stored*/
        int h_limit_fail        = 0;  /**< limit for the amount of backtracking allowed*/

        int h_base_max_diff     = 5; /**< only allow a base that is at most `h_base_max_diff` times larger than the
                                       *  previous base */

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
         *
         * @todo no-pruning DFS to fill up available leafs more efficiently
         */
        void automorphisms(sgraph* g, int* colmap = nullptr, dejavu_hook* hook = nullptr) {
            // first, we try to preprocess
            sassy::preprocessor m_prep; /*< initializes the preprocessor */

            // preprocess the graph using sassy
            PRINT("preprocessing...");
            sassy::sgraph gg;
            transfer_sgraph_to_sassy_sgraph(g, &gg);
            m_prep.reduce(&gg, colmap, hook); /*< reduces the graph */
            multiply_to_group_size(m_prep.base, m_prep.exp); /*< group size needed if the
                                                                                           *  early out below is used */

            // early-out if preprocessor finished solving the graph
            if(gg.v_size <= 1) return;

            // if the preprocessor changed the vertex set of the graph, need to use reverse translation
            dejavu_hook dhook = sassy::preprocessor::dejavu_hook;
            if(g->v_size != gg.v_size) hook = &dhook; /*< change hook to sassy hook*/

            // print that we are solving now...
            PRINT(std::endl << "solving..." << std::endl);
            progress_print_header();

            // we first need to set up quite a few datastructures and modules
            // we start by transferring the graph
            transfer_sassy_sgraph_to_sgraph(g, &gg);
            g->dense = !(g->e_size < g->v_size || g->e_size / g->v_size < g->v_size / (g->e_size / g->v_size));

            // settings of heuristics
            int  h_budget           = 1;  /*< budget of current restart iteration    */
            int  h_budget_inc_fact  = 5;  /*< factor used when increasing the budget */

            // statistics used to steer heuristics
            int  s_restarts         = -1;    /*< number of restarts performed                  */
            int  s_cost             = 0;     /*< cost induced by current restart iteration     */
            bool s_long_base        = false; /*< flag for bases that are considered very long  */
            bool s_short_base       = false; /*< flag for bases that are considered very short */
            bool s_few_cells        = false;
            int  s_leaves_added_this_restart = 0;
            bool s_inprocessed = false;     /*< flag which says that we've inprocessed the graph in the last restart */
            int  s_consecutive_discard = 0; /*< count how many bases have been discarded in a row -- prevents edge
                                             *  case where all searchable bases became much larger due to BFS or other
                                             *  changes*/
            bool s_last_bfs_pruned = false; /*< keep track of whether last BFS calculation pruned some nodes */

            // some data we want to keep track of
            std::vector<int> base;       /*< current base of vertices       */
            std::vector<int> base_sizes; /*< current size of colors on base */

            // local modules and workspace, to be used by other modules
            ir::refinement       m_refinement; /*< workspace for color refinement and other utilities */
            ir::selector_factory m_selectors;  /*< cell selector creation                             */
            groups::automorphism_workspace automorphism(g->v_size); /*< workspace to keep an automorphism */
            groups::schreier_workspace     schreierw(g->v_size);    /*< workspace for Schreier-Sims       */

            // shared, global modules
            sh_tree     = new ir::shared_tree(g->v_size); /*< BFS levels, shared leaves, ...           */
            sh_schreier = new groups::shared_schreier();             /*< Schreier structure to sift automorphisms */
            sh_schreier->h_error_bound = h_error_bound;              /*< pass the error bound parameter           */

            // initialize modules for high-level search strategies
            search_strategy::dfs_ir      m_dfs;       /*< depth-first search   */
            search_strategy::bfs_ir      m_bfs;       /*< breadth-first search */
            search_strategy::random_ir   m_rand;      /*< randomized search    */
            search_strategy::inprocessor m_inprocess; /*< inprocessing         */

            // link high-level, local strategy modules to local workspace modules
            m_dfs.link_to_workspace(&automorphism);
            m_bfs.link_to_workspace(&automorphism);
            m_rand.link_to_workspace(&schreierw, &automorphism);

            // initialize a coloring using colors of preprocessed graph
            coloring local_coloring;
            g->initialize_coloring(&local_coloring, colmap);
            const bool s_regular = local_coloring.cells == 1; // Is the graph regular?

            // set up a local state for IR computations
            ir::controller local_state(&m_refinement, &local_coloring); /*< controls movement in IR tree*/

            // save root state for random and BFS search, as well as restarts
            ir::reduced_save root_save;
            local_state.save_reduced_state(root_save); /*< root of the IR tree */
            int s_last_base_size = g->v_size + 1;         /*< v_size + 1 is larger than any actual base*/

            // now that we are set up, let's start solving the graph
            // loop for restarts
            while (true) {
                // "Dry land is not a myth, I've seen it!"
                ++s_restarts;
                if (s_restarts >= 1) { /* no need for this on first iteration of loop */
                    local_state.load_reduced_state(root_save); /*< start over from root */
                    progress_print_split();
                    const int increase_fac = (s_restarts >= 3) ? h_budget_inc_fact : 2;
                    h_budget *= increase_fac;
                    s_cost = 0;

                    if(s_inprocessed && sh_tree) sh_tree->clear_leaves(); // TODO if I link up the corresponding root to each leaf, I can keep using them... also "from BFS" would be more efficient in leaf recovery?
                    m_rand.reset_statistics();

                    s_leaves_added_this_restart = 0;
                }

                // keep vertices which are save to individualize for in-processing
                m_inprocess.inproc_can_individualize.clear();

                // find a selector, moves local_state to a leaf of IR shared_tree
                m_selectors.find_base(g, &local_state, s_restarts);
                auto selector = m_selectors.get_selector_hook();
                progress_print("selector", local_state.s_base_pos,local_state.T->get_position());
                int base_size = local_state.s_base_pos;
                // immediately discard this base if deemed too large, unless we are discarding too often
                if (base_size > h_base_max_diff * s_last_base_size && s_consecutive_discard < 3) {
                    progress_print("skip", local_state.s_base_pos,s_last_base_size);
                    ++s_consecutive_discard;
                    continue;
                }
                s_consecutive_discard = 0;    /*< reset counter since we are not discarding this one   */
                s_last_base_size = base_size; /*< also reset `s_last_base_size` to current `base_size` */

                // determine whether base is "large", or "short"
                s_long_base  = base_size >  sqrt(g->v_size);  /*< a "long"  base  */
                s_short_base = base_size <= 2;                   /*< a "short" base  */

                // make snapshot of trace and leaf, used by following search strategies
                local_state.compare_to_this();
                base       = local_state.base_vertex;   // we want to keep this for later
                base_sizes = local_state.base_color_sz; // used for upper bounds in Schreier structure

                // we first perform a depth-first search, starting from the computed leaf in local_state
                m_dfs.h_recent_cost_snapshot_limit = s_long_base ? 0.33 : 0.25; // set up DFS heuristic
                const int dfs_level =
                        m_dfs.do_dfs(hook, g, root_save.get_coloring(), local_state,
                                     &m_inprocess.inproc_can_individualize);
                progress_print("dfs", std::to_string(base_size) + "-" + std::to_string(dfs_level),
                               "~"+std::to_string((int)m_dfs.grp_sz_man) + "*10^" + std::to_string(m_dfs.grp_sz_exp));
                if (dfs_level == 0) break; // DFS finished the graph -- we are done!

                // next, we go into the random path + BFS algorithm, so we set up a Schreier structure and IR tree for
                // the given base

                // first, the schreier structure
                const bool can_keep_previous = (s_restarts >= 3) && !s_inprocessed;
                const bool reset_prob =
                        sh_schreier->reset(g->v_size, schreierw, base, base_sizes, dfs_level, can_keep_previous,
                                           m_inprocess.inproc_fixed_points);
                if(reset_prob) sh_schreier->reset_probabilistic_criterion(); /*< need to reset probabilistic abort
                                                                             *  criterion*/

                // then, we set up the shared IR tree, which contains computed BFS levels, stored leaves, as well as
                // further gathered data for heuristics
                sh_tree->reset(base, &root_save, can_keep_previous);

                s_inprocessed = false; /*< tracks whether we changed the root of IR tree since we've initialized BFS and
                                        *   Schreier last time */

                // first, we try to add the leaf of the base to the IR tree
                m_rand.specific_walk(g, *sh_tree, local_state, base);
                // TODO find a better way to add canonical leaf to leaf store

                // TODO continue refactor here

                int bfs_cost_estimate = m_bfs.next_level_estimate(*sh_tree, selector);

                const int flat_leaf_store_inc = sh_tree->stat_leaves();
                //leaf_store_limit = std::max(std::min(leaf_store_limit, h_budget/2), 2);

                int s_iterations = 0;
                int h_rand_allowed_fails = 0;
                s_few_cells = root_save.get_coloring()->cells <= 2; /*< flag to tune heuristics */
                if(s_inprocessed) m_rand.reset_statistics();

                // in the following, we will always decide between 3 different strategies: random paths, BFS, or
                // inprocessing followed by a restart
                enum decision {random_ir, bfs_ir, restart}; /*< this enum tracks the decision */
                decision h_next_routine = random_ir;        /*< here, we store the decision   */

                // we perform these strategies in a loop, with the following abort criteria
                bool do_a_restart        = false; /*< we want to do a full restart */
                bool finished_symmetries = false; /*< we found all the symmetries  */

                while (!do_a_restart && !finished_symmetries) {
                    // TODO decision heuristic goes here!
                    // update some estimations
                    bfs_cost_estimate = m_bfs.next_level_estimate(*sh_tree, selector);

                    h_next_routine = restart;

                    bool   s_have_estimate    = (m_rand.s_paths >= 5);
                    bool no_first_level_success = false;
                    double first_level_success = 1;
                    if(s_have_estimate) {
                        no_first_level_success = (m_rand.s_paths_fail1 == m_rand.s_paths);
                        first_level_success = (m_rand.s_paths_fail1*1.0) / (m_rand.s_paths*1.0);
                    }

                    // main decision procedure
                    const bool h_look_close = (m_rand.s_rolling_first_level_success > 0.5) && !s_short_base;
                    if(h_look_close) h_next_routine = bfs_ir;
                    if(first_level_success < 0.1) h_next_routine = bfs_ir;


                    // decision overrides
                    if(s_iterations == 0) {
                        h_next_routine = random_ir; /*< need to gather more information */
                        h_rand_allowed_fails += 5;
                    }
                    if(s_cost > h_budget) h_next_routine = restart;

                    /*
                     *
                     *                         if (m_rand.s_rolling_first_level_success > 0.99) s_few_cells = false;
                        if (s_short_base && s_regular && m_rand.s_rolling_success < 0.01) s_few_cells = true;
                        else if (!(((bfs_cost_estimate <= g->v_size ||
                                     bfs_cost_estimate * m_rand.s_rolling_first_level_success < 12)) &&
                                   m_rand.s_rolling_first_level_success < 0.99) &&
                                 leaf_store_limit <= (bfs_cost_estimate * m_rand.s_rolling_first_level_success) / 2) {
                            if (m_rand.s_rolling_success > 0.1 && !s_few_cells) {
                                leaf_store_limit *= 2;
                                continue;
                            }
                            if (m_rand.s_rolling_first_level_success * bfs_cost_estimate > leaf_store_limit &&
                                !s_few_cells) {
                                leaf_store_limit = std::min(std::max(
                                        (int) m_rand.s_rolling_first_level_success * bfs_cost_estimate,
                                        ((int) (1.5 * (leaf_store_limit + 1.0)))), h_budget + 1);
                                continue;
                            }
                        }

                        if (s_restarts < 2 && s_cost + bfs_cost_estimate > h_budget && !s_few_cells) {
                            do_a_restart = true;
                            progress_print("restart", std::to_string(s_cost), std::to_string(h_budget));
                            break;
                        }
                        s_few_cells = false;
                     */

                    switch(h_next_routine) {
                        case random_ir: {
                            m_rand.setup(h_look_close); // TODO: look_close does not seem worth it

                            const int leaves_pre = sh_tree->stat_leaves();
                            if (sh_tree->get_finished_up_to() == 0) {
                                // random automorphisms, single origin
                                m_rand.random_walks(g, hook, selector, *sh_tree, *sh_schreier, local_state, h_rand_allowed_fails);
                            } else {
                                // random automorphisms, sampled from shared_tree
                                m_rand.random_walks_from_tree(g, hook, selector, *sh_tree, *sh_schreier, local_state, h_rand_allowed_fails);
                            }
                            s_leaves_added_this_restart = sh_tree->stat_leaves() - leaves_pre;

                            progress_print("urandom", sh_tree->stat_leaves(), m_rand.s_rolling_success);
                            progress_print("schreier", "s" + std::to_string(sh_schreier->s_sparsegen()) + "/d" +
                                                       std::to_string(sh_schreier->s_densegen()), "_");

                            if (sh_schreier->deterministic_abort_criterion() ||
                                sh_schreier->probabilistic_abort_criterion()) {
                                finished_symmetries = true;
                            }
                            s_cost += s_leaves_added_this_restart; // TODO: should only add diff
                        }
                        break;
                        case bfs_ir: {
                            bfs_cost_estimate = sh_tree->get_finished_up_to() < base_size?
                                                m_bfs.next_level_estimate(*sh_tree,selector) : -1;
                            m_bfs.do_a_level(g, *sh_tree, local_state, selector);
                            progress_print("bfs", "0-" + std::to_string(sh_tree->get_finished_up_to()) + "(" +
                                                  std::to_string(bfs_cost_estimate) + ")",
                                           std::to_string(sh_tree->get_level_size(sh_tree->get_finished_up_to())));

                            // A bit of a mess here! Manage correct group size whenever BFS finishes the graph. A bit complicated
                            // because of the special code for `base_size == 2`, which can perform automorphism pruning. The
                            // removed automorphisms must be accounted for when calculating the group size.
                            if (sh_tree->get_finished_up_to() == base_size) {
                                finished_symmetries = true;
                                if (base_size == 2) { /*< case for the special code */
                                    m_prep.multiply_to_group_size(
                                            (double) sh_tree->h_bfs_top_level_orbit.orbit_size(base[0]) *
                                            sh_tree->h_bfs_automorphism_pw, 0);
                                } else { /*< base case which just multiplies the remaining elements of the level */
                                    m_prep.multiply_to_group_size(
                                            (double) sh_tree->get_level_size(sh_tree->get_finished_up_to()), 0);
                                }
                            }

                            // if there are less remaining nodes than we expected, we pruned some nodes
                            s_last_bfs_pruned =
                                    sh_tree->get_level_size(sh_tree->get_finished_up_to()) < bfs_cost_estimate;

                            // TODO this should go into decision code?
                            // probably want to inprocess if BFS was successful in pruning on the first level
                            if (sh_tree->get_finished_up_to() == 1 && s_last_bfs_pruned) {
                                do_a_restart = true;
                            }
                        }
                        break;
                        case restart:
                            do_a_restart = true;
                        break;
                    }

                    ++s_iterations;
                }

                // Are we done or just restarting?
                if (finished_symmetries) break; // we are done

                // If not done, we are restarting -- so we try to inprocess using the gathered data
                s_inprocessed = m_inprocess.inprocess(g, sh_tree, sh_schreier, local_state, root_save);
            }

            // We are done! Let's add up the total group size from all the different modules.
            s_grp_sz_man = 1.0;
            s_grp_sz_exp = 0;
            multiply_to_group_size(m_inprocess.s_grp_sz_man, m_inprocess.s_grp_sz_exp);
            multiply_to_group_size(m_prep.base, m_prep.exp);
            multiply_to_group_size(m_dfs.grp_sz_man, m_dfs.grp_sz_exp);
            if(sh_schreier) {
                multiply_to_group_size(sh_schreier->compute_group_size());
                delete sh_schreier; /*< also, clean up allocated schreier structure*/
            }
            delete sh_tree;  /*< ... and also the shared IR tree */
        }
    };
}

#endif //DEJAVU_DEJAVU_H
