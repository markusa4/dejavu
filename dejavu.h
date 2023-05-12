// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_DEJAVU_H
#define DEJAVU_DEJAVU_H

#include "dfs.h"
#include "bfs.h"
#include "rand.h"
#include "preprocess.h"
#include "inprocess.h"

// structures for testing
extern dejavu::ir::refinement test_r;
extern dejavu::sgraph dej_test_graph;
extern int*           dej_test_col;


namespace dejavu {
    void test_hook([[maybe_unused]] int n, [[maybe_unused]] const int *p, [[maybe_unused]] int nsupp,
                   [[maybe_unused]] const int *supp) {
        std::cout << "certifying..." << std::endl;
        assert(test_r.certify_automorphism_sparse(&dej_test_graph, p, nsupp, supp));
        assert(test_r.certify_automorphism(&dej_test_graph, dej_test_col, p));
        assert(test_r.certify_automorphism_sparse(&dej_test_graph, dej_test_col, p, nsupp, supp));
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

    public:
        // settings
        int h_error_bound       = 10; /**< error probability is below `(1/2)^h_error_bound`, default value of 10 thus
                                       * misses an automorphism with probabiltiy at most (1/2)^10 < 0.098% */
        //int h_limit_fail        = 0; /**< limit for the amount of backtracking allowed*/
        int h_base_max_diff     = 5; /**< only allow a base that is at most `h_base_max_diff` times larger than the
                                       *  previous base */

        // statistics
        big_number s_grp_sz; /**< size of the automorphism group computed */
        bool s_deterministic_termination = true; /**< did the last run terminate deterministically? */

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
            m_prep.reduce(g, colmap, hook); /*< reduces the graph */
            s_grp_sz.multiply(m_prep.base, m_prep.exp); /*< group size needed if the
                                                                                    *  early out below is used */

            // early-out if preprocessor finished solving the graph
            if(g->v_size <= 1) return;

            // if the preprocessor changed the vertex set of the graph, need to use reverse translation
            dejavu_hook dhook = sassy::preprocessor::_dejavu_hook;
            hook = &dhook; /*< change hook to sassy hook*/

            // print that we are solving now...
            PRINT("\nsolving...");
            progress_print_header();

            // flag to denote which color refinement version is used
            g->dense = !(g->e_size < g->v_size || g->e_size / g->v_size < g->v_size / (g->e_size / g->v_size));

            // settings of heuristics
            int  h_budget           = 1;  /*< budget of current restart iteration    */
            int  h_budget_inc_fact  = 5;  /*< factor used when increasing the budget */

            // statistics used to steer heuristics
            int  s_restarts         = -1;    /*< number of restarts performed                     */
            int  s_inproc_success   = 0;     /*< how many times inprocessing succeeded            */
            int  s_cost             = 0;     /*< cost induced so far by current restart iteration */
            bool s_long_base        = false; /*< flag for bases that are considered very long     */
            bool s_short_base       = false; /*< flag for bases that are considered very short    */
            bool s_few_cells        = false;
            [[maybe_unused]] bool s_many_cells       = false;
            bool s_inprocessed = false;     /*< flag which says that we've inprocessed the graph since last restart */
            int  s_consecutive_discard = 0; /*< count how many bases have been discarded in a row -- prevents edge
                                             *  case where all searchable bases became much larger due to BFS or other
                                             *  changes */
            bool s_last_bfs_pruned = false; /*< keep track of whether last BFS calculation pruned some nodes */
            big_number s_last_tree_sz;

            // some data we want to keep track of
            std::vector<int> base;       /*< current base of vertices       */
            std::vector<int> base_sizes; /*< current size of colors on base */

            // local modules and workspace, to be used by other modules
            ir::refinement    m_refinement; /*< workspace for color refinement and other utilities */
            ir::base_selector m_selectors;  /*< cell selector creation                             */
            groups::automorphism_workspace automorphism(g->v_size); /*< workspace to keep an automorphism */
            groups::schreier_workspace     schreierw(g->v_size);    /*< workspace for Schreier-Sims   */

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
            local_state.h_deviation_inc = std::min(static_cast<int>(floor(2*sqrt(g->v_size))), 128);

            // save root state for random and BFS search, as well as restarts
            ir::limited_save root_save;
            local_state.save_reduced_state(root_save); /*< root of the IR tree */
            int s_last_base_size = g->v_size + 1;         /*< v_size + 1 is larger than any actual base*/

            // now that we are set up, let's start solving the graph
            // loop for restarts
            while (true) {
                // "Dry land is not a myth, I've seen it!"
                const bool s_hard = h_budget > 10000; /* graph is "hard" */
                const bool s_easy = s_restarts == -1; /* graph is "easy" */

                ++s_restarts; /*< increase the restart counter */
                if(s_restarts > 0) {
                    // now, we manage the restart....
                    local_state.load_reduced_state(root_save); /*< start over from root */
                    progress_print_split();
                    const int increase_fac = (s_restarts >= 3)? h_budget_inc_fact :2;/*<factor to increase the budget */
                    if(s_inprocessed) h_budget = 1;       /*< if we inprocessed, we hope that the graph is easy again */
                    h_budget *= increase_fac;                                                /*< increase the budget! */
                    s_cost = 0;                                                                     /* reset the cost */
                }

                // keep vertices which are save to individualize for in-processing
                m_inprocess.inproc_can_individualize.clear();


                // find a selector, moves local_state to a leaf of IR shared_tree
                m_selectors.find_base(g, &local_state, s_restarts);
                auto selector = m_selectors.get_selector_hook();
                progress_print("sel", local_state.s_base_pos,local_state.T->get_position());
                int base_size = local_state.s_base_pos;

                // TODO if base sizes are equal, we should always prefer the one with smallest number of IR leaves

                // determine whether base is "large", or "short"
                s_long_base  = base_size >  sqrt(g->v_size);  /*< a "long"  base  */
                s_short_base = base_size <= 2;                   /*< a "short" base  */
                s_few_cells  = root_save.get_coloring()->cells <= 2;                  /*< "few"  cells in initial  */
                s_many_cells = root_save.get_coloring()->cells >  sqrt(g->v_size); /*< "many" cells in initial  */

                const bool s_same_length     = base_size == s_last_base_size;

                const bool s_too_long_anyway = base_size > h_base_max_diff * s_last_base_size;
                const bool s_too_long_long   = s_long_base && (base_size >  1.25 * s_last_base_size);
                const bool s_too_long_hard   = s_hard && !s_inprocessed && (base_size > s_last_base_size);
                const bool s_too_long        = s_too_long_anyway || s_too_long_long || s_too_long_hard;

                const bool s_too_big         = s_same_length && (s_last_tree_sz < m_selectors.get_ir_size_estimate());

                // immediately discard this base if deemed too large, unless we are discarding too often
                if ((s_too_big || s_too_long) && s_inproc_success < 2 && s_consecutive_discard < 3) {
                    progress_print("skip", local_state.s_base_pos,s_last_base_size);
                    ++s_consecutive_discard;
                    continue;
                }
                s_consecutive_discard = 0;    /*< reset counter since we are not discarding this one   */
                s_last_base_size = base_size; /*< also reset `s_last_base_size` to current `base_size` */
                s_last_tree_sz = m_selectors.get_ir_size_estimate();

                // make snapshot of trace and leaf, used by following search strategies
                local_state.compare_to_this();
                base       = local_state.base_vertex;   // we want to keep this for later
                base_sizes = local_state.base_color_sz; // used for upper bounds in Schreier structure
                const int s_trace_full_cost = local_state.T->get_position(); /*< total trace cost of this base */

                // we first perform a depth-first search, starting from the computed leaf in local_state
                m_dfs.h_recent_cost_snapshot_limit = s_long_base ? 0.33 : 0.25; // set up DFS heuristic
                const int dfs_level =
                        m_dfs.do_dfs(hook, g, root_save.get_coloring(), local_state,
                                     &m_inprocess.inproc_can_individualize);
                progress_print("dfs", std::to_string(base_size) + "-" + std::to_string(dfs_level),
                               "~"+std::to_string((int)m_dfs.s_grp_sz.mantissa) + "*10^" + std::to_string(m_dfs.s_grp_sz.exponent));
                if (dfs_level == 0) break; // DFS finished the graph -- we are done!
                const bool s_dfs_backtrack = m_dfs.s_termination == search_strategy::dfs_ir::termination_reason::r_fail;

                // next, we go into the random path + BFS algorithm, so we set up a Schreier structure and IR tree for
                // the given base

                // first, the schreier structure
                const bool can_keep_previous = (s_restarts >= 3) && !s_inprocessed;
                const bool reset_prob =
                        sh_schreier->reset(g->v_size, schreierw, base, base_sizes, dfs_level,
                                           can_keep_previous,m_inprocess.inproc_fixed_points);
                if(reset_prob) sh_schreier->reset_probabilistic_criterion(); /*< need to reset probabilistic abort
                                                                             *  criterion after inprocessing          */
                m_rand.reset_statistics();

                // then, we set up the shared IR tree, which contains computed BFS levels, stored leaves, as well as
                // further gathered data for heuristics
                sh_tree->reset(base, &root_save, can_keep_previous);
                if(s_inprocessed) sh_tree->clear_leaves();

                s_inprocessed = false; /*< tracks whether we changed the root of IR tree since we've initialized BFS and
                                        *   Schreier last time */

                // first, we try to add the leaf of the base to the IR tree
                m_rand.specific_walk(g, *sh_tree, local_state, base);
                // TODO find a better way to add canonical leaf to leaf store

                s_last_bfs_pruned = false;
                int s_bfs_next_level_nodes = m_bfs.next_level_estimate(*sh_tree, selector);
                int h_rand_fail_lim_total  = 0;
                int h_rand_fail_lim_now    = 0;

                // in the following, we will always decide between 3 different strategies: random paths, BFS, or
                // inprocessing followed by a restart
                enum decision {random_ir, bfs_ir, restart}; /*< this enum tracks the decision */
                decision h_next_routine = random_ir;        /*< here, we store the decision   */

                // we perform these strategies in a loop, with the following abort criteria
                bool do_a_restart        = false; /*< we want to do a full restart */
                bool finished_symmetries = false; /*< we found all the symmetries  */

                h_rand_fail_lim_now = 4;

                while (!do_a_restart && !finished_symmetries) {
                    // What do we do next? Random search, BFS, or a restart?
                    // here are our decision heuristics (AKA dark magic)

                    // first, we update our estimations / statistics
                    s_bfs_next_level_nodes = m_bfs.next_level_estimate(*sh_tree, selector);

                    const bool s_have_rand_estimate   = (m_rand.s_paths >= 4);  /*< for some estimations, we need a few
                                                                                 *  random paths */
                    double s_path_fail1_avg  = 0;
                    double s_trace_cost1_avg = s_trace_full_cost;

                    int s_random_path_trace_cost = s_trace_full_cost - sh_tree->get_current_level_tracepos();
                    if(s_have_rand_estimate) {
                        s_path_fail1_avg    = (double) m_rand.s_paths_fail1 / (double) m_rand.s_paths;
                        s_trace_cost1_avg   = (double) m_rand.s_trace_cost1 / (double) m_rand.s_paths;
                    }
                    const double reset_cost_rand = g->v_size;
                    const double reset_cost_bfs  = std::min(s_trace_cost1_avg, (double) g->v_size);
                    double s_bfs_cost_estimate  = (s_trace_cost1_avg        + reset_cost_bfs) * (s_bfs_next_level_nodes);
                    double s_rand_cost_estimate = (s_random_path_trace_cost + reset_cost_rand) * h_rand_fail_lim_now ;

                    const bool h_look_close = (m_rand.s_rolling_first_level_success > 0.5) && !s_short_base;

                    // using this data, we now make a decision
                    h_next_routine = restart; /*< undecided? do a restart */

                    if(!s_have_rand_estimate) { /*< don't have an estimate, so let's poke a bit with random! */
                        h_next_routine = random_ir; /*< need to gather more information */
                    } else {
                        // now that we have some data, we attempt to model how effective and costly random search and
                        // BFS is, to then make an informed decision of what to do next

                        // we do so by negatively scoring each method: higher score, worse technique
                        //double score_rand = s_trace_full_cost * h_rand_fail_lim_now * (1-m_rand.s_rolling_success);
                        double score_rand = s_rand_cost_estimate * (1-m_rand.s_rolling_success);
                        double score_bfs  = s_bfs_cost_estimate * (0.1 + 1-s_path_fail1_avg);

                        // we make some adjustments to try to model effect of techniques better
                        // increase BFS score if BFS does not prune nodes on the next level -- we want to be somewhat
                        // more reluctant to perform BFS in this case
                        if(s_path_fail1_avg < 0.01) score_bfs *= 2;

                        // we decrease the BFS score if we are beyond the first level, in the hope that this models
                        // the effect of trace deviation maps
                        if(sh_tree->get_finished_up_to() >= 1) score_bfs *= (1-s_path_fail1_avg);

                        // we make a decision...
                        h_next_routine    = (score_rand < score_bfs)? random_ir : bfs_ir;

                        // if we do random_ir next, increase its budget
                        h_rand_fail_lim_now = h_next_routine == random_ir? 2*h_rand_fail_lim_now : h_rand_fail_lim_now;
                    }

                    // we override the above decisions in specific cases...
                    if(h_next_routine == bfs_ir && s_bfs_next_level_nodes * (1-s_path_fail1_avg) > 2*h_budget)
                            h_next_routine = restart; /* best decision would be BFS, but it would exceed the budget a
                                                       * lot! */
                    if(s_cost > h_budget) h_next_routine = restart; /*< we exceeded our budget, restart */
                    if(s_dfs_backtrack && s_regular && s_few_cells && s_restarts == 0 &&
                       s_path_fail1_avg > 0.01 && sh_tree->get_finished_up_to() == 0)
                        h_next_routine = bfs_ir; /*< surely BFS will help in this case, so let's fast-track */

                    // ... but if we are "almost done" with random search... we stretch the budget and continue!
                    if(search_strategy::random_ir::h_almost_done(*sh_schreier)) h_next_routine = random_ir;
                    if(m_rand.s_rolling_success > 0.1 && s_cost <= h_budget * 4)   h_next_routine = random_ir;
                    if(s_hard && m_rand.s_succeed >= 1 && s_cost <= h_budget * 10) h_next_routine = random_ir;

                    // immediately inprocess if BFS was successful in pruning on the first level
                    if (sh_tree->get_finished_up_to() == 1 && s_last_bfs_pruned) h_next_routine = restart;

                    switch(h_next_routine) {
                        case random_ir: {
                            h_rand_fail_lim_total += h_rand_fail_lim_now;
                            m_rand.setup(h_look_close); // TODO: look_close does not seem worth it
                            m_rand.h_sift_random = !s_easy;

                            if (sh_tree->get_finished_up_to() == 0) {
                                // random automorphisms, single origin
                                m_rand.random_walks(g, hook, selector, *sh_tree, *sh_schreier, local_state,
                                                    h_rand_fail_lim_total);
                            } else {
                                // random automorphisms, sampled from shared_tree
                                m_rand.random_walks_from_tree(g, hook, selector, *sh_tree, *sh_schreier,
                                                              local_state, h_rand_fail_lim_total);
                            }

                            progress_print("urandom", sh_tree->stat_leaves(), m_rand.s_rolling_success);
                            if(sh_schreier->s_densegen() + sh_schreier->s_sparsegen() > 0) {
                                progress_print("schreier", "s" + std::to_string(sh_schreier->s_sparsegen()) + "/d" +
                                                           std::to_string(sh_schreier->s_densegen()), "_");
                            }

                            finished_symmetries = sh_schreier->deterministic_abort_criterion() ||
                                                  sh_schreier->probabilistic_abort_criterion();
                            s_deterministic_termination = !(sh_schreier->probabilistic_abort_criterion() &&
                                                            !sh_schreier->deterministic_abort_criterion());
                            s_cost += h_rand_fail_lim_now;
                        }
                        break;
                        case bfs_ir: {
                            // perform one level of BFS
                            m_bfs.h_use_deviation_pruning = !((s_inproc_success >= 2) && s_path_fail1_avg > 0.1);
                            m_bfs.do_a_level(g, *sh_tree, local_state, selector);
                            progress_print("bfs", "0-" + std::to_string(sh_tree->get_finished_up_to()) + "(" +
                                                  std::to_string(s_bfs_next_level_nodes) + ")",
                                           std::to_string(sh_tree->get_level_size(sh_tree->get_finished_up_to())));

                            // A bit of a mess here! Manage correct group size whenever BFS finishes the graph. A bit
                            // complicated because of the special code for `base_size == 2`, which can perform
                            // automorphism pruning. (The removed automorphisms must be accounted for when calculating
                            // the group size.)
                            if (sh_tree->get_finished_up_to() == base_size) {
                                finished_symmetries = true;
                                if (base_size == 2) { /*< case for the special code */
                                    m_prep.multiply_to_group_size(
                                            (double) sh_tree->h_bfs_top_level_orbit.orbit_size(base[0]) *
                                            sh_tree->h_bfs_automorphism_pw, 0);
                                } else if (base_size == 1) {
                                    m_prep.multiply_to_group_size(
                                            (double) sh_tree->h_bfs_top_level_orbit.orbit_size(base[0]), 0);
                                } else { /*< base case which just multiplies the remaining elements of the level */
                                    m_prep.multiply_to_group_size(
                                            (double) sh_tree->get_level_size(sh_tree->get_finished_up_to()), 0);
                                }
                            }

                            // if there are less remaining nodes than we expected, we pruned some nodes
                            s_last_bfs_pruned = sh_tree->get_current_level_size() < s_bfs_next_level_nodes;
                            m_rand.reset_statistics();
                            s_cost += sh_tree->get_current_level_size();
                        }
                        break;
                        case restart:
                            do_a_restart = true;
                        break;
                    }
                }

                // Are we done or just restarting?
                if (finished_symmetries) {
                    sh_schreier->compute_group_size(); // need to compute the group size now
                    break; // we are done
                }

                // we are restarting -- so we try to inprocess using the gathered data
                s_inprocessed = m_inprocess.inprocess(g, sh_tree, sh_schreier, local_state, root_save);
                s_inproc_success += s_inprocessed;

                // edge case where inprocessing might finish the graph
                if(root_save.get_coloring()->cells == g->v_size) {
                    finished_symmetries = true;
                    std::cout << "this" << std::endl;
                    break;
                }
            }

            // We are done! Let's add up the total group size from all the different modules.
            s_grp_sz.mantissa = 1.0;
            s_grp_sz.exponent = 0;
            s_grp_sz.multiply(m_inprocess.s_grp_sz);
            s_grp_sz.multiply(m_prep.base, m_prep.exp);
            s_grp_sz.multiply(m_dfs.s_grp_sz);
            s_grp_sz.multiply(sh_schreier->s_grp_sz);

            delete sh_schreier; /*< clean up allocated schreier structure */
            delete sh_tree;     /*< ... and also the shared IR tree       */
        }
    };
}

#endif //DEJAVU_DEJAVU_H
