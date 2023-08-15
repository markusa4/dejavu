// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_DFS_H
#define DEJAVU_DFS_H

#include <random>
#include <chrono>
#include "refinement.h"
#include "coloring.h"
#include "sgraph.h"
#include "trace.h"
#include "groups.h"
#include "ir.h"

namespace dejavu {

    /**
     * \brief High-level IR search strategies.
     *
     * Contains different high-level individualization-refinement search strategies. In particular, depth-first search,
     * breadth-first search as well as different methods of random search.
     */
    namespace search_strategy {

        /**
         * \brief Depth-first search without backtracking.
         *
         * Depth-first IR search which does not backtrack. Can parallelize along the base.
         * Due to the search not back-tracking, this module can not deal with difficult parts of combinatorial graphs.
         */
        class dfs_ir {
            mark_set singleton_diffs;
            int      singleton_diff = 0;
            int cost_snapshot = 0; /**< used to track cost-based abort criterion */
            groups::automorphism_workspace* ws_automorphism = nullptr;

        public:
            enum termination_reason {r_none, r_fail, r_cost};
            termination_reason s_termination = r_none; /**< Why did we stop performing DFS? Cost too high, or did we
                                                         *  have to backtrack? */

            double h_recent_cost_snapshot_limit = 0.25; /**< A float in the range [0-1]. Limits continuation of DFS
                                                          * search to whethercomputing recent elements only cost this
                                                          * fraction of the cost of an entire root-to-leaf walk. */
            big_number s_grp_sz; /**< group size */

            work_list check_color_deg1;
            work_list revert1;
            work_list check_color_deg2;
            work_list revert2;

            bool assert_equitable(sgraph *g, coloring* c) {
                check_color_deg1.resize(g->v_size);
                check_color_deg2.resize(g->v_size);
                revert1.resize(g->v_size);
                revert2.resize(g->v_size);
                for(int i = 0; i < c->domain_size; ++i) {
                    check_color_deg1[i] = 0;
                    check_color_deg2[i] = 0;
                }

                revert1.reset();
                revert2.reset();

                for(int i = 0; i < c->domain_size;) {
                    const int col = i;
                    const int col_sz = c->ptn[i] + 1;
                    if(col_sz > 1) {
                        const int vert = c->lab[col + 0];
                        assert(revert1.cur_pos == 0);
                        for(int k = 0; k < g->d[vert]; ++k) {
                            const int neighbour = g->e[g->v[vert] + k];
                            const int neighbour_col = c->vertex_to_col[neighbour];
                            if(check_color_deg1[neighbour_col] == 0) revert1.push_back(neighbour_col);
                            ++check_color_deg1[neighbour_col];
                        }

                        for (int j = 1; j < col_sz; ++j) {
                            const int vert =c->lab[col + j];
                            assert(revert2.cur_pos == 0);
                            for(int k = 0; k < g->d[vert]; ++k) {
                                const int neighbour     = g->e[g->v[vert] + k];
                                const int neighbour_col = c->vertex_to_col[neighbour];
                                if(check_color_deg2[neighbour_col] == 0) revert2.push_back(neighbour_col);
                                ++check_color_deg2[neighbour_col];
                            }
                            for(int k = 0; k < revert2.cur_pos; ++k) {
                                const int neighbour_col = revert2[k];
                                assert(check_color_deg1[neighbour_col] == check_color_deg2[neighbour_col]);
                            }
                            for(int k = 0; k < revert2.cur_pos; ++k) {
                                const int neighbour_col = revert2[k];
                                check_color_deg2[neighbour_col] = 0;
                            }
                            revert2.reset();
                        }
                        for(int k = 0; k < revert1.cur_pos; ++k) {
                            const int neighbour_col = revert1[k];
                            check_color_deg1[neighbour_col] = 0;
                        }
                        revert1.reset();
                    }
                    i += col_sz;
                }
                return true;
            }


            void link_to_workspace(groups::automorphism_workspace* automorphism) {
                ws_automorphism = automorphism;
            }
            /**
             * Recurses to an equal leaf (according to trace compared to within \p state), unless backtracking is
             * necessary.
             *
             * @param g The graph.
             * @param initial_colors The initial coloring of the graph \p g.
             * @param state IR controller used for computations. On successful return, \p state will be in a leaf of the
             * IR tree equal to the stored trace.
             * @param automorphism If an equal leaf is found, \p automorphism might contain an automorphism mapping the
             * found leaf to the stored leaf.
             * @return
             *
             * @todo should just always certify the leaf...
             */
            std::pair<bool, bool> recurse_to_equal_leaf(sgraph *g, coloring *initial_colors, ir::controller *state) {
                bool prev_fail = false;
                int prev_fail_pos = -1;
                int cert_pos = 0;

                while ((size_t) state->s_base_pos < state->compare_base_color.size()) {
                    const int col    = state->compare_base_color[state->s_base_pos]; // the base color
                    const int col_sz = state->c->ptn[col] + 1;

                    if (col_sz < 2) return {false, false}; // color of original base is trivial here, abort
                    int ind_v = state->c->lab[col];
                    if(state->c->vertex_to_col[state->base_vertex[state->s_base_pos]] == col)
                        ind_v = state->base_vertex[state->s_base_pos];

                    // individualize a vertex of base color
                    state->move_to_child(g, ind_v);

                    const int pos_start = state->base_singleton_pt[state->base_singleton_pt.size() - 1];
                    const int pos_end   = (int) state->singletons.size();

                    // write sparse automorphism with singletons derived so far
                    ws_automorphism->write_singleton(&state->compare_singletons, &state->singletons,
                                                 pos_start, pos_end);

                    singleton_diff += write_singleton_diffs(&singleton_diffs, &state->compare_singletons,
                                                            &state->singletons, pos_start,
                                                            pos_end);

                    if(singleton_diff == 0) {
                        const bool found_auto = state->certify(g, *ws_automorphism);
                        if(found_auto) return {true, true};
                        else {
                            //std::cout << "unnecessary check" << std::endl;
                        }
                    }
                }
                if (state->c->cells == g->v_size) {
                    return {true, false}; // we reached a leaf that looks equivalent, but could not certify
                                                // automorphism
                } else {
                    return {false, false}; // failed to reach equivalent leaf
                }
            }

            int __attribute__((noinline))  paired_recurse_to_equal_leaf(sgraph *g, ir::controller *state_left, ir::controller *state_right) {
                //const bool is_diffed_pre = state_left->update_diff_last_individualization(state_right);
                const bool is_diffed_pre = state_right->update_diff_vertices_last_individualization(state_left);
                if(!is_diffed_pre) {
                    std::cout << "diffed pre" << std::endl;
                    return 2;
                }
                if(state_left->get_diff_diverge()) {
                    std::cout << "diff diverge" << std::endl;
                    return 0;
                }

                //assert_equitable(g, state_right->c);
                //assert_equitable(g, state_left->c);

                int depth = 0;
                while ((size_t) state_right->c->cells < g->v_size) {
                    // pick a color that prevents the "matching OPP"

                    // pick vertices that differ in the color
                    int ind_v_right, ind_v_left;
                    //std::tie(ind_v_left, ind_v_right) = state_right->color_diff_pair(state_left, col);
                    std::tie(ind_v_left, ind_v_right) = state_right->diff_pair(state_left);

                    const int col = state_left->c->vertex_to_col[ind_v_left];
                    //const int col = state_left->get_diff_color();
                    ++depth;

                    assert(col >= 0);
                    assert(col < g->v_size);
                    const int col_sz_right = state_right->c->ptn[col] + 1;
                    const int col_sz_left  = state_left->c->ptn[col]  + 1;

                    assert(col_sz_right > 1);
                    assert(col_sz_left  > 1);

                    if (col_sz_right < 2 || col_sz_left != col_sz_right) return 0;

                    assert(ind_v_left >= -1);
                    assert(ind_v_right >= -1);

                    assert(state_left->c->vertex_to_col[ind_v_left] == col);
                    assert(state_right->c->vertex_to_col[ind_v_right] == col);
                    assert(state_left->c->vertex_to_col[ind_v_right] != col);
                    //assert(state_right->c->vertex_to_col[ind_v_left] != col);

                    state_right->use_trace_early_out(false);
                    state_left->use_trace_early_out(false);
                    state_right->T->set_hash(0);
                    state_left->T->set_hash(0);
                    state_right->T->reset_trace_unequal();
                    state_left->T->reset_trace_unequal();

                    // individualize a vertex of base color
                    state_right->move_to_child(g, ind_v_right);
                    state_left->move_to_child(g,  ind_v_left);

                    // invariant differs
                    if(state_right->T->get_hash() != state_left->T->get_hash()) return 0;

                    // can not perform diff udpates on this (AKA invariant should have differed, but is not
                    // guaranteed to differ always)
                    //const bool is_diffed = state_left->update_diff_last_individualization(state_right);
                    const bool is_diffed = state_right->update_diff_vertices_last_individualization(state_left);
                    if(state_left->get_diff_diverge()) return 0;

                    // no difference? check for automorphism now -- if it doesn't succeed there is no hope on this path
                    if(!is_diffed) {
                        ws_automorphism->reset();
                        state_right->singleton_automorphism(state_left, ws_automorphism);
                        const bool found_auto = state_right->certify(g, *ws_automorphism);
                        return found_auto;
                    }

                    assert(state_right->c->cells < g->v_size); // there can not be diff on a leaf
                }

                // unreachable
                assert(false);
                return false;
            }

            int write_singleton_diffs(mark_set* singleton_diffs, const std::vector<int> *singletons1, const std::vector<int> *singletons2,
                                 const int pos_start, const int pos_end) {
                int diff = 0;
                for (int i = pos_start; i < pos_end; ++i) {
                    if(!singleton_diffs->get((*singletons1)[i])) {
                        singleton_diffs->set((*singletons1)[i]);
                        ++diff;
                    } else {
                        singleton_diffs->unset((*singletons1)[i]);
                        --diff;
                    }
                    if(!singleton_diffs->get((*singletons2)[i])) {
                        singleton_diffs->set((*singletons2)[i]);
                        ++diff;
                    } else {
                        singleton_diffs->unset((*singletons2)[i]);
                        --diff;
                    }
                }
                return diff;
            }

            int do_paired_dfs(dejavu_hook* hook, sgraph *g, coloring *initial_colors, ir::controller &state_left,
                              ir::controller &state_right, std::vector<std::pair<int, int>>* save_to_individualize) {
                s_grp_sz.mantissa = 1.0;
                s_grp_sz.exponent = 0;

                // orbit algorithm structure
                groups::orbit orbs(g->v_size);

                // automorphism workspace
                ws_automorphism->reset();

                // saucy strategy of paired states
                state_left.copy(&state_right);

                // we want to terminate if things become to costly, we save the current trace position to track cost
                cost_snapshot = state_right.T->get_position();

                // tell the controller we are performing DFS now
                state_right.mode_compare_base();

                // abort criteria
                double recent_cost_snapshot = 0;
                bool   fail = false;
                int    trace_pos_reset = 0;

                // loop that serves to optimize Tinhofer graphs
                while (recent_cost_snapshot < h_recent_cost_snapshot_limit && state_right.s_base_pos > 0 && !fail) {
                    // backtrack one level
                    state_right.move_to_parent();
                    if((state_right.s_base_pos & 0x00000FFF) == 0)
                        std::cout << "                  " << ">partial-cost " << recent_cost_snapshot << "/" << h_recent_cost_snapshot_limit << std::endl;

                    // remember which color we are individualizing
                    const int col         = state_right.base_color[state_right.s_base_pos];
                    const int col_sz      = state_right.c->ptn[col] + 1;
                    const int base_vertex = state_right.base_vertex[state_right.s_base_pos]; // base vertex

                    int prune_cost_snapshot = 0;

                    // iterate over current color class
                    //for (int i = col_sz - 1; i >= 0; --i) {
                    for (int i = 0; i < col_sz; ++i) {
                        const int ind_v = state_right.leaf_color.lab[col + i];
                        assert(state_right.c->vertex_to_col[base_vertex] == state_right.c->vertex_to_col[ind_v]);

                        // only consider every orbit once
                        if (ind_v == base_vertex || !orbs.represents_orbit(ind_v)) continue;
                        if (orbs.are_in_same_orbit(ind_v, base_vertex))            continue;

                        // track cost of this refinement for whatever is to come
                        const int cost_start = state_right.T->get_position();

                        // put the "left" state in the correct base position
                        while (state_left.s_base_pos > state_right.s_base_pos) state_left.move_to_parent();
                        state_left.T->set_position(state_right.T->get_position());
                        assert(state_left.c->vertex_to_col[base_vertex] == state_left.c->vertex_to_col[ind_v]);
                        assert(state_left.s_base_pos == state_right.s_base_pos);
                        assert(state_left.T->get_position() == state_right.T->get_position());

                        // reset difference, since state_left and state_right are in the same node of the tree
                        // now -- so there is no difference
                        state_right.reset_diff();

                        // perform individualization-refinement in state_right
                        //TODO make a mode for this in the IR controller module
                        const int prev_base_pos = state_right.s_base_pos;
                        trace_pos_reset = state_right.T->get_position(); // TODO is there an elegant solution to this?
                        state_right.T->reset_trace_equal();
                        state_right.T->set_hash(0);
                        state_right.use_trace_early_out(true);
                        state_right.move_to_child(g, ind_v);
                        bool pruned     = !state_right.T->trace_equal(); // TODO keep track of cost here as well!

                        bool found_auto = false;
                        assert(state_right.c->vertex_to_col[ind_v] == state_right.leaf_color.vertex_to_col[base_vertex]);

                        if(!pruned) {
                            // write singleton diff into automorphism...
                            const int wr_pos_st  = state_right.base_singleton_pt[state_right.base_singleton_pt.size() - 1];
                            const int wr_pos_end = (int) state_right.singletons.size();

                            ws_automorphism->reset();
                            // ... and then check whether this implies a (sparse) automorphism
                            ws_automorphism->write_singleton(&state_right.compare_singletons,
                                                             &state_right.singletons, wr_pos_st,
                                                             wr_pos_end);
                            found_auto = state_right.certify(g, *ws_automorphism);

                            assert(ws_automorphism->perm()[base_vertex] == ind_v);
                            // if no luck with sparse automorphism, try more proper walk to leaf node
                            if (!found_auto) {
                                ws_automorphism->reset();

                                //TODO make a mode for this in the IR controller module
                                state_left.T->reset_trace_equal();
                                state_left.T->set_hash(0);
                                state_left.T->set_position(trace_pos_reset);
                                state_left.use_trace_early_out(true);
                                state_left.move_to_child(g, base_vertex); // need to move left to base vertex
                                assert(state_left.T->trace_equal());
                                //assert_equitable(g, state_left.c);
                                auto return_code = paired_recurse_to_equal_leaf(g, &state_left, &state_right);
                                found_auto = (return_code == 1);
                                pruned     = (return_code == 2);
                            }
                        }

                        // track cost-based abort criterion
                        const int cost_end   = state_right.T->get_position();
                        double cost_partial  = (cost_end - cost_start) / (cost_snapshot*1.0);
                        recent_cost_snapshot = (cost_partial + recent_cost_snapshot * 3) / 4;
                        prune_cost_snapshot += pruned?(cost_end - cost_start):0;

                        // if we found automorphism, add to orbit and call hook
                        if (found_auto) {
                            assert(ws_automorphism->nsupport() > 0);
                            assert(ws_automorphism->perm()[ind_v] == base_vertex || ws_automorphism->perm()[base_vertex] == ind_v);
                            orbs.add_automorphism_to_orbit(*ws_automorphism);
                            if(hook) (*hook)(g->v_size, ws_automorphism->perm(), ws_automorphism->nsupport(),
                                             ws_automorphism->support());
                        }
                        ws_automorphism->reset();

                        // move state back up where we started in this iteration
                        while (prev_base_pos < state_right.s_base_pos) {
                            state_right.move_to_parent();
                        }
                        state_right.T->set_position(trace_pos_reset);

                        // if no automorphism could be determined we would have to backtrack -- so stop!
                        if ((!found_auto && !pruned) || ((4*prune_cost_snapshot > cost_snapshot) &&
                            (prev_base_pos > 0))) {
                            fail = true;
                            break;
                        }

                        // if orbit size equals color size now, we are done on this DFS level
                        if (orbs.orbit_size(base_vertex) == col_sz) break;
                    }

                    // if we did not fail, accumulate size of current level to group size
                    if (!fail) {
                        save_to_individualize->emplace_back(state_right.base_vertex[state_right.s_base_pos],
                                                            orbs.orbit_size(base_vertex));
                        s_grp_sz.multiply(orbs.orbit_size(base_vertex));
                    }
                }

                // set reason for termination
                s_termination = fail? r_fail : r_cost;

                // if DFS failed on current level, we did not finish the current level -- has to be accounted for
                return state_right.s_base_pos + (fail);
            }
        };
    }
}

#endif //DEJAVU_DFS_H
