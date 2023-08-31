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
            int cost_snapshot = 0; /**< used to track cost-based abort criterion */
            groups::automorphism_workspace& ws_automorphism;

        public:
            enum termination_reason {r_none, r_fail, r_cost};
            termination_reason s_termination = r_none; /**< Why did we stop performing DFS? Cost too high, or did we
                                                         *  have to backtrack? */

            double h_recent_cost_snapshot_limit = 0.25; /**< A float in the range [0-1]. Limits continuation of DFS
                                                          * search to whethercomputing recent elements only cost this
                                                          * fraction of the cost of an entire root-to-leaf walk. */
            big_number s_grp_sz; /**< group size */

            explicit dfs_ir(groups::automorphism_workspace& automorphism) : ws_automorphism(automorphism) {}

            int paired_recurse_to_equal_leaf(sgraph* g, ir::controller& state_left, ir::controller& state_right,
                                             int max_depth, bool recurse=false) {
                if(max_depth == 0) return paired_recurse_to_equal_leaf(g, state_left, state_right, recurse);

                if(!recurse) {
                    const bool is_diffed_pre = state_right.update_diff_vertices_last_individualization(state_left);
                    if (!is_diffed_pre || state_right.get_diff_diverge()) return 2;
                }

                // pick vertices that differ in the color
                int ind_v_right, ind_v_left;
                std::tie(ind_v_left, ind_v_right) = state_right.diff_pair(state_left);

                // pick a color that prevents the "matching OPP"
                const int col          = state_left.c->vertex_to_col[ind_v_left];
                const int col_sz_right = state_right.c->ptn[col] + 1;
                const int col_sz_left  = state_left.c->ptn[col]  + 1;

                assert(ind_v_left >= -1);
                assert(ind_v_right >= -1);

                assert(state_left.c->vertex_to_col[ind_v_left] == col);
                assert(state_right.c->vertex_to_col[ind_v_right] == col);
                assert(state_left.c->vertex_to_col[ind_v_right] != col);

                if (col_sz_right < 2 || col_sz_left != col_sz_right) return 2;

                const int base_pos = state_left.s_base_pos;

                std::vector<int> color;
                bool first_ind = true;

                int diff_verts_initial = state_right.diff_num();
                int cells_initial = state_right.c->cells;
                //std::cout << "initial:" << diff_verts_initial << std::endl;

                state_left.use_trace_early_out(false);
                state_left.T->reset_trace_unequal();
                state_left.T->set_hash(0);
                state_left.move_to_child(g,  ind_v_left);
                unsigned long saved_hash = state_left.T->get_hash();

                for(int i = 0; i < col_sz_left; ++i) {
                    while(state_right.s_base_pos > base_pos) {
                        state_right.revert_diff_last_individualization();
                        state_right.move_to_parent();
                    }
                    if(diff_verts_initial != state_right.diff_num()){
                        std::cout << diff_verts_initial << "/" << state_right.diff_num() << std::endl;
                    }
                    assert(diff_verts_initial == state_right.diff_num());
                    assert(cells_initial == state_right.c->cells);
                    while(state_left.s_base_pos > base_pos + 1) {
                        state_left.move_to_parent();
                    }

                    state_right.use_trace_early_out(false);
                    state_right.T->set_hash(0);
                    state_right.T->reset_trace_unequal();

                    int ind_v_right_now;
                    if(first_ind) {
                        ind_v_right_now = ind_v_right;
                        first_ind = false;
                    } else if(color.size() == 0) {
                        color.reserve(col_sz_left);
                        color.push_back(ind_v_right);
                        for(int j = 0; j < col_sz_left; ++j) {
                            const int v = state_right.c->lab[col + j];
                            if(v != ind_v_right) color.push_back(v);
                        }
                        ind_v_right_now = color[i];
                    } else {
                        ind_v_right_now = color[i];
                    }

                    // individualize a vertex
                    state_right.move_to_child(g, ind_v_right_now);
                    assert(state_left.s_base_pos == state_right.s_base_pos);

                    // can not perform diff udpates on this (AKA invariant should have differed, but is not
                    // guaranteed to differ always)

                    // invariant differs
                    if(state_right.T->get_hash() != saved_hash) {
                        state_right.move_to_parent();
                        continue;
                    }

                    const bool is_diffed = state_right.update_diff_vertices_last_individualization(state_left);
                    //std::cout << "->" << state_right.diff_num() << std::endl;

                    if(state_right.get_diff_diverge()) {
                        //std::cout << "diff diverge" << std::endl;
                        continue;
                    }

                    // no difference? check for automorphism now -- if it doesn't succeed there is no hope on this path
                    if(!is_diffed) {
                        //std::cout << state_right.s_base_pos << std::endl;
                        ws_automorphism.reset();
                        state_right.singleton_automorphism(state_left, ws_automorphism);
                        const bool found_auto = state_right.certify(g, ws_automorphism);
                        if(found_auto) {
                            return true;
                        } else {
                            continue;
                        }
                    }

                    if(state_right.diff_num() <= 0) std::cout << state_right.diff_num() << std::endl;
                    assert(state_right.diff_num() > 0);
                    //std::cout << "recursive call" << max_depth - 1 << std::endl;
                    int res = paired_recurse_to_equal_leaf(g, state_left, state_right, max_depth - 1, true);
                    //std::cout << "comeback from recursion" << max_depth - 1 << std::endl;
                    if(res == 0) {
                        //std::cout << "failed at " << max_depth << std::endl;
                        return 0;
                    } else if (res == 1) {
                        //std::cout << "fail" << std::endl;
                        return 1;
                    }
                }

                //std::cout << "end loop"<< std::endl;

                // TODO something here
                return 2;
            }

            int paired_recurse_to_equal_leaf(sgraph* g, ir::controller& state_left, ir::controller& state_right,
                                             bool recurse=false) {
                if(!recurse) {
                    const bool is_diffed_pre = state_right.update_diff_vertices_last_individualization(state_left);
                    if (!is_diffed_pre || state_right.get_diff_diverge()) return 2;
                }

                while (state_right.c->cells < g->v_size) {
                    // pick a color that prevents the "matching OPP"

                    // pick vertices that differ in the color
                    int ind_v_right, ind_v_left;
                    std::tie(ind_v_left, ind_v_right) = state_right.diff_pair(state_left);

                    const int col = state_left.c->vertex_to_col[ind_v_left];

                    assert(col >= 0);
                    assert(col < g->v_size);
                    const int col_sz_right = state_right.c->ptn[col] + 1;
                    const int col_sz_left  = state_left.c->ptn[col]  + 1;

                    assert(col_sz_right > 1);
                    assert(col_sz_left  > 1);

                    if (col_sz_right < 2 || col_sz_left != col_sz_right) return 0;

                    assert(ind_v_left >= -1);
                    assert(ind_v_right >= -1);

                    assert(state_left.c->vertex_to_col[ind_v_left] == col);
                    assert(state_right.c->vertex_to_col[ind_v_right] == col);
                    assert(state_left.c->vertex_to_col[ind_v_right] != col);

                    state_right.use_trace_early_out(false);
                    state_left.use_trace_early_out(false);
                    state_right.T->reset_trace_unequal();
                    state_left.T->reset_trace_unequal();

                    // individualize a vertex of base color
                    state_right.move_to_child(g, ind_v_right);
                    state_left.move_to_child(g,  ind_v_left);

                    // invariant differs
                    if(state_right.T->get_hash() != state_left.T->get_hash()) return 0;

                    // can not perform diff udpates on this (AKA invariant should have differed, but is not
                    // guaranteed to differ always)
                    const bool is_diffed = state_right.update_diff_vertices_last_individualization(state_left);
                    if(state_left.get_diff_diverge()) return 0;

                    // no difference? check for automorphism now -- if it doesn't succeed there is no hope on this path
                    if(!is_diffed) {
                        ws_automorphism.reset();
                        state_right.singleton_automorphism(state_left, ws_automorphism);
                        const bool found_auto = state_right.certify(g, ws_automorphism);
                        return found_auto;
                    }

                    assert(state_right.c->cells < g->v_size); // there can not be diff on a leaf
                }

                // unreachable
                assert(false);
                return false;
            }

            int do_paired_dfs(dejavu_hook* hook, sgraph *g, ir::controller &state_left, ir::controller& state_right,
                              std::vector<std::pair<int, int>>& computed_orbits, int max_depth = 0) {
                if(h_recent_cost_snapshot_limit < 0) return state_right.s_base_pos;
                s_grp_sz.mantissa = 1.0;
                s_grp_sz.exponent = 0;

                // orbit algorithm structure
                groups::orbit orbs(g->v_size);
                mark_set orbit_handled(g->v_size);

                // automorphism workspace
                ws_automorphism.reset();

                // saucy strategy of paired states
                state_left.link_compare(&state_right);

                // we want to terminate if things become to costly, we save the current trace position to track cost
                cost_snapshot = state_right.T->get_position();

                // tell the controller we are performing DFS now
                state_right.mode_compare_base();

                int failed_first_level = -1;

                // abort criteria
                double recent_cost_snapshot = 0;
                bool   fail = false;
                int    trace_pos_reset = 0;

                // loop that serves to optimize Tinhofer graphs
                while ((recent_cost_snapshot < h_recent_cost_snapshot_limit || state_right.s_base_pos <= 1) &&
                        state_right.s_base_pos > 0 && !fail) {
                //while (state_right.s_base_pos > 0) {
                    // backtrack one level
                    state_right.move_to_parent();
                    if((state_right.s_base_pos & 0x00000FFF) == 0x00000FFD)
                        progress_current_method("dfs", "base_pos", static_cast<double>(state_right.compare_base->size())
                                                -state_right.s_base_pos, "cost_snapshot", recent_cost_snapshot);
                    // remember which color we are individualizing
                    const int col         = state_right.base[state_right.s_base_pos].color;
                    const int col_sz      = state_right.c->ptn[col] + 1;
                    const int base_vertex = state_right.base_vertex[state_right.s_base_pos]; // base vertex
                              assert(col_sz >= 2);
                    assert(!state_right.there_is_difference_to_base_including_singles(g->v_size));

                    int prune_cost_snapshot = 0; /*< if we prune, keep track of how costly it is */
                    orbit_handled.reset();

                    // iterate over current color class
                    for (int i = 0; i < col_sz; ++i) {
                        int ind_v = state_right.leaf_color.lab[col + i];
                        assert(state_right.c->vertex_to_col[base_vertex] == state_right.c->vertex_to_col[ind_v]);

                        // only consider every orbit once
                        if(orbs.are_in_same_orbit(ind_v, base_vertex)) continue;
                        //if(!orbs.represents_orbit(ind_v)) continue;
                        if(!orbs.represents_orbit(ind_v)) ind_v = orbs.find_orbit(ind_v);
                        if(orbit_handled.get(ind_v)) continue;
                        orbit_handled.set(ind_v);


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
                        state_right.use_trace_early_out(true);
                        state_right.move_to_child(g, ind_v);

                        bool pruned = !state_right.T->trace_equal();
                        bool found_auto = false;

                        assert(state_right.c->vertex_to_col[ind_v] ==
                               state_right.leaf_color.vertex_to_col[base_vertex]);

                        if(!pruned) {
                            // write singleton diff into automorphism...
                            const int wr_pos_st  = state_right.base[state_right.base.size() - 1].singleton_pt;
                            const int wr_pos_end = (int) state_right.singletons.size();

                            ws_automorphism.reset();
                            // ... and then check whether this implies a (sparse) automorphism
                            ws_automorphism.write_singleton(state_right.compare_singletons,
                                                             &state_right.singletons, wr_pos_st,
                                                             wr_pos_end);
                            found_auto = state_right.certify(g, ws_automorphism);

                            assert(ws_automorphism.perm()[base_vertex] == ind_v);
                            // if no luck with sparse automorphism, try more proper walk to leaf node
                            if (!found_auto) {
                                ws_automorphism.reset();

                                state_left.T->reset_trace_equal();
                                state_left.T->set_position(trace_pos_reset);
                                state_left.use_trace_early_out(true);
                                state_left.move_to_child(g, base_vertex); // need to move left to base vertex

                                assert(state_left.T->trace_equal());
                                assert(state_left.T->get_position() == state_right.T->get_position());

                                auto return_code = paired_recurse_to_equal_leaf(g, state_left, state_right);
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
                            assert(ws_automorphism.nsupport() > 0);
                            assert(ws_automorphism.perm()[ind_v] == base_vertex ||
                                   ws_automorphism.perm()[base_vertex] == ind_v);
                            orbs.add_automorphism_to_orbit(ws_automorphism);
                            if(hook) (*hook)(g->v_size, ws_automorphism.perm(), ws_automorphism.nsupport(),
                                             ws_automorphism.support());
                        }
                        ws_automorphism.reset();

                        // move state back up where we started in this iteration
                        while (prev_base_pos < state_right.s_base_pos) {
                            state_right.move_to_parent();
                        }
                        state_right.T->set_position(trace_pos_reset);

                        // if no automorphism could be determined we would have to backtrack -- so stop!

                        // ((4*prune_cost_snapshot > cost_snapshot) &&
                        //                                                         (prev_base_pos > 0))
                        if ((!found_auto && !pruned) ) {
                            //if(((4*prune_cost_snapshot > cost_snapshot) &&
                            //    (prev_base_pos > 0))) std::cout << "cost fail" << std::endl;
                            fail = true;
                            if(failed_first_level == -1)
                                failed_first_level = state_right.s_base_pos;
                            break;
                        }
                        // if orbit size equals color size now, we are done on this DFS level
                        if (orbs.orbit_size(base_vertex) == col_sz) break;
                    }

                    // if we did not fail, accumulate size of current level to group size
                    if (!fail) {
                        computed_orbits.emplace_back(state_right.base_vertex[state_right.s_base_pos],
                                                     orbs.orbit_size(base_vertex));
                        s_grp_sz.multiply(orbs.orbit_size(base_vertex));
                    }
                }

                // set reason for termination
                s_termination = fail? r_fail : r_cost;

                // if DFS failed on current level, we did not finish the current level -- has to be accounted for
                return fail?failed_first_level + (fail):state_right.s_base_pos;
            }
        };
    }
}

#endif //DEJAVU_DFS_H
