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
            groups::automorphism_workspace* ws_automorphism = nullptr;

        public:
            enum termination_reason {r_none, r_fail, r_cost};
            termination_reason s_termination = r_none; /**< Why did we stop performing DFS? Cost too high, or did we
                                                         *  have to backtrack? */

            double h_recent_cost_snapshot_limit = 0.25; /**< A float in the range [0-1]. Limits continuation of DFS
                                                          * search to whethercomputing recent elements only cost this
                                                          * fraction of the cost of an entire root-to-leaf walk. */
            big_number s_grp_sz; /**< group size */


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
                    const int ind_v = state->c->lab[col];

                    // individualize a vertex of base color
                    state->move_to_child(g, ind_v);
                    const int pos_start = state->base_singleton_pt[state->base_singleton_pt.size() - 1];
                    const int pos_end   = (int) state->singletons.size();

                    // write sparse automorphism with singletons derived so far
                    ws_automorphism->write_singleton(&state->compare_singletons, &state->singletons,
                                                 pos_start, pos_end);

                    // try to certify the automorphism, but if it failed previously, first check whether that specific
                    // failure was even rectified
                    bool prev_cert = true;
                    if (prev_fail && state->s_last_refinement_singleton_only && state->s_hint_color_is_singleton_now) {
                        prev_cert = state->check_single_failure(g, initial_colors->vertex_to_col,
                                                                *ws_automorphism,prev_fail_pos);
                    }

                    if (prev_cert && state->s_last_refinement_singleton_only && state->s_hint_color_is_singleton_now) {
                        auto [certified, fail_pos, new_cert_pos] =
                                state->certify_sparse_report_fail_resume(g,
                                                                         initial_colors->vertex_to_col,
                                                                         *ws_automorphism, cert_pos);
                        cert_pos = new_cert_pos;
                        if (certified) {
                            std::tie(certified, fail_pos, new_cert_pos) =
                                    state->certify_sparse_report_fail_resume(g,initial_colors->vertex_to_col,
                                                                             *ws_automorphism, 0);
                        }

                        if (certified) {
                            return {true, true}; // we found an automorphism
                        } else {
                            prev_fail = true;
                            prev_fail_pos = fail_pos;
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

            /**
             * Performs DFS from a given leaf node. Does not backtrack and returns the level up to which DFS succeeded.
             *
             * @param g The graph.
             * @param initial_colors The initial coloring of the graph \p g.
             * @param local_state The state from which DFS will be performed. Must be a leaf node of the IR shared_tree.
             * @return The level up to which DFS succeeded.
             */
            int do_dfs(dejavu_hook* hook, sgraph *g, coloring *initial_colors, ir::controller &local_state,
                       std::vector<std::pair<int, int>>* save_to_individualize) {
                s_grp_sz.mantissa = 1.0;
                s_grp_sz.exponent = 0;

                // orbit algorithm structure
                groups::orbit orbs(g->v_size);

                // automorphism workspace
                ws_automorphism->reset();

                // we want to terminate if things become to costly, we save the current trace position to track cost
                cost_snapshot = local_state.T->get_position();

                // tell the controller we are performing DFS now
                local_state.mode_compare_base();

                // abort criteria
                double recent_cost_snapshot = 0;
                bool   fail = false;

                // loop that serves to optimize Tinhofer graphs
                while (recent_cost_snapshot < h_recent_cost_snapshot_limit && local_state.s_base_pos > 0 && !fail) {
                    // backtrack one level
                    local_state.move_to_parent();

                    // remember which color we are individualizing
                    const int col = local_state.base_color[local_state.s_base_pos];
                    const int col_sz = local_state.c->ptn[col] + 1;
                    const int vert = local_state.base_vertex[local_state.s_base_pos]; // base vertex

                    // iterate over current color class
                    for (int i = col_sz - 1; i >= 0; --i) {
                        const int ind_v = local_state.leaf_color.lab[col + i];
                        assert(local_state.c->vertex_to_col[vert] == local_state.c->vertex_to_col[ind_v]);

                        // only consider every orbit once
                        if (ind_v == vert || !orbs.represents_orbit(ind_v))  continue;
                        if (orbs.are_in_same_orbit(ind_v, vert))        continue;

                        // track cost of this refinement for whatever is to come
                        const int cost_start = local_state.T->get_position();

                        // perform individualization-refinement
                        const int prev_base_pos = local_state.s_base_pos;
                        local_state.T->reset_trace_equal();
                        local_state.move_to_child(g, ind_v);
                        assert(local_state.c->vertex_to_col[ind_v] == local_state.leaf_color.vertex_to_col[vert]);

                        // write singleton diff into automorphism...
                        const int wr_pos_st = local_state.base_singleton_pt[local_state.base_singleton_pt.size() - 1];
                        const int wr_pos_end= (int) local_state.singletons.size();

                        // ... and then check whether this implies a (sparse) automorphism
                        ws_automorphism->write_singleton(&local_state.compare_singletons,
                                                         &local_state.singletons, wr_pos_st,
                                                         wr_pos_end);
                        bool found_auto = local_state.certify(g, *ws_automorphism);
                        assert(ws_automorphism->perm()[vert] == ind_v);
                        // if no luck with sparse automorphism, try more proper walk to leaf node
                        if (!found_auto) {
                            auto [success, certified] = recurse_to_equal_leaf(g, initial_colors,
                                                                              &local_state);
                            found_auto = (success && certified);
                            if (success && !certified) {
                                ws_automorphism->reset();
                                assert(ws_automorphism->perm()[vert] == vert);
                                assert(ws_automorphism->perm()[ind_v] == ind_v);
                                assert(local_state.c->vertex_to_col[ind_v] == local_state.leaf_color.vertex_to_col[vert]);
                                ws_automorphism->write_color_diff(local_state.c->vertex_to_col,
                                                                  local_state.leaf_color.lab);
                                assert(ws_automorphism->perm()[ind_v] == vert || ws_automorphism->perm()[vert] == ind_v);
                                found_auto = local_state.certify(g, *ws_automorphism);
                            }
                        }

                        // track cost-based abort criterion
                        const int cost_end   = local_state.T->get_position();
                        double cost_partial  = (cost_end - cost_start) / (cost_snapshot*1.0);
                        recent_cost_snapshot = (cost_partial + recent_cost_snapshot * 3) / 4;

                        // if we found automorphism, add to orbit and call hook
                        if (found_auto) {
                            assert(ws_automorphism->nsupport() > 0);
                            assert(ws_automorphism->perm()[ind_v] == vert || ws_automorphism->perm()[vert] == ind_v);
                            if(hook) (*hook)(g->v_size, ws_automorphism->perm(), ws_automorphism->nsupport(),
                                             ws_automorphism->support());
                            orbs.add_automorphism_to_orbit(*ws_automorphism);
                        }
                        ws_automorphism->reset();

                        // move state back up where we started in this iteration
                        while (prev_base_pos < local_state.s_base_pos) {
                            local_state.move_to_parent();
                        }

                        // if no automorphism could be determined we would have to backtrack -- so stop!
                        if (!found_auto) {
                            fail = true;
                            break;
                        }

                        // if orbit size equals color size now, we are done on this DFS level
                        if (orbs.orbit_size(vert) == col_sz) break;
                    }

                    // if we did not fail, accumulate size of current level to group size
                    if (!fail) {
                        if(initial_colors->vertex_to_col[initial_colors->lab[col]] == col &&
                           initial_colors->ptn[col] + 1 == col_sz) {
                            save_to_individualize->emplace_back(local_state.base_vertex[local_state.s_base_pos], col_sz);
                        }
                        s_grp_sz.multiply(col_sz);
                    }
                }

                // set reason for termination
                s_termination = fail? r_fail : r_cost;

                // if DFS failed on current level, we did not finish the current level -- has to be accounted for
                return local_state.s_base_pos + (fail);
            }
        };
    }
}

#endif //DEJAVU_DFS_H
