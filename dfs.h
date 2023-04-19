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
            int fail_cnt = 0;
            int threads = 1;
            ir::refinement   *R = nullptr;
            int cost_snapshot = 0; /**< used to track cost-based abort criterion */
            ir::trace compare_T; // TODO should not be inside dfs_ir
            double h_recent_cost_snapshot_limit = 0.25;

        public:
            long double grp_sz_man = 1.0; /**< group size mantissa */
            int         grp_sz_exp = 0;   /**< group size exponent */

            /**
             * Setup the DFS module.
             *
             * @param threads number of threads we are allowed to dispatch
             * @param R refinement workspace
             */
            void setup(int threads, ir::refinement *R, double recent_cost_snapshot_limit = 0.25) {
                this->R = R;
                this->threads = threads;
                grp_sz_man = 1.0;
                grp_sz_exp = 0;
                h_recent_cost_snapshot_limit = recent_cost_snapshot_limit;
            }

            std::pair<bool, bool> recurse_to_equal_leaf(sgraph *g, coloring *initial_colors, ir::controller *state,
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

                    //assert(state->h_hint_color_is_singleton_now ? state->h_last_refinement_singleton_only : true);

                    if (prev_fail && state->h_last_refinement_singleton_only && state->h_hint_color_is_singleton_now) {
                        prev_cert = R->check_single_failure(g, initial_colors->vertex_to_col, automorphism.perm(), prev_fail_pos);
                    }

                    //if(state->c->cells == g->v_size) {
                    if (prev_cert && state->h_last_refinement_singleton_only && state->h_hint_color_is_singleton_now) {
                        // TODO: add better heuristic to not always do this check, too expensive!
                        auto cert_res = R->certify_automorphism_sparse_report_fail_resume(g, initial_colors->vertex_to_col,
                                                                                          automorphism.perm(),
                                                                                          automorphism.nsupport(),
                                                                                          automorphism.support(),
                                                                                          cert_pos);
                        cert_pos = std::get<2>(cert_res);
                        if (std::get<0>(cert_res)) {
                            cert_res = R->certify_automorphism_sparse_report_fail_resume(g, initial_colors->vertex_to_col,
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

            /**
             * Performs DFS from a given leaf node. Does not backtrack and returns the level up to which DFS succeeded.
             *
             * @param g The graph.
             * @param initial_colors The initial coloring of the graph \p g.
             * @param local_state The state from which DFS will be performed. Must be a leaf node of the IR shared_tree.
             * @return The level up to which DFS succeeded.
             */
            int do_dfs(sgraph *g, coloring *initial_colors, ir::controller &local_state, std::vector<std::pair<int, int>>* save_to_individualize) {
                // orbit algorithm structure
                groups::orbit orbs;
                orbs.initialize(g->v_size);

                // automorphism workspace
                groups::automorphism_workspace pautomorphism(g->v_size);

                // we want to terminate if things become to costly, we save the current trace position to track cost
                cost_snapshot = local_state.T->get_position();

                // tell the controller we are performing DFS now
                local_state.mode_dfs();

                // abort criteria
                double recent_cost_snapshot = 0;
                bool   fail = false;

                // loop that serves to optimize Tinhofer graphs
                while (recent_cost_snapshot < h_recent_cost_snapshot_limit && local_state.base_pos > 0 && !fail) {
                    // backtrack one level
                    local_state.move_to_parent();

                    // remember which color we are individualizing
                    const int col = local_state.base_color[local_state.base_pos];
                    const int col_sz = local_state.c->ptn[col] + 1;
                    const int vert = local_state.base_vertex[local_state.base_pos]; // base vertex

                    // iterate over current color class
                    for (int i = col_sz - 1; i >= 0; --i) {
                        const int ind_v = local_state.leaf_color.lab[col + i];

                        // only consider every orbit once
                        if (ind_v == vert || !orbs.represents_orbit(ind_v)) continue;
                        if (orbs.are_in_same_orbit(ind_v, vert))        continue;

                        // track cost of this refinement for whatever is to come
                        const int cost_start = local_state.T->get_position();

                        // perform individualization-refinement
                        const int prev_base_pos = local_state.base_pos;
                        local_state.T->reset_trace_equal();
                        local_state.move_to_child(R, g, ind_v);

                        // write singleton diff into automorphism...
                        const int wr_pos_st = local_state.base_singleton_pt[local_state.base_singleton_pt.size() - 1];
                        const int wr_pos_end= (int) local_state.singletons.size();

                        // ... and then check whether this implies a (sparse) automorphism
                        pautomorphism.write_singleton(&local_state.compare_singletons, &local_state.singletons,
                                                      wr_pos_st,wr_pos_end);
                        bool found_auto = R->certify_automorphism_sparse(g, initial_colors->vertex_to_col, pautomorphism.perm(),
                                                                    pautomorphism.nsupport(),
                                                                    pautomorphism.support());
                        assert(pautomorphism.perm()[vert] == ind_v);

                        // if no luck with sparse automorphism, try more proper walk to leaf node
                        if (!found_auto) {
                            auto rec_succeeded = recurse_to_equal_leaf(g, initial_colors, &local_state, pautomorphism);
                            found_auto = (rec_succeeded.first && rec_succeeded.second);
                            if (rec_succeeded.first && !rec_succeeded.second) {
                                pautomorphism.reset();
                                pautomorphism.write_color_diff(local_state.c->vertex_to_col, local_state.leaf_color.lab);
                                found_auto = R->certify_automorphism_sparse(g, initial_colors->vertex_to_col,
                                                                            pautomorphism.perm(),
                                                                            pautomorphism.nsupport(),
                                                                            pautomorphism.support());
                            }
                        }

                        // track cost-based abort criterion
                        const int cost_end   = local_state.T->get_position();
                        double cost_partial  = (cost_end - cost_start) / (cost_snapshot*1.0);
                        recent_cost_snapshot = (cost_partial + recent_cost_snapshot * 3) / 4;

                        // if we found automorphism, add to orbit, (and TODO: call hook)
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
                        if (orbs.orbit_size(vert) == col_sz) break;
                    }

                    // if we did not fail, accumulate size of current level to group size
                    if (!fail) {

                        if(initial_colors->vertex_to_col[initial_colors->lab[col]] == col && initial_colors->ptn[col] + 1 == col_sz) {
                            save_to_individualize->push_back({local_state.leaf_color.lab[col], col_sz});
                        } else {
                            //std::cout << (initial_colors->vertex_to_col[initial_colors->lab[col]]  == col) << ", " << (initial_colors->ptn[col] + 1 == col_sz) << std::endl;
                        }

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
}

#endif //DEJAVU_DFS_H
