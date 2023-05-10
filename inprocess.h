// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_INPROCESS_H
#define DEJAVU_INPROCESS_H

#include "dfs.h"
#include "bfs.h"
#include "rand.h"

namespace dejavu::search_strategy {
    class inprocessor {
    public:
        // statistics
        big_number s_grp_sz; /**< group size */

        // TODO: option to compute add stronger invariants on the furthest BFS level

        std::vector<std::pair<int, int>> inproc_can_individualize; /**< vertices that can be individualized */
        std::vector<int>                 inproc_fixed_points;      /**< vertices fixed by inprocessing      */

        bool inprocess(sgraph *g, ir::shared_tree *tree, groups::shared_schreier *group, ir::controller &local_state,
                       ir::limited_save &root_save) {
            local_state.load_reduced_state(root_save);

            const int cell_prev = root_save.get_coloring()->cells; /*< keep track how many cells we have initially*/
            bool touched_coloring = false; /*< whether we change the root_save or not, i.e., whether we change
                                            *  anything */

            if (tree->get_finished_up_to() >= 1) { // did we do BFS?
                // TODO: hashing actually makes this worse!
                work_list hash(g->v_size);
                mark_set is_pruned(g->v_size);
                tree->mark_first_level(is_pruned);

                tree->make_node_invariant(); // "compresses" node invariant from all levels into first level


                for (int i = 0; i < g->v_size; ++i) {
                    hash[i]  = is_pruned.get(i);
                    hash[i] += (int) (*tree->get_node_invariant())[i] % 256;
                }
                for (int i = 0; i < g->v_size; ++i) {
                    hash[i] = 257 * local_state.c->vertex_to_col[i] + hash[i];
                }
                g->initialize_coloring(local_state.c, hash.get_array());
                const int cell_after = local_state.c->cells;
                progress_print("inpr_bfs", std::to_string(cell_prev),
                               std::to_string(cell_after));
                if (cell_after != cell_prev) {
                    local_state.refine(g);
                    progress_print("inpr_ref", std::to_string(cell_after),
                                   std::to_string(local_state.c->cells));
                } else {
                    //assert(!something_was_pruned);
                }
            }
            group->determine_potential_individualization(&inproc_can_individualize,
                                                         local_state.get_coloring());

            if (!inproc_can_individualize.empty()) {
                for (auto &i: inproc_can_individualize) {
                    const int ind_v = i.first;
                    const int ind_col = local_state.c->vertex_to_col[ind_v];
                    assert(i.second == local_state.c->ptn[ind_col] + 1);
                    s_grp_sz.multiply(local_state.c->ptn[ind_col] + 1);
                    local_state.move_to_child_no_trace(g, ind_v);
                    inproc_fixed_points.push_back(ind_v);
                }
                progress_print("inpr_ind", "_",
                               std::to_string(local_state.c->cells));
                inproc_can_individualize.clear();
            }

            inproc_can_individualize.clear();
            if(cell_prev != local_state.c->cells) {
                local_state.save_reduced_state(root_save);
                touched_coloring = true;
            }
            return touched_coloring;
        }
    };
}

#endif //DEJAVU_INPROCESS_H
