// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_INPROCESS_H
#define DEJAVU_INPROCESS_H

#include "dfs.h"
#include "bfs.h"
#include "rand.h"

namespace dejavu {
    namespace search_strategy {
        class inprocessor {
        private:
            void multiply_to_group_size(long double add_grp_sz_man, int add_grp_sz_exp) {
                while (add_grp_sz_man >= 10.0) {
                    s_grp_sz_exp += 1;
                    add_grp_sz_man = add_grp_sz_man / 10;
                }
                s_grp_sz_exp += add_grp_sz_exp;
                s_grp_sz_man *= add_grp_sz_man;
                while (s_grp_sz_man >= 10.0) {
                    s_grp_sz_exp += 1;
                    s_grp_sz_man = s_grp_sz_man / 10;
                }
            }

        public:
            // statistics
            long double s_grp_sz_man = 1.0; /**< group size mantissa, group size is `s_grp_sz_man^s_grp_sz_exp`
                                          * \sa s_grp_sz_exp */
            int s_grp_sz_exp = 0;   /**< group size exponent, group size is `s_grp_sz_man^s_grp_sz_exp`
                                          * \sa s_grp_sz_man  */

            std::vector<std::pair<int, int>> inproc_can_individualize; /**< vertices that can be individualized */
            std::vector<int>                 inproc_fixed_points;      /**< vertices fixed by inprocessing      */

            bool
            inprocess(sgraph *g, ir::shared_tree *tree, groups::shared_schreier *group, ir::controller &local_state,
                      ir::reduced_save &root_save) {
                local_state.load_reduced_state(root_save);
                const int cell_prev = root_save.get_coloring()->cells;

                if (tree->get_finished_up_to() >= 1) {
                    // TODO: improve saving hash for pruned nodes, then propagate this to first level
                    // TODO: could fix "usr" maybe?
                    // TODO: use "pruned" flag, etc.
                    // TODO: hashing actually makes this worse!
                    work_list hash(g->v_size);
                    mark_set is_pruned(g->v_size);
                    tree->mark_first_level(is_pruned);

                    for (int i = 0; i < g->v_size; ++i) {
                        if (is_pruned.get(i)) hash[i] = 1;
                        else hash[i] = 0;
                        hash[i] += (int) (*tree->get_node_invariant())[i] % 16;
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

                group->determine_save_to_individualize(&inproc_can_individualize, local_state.get_coloring());

                if (!inproc_can_individualize.empty()) {
                    for (auto &i: inproc_can_individualize) {
                        const int ind_v = i.first;
                        const int test_col_sz = i.second;
                        const int ind_col = local_state.c->vertex_to_col[ind_v];
                        assert(test_col_sz == local_state.c->ptn[ind_col] + 1);
                        multiply_to_group_size(local_state.c->ptn[ind_col] + 1, 0);
                        local_state.move_to_child_no_trace(g, ind_v);
                        inproc_fixed_points.push_back(ind_v);
                    }
                    progress_print("inproc_ind", "_",
                                   std::to_string(local_state.c->cells));
                    inproc_can_individualize.clear();
                }

                inproc_can_individualize.clear();
                local_state.save_reduced_state(root_save);

                return cell_prev != local_state.c->cells;
            }
        };
    }
}

#endif //DEJAVU_INPROCESS_H
