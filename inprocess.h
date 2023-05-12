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


        static unsigned long invariant_path2(sgraph *g, coloring* c, int v) {
            unsigned long hash_inv = 0;
            const int ptn_st = g->v[v];
            const int ptn_en = ptn_st + g->d[v];
            for(int i = ptn_st; i < ptn_en; ++i) {
                const int neighb_v = g->e[i];
                const int ptn_st_i = g->v[neighb_v];
                const int ptn_en_i = ptn_st_i + g->d[neighb_v];
                for(int j = ptn_st_i; j < ptn_en_i; ++j) {
                    const int neighb2_v = g->e[j];
                    const unsigned int add_hash = hash(c->vertex_to_col[neighb2_v]);
                    hash_inv += add_hash;
                }
            }
            return hash_inv;
        }

        std::vector<std::pair<int, int>> inproc_can_individualize; /**< vertices that can be individualized */
        std::vector<int>                 inproc_fixed_points;      /**< vertices fixed by inprocessing      */

        bool inprocess(sgraph *g, ir::shared_tree *tree, groups::shared_schreier *group, ir::controller &local_state,
                       ir::limited_save &root_save, int budget) {
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

            // A stronger invariant from the BFS tree
            /*if (tree->get_finished_up_to() >= 2 || cell_prev != local_state.c->cells) {
                const int cell_prev_inv = local_state.get_coloring()->cells;
                work_list hash(g->v_size);
                std::vector<unsigned long> new_invariant;
                new_invariant.reserve(g->v_size);
                for(int i = 0; i < g->v_size; ++i) {
                    new_invariant.push_back(local_state.get_coloring()->vertex_to_col[i]);
                }
                //const int level = tree->get_finished_up_to();
                //for(int i = 0; i < tree->get_level_size(level); ++i) {
                //    auto node = tree->pick_node_from_level(level, i);
                //    auto node_save = node->get_save();
                //    const int first_level_node = node_save->get_base()[0];
                //    auto inv = invariant_path2(g, node_save->get_coloring(), first_level_node);
                //    new_invariant[first_level_node] += inv;
                //}
                //for(int i = 0; i < g->v_size; ++i) {
                //    const int v = i;
                //    auto inv = invariant_path2(g, local_state.get_coloring(), v);
                //    new_invariant[v] += inv;
                //}

                const int level = tree->get_finished_up_to();
                const int level_sz = tree->get_level_size(level);
                const int degree = g->d[0];
                //auto test_node = tree->pick_node_from_level(level, 0);
                //auto test_node_save = test_node->get_save();
                auto test_node_save = &root_save;
                int smallest_non_trivial_color    = -1;
                int smallest_non_trivial_color_sz = INT32_MAX;
                for (int i = 0; i < test_node_save->get_coloring()->ptn_sz;) {
                    const int col_sz = test_node_save->get_coloring()->ptn[i];
                    if (col_sz > 0 && col_sz < smallest_non_trivial_color_sz) {
                        smallest_non_trivial_color    = i;
                        smallest_non_trivial_color_sz = col_sz;
                    }
                    i += col_sz + 1;
                }

                if(smallest_non_trivial_color >= 0 && level_sz * degree < budget && smallest_non_trivial_color_sz < g->v_size) {
                    std::cout << smallest_non_trivial_color_sz << std::endl;
                    for (int i = 0; i < level_sz; ++i) {
                        auto node = tree->pick_node_from_level(level, i);
                        auto node_save = node->get_save();
                        //for (int j = 0; j < smallest_non_trivial_color_sz + 1; ++j) {
                            //const int v = root_save.get_coloring()->lab[smallest_non_trivial_color + j];
                        for (int j = 0; j < g->v_size; ++j) {
                            const int v = j;
                            auto inv = invariant_path2(g, node_save->get_coloring(), v);
                            new_invariant[v] += inv;
                        }
                    }
                }

                for (int i = 0; i < g->v_size; ++i) {
                    hash[i] = (int) new_invariant[i] % 256;
                }
                for (int i = 0; i < g->v_size; ++i) {
                    hash[i] = 256 * local_state.c->vertex_to_col[i] + hash[i];
                }
                g->initialize_coloring(local_state.c, hash.get_array());
                const int cell_after = local_state.c->cells;
                progress_print("inpr_inv", std::to_string(cell_prev_inv),
                               std::to_string(cell_after));
                if (cell_after != cell_prev_inv) {
                    local_state.refine(g);
                    progress_print("inpr_ref", std::to_string(cell_after),
                                   std::to_string(local_state.c->cells));
                }
            }*/

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
