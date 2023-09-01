// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_INPROCESS_H
#define DEJAVU_INPROCESS_H

#include "dfs.h"
#include "bfs.h"
#include "rand.h"
#include "components.h"

namespace dejavu::search_strategy {
    class inprocessor {
    public:
        // statistics
        big_number s_grp_sz; /**< group size */

        std::vector<std::pair<int, int>> inproc_can_individualize; /**< vertices that can be individualized           */
        std::vector<std::pair<int, int>> inproc_maybe_individualize; /**< vertices that can be individualized         */
        std::vector<int>                 inproc_fixed_points;      /**< vertices fixed by inprocessing                */

        std::vector<int> nodes;

        void sort_nodes_map(std::vector<unsigned long>* map, int* colmap) {
            struct comparator_map {
                std::vector<unsigned long> *map;
                int* colmap;

                explicit comparator_map(std::vector<unsigned long> *map, int* colmap) {
                    this->map = map;
                    this->colmap = colmap;
                }

                bool operator()(const int &a, const int &b) const {
                    return (colmap[a] < colmap[b]) || ((colmap[a] == colmap[b]) && ((*map)[a] < (*map)[b]));
                }
            };
            auto c = comparator_map(map, colmap);
            std::sort(nodes.begin(), nodes.end(), c);
        }

        void sort_nodes_map(unsigned long* map, int* colmap) {
            struct comparator_map {
                unsigned long *map;
                int* colmap;

                explicit comparator_map(unsigned long *map, int* colmap) {
                    this->map = map;
                    this->colmap = colmap;
                }

                bool operator()(const int &a, const int &b) const {
                    return (colmap[a] < colmap[b]) || ((colmap[a] == colmap[b]) && (map[a] < map[b]));
                }
            };
            auto c = comparator_map(map, colmap);
            std::sort(nodes.begin(), nodes.end(), c);
        }


        void shallow_split_invariant(sgraph* g, ir::controller &local_state, worklist_t<unsigned long>& inv) {
            local_state.use_reversible(true);
            local_state.use_trace_early_out(true);
            local_state.use_increase_deviation(true);
            local_state.use_split_limit(true, 20); // 20

            for(int i = 0; i < g->v_size; ++i) {
                local_state.T->set_hash(0);
                local_state.reset_trace_equal();
                const int col = local_state.c->vertex_to_col[i];
                const int col_sz = local_state.c->ptn[col] + 1;
                if(col_sz >= 2) {
                    local_state.move_to_child(g, i);
                    inv[i] += local_state.T->get_hash();
                    local_state.move_to_parent();
                } else {
                    inv[i] = 0;
                }
            }

            local_state.use_trace_early_out(false);
            local_state.use_increase_deviation(false);
            local_state.use_split_limit(false);
        }

        void shallow_split_invariant2(sgraph* g, ir::controller &local_state, worklist_t<unsigned long>& inv) {
            local_state.use_reversible(true);
            local_state.use_trace_early_out(false);
            local_state.use_increase_deviation(true);
            local_state.use_split_limit(true, 8);

            mark_set original_colors(g->v_size);
            for(int _col = 0; _col < g->v_size;) {
                const int _col_sz = local_state.c->ptn[_col] + 1;
                original_colors.set(_col);
                _col += _col_sz;
            }

            for(int i = 0; i < g->v_size; ++i) {
                local_state.T->set_hash(0);
                local_state.reset_trace_equal();
                const int col = local_state.c->vertex_to_col[i];
                const int col_sz = local_state.c->ptn[col] + 1;
                if(col_sz >= 2) {
                    local_state.move_to_child(g, i);
                    inv[i] += local_state.T->get_hash();
                    for(int _col = 0; _col < g->v_size;) {
                        const int _col_sz = local_state.c->ptn[_col] + 1;
                        if (_col_sz >= 2 && _col_sz <= 16 && !original_colors.get(_col)) {
                            for (int jj = _col; jj < _col + _col_sz; ++jj) {
                                const int j = local_state.c->lab[jj];
                                local_state.move_to_child(g, j);
                                inv[i] += local_state.T->get_hash();
                                local_state.move_to_parent();
                            }
                        }
                        _col += _col_sz;
                    }
                    local_state.move_to_parent();
                } else {
                    inv[i] = 0;
                }
            }

            local_state.use_trace_early_out(false);
            local_state.use_increase_deviation(false);
            local_state.use_split_limit(false);
        }

        unsigned long shallow_invariant_comb(sgraph* g, ir::controller &local_state) {
            local_state.use_reversible(true);
            local_state.use_trace_early_out(true);
            local_state.use_increase_deviation(true);
            local_state.use_split_limit(true, 64); // 20

            unsigned long acc = 0;
            int pick_col = -1;

            for(int i = 0; i < g->v_size;) {
                const int col = i;
                const int col_sz = local_state.c->ptn[col] + 1;
                if (col_sz >= 2) {
                    pick_col = col;
                    break;
                }
                i += col_sz;
            }
            if(pick_col != -1) {
                const int pick_col_sz = local_state.c->ptn[pick_col] + 1;
                std::vector<int> verts;
                for (int i = 0; i < pick_col_sz; ++i) verts.push_back(local_state.c->lab[pick_col + i]);
                for (int i = 0; i < verts.size(); ++i) {
                    local_state.T->set_hash(0);
                    local_state.reset_trace_equal();
                    local_state.move_to_child(g, verts[i]);
                    acc += local_state.T->get_hash();
                    local_state.move_to_parent();
                }
            }

            local_state.use_trace_early_out(false);
            local_state.use_increase_deviation(false);
            local_state.use_split_limit(false);

            return acc;
        }

        static int edgeflip(sgraph *g, coloring& c, ds::mark_set& del_e) {
            worklist worklist_deg0(g->v_size); // serves as the "test vector"

            int flipped_edge = 0;
            int potential = 0;

            for (int y = 0; y < g->v_size; ++y) worklist_deg0[y] = 0;

            for (int i = 0; i < g->v_size;) {
                int v = c.lab[i];
                for (int f = g->v[v]; f < g->v[v] + g->d[v]; ++f) {
                    const int v_neigh = g->e[f];
                    worklist_deg0[c.vertex_to_col[v_neigh]] += 1;
                }

                for (int ii = 0; ii < c.ptn[i] + 1; ++ii) {
                    const int vx = c.lab[i + ii];
                    bool skipped_none = true;
                    for (int f = g->v[vx]; f < g->v[vx] + g->d[vx]; ++f) {
                        const int v_neigh = g->e[f];
                        if (worklist_deg0[c.vertex_to_col[v_neigh]] ==
                            c.ptn[c.vertex_to_col[v_neigh]] + 1) {
                            del_e.set(f); // mark edge for deletion (reverse edge is handled later automatically)
                            skipped_none = false;
                            ++flipped_edge;
                        }
                    }
                    if(skipped_none) break;
                }

                for (int f = g->v[v]; f < g->v[v] + g->d[v]; ++f) {
                    const int v_neigh = g->e[f];
                    worklist_deg0[c.vertex_to_col[v_neigh]] = 0;
                }

                i += c.ptn[i] + 1;
            }

            return flipped_edge;
        }

        void perform_del_edge(dejavu::sgraph *g, mark_set& del_e) {
            if (g->v_size <= 1)
                return;

            // int pre_esize = g->e_size;
            // copy some stuff
            std::vector<int> g_old_v;
            g_old_v.clear();
            g_old_v.reserve(g->v_size);

            for (int i = 0; i < g->v_size; ++i) {
                g_old_v.push_back(g->v[i]);
            }

            // make graph smaller using the translation array
            int epos = 0;
            for (int i = 0; i < g->v_size; ++i) {
                const int old_v = i;
                const int new_v = old_v;

                if (new_v >= 0) {
                    int new_d = 0;
                    g->v[new_v] = epos;
                    for (int j = g_old_v[old_v]; j < g_old_v[old_v] + g->d[old_v]; ++j) {
                        const int ve = g->e[j];
                        const int new_ve = ve;
                        if (!del_e.get(j)) {
                            assert(new_ve < new_vsize);
                            assert(new_ve >= 0);
                            ++new_d;
                            g->e[epos] = new_ve;
                            ++epos;
                        }
                    }

                    g->d[new_v] = new_d;
                }
            }

            g->e_size = epos;
        }

        /**
         * Inprocess the (colored) graph using all the available solver data.
         *
         * @param g graph
         * @param tree currently available  ir tree
         * @param group currently available group of symmetries
         * @param local_state workspace to perform individualization&refinement in
         * @param root_save the current coloring of the IR tree root
         * @param budget
         * @return whether any preprocessing was performed
         */
        bool inprocess(sgraph *g, ir::shared_tree &tree, groups::compressed_schreier &group, ir::controller &local_state,
                       ir::limited_save &root_save, [[maybe_unused]] int budget, bool use_bfs_inprocess,
                       bool use_shallow_inprocess, bool use_shallow_quadratic_inprocess) {
            local_state.load_reduced_state(root_save);

            const int cell_prev = root_save.get_coloring()->cells; /*< keep track how many cells we have initially*/
            bool touched_coloring = false; /*< whether we change the root_save or not, i.e., whether we change
                                            *  anything */

            // computes a shallow breadth-first invariant
            if (use_shallow_inprocess && !(tree.get_finished_up_to() >= 1 && use_bfs_inprocess)) {
                worklist hash(g->v_size);
                worklist_t<unsigned long> inv(g->v_size);

                nodes.reserve(g->v_size);
                nodes.clear();

                for(int i = 0; i < g->v_size; ++i) inv[i] = 0;
                shallow_split_invariant(g, local_state, inv);

                int num_of_hashs = 0;
                for(int i = 0; i < g->v_size; ++i) nodes.push_back(i);
                sort_nodes_map(inv.get_array(), local_state.c->vertex_to_col);
                unsigned long last_inv = -1;
                int last_col  = local_state.c->vertex_to_col[0];
                for(int i = 0; i < g->v_size; ++i) {
                    const int v      = nodes[i];
                    const unsigned long v_inv = inv[v];
                    const int  v_col = local_state.c->vertex_to_col[v];
                    if(last_col != v_col) {
                        last_col = v_col;
                        last_inv = v_inv;
                        ++num_of_hashs;
                    }
                    if(v_inv != last_inv) {
                        last_inv = v_inv;
                        ++num_of_hashs;
                    }
                    hash[v] = num_of_hashs;
                }

                g->initialize_coloring(local_state.c, hash.get_array());
                const int cell_after = local_state.c->cells;
                if (cell_after != cell_prev) {
                    local_state.refine(g);
                }
            }

            if (use_shallow_quadratic_inprocess) {
                worklist hash(g->v_size);
                worklist_t<unsigned long> inv(g->v_size);

                nodes.reserve(g->v_size);
                nodes.clear();

                for(int i = 0; i < g->v_size; ++i) inv[i] = 0;
                shallow_split_invariant2(g, local_state, inv);

                int num_of_hashs = 0;
                for(int i = 0; i < g->v_size; ++i) nodes.push_back(i);
                sort_nodes_map(inv.get_array(), local_state.c->vertex_to_col);
                unsigned long last_inv = -1;
                int last_col  = local_state.c->vertex_to_col[0];
                for(int i = 0; i < g->v_size; ++i) {
                    const int v      = nodes[i];
                    const unsigned long v_inv = inv[v];
                    const int  v_col = local_state.c->vertex_to_col[v];
                    if(last_col != v_col) {
                        last_col = v_col;
                        last_inv = v_inv;
                        ++num_of_hashs;
                    }
                    if(v_inv != last_inv) {
                        last_inv = v_inv;
                        ++num_of_hashs;
                    }
                    hash[v] = num_of_hashs;
                }

                g->initialize_coloring(local_state.c, hash.get_array());
                const int cell_after = local_state.c->cells;
                if (cell_after != cell_prev) {
                    local_state.refine(g);
                }
            }

            if (tree.get_finished_up_to() >= 1 && use_bfs_inprocess) {
                worklist hash(g->v_size);
                mark_set is_pruned(g->v_size);
                tree.mark_first_level(is_pruned);

                tree.make_node_invariant(); // "compresses" node invariant from all levels into first level

                nodes.reserve(g->v_size);
                nodes.clear();

                int num_of_hashs = 0;

                for(int i = 0; i < g->v_size; ++i) nodes.push_back(i);
                sort_nodes_map(tree.get_node_invariant(), local_state.c->vertex_to_col);
                unsigned long last_inv = (*tree.get_node_invariant())[0];
                int last_col           = local_state.c->vertex_to_col[0];
                for(int i = 0; i < g->v_size; ++i) {
                    const int v               = nodes[i];
                    const unsigned long v_inv = (*tree.get_node_invariant())[v];
                    const int  v_col = local_state.c->vertex_to_col[v];
                    if(last_col != v_col) {
                        last_col = v_col;
                        last_inv = v_inv;
                        ++num_of_hashs;
                    }
                    if(v_inv != last_inv) {
                        last_inv = v_inv;
                        ++num_of_hashs;
                    }
                    hash[v] = num_of_hashs;
                }

                g->initialize_coloring(local_state.c, hash.get_array());
                const int cell_after = local_state.c->cells;
                if (cell_after != cell_prev) {
                    local_state.refine(g);
                }
            }

            group.determine_potential_individualization(&inproc_can_individualize, local_state.get_coloring());
            if (!inproc_can_individualize.empty() || !inproc_maybe_individualize.empty()) {
                int num_inds = 0;
                for (auto &i: inproc_can_individualize) {
                    const int ind_v = i.first;
                    const int ind_col = local_state.c->vertex_to_col[ind_v];
                    assert(i.second == local_state.c->ptn[ind_col] + 1);
                    s_grp_sz.multiply(local_state.c->ptn[ind_col] + 1);
                    local_state.move_to_child_no_trace(g, ind_v);
                    inproc_fixed_points.push_back(ind_v);
                    ++num_inds;
                }
                for (auto &i: inproc_maybe_individualize) {
                    const int ind_v      = i.first;
                    const int ind_col    = local_state.c->vertex_to_col[ind_v];
                    const int ind_col_sz = local_state.c->ptn[ind_col] + 1;
                    const int orb_sz_det = i.second;
                    if(ind_col_sz > 1 && ind_col_sz == orb_sz_det) {
                        s_grp_sz.multiply(ind_col_sz);
                        local_state.move_to_child_no_trace(g, ind_v);
                        inproc_fixed_points.push_back(ind_v);
                        ++num_inds;
                    }
                }
                if(num_inds > 0) {
                    /*progress_print("inpr_ind", "_",
                                   std::to_string(local_state.c->cells));*/
                    inproc_can_individualize.clear();
                }
            }

            inproc_can_individualize.clear();
            inproc_maybe_individualize.clear();

            if(cell_prev != local_state.c->cells) {
                local_state.save_reduced_state(root_save);
                touched_coloring = true;
            }

            return touched_coloring;
        }
    };
}

#endif //DEJAVU_INPROCESS_H
