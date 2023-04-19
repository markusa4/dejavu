// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_GRAPH_H
#define DEJAVU_GRAPH_H

#include <vector>
#include "coloring.h"
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <set>
#include <cstring>
#include "sgraph.h"

namespace dejavu {
    class graph {
        struct vertexComparator {
            vertexComparator(const graph &g) : g(g) {}
            const graph &g;
            bool operator()(const int &v1, const int &v2) {
                return g.d(v1) < g.d(v2);
            }
        };

        struct vertexComparatorColor {
            vertexComparatorColor(const graph &g, const int *vertex_to_col) : g(g), vertex_to_col(vertex_to_col) {}
            const graph &g;
            const int *vertex_to_col;

            bool operator()(const int &v1, const int &v2) {
                //return (g.d[v1] < g.d[v2]) || ((g.d[v1] == g.d[v2]) && (vertex_to_col[v1] < vertex_to_col[v2]));
                return (vertex_to_col[v1] < vertex_to_col[v2]);
            }
        };

    public:
        bool initialized = false;
        int *data;
        int *vd;
        int *e;

        int v_size;
        int e_size;

        bool dense = false;

        int max_degree;

        void initialize(int nv, int ne) {
            initialized = true;
            data = new int[2*nv + ne];
            vd   = data;
            e    = data + 2*nv;
        }

        [[nodiscard]] inline int v(const int v) const {
            return vd[v+v];
        }

        [[nodiscard]] inline int d(const int v) const {
            return vd[v+v+1];
        }

        // initialize a coloring of this sgraph, partitioning degrees of vertices
        void initialize_coloring(coloring *c, int *vertex_to_col) {
            c->alloc(this->v_size);
            std::memset(c->ptn, 1, sizeof(int) * v_size);

            if (this->v_size < c->lab_sz) {
                c->lab_sz = this->v_size;
                c->ptn_sz = this->v_size;
            }

            if (v_size == 0)
                return;

            int cells = 0;
            int last_new_cell = 0;

            if (vertex_to_col == nullptr) {
                for (int i = 0; i < v_size; i++) {
                    c->lab[i] = i;
                }
                std::sort(c->lab, c->lab + c->lab_sz, vertexComparator(*this));
                for (int i = 0; i < c->lab_sz; i++) {
                    c->vertex_to_col[c->lab[i]] = last_new_cell;
                    c->vertex_to_lab[c->lab[i]] = i;
                    if (i + 1 == c->lab_sz) {
                        cells += 1;
                        c->ptn[last_new_cell] = i - last_new_cell;
                        c->ptn[i] = 0;
                        break;
                    }
                    assert(this->d(c->lab[i]) <= this->d(c->lab[i + 1]));
                    if (this->d(c->lab[i]) != this->d(c->lab[i + 1])) {
                        c->ptn[i] = 0;
                        cells += 1;
                        c->ptn[last_new_cell] = i - last_new_cell;
                        last_new_cell = i + 1;
                        continue;
                    }
                }
            } else {
                int col = 0;

                int min_col = INT32_MAX;
                int max_col = INT32_MIN;
                for (int i = 0; i < v_size; i++) {
                    if (vertex_to_col[i] < min_col)
                        min_col = vertex_to_col[i];
                    if (vertex_to_col[i] > max_col)
                        max_col = vertex_to_col[i];
                }

                std::vector<int> colsize;
                colsize.reserve(std::min(this->v_size, (max_col - min_col) + 1));

                if (min_col < 0 || max_col > 4 * this->v_size) {
                    std::unordered_map<int, int> colors; // TODO: should not use unordered_map!
                    colors.reserve(this->v_size);
                    for (int i = 0; i < v_size; i++) {
                        auto it = colors.find(vertex_to_col[i]);
                        if (it == colors.end()) {
                            colors.insert(std::pair<int, int>(vertex_to_col[i], col));
                            colsize.push_back(1);
                            assert(col < this->v_size);
                            vertex_to_col[i] = col;
                            ++col;
                        } else {
                            const int found_col = it->second;
                            assert(found_col < this->v_size);
                            vertex_to_col[i] = found_col;
                            ++colsize[found_col];
                        }
                    }
                } else {
                    std::vector<int> colors;
                    colors.reserve(max_col + 1);
                    for (int i = 0; i < max_col + 1; ++i)
                        colors.push_back(-1);
                    for (int i = 0; i < v_size; i++) {
                        if (colors[vertex_to_col[i]] == -1) {
                            colors[vertex_to_col[i]] = col;
                            colsize.push_back(1);
                            assert(col < this->v_size);
                            vertex_to_col[i] = col;
                            ++col;
                        } else {
                            const int found_col = colors[vertex_to_col[i]];
                            assert(found_col < this->v_size);
                            vertex_to_col[i] = found_col;
                            ++colsize[found_col];
                        }
                    }
                }

                int increment = 0;
                for (int i = 0; i < colsize.size(); i++) {
                    const int col_sz = colsize[i];
                    colsize[i] += increment;
                    increment += col_sz;
                }
                assert(increment == v_size);

                for (int i = 0; i < v_size; i++) {
                    const int v_col = vertex_to_col[i];
                    --colsize[v_col];
                    const int v_lab_pos = colsize[v_col];
                    c->lab[v_lab_pos] = i;
                }

                /*for(int i = 0; i < v_size; i++) {
                    c->lab[i] = i;
                }
                std::sort(c->lab, c->lab + c->lab_sz, vertexComparatorColor(*this, vertex_to_col));*/
                for (int i = 0; i < c->lab_sz; i++) {
                    c->vertex_to_col[c->lab[i]] = last_new_cell;
                    c->vertex_to_lab[c->lab[i]] = i;
                    if (i + 1 == c->lab_sz) {
                        cells += 1;
                        c->ptn[last_new_cell] = i - last_new_cell;
                        c->ptn[i] = 0;
                        break;
                    }
                    //assert(this->d[c->lab[i]] <= this->d[c->lab[i + 1]]);
                    //if(this->d[c->lab[i]] < this->d[c->lab[i + 1]]  || (this->d[c->lab[i]] == this->d[c->lab[i + 1]]
                    //&& (vertex_to_col[c->lab[i]] < vertex_to_col[c->lab[i + 1]]))) {
                    if (vertex_to_col[c->lab[i]] != vertex_to_col[c->lab[i + 1]]) {
                        c->ptn[i] = 0;
                        cells += 1;
                        c->ptn[last_new_cell] = i - last_new_cell;
                        last_new_cell = i + 1;
                        continue;
                    }
                }
            }

            c->cells = cells;
        }

        void initialize_coloring_raw(coloring *c) {
            c->alloc(this->v_size);

            for (int i = 0; i < v_size; i++) {
                c->lab[i] = i;
                c->vertex_to_lab[i] = i;
            }

            std::memset(c->vertex_to_col, 0, sizeof(int) * v_size);
            c->ptn[0] = v_size - 1;
            c->ptn[v_size - 1] = 0;
            c->cells = 1;
        }

        void sanity_check() {
#ifndef NDEBUG
            for (int i = 0; i < v_size; ++i) {
                assert(d(i) > 0 ? v(i) < e_size : true);
                assert(d(i) > 0 ? v(i) >= 0 : true);
                assert(d(i) >= 0);
                assert(d(i) < v_size);
            }
            for (int i = 0; i < e_size; ++i) {
                assert(e[i] < v_size);
                assert(e[i] >= 0);
            }

            // multiedge test
            mark_set multiedge_test;
            multiedge_test.initialize(v_size);
            for (int i = 0; i < v_size; ++i) {
                multiedge_test.reset();
                for (int j = 0; j < d(i); ++j) {
                    const int neigh = e[v(i) + j];
                    assert(!multiedge_test.get(neigh));
                    multiedge_test.set(neigh);
                }
            }

            // fwd - bwd test
            multiedge_test.initialize(v_size);
            for (int i = 0; i < v_size; ++i) {
                multiedge_test.reset();
                for (int j = 0; j < d(i); ++j) {
                    const int neigh = e[v(i) + j];
                    bool found = false;
                    for (int k = 0; k < d(neigh); ++k) {
                        const int neigh_neigh = e[v(neigh) + k];
                        if (neigh_neigh == i) {
                            found = true;
                            break;
                        }
                    }
                    assert(found);
                }
            }
#endif
        }

        void copy_graph(graph *g) {
            if (initialized) {
                delete[] data;
            }
            initialize(g->v_size, g->e_size);

            memcpy(vd, g->vd, 2*g->v_size * sizeof(int));
            memcpy(e, g->e, g->e_size * sizeof(int));
            v_size = g->v_size;
            e_size = g->e_size;
            max_degree = g->max_degree;
        }

        void copy_graph(sgraph *g) {
            if (initialized) {
                delete[] data;
            }
            initialize(g->v_size, g->e_size);

            for (int i = 0; i < g->v_size; ++i) {
                vd[2*i]   = g->v[i];
                vd[2*i+1] = g->d[i];
            }
            memcpy(e, g->e, g->e_size * sizeof(int));
            v_size = g->v_size;
            e_size = g->e_size;
            max_degree = g->max_degree;
        }

        void sort_edgelist() {
            for (int i = 0; i < v_size; ++i) {
                const int estart = v(i);
                const int eend = estart + d(i);
                std::sort(e + estart, e + eend);
            }
        }

        void print() {
            PRINT("nv: " << v_size);
            PRINT("ne: " << e_size);
            /*for(int i = 0; i < v_size; ++i) {
                const int estart = v[i];
                const int eend   = estart + d[i];
                for(int j = estart; j < eend; ++j) {
                    PRINT("[api] (" << i << ", " << e[j] << ")")
                }
            }*/
        }

        ~graph() {
            if (initialized) {
                delete[] data;
            }
        }
    };


    static void copy_graph(graph *g1,
                           graph *g2) {
        g2->v_size = g1->v_size;
        g2->e_size = g1->e_size;
        g2->max_degree = g1->max_degree;

        g2->initialize(g2->v_size, g2->e_size);


        memcpy(g2->vd, g1->vd, 2*g1->v_size * sizeof(int));
        for (int i = 0; i < g1->e_size; ++i) {
            g2->e[i] = g1->e[i];
        }
    }

    void permute_colmap(int **colmap, int colmap_sz, int *p) {
        int *new_colmap = new int[colmap_sz];
        for (int i = 0; i < colmap_sz; ++i) {
            new_colmap[i] = (*colmap)[p[i]];
        }
        int *old_colmap = *colmap;
        *colmap = new_colmap;
        delete[] old_colmap;
    }
}
#endif //DEJAVU_GRAPH_H
