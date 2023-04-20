// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_SGRAPH_H
#define DEJAVU_SGRAPH_H

#include <vector>
#include "coloring.h"
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <set>
#include <cstring>
#include "ds.h"

class sgraph {
    struct vertexComparator {
        vertexComparator(const sgraph& g) : g(g) {}
        const sgraph& g;

         bool operator()(const int & v1, const int & v2) {
            return g.d[v1] < g.d[v2];
        }
    };

    struct vertexComparatorColor {
        vertexComparatorColor(const sgraph& g, const int* vertex_to_col) : g(g), vertex_to_col(vertex_to_col) {}
        const sgraph& g;
        const int* vertex_to_col;

        bool operator()(const int & v1, const int & v2) {
            //return (g.d[v1] < g.d[v2]) || ((g.d[v1] == g.d[v2]) && (vertex_to_col[v1] < vertex_to_col[v2]));
            return (vertex_to_col[v1] < vertex_to_col[v2]);
        }
    };
public:
    bool initialized = false;
    int* v;
    int* d;
    int* e;

    int v_size;
    int e_size;

    bool dense = false;

    int max_degree;

    void initialize(int nv, int ne) {
        initialized = true;
        v = new int[nv];
        d = new int[nv];
        e = new int[ne];
    }

    // initialize a coloring of this sgraph, partitioning degrees of vertices
    void initialize_coloring(coloring *c, int* vertex_to_col) {
        c->alloc(this->v_size);
        std::memset(c->ptn, 1, sizeof(int) * v_size);

        if(this->v_size < c->lab_sz) {
            c->lab_sz = this->v_size;
            c->ptn_sz = this->v_size;
        }

        if(v_size == 0)
            return;

        int cells = 0;
        int last_new_cell   = 0;

        if(vertex_to_col == nullptr) {
            for(int i = 0; i < v_size; i++) {
                c->lab[i] = i;
            }
            std::sort(c->lab, c->lab + c->lab_sz, vertexComparator(*this));
            for(int i = 0; i < c->lab_sz; i++) {
                c->vertex_to_col[c->lab[i]] = last_new_cell;
                c->vertex_to_lab[c->lab[i]] = i;
                if(i + 1 == c->lab_sz) {
                    cells += 1;
                    c->ptn[last_new_cell] = i - last_new_cell;
                    c->ptn[i] = 0;
                    break;
                }
                assert(this->d[c->lab[i]] <= this->d[c->lab[i + 1]]);
                if(this->d[c->lab[i]] != this->d[c->lab[i + 1]]) {
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
            for(int i = 0; i < v_size; i++) {
                if(vertex_to_col[i] < min_col)
                    min_col = vertex_to_col[i];
                if(vertex_to_col[i] > max_col)
                    max_col = vertex_to_col[i];
            }

            std::vector<int> colsize;
            colsize.reserve(std::min(this->v_size, (max_col - min_col)  + 1));

            if(min_col < 0 || max_col > 4*this->v_size) {
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
                for(int i = 0; i < max_col + 1; ++i)
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
            for(int i = 0; i < colsize.size(); i++) {
                const int col_sz = colsize[i];
                colsize[i] += increment;
                increment += col_sz;
            }
            assert(increment == v_size);

            for(int i = 0; i < v_size; i++) {
                const int v_col     = vertex_to_col[i];
                --colsize[v_col];
                const int v_lab_pos = colsize[v_col];
                c->lab[v_lab_pos] = i;
            }

            /*for(int i = 0; i < v_size; i++) {
                c->lab[i] = i;
            }
            std::sort(c->lab, c->lab + c->lab_sz, vertexComparatorColor(*this, vertex_to_col));*/
            for(int i = 0; i < c->lab_sz; i++) {
                c->vertex_to_col[c->lab[i]] = last_new_cell;
                c->vertex_to_lab[c->lab[i]] = i;
                if(i + 1 == c->lab_sz) {
                    cells += 1;
                    c->ptn[last_new_cell] = i - last_new_cell;
                    c->ptn[i] = 0;
                    break;
                }
                //assert(this->d[c->lab[i]] <= this->d[c->lab[i + 1]]);
                //if(this->d[c->lab[i]] < this->d[c->lab[i + 1]]  || (this->d[c->lab[i]] == this->d[c->lab[i + 1]]
                //&& (vertex_to_col[c->lab[i]] < vertex_to_col[c->lab[i + 1]]))) {
                if(vertex_to_col[c->lab[i]] != vertex_to_col[c->lab[i + 1]]) {
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

        for(int i = 0; i < v_size; i++) {
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
        for(int i = 0; i < v_size; ++i) {
            assert(d[i]>0?v[i] < e_size:true);
            assert(d[i]>0?v[i] >= 0:true);
            assert(d[i] >= 0);
            assert(d[i] < v_size);
        }
        for(int i = 0; i < e_size; ++i) {
            assert(e[i] < v_size);
            assert(e[i] >= 0);
        }

        // multiedge test
        dejavu::ds::mark_set multiedge_test;
        multiedge_test.initialize(v_size);
        for(int i = 0; i < v_size; ++i) {
            multiedge_test.reset();
            for(int j = 0; j < d[i]; ++j) {
                const int neigh = e[v[i] + j];
                assert(!multiedge_test.get(neigh));
                multiedge_test.set(neigh);
            }
        }

        // fwd - bwd test
        multiedge_test.initialize(v_size);
        for(int i = 0; i < v_size; ++i) {
            multiedge_test.reset();
            for(int j = 0; j < d[i]; ++j) {
                const int neigh = e[v[i] + j];
                bool found = false;
                for(int k = 0; k < d[neigh]; ++k) {
                    const int neigh_neigh = e[v[neigh] + k];
                    if(neigh_neigh == i) {
                        found = true;
                        break;
                    }
                }
                assert(found);
            }
        }
#endif
    }

    void copy_graph(sgraph* g) {
        if(initialized) {
            delete[] v;
            delete[] d;
            delete[] e;
        }
        initialize(g->v_size, g->e_size);

        memcpy(v, g->v, g->v_size*sizeof(int));
        memcpy(d, g->d, g->v_size*sizeof(int));
        memcpy(e, g->e, g->e_size*sizeof(int));
        v_size = g->v_size;
        e_size = g->e_size;
        max_degree = g->max_degree;
    }

    void sort_edgelist() {
        for(int i = 0; i < v_size; ++i) {
            const int estart = v[i];
            const int eend   = estart + d[i];
            std::sort(e + estart, e + eend);
        }
    }

    void print() {
        PRINT("[api] v_size: " << v_size);
        PRINT("[api] e_size: " << e_size);
        /*for(int i = 0; i < v_size; ++i) {
            const int estart = v[i];
            const int eend   = estart + d[i];
            for(int j = estart; j < eend; ++j) {
                PRINT("[api] (" << i << ", " << e[j] << ")")
            }
        }*/
    }

    ~sgraph() {
        if(initialized) {
            delete[] v;
            delete[] d;
            delete[] e;
        }
    }
};

void permute_colmap(int** colmap, int colmap_sz, int* p) {
    int* new_colmap = new int[colmap_sz];
    for(int i = 0; i < colmap_sz; ++i) {
        new_colmap[i] = (*colmap)[p[i]];
    }
    int* old_colmap = *colmap;
    *colmap = new_colmap;
    delete[] old_colmap;
}

static sgraph test_graph;

#endif //DEJAVU_SGRAPH_H
