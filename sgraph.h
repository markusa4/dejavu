#ifndef DEJAVU_GRAPH_H
#define DEJAVU_GRAPH_H

#include <vector>
#include "bijection.h"
#include "coloring.h"
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <set>
#include <cstring>

template<class vertex_t, class degree_t, class edge_t>
class alignas(16) sgraph_t {
    struct vertexComparator {
        vertexComparator(const sgraph_t<vertex_t, degree_t, edge_t>& g) : g(g) {}
        const sgraph_t& g;

         bool operator()(const vertex_t & v1, const vertex_t & v2) {
            return g.d[v1] < g.d[v2];
        }
    };

    struct vertexComparatorColor {
        vertexComparatorColor(const sgraph_t<vertex_t, degree_t, edge_t>& g, const vertex_t* vertex_to_col) : g(g), vertex_to_col(vertex_to_col) {}
        const sgraph_t& g;
        const vertex_t* vertex_to_col;

        bool operator()(const vertex_t & v1, const vertex_t & v2) {
            //return (g.d[v1] < g.d[v2]) || ((g.d[v1] == g.d[v2]) && (vertex_to_col[v1] < vertex_to_col[v2]));
            return (vertex_to_col[v1] < vertex_to_col[v2]);
        }
    };
public:
    bool initialized = false;
    edge_t*   v;
    degree_t* d;
    vertex_t* e;

    int v_size;
    int d_size;
    int e_size;

    int max_degree;

    void initialize(int nv, int ne) {
        initialized = true;
        v = new edge_t[nv];
        d = new degree_t[nv];
        e = new vertex_t[ne];
    }

    // initialize a coloring of this sgraph, partitioning degrees of vertices
    void initialize_coloring(coloring<vertex_t> *c, vertex_t* vertex_to_col) {
        c->alloc(this->v_size);
        for(int i = 0; i < v_size; i++) {
            c->lab[i] = i;
        }

        std::memset(c->ptn, 1, sizeof(int) * v_size);

        int cells = 0;
        int last_new_cell   = 0;

        if(vertex_to_col == nullptr) {
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
                if(this->d[c->lab[i]] < this->d[c->lab[i + 1]]) {
                    c->ptn[i] = 0;
                    cells += 1;
                    c->ptn[last_new_cell] = i - last_new_cell;
                    last_new_cell = i + 1;
                    continue;
                }
            }
        } else {
            std::sort(c->lab, c->lab + c->lab_sz, vertexComparatorColor(*this, vertex_to_col));
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
                if(vertex_to_col[c->lab[i]] < vertex_to_col[c->lab[i + 1]]) {
                    c->ptn[i] = 0;
                    cells += 1;
                    c->ptn[last_new_cell] = i - last_new_cell;
                    last_new_cell = i + 1;
                    continue;
                }
            }
        }

        c->cells = cells;
        //std::cout << "Cells after initialize: " << cells << std::endl;
    }

    void initialize_coloring_raw(coloring<vertex_t> *c) {
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

    // certify that a permutation is an automorphism of the sgraph
    bool certify_automorphism(bijection<vertex_t> p) {
        assert(p.map_sz == v_size);

        std::set<int> image_neighbours_of_i;
        for(int i = 0; i < v_size; ++i) {
            int image_i = p.map_vertex(i);
            if(d[i] != d[image_i]) // degrees must be equal
                return false;

            image_neighbours_of_i.clear();
            // automorphism must preserve neighbours
            for(int j = v[i]; j < v[i] + d[i]; ++j) {
                int vertex_j = e[j];
                int image_j  = p.map_vertex(vertex_j);
                image_neighbours_of_i.insert(image_j);
            }
            for(int j = v[image_i]; j < v[image_i] + d[image_i]; ++j) {
                int vertex_j = e[j];
                if(image_neighbours_of_i.find(vertex_j) == image_neighbours_of_i.end()) {
                    return false;
                }
                image_neighbours_of_i.erase(image_neighbours_of_i.find(vertex_j));
            }
            if(!image_neighbours_of_i.empty()) {
                return false;
            }
        }

        return true;
    }

    void permute_graph(sgraph_t<vertex_t, degree_t, edge_t>* ng, bijection<vertex_t>* p) {
        ng->initialize(v_size, e_size);
        ng->v_size = v_size;
        ng->d_size = d_size;
        ng->e_size = e_size;
        ng->max_degree = max_degree;

        bijection<vertex_t> p_inv;
        p_inv.initialize_empty(p->map_sz);
        p_inv.copy(p);
        p_inv.inverse();

        std::set<int> vertices_hit;

        int epos = 0;
        for(int i = 0; i < v_size; ++i) {
            int mapped_v = p->map_vertex(i);
            assert(p_inv.map_vertex(mapped_v) == i);
            assert(mapped_v < v_size);
            vertices_hit.insert(mapped_v);
            ng->d[i] = d[mapped_v];
            ng->v[i] = epos;
            for(int j = v[mapped_v]; j < v[mapped_v] + d[mapped_v]; j++) {
                assert(j < e_size);
                ng->e[epos] = p_inv.map_vertex(e[j]);
                epos += 1;
            }
            //epos += ng->d[i];
        }

        //assert(v_size == vertices_hit.size());

        assert(ng->v_size == v_size);
        assert(ng->e_size == e_size);
        assert(ng->d_size == d_size);
        assert(epos == ng->e_size);

        return;
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
        mark_set multiedge_test;
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

    void copy_graph(sgraph_t<vertex_t, degree_t, edge_t>* g) {
        if(initialized) {
            delete[] v;
            delete[] d;
            delete[] e;
        }
        initialize(g->v_size, g->e_size);

        memcpy(v, g->v, g->v_size*sizeof(edge_t));
        memcpy(d, g->d, g->d_size*sizeof(degree_t));
        memcpy(e, g->e, g->e_size*sizeof(vertex_t));
        v_size = g->v_size;
        d_size = g->d_size;
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
        PRINT("[api] d_size: " << d_size);
        PRINT("[api] e_size: " << e_size);
        /*for(int i = 0; i < v_size; ++i) {
            const int estart = v[i];
            const int eend   = estart + d[i];
            for(int j = estart; j < eend; ++j) {
                PRINT("[api] (" << i << ", " << e[j] << ")")
            }
        }*/
    }

    ~sgraph_t() {
        if(initialized) {
            delete[] v;
            delete[] d;
            delete[] e;
        }
    }
};


template<class vertex_type_src, class degree_type_src, class edge_type_src,
         class vertex_type_tgt, class degree_type_tgt, class edge_type_tgt>
static void copy_graph(sgraph_t<vertex_type_src, degree_type_src, edge_type_src>* g1,
                       sgraph_t<vertex_type_tgt, degree_type_tgt, edge_type_tgt>* g2) {
    g2->v_size = g1->v_size;
    g2->d_size = g1->d_size;
    g2->e_size = g1->e_size;
    g2->max_degree = g1->max_degree;

    g2->initialize(g2->v_size, g2->e_size);

    for(int i = 0; i < g1->v_size; ++i) {
        g2->v[i] = static_cast<edge_type_tgt>(g1->v[i]);
    }
    for(int i = 0; i < g1->d_size; ++i) {
        g2->d[i] = static_cast<degree_type_tgt>(g1->d[i]);
    }
    for(int i = 0; i < g1->e_size; ++i) {
        g2->e[i] = static_cast<vertex_type_tgt>(g1->e[i]);
    }
}

enum sgraph_type {DSG_INT_INT_INT, DSG_SHORT_SHORT_INT, DSG_SHORT_SHORT_SHORT, DSG_CHAR_CHAR_SHORT, DSG_CHAR_CHAR_CHAR};
struct dynamic_sgraph {
    sgraph_type type;
    sgraph_t<int, int, int>*             sgraph_0 = nullptr;
    sgraph_t<int16_t, int16_t, int>*     sgraph_1 = nullptr;
    sgraph_t<int16_t, int16_t, int16_t>* sgraph_2 = nullptr;
    sgraph_t<int8_t,  int8_t,  int16_t>* sgraph_3 = nullptr;
    sgraph_t<int8_t,  int8_t,  int8_t>*  sgraph_4 = nullptr;

    static void read(sgraph_t<int, int, int>* g, dynamic_sgraph* sg) {
        bool short_v = (g->v_size <= 32767);
        bool short_e = (g->e_size <= 32767);
        bool char_v  = (g->v_size <= 127);
        bool char_e  = (g->e_size <= 127);

        if(char_v && char_e) {
            sg->sgraph_4 = new sgraph_t<int8_t,  int8_t,  int8_t>;
            copy_graph<int, int, int, int8_t,  int8_t,  int8_t>(g, sg->sgraph_4);
            sg->type = DSG_CHAR_CHAR_CHAR;
            return;
        }

        if(char_v && short_e) {
            sg->sgraph_3 = new sgraph_t<int8_t,  int8_t,  int16_t>;
            copy_graph<int, int, int, int8_t,  int8_t,  int16_t>(g, sg->sgraph_3);
            sg->type = DSG_CHAR_CHAR_SHORT;
            return;
        }

        if(short_v && short_e) {
            sg->sgraph_2 = new sgraph_t<int16_t,  int16_t,  int16_t>;
            copy_graph<int, int, int, int16_t,  int16_t,  int16_t>(g, sg->sgraph_2);
            sg->type = DSG_SHORT_SHORT_SHORT;
            return;
        }

        if(short_v && !short_e) {
            sg->sgraph_1 = new sgraph_t<int16_t,  int16_t,  int>;
            copy_graph<int, int, int, int16_t,  int16_t,  int>(g, sg->sgraph_1);
            sg->type = DSG_SHORT_SHORT_INT;
            return;
        }

        sg->sgraph_0 = g;
        sg->type = DSG_INT_INT_INT;
        return;
    }
};


template<class vertex_t, class degree_t, class edge_t>
sgraph_t<vertex_t, degree_t, edge_t>* disjoint_union(sgraph_t<vertex_t, degree_t, edge_t>* g1,
                                                     sgraph_t<vertex_t, degree_t, edge_t>* g2) {
    sgraph_t<vertex_t, degree_t, edge_t>* union_g = new  sgraph_t<vertex_t, degree_t, edge_t>();
    union_g->v_size = g1->v_size + g2->v_size;
    union_g->e_size = g1->e_size + g2->e_size;
    union_g->d_size = g1->d_size + g2->d_size;

    int g2_vshift= g1->v_size;
    int g2_eshift= g1->e_size;
    int g2_dshift= g1->d_size;

    union_g->initialize(union_g->v_size, union_g->e_size);

    for(int i = 0; i < g1->v_size; ++i)
        union_g->v[i] = g1->v[i];
    for(int i = 0; i < g2->v_size; ++i)
        union_g->v[i + g2_vshift] = g2->v[i] + g2_eshift;

    for(int i = 0; i < g1->d_size; ++i)
        union_g->d[i] = g1->d[i];
    for(int i = 0; i < g2->d_size; ++i)
        union_g->d[i + g2_dshift] = g2->d[i];

    for(int i = 0; i < g1->e_size; ++i)
        union_g->e[i] = g1->e[i];
    for(int i = 0; i < g2->e_size; ++i)
        union_g->e[i + g2_eshift] = g2->e[i] + g2_vshift;

    union_g->max_degree = std::max(g1->max_degree, g2->max_degree);

    return union_g;
}


typedef sgraph_t<int, int, int> sgraph;

template<class vertex_t>
void permute_colmap(int** colmap, int colmap_sz, vertex_t* p) {
    int* new_colmap = new vertex_t[colmap_sz];
    for(int i = 0; i < colmap_sz; ++i) {
        new_colmap[i] = (*colmap)[p[i]];
    }
    int* old_colmap = *colmap;
    *colmap = new_colmap;
    delete[] old_colmap;
}

static sgraph test_graph;

#endif //DEJAVU_GRAPH_H
