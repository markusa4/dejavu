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

template<class vertex_type, class degree_type, class edge_type>
class alignas(16) sgraph_temp {
    struct vertexComparator {
        vertexComparator(const sgraph_temp<vertex_type, degree_type, edge_type>& g) : g(g) {}
        const sgraph_temp& g;

        bool operator()( const vertex_type & v1, const vertex_type & v2) {
            return g.d[v1] < g.d[v2];
        }
    };
public:
    edge_type*   v;
    degree_type* d;
    vertex_type* e;

    int v_size;
    int d_size;
    int e_size;

    int max_degree;

    // initialize a coloring of this sgraph, partitioning degrees of vertices
    void initialize_coloring(coloring_temp<vertex_type> *c) {
        c->lab = new vertex_type[this->v_size];
        c->ptn = new vertex_type[this->v_size];
        c->vertex_to_col = new vertex_type[this->v_size];
        c->vertex_to_lab = new vertex_type[this->v_size];
        c->lab_sz = this->v_size;
        c->ptn_sz = this->v_size;
        c->init = true;

        for(int i = 0; i < v_size; i++) {
            c->vertex_to_col[i] = -1;
            c->vertex_to_lab[i] = -1;
            c->lab[i] = i;
            c->ptn[i] = 1;
        }

        std::sort(c->lab, c->lab + c->lab_sz, vertexComparator(*this));

        int cells = 0;
        int last_new_cell   = 0;
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

        c->cells = cells;
    }

    // certify that a permutation is an automorphism of the sgraph
    bool certify_automorphism(bijection_temp<vertex_type> p) {
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

    void permute_graph(sgraph_temp<vertex_type, degree_type, edge_type>* ng, bijection_temp<vertex_type>* p) {
        ng->v = new edge_type[v_size];
        ng->d = new degree_type[d_size];
        ng->e = new vertex_type[e_size];
        ng->v_size = v_size;
        ng->d_size = d_size;
        ng->e_size = e_size;
        ng->max_degree = max_degree;

        bijection_temp<vertex_type> p_inv;
        p_inv.map = new vertex_type[p->map_sz];
        p_inv.map_sz = p->map_sz;
        p_inv.deletable();
        memcpy(p_inv.map, p->map, p->map_sz*sizeof(vertex_type));
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

    void copy_graph(sgraph_temp<vertex_type, degree_type, edge_type>* g) {
        v = new edge_type[g->v_size];
        d = new degree_type[g->d_size];
        e = new vertex_type[g->e_size];

        memcpy(v, g->v, g->v_size*sizeof(edge_type));
        memcpy(d, g->d, g->d_size*sizeof(degree_type));
        memcpy(e, g->e, g->e_size*sizeof(vertex_type));
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
};


template<class vertex_type_src, class degree_type_src, class edge_type_src,
         class vertex_type_tgt, class degree_type_tgt, class edge_type_tgt>
static void copy_graph(sgraph_temp<vertex_type_src, degree_type_src, edge_type_src>* g1,
                       sgraph_temp<vertex_type_tgt, degree_type_tgt, edge_type_tgt>* g2) {
    g2->v_size = g1->v_size;
    g2->d_size = g1->d_size;
    g2->e_size = g1->e_size;
    g2->max_degree = g1->max_degree;

    g2->v = new edge_type_tgt[g2->v_size];
    g2->d = new degree_type_tgt[g2->d_size];
    g2->e = new vertex_type_tgt[g2->e_size];

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
    sgraph_temp<int, int, int>*             sgraph_0 = nullptr;
    sgraph_temp<int16_t, int16_t, int>*     sgraph_1 = nullptr;
    sgraph_temp<int16_t, int16_t, int16_t>* sgraph_2 = nullptr;
    sgraph_temp<int8_t,  int8_t,  int16_t>* sgraph_3 = nullptr;
    sgraph_temp<int8_t,  int8_t,  int8_t>*  sgraph_4 = nullptr;

    static void read(sgraph_temp<int, int, int>* g, dynamic_sgraph* sg) {
        bool short_v = (g->v_size <= 32767);
        bool short_e = (g->e_size <= 32767);
        bool char_v  = (g->v_size <= 127);
        bool char_e  = (g->e_size <= 127);

        if(char_v && char_e) {
            sg->sgraph_4 = new sgraph_temp<int8_t,  int8_t,  int8_t>;
            copy_graph<int, int, int, int8_t,  int8_t,  int8_t>(g, sg->sgraph_4);
            sg->type = DSG_CHAR_CHAR_CHAR;
            return;
        }

        if(char_v && short_e) {
            sg->sgraph_3 = new sgraph_temp<int8_t,  int8_t,  int16_t>;
            copy_graph<int, int, int, int8_t,  int8_t,  int16_t>(g, sg->sgraph_3);
            sg->type = DSG_CHAR_CHAR_SHORT;
            return;
        }

        if(short_v && short_e) {
            sg->sgraph_2 = new sgraph_temp<int16_t,  int16_t,  int16_t>;
            copy_graph<int, int, int, int16_t,  int16_t,  int16_t>(g, sg->sgraph_2);
            sg->type = DSG_SHORT_SHORT_SHORT;
            return;
        }

        if(short_v && !short_e) {
            sg->sgraph_1 = new sgraph_temp<int16_t,  int16_t,  int>;
            copy_graph<int, int, int, int16_t,  int16_t,  int>(g, sg->sgraph_1);
            sg->type = DSG_SHORT_SHORT_INT;
            return;
        }

        sg->sgraph_0 = g;
        sg->type = DSG_INT_INT_INT;
        return;
    }

};

typedef sgraph_temp<int, int, int> sgraph;

#endif //DEJAVU_GRAPH_H
