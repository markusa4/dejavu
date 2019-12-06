#ifndef DEJAVU_GRAPH_H
#define DEJAVU_GRAPH_H

#include <vector>
#include "bijection.h"
#include "coloring.h"

class alignas(16) sgraph {
    struct vertexComparator {
        vertexComparator(const sgraph& g) : g(g) {}
        const sgraph& g;

        bool operator()( const int & v1, const int & v2) {
            return g.d[v1] < g.d[v2];
        }
    };
public:
    int* v;
    int* d;
    int* e;

    int v_size;
    int d_size;
    int e_size;

    int max_degree;

    void permute_graph(sgraph* ng, bijection* p);
    bool certify_automorphism(bijection p);
    void initialize_coloring(coloring* c);
    void copy_graph(sgraph *g);
};

#endif //DEJAVU_GRAPH_H
