#ifndef BRUTUS_GRAPH_H
#define BRUTUS_GRAPH_H

#include <vector>
#include "bijection.h"
#include "coloring.h"
#include "coloring_bucket.h"

class alignas(64) sgraph {
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

    sgraph permute_graph(bijection p);
    bool certify_automorphism(bijection p);
    void initialize_coloring(coloring* c);
    void initialize_coloring_bucket(coloring_bucket* c);

    void copy_graph(sgraph *g);
};

#endif //BRUTUS_GRAPH_H
