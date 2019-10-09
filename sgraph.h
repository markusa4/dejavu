#ifndef BRUTUS_GRAPH_H
#define BRUTUS_GRAPH_H

#include <vector>
#include "bijection.h"
#include "coloring.h"
#include "coloring_bucket.h"

class sgraph {
    struct vertexComparator {
        vertexComparator(const sgraph& g) : g(g) {}
        const sgraph& g;

        bool operator()( const int & v1, const int & v2) {
            return g.d[v1] < g.d[v2];
        }
    };
public:
    std::vector<int> v;
    std::vector<int> d;
    std::vector<int> e;

    sgraph permute_graph(bijection p);
    bool certify_automorphism(bijection p);
    void initialize_coloring(coloring* c);
    void initialize_coloring_bucket(coloring_bucket* c);
};

#endif //BRUTUS_GRAPH_H
