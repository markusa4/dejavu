#ifndef BRUTUS_GRAPH_H
#define BRUTUS_GRAPH_H

#include <vector>
#include "bijection.h"
#include "coloring.h"

class graph {
    struct vertexComparator {
        vertexComparator(const graph& g) : g(g) {}
        const graph& g;

        bool operator()( const int & v1, const int & v2) {
            return g.d[v1] < g.d[v2];
        }
    };
public:
    std::vector<int> v;
    std::vector<int> d;
    std::vector<int> e;

    graph permute_graph(bijection p);
    void initialize_coloring(coloring* c);
};

#endif //BRUTUS_GRAPH_H
