#include <iostream>
#include "parser.h"
#include "refinement.h"
#include "selector.h"
#include <assert.h>

void label_graph(graph* g, bijection* canon_p) {
    coloring c;
    g->initialize_coloring(&c);
    std::cout << "------------------" << std::endl;

    refinement R;
    selector S;

    std::set<std::pair<int, int>> changes;

    R.refine_coloring(g, &c, &changes);
    R.undo_changes(g, &c, &changes);

    int s = S.select_color(g, &c);

    if(s == - 1) {
        std::cout << "Discrete coloring found." << std::endl;
    }

    int v = c.lab[s];
    std::cout << "Individualizing " << v << ", color class " << s << std::endl;
    R.individualize_vertex(g, &c, v);

    R.refine_coloring(g, &c, &changes);
}

int main() {
    std::cout << "dejavu 0.1" << std::endl;
    std::cout << "------------------" << std::endl;

    // parse a graph
    parser p;
    graph g;
    p.parse_dimacs_file("/home/markus/Downloads/graphs/rantree/rantree/rantree-10.bliss", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/k/k/k-100", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/ag/ag/ag2-5", &g);
    // canonically label the graph
    bijection canon_p;
    label_graph(&g, &canon_p);
    return 0;
}