#include <iostream>
#include "parser.h"
#include "refinement.h"
#include "selector.h"
#include "ir_tools.h"
#include <assert.h>

void label_graph(graph* g, bijection* canon_p) {
    coloring c;
    g->initialize_coloring(&c);
    std::cout << "------------------" << std::endl;

    ir_tools IR;
    IR.label_graph(g, canon_p);

}

int main() {
    std::cout << "dejavu 0.1" << std::endl;
    std::cout << "------------------" << std::endl;

    // parse a graph
    parser p;
    graph g;
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/rantree/rantree/rantree-10.bliss", &g);
    p.parse_dimacs_file("/home/markus/Downloads/graphs/k/k/k-4", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/ag/ag/ag2-5", &g);
    // canonically label the graph
    bijection canon_p;
    label_graph(&g, &canon_p);
    return 0;
}