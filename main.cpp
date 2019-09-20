#include <iostream>
#include "parser.h"
#include "refinement.h"
#include "selector.h"
#include "ir_tools.h"
#include <assert.h>

void label_graph(sgraph* g, bijection* canon_p) {
    coloring c;
    g->initialize_coloring(&c);
    std::cout << "------------------" << std::endl;

    ir_tools IR;
    IR.label_graph(g, canon_p);
}

int main() {
    std::cout << "dejavu 0.1" << std::endl;
    std::cout << "------------------" << std::endl;

    // parse a sgraph
    parser p;
    sgraph g;
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/rantree/rantree/rantree-20.bliss", &g);
    p.parse_dimacs_file("/home/markus/Downloads/graphs/k/k/k-5", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/ag/ag/ag2-2", &g);
    // canonically label the sgraph


    //bijection test;
    //test.map = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    //std::cout << g.certify_automorphism(test) << std::endl;

    bijection canon_p;
    label_graph(&g, &canon_p);
    return 0;
}