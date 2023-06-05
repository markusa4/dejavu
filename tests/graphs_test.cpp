#include "gtest/gtest.h"
#include "../dejavu.h"
#include "../parser.h"
#include <filesystem>

void test_graph(std::string filename) {
    dejavu::sgraph *g = new dejavu::sgraph();
    std::cout << "Parsing " << filename << "..." << std::endl;
    int* colmap = nullptr;
    parse_dimacs_file_fast(filename, g, &colmap);

    auto test_hook_ = dejavu_hook(dejavu::test_hook);

    dejavu::dejavu2 d;
    d.automorphisms(g, colmap);
}

TEST(graphs_test, graph_suite) {
    std::string directory = TEST_RESOURCE_DIR;
    test_graph(directory + "X2.bliss");
    test_graph(directory + "rantree-100.bliss");
    test_graph(directory + "ransq_1000_a.bliss");
    test_graph(directory + "Ranreg1024.bliss");
    test_graph(directory + "pp-16-10.dimacs");
    test_graph(directory + "NCKL900K.bliss");
    test_graph(directory + "kef14.dimacs");
    test_graph(directory + "f-lex-reg-20-1.dimacs");
    test_graph(directory + "latin-20.dimacs");
    test_graph(directory + "latin-30.dimacs");
    test_graph(directory + "latin-sw-23-4.dimacs");
    test_graph(directory + "cfi-172.dimacs");
    test_graph(directory + "highschool1-aigio.dimacs");
    test_graph(directory + "AS.bliss");
    test_graph(directory + "usr4_116-1.dimacs");
}