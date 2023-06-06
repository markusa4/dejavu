// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

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

TEST(graphs_test, graph_suite_comb) {
    std::string directory = TEST_RESOURCE_DIR;
    test_graph(directory + "X2.bliss");
    test_graph(directory + "Ranreg1024.bliss");
    test_graph(directory + "pp-16-10.dimacs");
    test_graph(directory + "pp-25-100");
    test_graph(directory + "NCKL900K.bliss");
    test_graph(directory + "kef14.dimacs");
    test_graph(directory + "f-lex-reg-20-1.dimacs");
    test_graph(directory + "latin-20.dimacs");
    test_graph(directory + "latin-30.dimacs");
    test_graph(directory + "latin-sw-23-4.dimacs");
    test_graph(directory + "cfi-172.dimacs");
    test_graph(directory + "usr4_116-1.dimacs");
    test_graph(directory + "CHH_cc7-7_1078-1.dimacs");
    test_graph(directory + "sts-79.dimacs");
    test_graph(directory + "triang-20.dimacs");
}

TEST(graphs_test, graph_suite_prep) {
    std::string directory = TEST_RESOURCE_DIR;
    test_graph(directory + "rantree-100.bliss");
    test_graph(directory + "highschool1-aigio.dimacs");
    test_graph(directory + "AS.bliss");
}

TEST(graphs_test, graph_suite_ref) {
    std::string directory = TEST_RESOURCE_DIR;
    test_graph(directory + "ransq_1000_a.bliss");
    test_graph(directory + "ran10_500_a.bliss");
    test_graph(directory + "ran2_200_a.bliss");
}
