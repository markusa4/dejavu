// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../dejavu.h"

int count_auto = 0;
void count_hook(int n, const int *p, int nsupp, const int *supp) {
    ++count_auto;
}

TEST(simple_graphs_test, trivial_graphs1) {
    auto test_hook = dejavu_hook(count_hook);

    dejavu::static_graph g1;
    g1.initialize_graph(0, 0);
    dejavu::solver d;
    d.automorphisms(&g1, &test_hook);

    EXPECT_EQ(count_auto, 0);
    EXPECT_EQ(d.get_automorphism_group_size().exponent, 0);
    EXPECT_NEAR(d.get_automorphism_group_size().mantissa, 1.0, 0.001);
}

TEST(simple_graphs_test, trivial_graphs2) {
    auto test_hook = dejavu_hook(count_hook);

    dejavu::static_graph g1;
    g1.initialize_graph(1, 0);
    g1.add_vertex(0,0);
    dejavu::solver d;
    d.automorphisms(&g1, &test_hook);

    EXPECT_EQ(count_auto, 0);
    EXPECT_EQ(d.get_automorphism_group_size().exponent, 0);
    EXPECT_NEAR(d.get_automorphism_group_size().mantissa, 1.0, 0.001);
}

TEST(simple_graphs_test, trivial_graphs3) {
    auto test_hook = dejavu_hook(count_hook);

    dejavu::static_graph g1;
    g1.initialize_graph(2, 1);
    g1.add_vertex(0,1);
    g1.add_vertex(1,1);
    g1.add_edge(0, 1);
    dejavu::solver d;
    d.automorphisms(&g1, &test_hook);

    EXPECT_EQ(count_auto, 0);
    EXPECT_EQ(d.get_automorphism_group_size().exponent, 0);
    EXPECT_NEAR(d.get_automorphism_group_size().mantissa, 1.0, 0.001);
}

TEST(simple_graphs_test, trivial_graphs4) {
    auto test_hook = dejavu_hook(count_hook);

    dejavu::static_graph g1;
    g1.initialize_graph(2, 1);
    g1.add_vertex(0,1);
    g1.add_vertex(0,1);
    g1.add_edge(0, 1);
    dejavu::solver d;
    d.automorphisms(&g1, &test_hook);
    EXPECT_EQ(d.get_automorphism_group_size().exponent, 0);
    EXPECT_NEAR(d.get_automorphism_group_size().mantissa, 2.0, 0.001);
}

TEST(simple_graphs_test, trivial_graphs5) {
    auto test_hook = dejavu_hook(count_hook);

    dejavu::static_graph g1;
    g1.initialize_graph(3, 2);
    g1.add_vertex(0,1);
    g1.add_vertex(0,2);
    g1.add_vertex(0,1);
    g1.add_edge(0, 1);
    g1.add_edge(1, 2);
    dejavu::solver d;
    d.automorphisms(&g1, &test_hook);
    EXPECT_EQ(d.get_automorphism_group_size().exponent, 0);
    EXPECT_NEAR(d.get_automorphism_group_size().mantissa, 2.0, 0.001);
}

TEST(simple_graphs_test, option_test1) {
    auto test_hook = dejavu_hook(count_hook);

    dejavu::static_graph g1;
    g1.initialize_graph(10, 0);
    for(int i = 0; i < 10; ++i) g1.add_vertex(0,0);

    dejavu::solver d;
    d.set_disallow_alteration();
    d.set_strong_certification();
    d.automorphisms(&g1, &test_hook);
    EXPECT_EQ(g1.get_sgraph()->v_size, 10);
}

TEST(simple_graphs_test, option_test2) {
    auto test_hook = dejavu_hook(count_hook);

    dejavu::static_graph g1;
    g1.initialize_graph(10, 0);
    for(int i = 0; i < 10; ++i) g1.add_vertex(0,0);

    dejavu::solver d;
    d.set_print(false); // to make sure it does not crash...
    d.automorphisms(&g1, &test_hook);
}

// Two copies of a component with two hubs h1, h2: h2 hangs off a star s
// (three further leaves), and two paths leaf - a - p - b - h2 with a and b also
// adjacent to h1. Every generator reported for this graph must be an
// automorphism of it; with the default preprocessing, some generators move
// a, p, b of one path onto the other path but leave the degree-1 leaves fixed.
TEST(simple_graphs_test, degree1_leaves_lifted_with_their_path) {
    const int nv = 28;
    const std::vector<std::pair<int, int>> edges = {
        {1, 2}, {2, 3}, {2, 4}, {2, 5}, {6, 7}, {7, 8}, {8, 9}, {0, 7}, {0, 9}, {1, 9},
        {10, 11}, {11, 12}, {12, 13}, {0, 11}, {0, 13}, {1, 13}, {15, 16}, {16, 17}, {16, 18},
        {16, 19}, {20, 21}, {21, 22}, {22, 23}, {14, 21}, {14, 23}, {15, 23}, {24, 25}, {25, 26},
        {26, 27}, {14, 25}, {14, 27}, {15, 27}};
    std::vector<int> degree(nv, 0);
    for (const auto& e : edges) {
        ++degree[e.first];
        ++degree[e.second];
    }

    dejavu::static_graph g1;
    g1.initialize_graph(nv, static_cast<int>(edges.size()));
    for (int i = 0; i < nv; ++i) g1.add_vertex(0, degree[i]);
    for (const auto& e : edges) g1.add_edge(e.first, e.second);

    // keep a copy of the input to certify the reported generators against
    dejavu::sgraph original;
    original.copy_graph(g1.get_sgraph());
    dejavu::markset scratch_set(nv);

    int generators = 0;
    int not_automorphisms = 0;
    auto test_hook = dejavu_hook([&](int, const int *p, int nsupp, const int *supp) {
        ++generators;
        if (!dejavu::ir::certification::certify_automorphism_sparse(scratch_set, &original, p, nsupp, supp))
            ++not_automorphisms;
    });

    dejavu::solver d;
    d.set_print(false);
    d.automorphisms(&g1, &test_hook);

    EXPECT_GT(generators, 0);
    EXPECT_EQ(not_automorphisms, 0);
    EXPECT_EQ(d.get_automorphism_group_size().exponent, 2);
    EXPECT_NEAR(d.get_automorphism_group_size().mantissa, 2.88, 0.001);
}
