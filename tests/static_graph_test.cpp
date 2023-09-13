// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../dejavu.h"

TEST(static_graph_test, edge_case_test1) {
    dejavu::static_graph g1;
    g1.initialize_graph(0, 0);
    EXPECT_ANY_THROW(g1.add_vertex(0, 0));
}

TEST(static_graph_test, edge_case_test2) {
    dejavu::static_graph g1;
    g1.initialize_graph(2, 3);
    g1.add_vertex(0, 1);
    g1.add_vertex(0, 2);
    EXPECT_ANY_THROW(g1.add_vertex(0, 0));
}

TEST(static_graph_test, edge_case_test3) {
    dejavu::static_graph g1;
    g1.initialize_graph(2, 0);
    g1.add_vertex(0, 1);
    g1.add_vertex(0, 1);
    EXPECT_ANY_THROW(g1.add_edge(0, 1));
}

TEST(static_graph_test, edge_case_test4) {
    dejavu::static_graph g1;
    EXPECT_ANY_THROW(g1.add_edge(0, 1));
    EXPECT_ANY_THROW(g1.add_vertex(0, 1));
    EXPECT_ANY_THROW(g1.get_sgraph());
    EXPECT_ANY_THROW(g1.get_coloring());
}

TEST(static_graph_test, edge_case_test5) {
    dejavu::static_graph g1;
    g1.initialize_graph(2, 1);
    g1.add_vertex(0, 1);
    g1.add_vertex(0, 1);
    EXPECT_ANY_THROW(g1.add_edge(1, 0));
    EXPECT_ANY_THROW(g1.add_edge(0, 3));
    EXPECT_ANY_THROW(g1.add_edge(3, 4));
}

TEST(static_graph_test, edge_case_test6) {
    dejavu::static_graph g1;
    g1.initialize_graph(2, 1);
    EXPECT_ANY_THROW(g1.initialize_graph(0, 1));
    g1.add_vertex(0, 1);
    g1.add_vertex(0, 1);
    g1.add_edge(0, 1);
    g1.get_sgraph();

    EXPECT_ANY_THROW(g1.add_edge(0, 1));
    EXPECT_ANY_THROW(g1.add_vertex(0, 1));
    EXPECT_ANY_THROW(g1.initialize_graph(0, 1));
}
