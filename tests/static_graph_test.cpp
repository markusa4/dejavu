// Copyright 2024 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../dejavu.h"


TEST(static_graph_test, construct) {
    dejavu::static_graph g1;
    g1.initialize_graph(4, 1);
    g1.add_vertex(0, 0);
    g1.add_vertex(0, 0);
    g1.add_vertex(0, 1);
    g1.add_vertex(0, 1);
    g1.add_edge(2, 3);

    assert(g1.get_sgraph()->d[0] == 0);
    assert(g1.get_sgraph()->d[1] == 0);
    assert(g1.get_sgraph()->d[2] == 1);
    assert(g1.get_sgraph()->d[3] == 1);

    assert(g1.get_coloring()[0] == 0);
    assert(g1.get_coloring()[1] == 0);
    assert(g1.get_coloring()[2] == 0);
    assert(g1.get_coloring()[3] == 0);
}

TEST(static_graph_test, copy1) {
    dejavu::static_graph g1;
    g1.initialize_graph(4, 1);
    dejavu::static_graph g2;
    g2 = g1;

    g1.add_vertex(0, 0);
    g1.add_vertex(0, 0);
    g1.add_vertex(0, 1);
    g1.add_vertex(0, 1);
    g1.add_edge(2, 3);

    g2.add_vertex(1, 1);
    g2.add_vertex(1, 1);
    g2.add_vertex(1, 0);
    g2.add_vertex(1, 0);
    g2.add_edge(0, 1);

    assert(g1.get_sgraph()->d[0] == 0);
    assert(g1.get_sgraph()->d[1] == 0);
    assert(g1.get_sgraph()->d[2] == 1);
    assert(g1.get_sgraph()->d[3] == 1);

    assert(g1.get_coloring()[0] == 0);
    assert(g1.get_coloring()[1] == 0);
    assert(g1.get_coloring()[2] == 0);
    assert(g1.get_coloring()[3] == 0);

    assert(g2.get_sgraph()->d[0] == 1);
    assert(g2.get_sgraph()->d[1] == 1);
    assert(g2.get_sgraph()->d[2] == 0);
    assert(g2.get_sgraph()->d[3] == 0);

    assert(g2.get_coloring()[0] == 1);
    assert(g2.get_coloring()[1] == 1);
    assert(g2.get_coloring()[2] == 1);
    assert(g2.get_coloring()[3] == 1);
}

TEST(static_graph_test, copy2) {
    dejavu::static_graph g1;
    g1.initialize_graph(4, 1);
    dejavu::static_graph g2 = g1;

    g1.add_vertex(0, 0);
    g1.add_vertex(0, 0);
    g1.add_vertex(0, 1);
    g1.add_vertex(0, 1);
    g1.add_edge(2, 3);

    g2.add_vertex(1, 1);
    g2.add_vertex(1, 1);
    g2.add_vertex(1, 0);
    g2.add_vertex(1, 0);
    g2.add_edge(0, 1);

    assert(g1.get_sgraph()->d[0] == 0);
    assert(g1.get_sgraph()->d[1] == 0);
    assert(g1.get_sgraph()->d[2] == 1);
    assert(g1.get_sgraph()->d[3] == 1);

    assert(g1.get_coloring()[0] == 0);
    assert(g1.get_coloring()[1] == 0);
    assert(g1.get_coloring()[2] == 0);
    assert(g1.get_coloring()[3] == 0);

    assert(g2.get_sgraph()->d[0] == 1);
    assert(g2.get_sgraph()->d[1] == 1);
    assert(g2.get_sgraph()->d[2] == 0);
    assert(g2.get_sgraph()->d[3] == 0);

    assert(g2.get_coloring()[0] == 1);
    assert(g2.get_coloring()[1] == 1);
    assert(g2.get_coloring()[2] == 1);
    assert(g2.get_coloring()[3] == 1);
}

TEST(static_graph_test, copy3) {
    dejavu::static_graph g1;
    dejavu::static_graph g2 = g1;
    g1.initialize_graph(4, 1);
    g2.initialize_graph(4, 1);

    g1.add_vertex(0, 0);
    g1.add_vertex(0, 0);
    g1.add_vertex(0, 1);
    g1.add_vertex(0, 1);
    g1.add_edge(2, 3);

    g2.add_vertex(1, 1);
    g2.add_vertex(1, 1);
    g2.add_vertex(1, 0);
    g2.add_vertex(1, 0);
    g2.add_edge(0, 1);

    assert(g1.get_sgraph()->d[0] == 0);
    assert(g1.get_sgraph()->d[1] == 0);
    assert(g1.get_sgraph()->d[2] == 1);
    assert(g1.get_sgraph()->d[3] == 1);

    assert(g1.get_coloring()[0] == 0);
    assert(g1.get_coloring()[1] == 0);
    assert(g1.get_coloring()[2] == 0);
    assert(g1.get_coloring()[3] == 0);

    assert(g2.get_sgraph()->d[0] == 1);
    assert(g2.get_sgraph()->d[1] == 1);
    assert(g2.get_sgraph()->d[2] == 0);
    assert(g2.get_sgraph()->d[3] == 0);

    assert(g2.get_coloring()[0] == 1);
    assert(g2.get_coloring()[1] == 1);
    assert(g2.get_coloring()[2] == 1);
    assert(g2.get_coloring()[3] == 1);
}

TEST(static_graph_test, copy4) {
    dejavu::static_graph g1;
    dejavu::static_graph g2;

    g2.initialize_graph(7, 12);
    g1.initialize_graph(4, 1);

    g2 = g1;

    g1.add_vertex(0, 0);
    g1.add_vertex(0, 0);
    g1.add_vertex(0, 1);
    g1.add_vertex(0, 1);
    g1.add_edge(2, 3);

    g2.add_vertex(1, 1);
    g2.add_vertex(1, 1);
    g2.add_vertex(1, 0);
    g2.add_vertex(1, 0);
    g2.add_edge(0, 1);

    assert(g1.get_sgraph()->d[0] == 0);
    assert(g1.get_sgraph()->d[1] == 0);
    assert(g1.get_sgraph()->d[2] == 1);
    assert(g1.get_sgraph()->d[3] == 1);

    assert(g1.get_coloring()[0] == 0);
    assert(g1.get_coloring()[1] == 0);
    assert(g1.get_coloring()[2] == 0);
    assert(g1.get_coloring()[3] == 0);

    assert(g2.get_sgraph()->d[0] == 1);
    assert(g2.get_sgraph()->d[1] == 1);
    assert(g2.get_sgraph()->d[2] == 0);
    assert(g2.get_sgraph()->d[3] == 0);

    assert(g2.get_coloring()[0] == 1);
    assert(g2.get_coloring()[1] == 1);
    assert(g2.get_coloring()[2] == 1);
    assert(g2.get_coloring()[3] == 1);
}

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
