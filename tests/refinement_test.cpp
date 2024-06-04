// Copyright 2024 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../dejavu.h"

using dejavu::ir::refinement;
using dejavu::ds::coloring;

TEST(refinement_test, trivial_graphs_first1) {
    dejavu::static_graph g1;
    g1.initialize_graph(0, 0);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring_first(g1.get_sgraph(), &c);
}

TEST(refinement_test, trivial_graphs_first2) {
    dejavu::static_graph g1;
    g1.initialize_graph(1, 0);
    g1.add_vertex(0,0);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring_first(g1.get_sgraph(), &c);
    c.check();
}

TEST(refinement_test, trivial_graphs_first3) {
    dejavu::static_graph g1;
    g1.initialize_graph(2, 1);
    g1.add_vertex(0,1);
    g1.add_vertex(1,1);
    g1.add_edge(0, 1);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring_first(g1.get_sgraph(), &c);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[1]);
    c.check();
}

TEST(refinement_test, trivial_graphs_first4) {
    dejavu::static_graph g1;
    g1.initialize_graph(2, 1);
    g1.add_vertex(0,1);
    g1.add_vertex(0,1);
    g1.add_edge(0, 1);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring_first(g1.get_sgraph(), &c);
    EXPECT_EQ(c.vertex_to_col[0], c.vertex_to_col[1]);
    R.individualize_vertex(&c, 0);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[1]);
    c.check();
}

TEST(refinement_test, trivial_graphs_first5) {
    dejavu::static_graph g1;
    g1.initialize_graph(3, 2);
    g1.add_vertex(0,1);
    g1.add_vertex(0,2);
    g1.add_vertex(0,1);
    g1.add_edge(0, 1);
    g1.add_edge(1, 2);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring_first(g1.get_sgraph(), &c);
    EXPECT_EQ(c.vertex_to_col[0], c.vertex_to_col[2]);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[1]);
    R.individualize_vertex(&c, 0);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[2]);
    c.check();
}

TEST(refinement_test, trivial_graphs1) {
    dejavu::static_graph g1;
    g1.initialize_graph(0, 0);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring(g1.get_sgraph(), &c);
}

TEST(refinement_test, trivial_graphs2) {
    dejavu::static_graph g1;
    g1.initialize_graph(1, 0);
    g1.add_vertex(0,0);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring(g1.get_sgraph(), &c);
    c.check();
}

TEST(refinement_test, trivial_graphs3) {
    dejavu::static_graph g1;
    g1.initialize_graph(2, 1);
    g1.add_vertex(0,1);
    g1.add_vertex(1,1);
    g1.add_edge(0, 1);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring(g1.get_sgraph(), &c);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[1]);
    c.check();
}

TEST(refinement_test, trivial_graphs4) {
    dejavu::static_graph g1;
    g1.initialize_graph(2, 1);
    g1.add_vertex(0,1);
    g1.add_vertex(0,1);
    g1.add_edge(0, 1);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring(g1.get_sgraph(), &c);
    EXPECT_EQ(c.vertex_to_col[0], c.vertex_to_col[1]);
    R.individualize_vertex(&c, 0);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[1]);
    c.check();
}

TEST(refinement_test, trivial_graphs5) {
    dejavu::static_graph g1;
    g1.initialize_graph(3, 2);
    g1.add_vertex(0,1);
    g1.add_vertex(0,2);
    g1.add_vertex(0,1);
    g1.add_edge(0, 1);
    g1.add_edge(1, 2);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring(g1.get_sgraph(), &c);
    EXPECT_EQ(c.vertex_to_col[0], c.vertex_to_col[2]);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[1]);
    R.individualize_vertex(&c, 0);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[2]);
    c.check();
}


TEST(refinement_test, repeated_test1) {
    dejavu::static_graph g1;
    g1.initialize_graph(3, 2);
    g1.add_vertex(0,1);
    g1.add_vertex(0,2);
    g1.add_vertex(0,1);
    g1.add_edge(0, 1);
    g1.add_edge(1, 2);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring_first(g1.get_sgraph(), &c);
    EXPECT_EQ(c.vertex_to_col[0], c.vertex_to_col[2]);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[1]);
    R.individualize_vertex(&c, 0);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[2]);
    c.check();

    dejavu::static_graph g2;
    g2.initialize_graph(4, 2);
    g2.add_vertex(0,1);
    g2.add_vertex(0,2);
    g2.add_vertex(0,1);
    g2.add_vertex(0,0);
    g2.add_edge(0, 1);
    g2.add_edge(1, 2);

    g2.get_sgraph()->initialize_coloring(&c, g2.get_coloring());
    R.refine_coloring_first(g2.get_sgraph(), &c);
    EXPECT_EQ(c.vertex_to_col[0], c.vertex_to_col[2]);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[1]);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[3]);
    R.individualize_vertex(&c, 0);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[2]);
    c.check();
}

TEST(refinement_test, repeated_test2) {
    dejavu::static_graph g1;
    g1.initialize_graph(3, 2);
    g1.add_vertex(0,1);
    g1.add_vertex(0,2);
    g1.add_vertex(0,1);
    g1.add_edge(0, 1);
    g1.add_edge(1, 2);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring(g1.get_sgraph(), &c);
    EXPECT_EQ(c.vertex_to_col[0], c.vertex_to_col[2]);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[1]);
    R.individualize_vertex(&c, 0);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[2]);
    c.check();

    dejavu::static_graph g2;
    g2.initialize_graph(4, 2);
    g2.add_vertex(0,1);
    g2.add_vertex(0,2);
    g2.add_vertex(0,1);
    g2.add_vertex(0,0);
    g2.add_edge(0, 1);
    g2.add_edge(1, 2);

    g2.get_sgraph()->initialize_coloring(&c, g2.get_coloring());
    R.refine_coloring(g2.get_sgraph(), &c);
    EXPECT_EQ(c.vertex_to_col[0], c.vertex_to_col[2]);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[1]);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[3]);
    R.individualize_vertex(&c, 0);
    EXPECT_NE(c.vertex_to_col[0], c.vertex_to_col[2]);
    c.check();
}

TEST(refinement_test, certify_test1) {
    dejavu::static_graph g1;
    g1.initialize_graph(3, 2);
    g1.add_vertex(0,1);
    g1.add_vertex(0,2);
    g1.add_vertex(0,1);
    g1.add_edge(0, 1);
    g1.add_edge(1, 2);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring(g1.get_sgraph(), &c);
    int* p = new int[3];
    int* supp = new int[3];
    p[0] = 0;
    p[1] = 1;
    p[2] = 2;

    EXPECT_TRUE(R.certify_automorphism(g1.get_sgraph(), p));
    EXPECT_TRUE(R.certify_automorphism_sparse(g1.get_sgraph(), p, 0, nullptr));

    p[0] = 2;
    p[1] = 1;
    p[2] = 0;

    supp[0] = 0;
    supp[1] = 2;

    EXPECT_TRUE(R.certify_automorphism(g1.get_sgraph(), p));
    EXPECT_TRUE(R.certify_automorphism_sparse(g1.get_sgraph(), p, 2, supp));

    p[0] = 1;
    p[1] = 0;
    p[2] = 2;

    supp[0] = 0;
    supp[1] = 1;

    EXPECT_FALSE(R.certify_automorphism(g1.get_sgraph(), p));
    EXPECT_FALSE(R.certify_automorphism_sparse(g1.get_sgraph(), p, 2, supp));

    EXPECT_FALSE(R.certify_automorphism(g1.get_sgraph(), p));
    EXPECT_FALSE(R.certify_automorphism_sparse(g1.get_sgraph(), p, 2, supp));
}


TEST(refinement_test, certify_test2) {
    dejavu::static_graph g1;
    g1.initialize_graph(4, 3);
    g1.add_vertex(0,1);
    g1.add_vertex(1,2);
    g1.add_vertex(2,2);
    g1.add_vertex(0,1);
    g1.add_edge(0, 1);
    g1.add_edge(1, 2);
    g1.add_edge(2, 3);
    refinement R;
    coloring c;
    g1.get_sgraph()->initialize_coloring(&c, g1.get_coloring());
    R.refine_coloring(g1.get_sgraph(), &c);
    int* p = new int[4];
    int* supp = new int[4];
    p[0] = 0;
    p[1] = 1;
    p[2] = 2;
    p[3] = 3;

    EXPECT_TRUE(R.certify_automorphism(g1.get_sgraph(), p));
    EXPECT_TRUE(R.certify_automorphism_sparse(g1.get_sgraph(), p, 0, nullptr));

    p[0] = 3;
    p[1] = 1;
    p[2] = 2;
    p[3] = 0;

    supp[0] = 0;
    supp[1] = 3;

    ASSERT_EQ(g1.get_coloring()[0], g1.get_coloring()[3]);
    ASSERT_NE(g1.get_coloring()[1], g1.get_coloring()[2]);

    EXPECT_FALSE(R.certify_automorphism(g1.get_sgraph(), p));
    EXPECT_FALSE(R.certify_automorphism_sparse(g1.get_sgraph(), p, 2, supp));

    ASSERT_NE(c.vertex_to_col[0], c.vertex_to_col[3]);
    ASSERT_NE(c.vertex_to_col[1], c.vertex_to_col[2]);

    p[0] = 3;
    p[1] = 2;
    p[2] = 1;
    p[3] = 0;

    supp[0] = 0;
    supp[1] = 3;
    supp[2] = 2;
    supp[3] = 1;

    ASSERT_EQ(g1.get_coloring()[0], g1.get_coloring()[3]);
    ASSERT_NE(g1.get_coloring()[1], g1.get_coloring()[2]);

    EXPECT_TRUE(R.certify_automorphism(g1.get_sgraph(), p));
    EXPECT_TRUE(R.certify_automorphism_sparse(g1.get_sgraph(), p, 4, supp));

    ASSERT_NE(c.vertex_to_col[0], c.vertex_to_col[3]);
    ASSERT_NE(c.vertex_to_col[1], c.vertex_to_col[2]);
}