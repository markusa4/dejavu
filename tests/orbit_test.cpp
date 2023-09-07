// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../dejavu.h"

using dejavu::groups::orbit;
using dejavu::hooks::orbit_hook;

TEST(orbit_test, orbit_construction_test) {
    orbit o;
    o.initialize(16);
    ASSERT_FALSE(o.are_in_same_orbit(1, 7));
    ASSERT_FALSE(o.are_in_same_orbit(1, 8));
    o.combine_orbits(1, 7);
    ASSERT_TRUE(o.are_in_same_orbit(1, 7));
    ASSERT_FALSE(o.are_in_same_orbit(1, 8));

    o.initialize(26);
    ASSERT_FALSE(o.are_in_same_orbit(1, 7));
    ASSERT_FALSE(o.are_in_same_orbit(1, 8));
    o.combine_orbits(1, 7);
    ASSERT_TRUE(o.are_in_same_orbit(1, 7));
    ASSERT_FALSE(o.are_in_same_orbit(1, 8));

    o.initialize(15);
    ASSERT_FALSE(o.are_in_same_orbit(1, 7));
    ASSERT_FALSE(o.are_in_same_orbit(1, 8));
    o.combine_orbits(1, 7);
    ASSERT_TRUE(o.are_in_same_orbit(1, 7));
    ASSERT_FALSE(o.are_in_same_orbit(1, 8));
}

TEST(orbit_test, orbit_combine_test) {
    orbit o;
    o.initialize(127);
    ASSERT_FALSE(o.are_in_same_orbit(1, 7));
    ASSERT_FALSE(o.are_in_same_orbit(1, 8));
    o.combine_orbits(1, 7);
    ASSERT_TRUE(o.are_in_same_orbit(1, 7));
    ASSERT_FALSE(o.are_in_same_orbit(1, 8));
    ASSERT_TRUE(o.represents_orbit(1));
    ASSERT_FALSE(o.represents_orbit(7));
    o.combine_orbits(1, 8);
    ASSERT_TRUE(o.are_in_same_orbit(1, 7));
    ASSERT_TRUE(o.are_in_same_orbit(1, 8));
    ASSERT_TRUE(o.represents_orbit(1));
    ASSERT_FALSE(o.represents_orbit(7));
    ASSERT_FALSE(o.represents_orbit(8));
    ASSERT_FALSE(o.are_in_same_orbit(15, 8));
    o.combine_orbits(8, 15);
    ASSERT_TRUE(o.are_in_same_orbit(1, 15));
    ASSERT_TRUE(o.are_in_same_orbit(15, 7));
    ASSERT_EQ(o.find_orbit(1), o.find_orbit(15));
    ASSERT_EQ(o.find_orbit(7), o.find_orbit(15));
    ASSERT_EQ(o.find_orbit(8), o.find_orbit(15));
    ASSERT_NE(o.find_orbit(27), o.find_orbit(15));

    ASSERT_TRUE(o.are_in_same_orbit(1, 1));
    ASSERT_TRUE(o.are_in_same_orbit(15, 15));
    ASSERT_TRUE(o.are_in_same_orbit(28, 28));
}

TEST(orbit_test, orbit_hook_test) {
    orbit o(8);
    orbit_hook hook(o);

    dejavu::static_graph g1;
    g1.initialize_graph(3, 3);
    g1.add_vertex(0,2);
    g1.add_vertex(0,2);
    g1.add_vertex(0,2);
    g1.add_edge(0, 1);
    g1.add_edge(1, 2);
    g1.add_edge(0, 2);

    dejavu::solver d1;
    d1.automorphisms(&g1, hook.get_hook());

    ASSERT_TRUE(o.are_in_same_orbit(0, 1));
    ASSERT_TRUE(o.are_in_same_orbit(1, 2));
    ASSERT_FALSE(o.are_in_same_orbit(2, 3));
    ASSERT_FALSE(o.are_in_same_orbit(3, 4));
    ASSERT_FALSE(o.are_in_same_orbit(4, 5));

    ASSERT_TRUE(o.represents_orbit(0));
    ASSERT_FALSE(o.represents_orbit(1));
    ASSERT_FALSE(o.represents_orbit(2));
    ASSERT_TRUE(o.represents_orbit(3));
    ASSERT_TRUE(o.represents_orbit(4));
    ASSERT_TRUE(o.represents_orbit(5));
    ASSERT_TRUE(o.represents_orbit(6));
    ASSERT_TRUE(o.represents_orbit(7));

    o.reset();
    ASSERT_FALSE(o.are_in_same_orbit(0, 1));

    dejavu::static_graph g2;
    g2.initialize_graph(8, 3);
    g2.add_vertex(0,2);
    g2.add_vertex(0,2);
    g2.add_vertex(0,2);
    g2.add_vertex(0,0);
    g2.add_vertex(0,0);
    g2.add_vertex(0,0);
    g2.add_vertex(0,0);
    g2.add_vertex(0,0);
    g2.add_edge(0, 1);
    g2.add_edge(1, 2);
    g2.add_edge(0, 2);

    dejavu::solver d2;
    d2.automorphisms(&g2, hook.get_hook());

    ASSERT_TRUE(o.are_in_same_orbit(0, 1));
    ASSERT_TRUE(o.are_in_same_orbit(1, 2));
    ASSERT_FALSE(o.are_in_same_orbit(2, 3));
    ASSERT_TRUE(o.are_in_same_orbit(3, 4));
    ASSERT_TRUE(o.are_in_same_orbit(4, 5));
    ASSERT_TRUE(o.are_in_same_orbit(5, 6));
    ASSERT_TRUE(o.are_in_same_orbit(6, 7));

    ASSERT_TRUE(o.represents_orbit(0));
    ASSERT_FALSE(o.represents_orbit(1));
    ASSERT_FALSE(o.represents_orbit(2));
    ASSERT_TRUE(o.represents_orbit(3));
    ASSERT_FALSE(o.represents_orbit(4));
    ASSERT_FALSE(o.represents_orbit(5));
    ASSERT_FALSE(o.represents_orbit(6));
    ASSERT_FALSE(o.represents_orbit(7));
}