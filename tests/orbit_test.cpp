// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../dejavu.h"

using dejavu::groups::orbit;

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
    o.combine_orbits(1, 8);
    ASSERT_TRUE(o.are_in_same_orbit(1, 7));
    ASSERT_TRUE(o.are_in_same_orbit(1, 8));
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