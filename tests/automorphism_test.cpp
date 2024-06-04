// Copyright 2024 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../groups.h"

using dejavu::groups::automorphism_workspace;

TEST(automorphism, construction) {
    automorphism_workspace aut;
    aut.resize(128);
    EXPECT_EQ(aut.nsupp(), 0);
    for(int i = 0; i < 128; ++i) EXPECT_EQ(aut[i], i);

    aut.write_single_map(7, 8);
    aut.write_single_map(8, 7);
    for(int i = 0; i < 128; ++i) {
        if(i < 7 || i > 8) {
            EXPECT_EQ(aut[i], i);
        }
    }
    EXPECT_EQ(aut[7], 8);
    EXPECT_EQ(aut[8], 7);
    EXPECT_EQ(aut.nsupp(), 2);
}

TEST(automorphism, reset) {
    automorphism_workspace aut;
    aut.resize(128);
    aut.write_single_map(7, 8);
    aut.write_single_map(8, 7);
    for(int i = 0; i < 128; ++i) {
        if(i < 7 || i > 8) {
            EXPECT_EQ(aut[i], i);
        }
    }
    EXPECT_EQ(aut[7], 8);
    EXPECT_EQ(aut[8], 7);
    EXPECT_EQ(aut.nsupp(), 2);

    aut.reset();

    EXPECT_EQ(aut.nsupp(), 0);
    for(int i = 0; i < 128; ++i) EXPECT_EQ(aut[i], i);
}

TEST(automorphism, copy) {
    automorphism_workspace aut1;
    aut1.resize(128);
    aut1.write_single_map(7, 8);
    aut1.write_single_map(8, 7);
    for(int i = 0; i < 128; ++i) {
        if(i < 7 || i > 8) {
            EXPECT_EQ(aut1[i], i);
        }
    }
    EXPECT_EQ(aut1[7], 8);
    EXPECT_EQ(aut1[8], 7);
    EXPECT_EQ(aut1.nsupp(), 2);

    automorphism_workspace aut2 = aut1;

    aut1.reset();

    EXPECT_EQ(aut1.nsupp(), 0);
    for(int i = 0; i < 128; ++i) EXPECT_EQ(aut1[i], i);
    for(int i = 0; i < 128; ++i) {
        if(i < 7 || i > 8) {
            EXPECT_EQ(aut2[i], i);
        }
    }
    EXPECT_EQ(aut2[7], 8);
    EXPECT_EQ(aut2[8], 7);
    EXPECT_EQ(aut2.nsupp(), 2);
}

TEST(automorphism, apply) {
    automorphism_workspace aut1;
    aut1.resize(128);
    aut1.write_single_map(7, 8);
    aut1.write_single_map(8, 9);
    aut1.write_single_map(9, 7);
    automorphism_workspace aut2 = aut1;

    aut2.write_single_map(5, 6);
    aut2.write_single_map(6, 5);

    aut1.apply(aut2, 1);

    for(int i = 0; i < 128; ++i) {
        if(i < 5 || i > 9) {
            EXPECT_EQ(aut1[i], i);
        }
    }
    EXPECT_EQ(aut1[7], 9);
    EXPECT_EQ(aut1[8], 7);
    EXPECT_EQ(aut1[9], 8);
    EXPECT_EQ(aut1[5], 6);
    EXPECT_EQ(aut1[6], 5);
    EXPECT_EQ(aut1.nsupp(), 5);
}
