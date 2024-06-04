// Copyright 2024 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../ds.h"

using dejavu::ds::workspace;

TEST(workspace_test, construction) {
    workspace m;
    m.allocate(128);
    m[127] = 1;
    m[42] = 8;
    EXPECT_EQ(m[127], 1);
    EXPECT_EQ(m[42], 8);
}

TEST(workspace_test, reset) {
    workspace m;
    m.allocate(128);
    m.reset();
    for(int i = 0; i < 128; ++i) EXPECT_EQ(m[i], 0);
}

TEST(workspace_test, resize) {
    workspace m;
    m.allocate(128);
    for(int i = 0; i < 128; ++i) m[i] = i;
    m.resize(256);
    for(int i = 0; i < 128; ++i) EXPECT_EQ(m[i], i);
}

TEST(workspace_test, copy) {
    workspace m1;
    m1.allocate(128);
    m1[127] = 1;
    m1[42] = 8;
    EXPECT_EQ(m1[127], 1);
    EXPECT_EQ(m1[42], 8);

    workspace m2 = m1;
    EXPECT_EQ(m2[127], 1);
    EXPECT_EQ(m2[42], 8);

    m2[42] = 7;
    EXPECT_EQ(m1[42], 8);
    EXPECT_EQ(m2[42], 7);
}