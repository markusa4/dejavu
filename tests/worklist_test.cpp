// Copyright 2024 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../ds.h"

using dejavu::ds::worklist;

TEST(worklist_test, basic_construction) {
    worklist m;
    m.allocate(128);
    m.push_back(0);
    m.push_back(1);
    m.push_back(2);
    m.push_back(3);
    EXPECT_EQ(m.size(), 4);
    for(int i = 0; i < m.size(); ++i) EXPECT_EQ(m[i], i);
}

TEST(worklist_test, basic_reset) {
    worklist m;
    m.allocate(128);
    m.push_back(0);
    m.push_back(1);
    m.push_back(2);
    m.push_back(3);
    EXPECT_EQ(m.size(), 4);
    m.reset();
    EXPECT_EQ(m.size(), 0);
}

TEST(worklist_test, copy) {
    worklist m1;
    m1.allocate(128);
    m1.push_back(0);
    m1.push_back(1);
    m1.push_back(2);
    m1.push_back(3);
    worklist m2 = m1;
    EXPECT_EQ(m1.size(), m2.size());
    for(int i = 0; i < m1.size(); ++i) EXPECT_EQ(m1[i], m2[i]);
    m1.reset();
    EXPECT_EQ(m2.size(), 4);
    EXPECT_EQ(m1.size(), 0);
    m1.push_back(1);
    EXPECT_EQ(m2[0], 0);
    EXPECT_EQ(m1[0], 1);
}


TEST(worklist_test, pop) {
    worklist m;
    m.allocate(128);
    m.push_back(0);
    m.push_back(1);
    m.push_back(2);
    m.push_back(3);
    EXPECT_EQ(m.size(), 4);
    for(int i = 0; i < m.size(); ++i) EXPECT_EQ(m[i], i);
    for(int i = 3; i >= 0; --i) EXPECT_EQ(m.pop_back(), i);
}

TEST(worklist_test, resize) {
    worklist m;
    m.allocate(5);
    m.push_back(0);
    m.push_back(1);
    m.push_back(2);
    m.push_back(3);
    EXPECT_EQ(m.size(), 4);
    for(int i = 0; i < m.size(); ++i) EXPECT_EQ(m[i], i);
    m.resize(256);
    m.push_back(4);
    m.push_back(5);
    m.push_back(6);
    m.push_back(7);
    m.push_back(8);
    for(int i = 0; i < m.size(); ++i) EXPECT_EQ(m[i], i);
}

TEST(worklist_test, sort) {
    worklist m;
    m.allocate(128);
    m.push_back(3);
    m.push_back(2);
    m.push_back(1);
    m.push_back(0);
    EXPECT_EQ(m.size(), 4);
    m.sort();
    for(int i = 0; i < m.size(); ++i) EXPECT_EQ(m[i], i);
}