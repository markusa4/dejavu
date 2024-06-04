// Copyright 2024 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../ds.h"

using dejavu::ds::markset;

TEST(markset_test, basic_construction) {
    markset m;
    m.initialize(128);
    for(int i = 0; i < 128; ++i) {
        EXPECT_EQ(m.get(i), false);
    }
    m.set(0);
    EXPECT_EQ(m.get(0), true);
    m.initialize(128);
    EXPECT_EQ(m.get(0), false);
    m.set(0);
    EXPECT_EQ(m.get(0), true);
    m.reset();
    m.initialize(127);
    for(int i = 0; i < 127; ++i) {
        EXPECT_EQ(m.get(i), false);
    }

    markset m2(128);
    for(int i = 0; i < 128; ++i) {
        EXPECT_EQ(m2.get(i), false);
    }
}

TEST(markset_test, basic_reset) {
    markset m;
    m.initialize(128);
    for(int i = 0; i < 128; ++i) {
        if(i % 2 == 0) m.set(i);
    }
    for(int i = 0; i < 128; ++i) {
        if(i % 2 == 0) EXPECT_EQ(m.get(i), true);
        else EXPECT_EQ(m.get(i), false);
    }
    m.reset();
    for(int i = 0; i < 128; ++i) {
        EXPECT_EQ(m.get(i), false);
    }
}

TEST(markset_test, unset) {
    markset m;
    m.initialize(128);
    EXPECT_EQ(m.get(47), false);
    m.set(47);
    EXPECT_EQ(m.get(47), true);
    m.unset(47);
    EXPECT_EQ(m.get(47), false);
}

TEST(markset_test, copy) {
    markset m1;
    m1.initialize(128);
    EXPECT_EQ(m1.get(47), false);
    m1.set(47);
    EXPECT_EQ(m1.get(47), true);

    markset m2 = m1;
    EXPECT_EQ(m2.get(47), true);
    EXPECT_EQ(m2.get(48), false);
    m1.reset();
    EXPECT_EQ(m1.get(47), false);
    EXPECT_EQ(m1.get(48), false);
    EXPECT_EQ(m2.get(47), true);
    EXPECT_EQ(m2.get(48), false);
}

TEST(markset_test, assignment) {
    markset m1;
    m1.initialize(128);
    EXPECT_EQ(m1.get(47), false);
    m1.set(47);
    EXPECT_EQ(m1.get(47), true);

    markset m2(12);
    m2.set(5);
    EXPECT_EQ(m2.get(5), true);

    m2 = m1;
    EXPECT_EQ(m2.get(47), true);
    EXPECT_EQ(m2.get(48), false);
    EXPECT_EQ(m2.get(5), false);
    m1.reset();
    EXPECT_EQ(m1.get(47), false);
    EXPECT_EQ(m1.get(48), false);
    EXPECT_EQ(m2.get(47), true);
    EXPECT_EQ(m2.get(48), false);
}

TEST(markset_test, many_resets) {
    markset m;
    m.initialize(4);
    for(long i = 0; i < INT16_MAX; ++i) {
        EXPECT_EQ(m.get(1), false);
        m.set(1);
        EXPECT_EQ(m.get(1), true);
        m.reset();
    }
}