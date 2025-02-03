// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../ds.h"

using dejavu::ds::workset_t;

TEST(workset_test, basic_construction) {
    workset_t<int> m;
    m.initialize(128);
    for(int i = 0; i < 128; ++i) {
        EXPECT_EQ(m.get(i), -1);
    }
    m.set(0, 1);
    EXPECT_EQ(m.get(0), 1);
    m.initialize(128);
    EXPECT_EQ(m.get(0), -1);
    m.set(0, 1);
    EXPECT_EQ(m.get(0), 1);
    m.reset();
    m.initialize(127);
    for(int i = 0; i < 127; ++i) {
        EXPECT_EQ(m.get(i), -1);
    }

    workset_t<int> m2(128);
    for(int i = 0; i < 128; ++i) {
        EXPECT_EQ(m2.get(i), -1);
    }
}

TEST(workset_test, basic_reset) {
    workset_t<int> m;
    m.initialize(128);
    for(int i = 0; i < 128; ++i) {
        if(i % 2 == 0) m.set(i, 1);
    }
    for(int i = 0; i < 128; ++i) {
        if(i % 2 == 0) EXPECT_EQ(m.get(i), 1);
        else EXPECT_EQ(m.get(i), -1);
    }
    m.reset_hard();
    for(int i = 0; i < 128; ++i) {
        EXPECT_EQ(m.get(i), -1);
    }
}

TEST(workset_test, copy) {
    workset_t<int> m1;
    m1.initialize(128);
    EXPECT_EQ(m1.get(47), -1);
    m1.set(47, 1);
    EXPECT_EQ(m1.get(47), 1);

    workset_t<int> m2 = m1;
    EXPECT_EQ(m2.get(47), 1);
    EXPECT_EQ(m2.get(48), -1);
    m1.reset_hard();
    EXPECT_EQ(m1.get(47), -1);
    EXPECT_EQ(m1.get(48), -1);
    EXPECT_EQ(m2.get(47), 1);
    EXPECT_EQ(m2.get(48), -1);
}

TEST(workset_test, assignment) {
    workset_t<int> m1;
    m1.initialize(128);
    EXPECT_EQ(m1.get(47), -1);
    m1.set(47, 1);
    EXPECT_EQ(m1.get(47), 1);

    workset_t<int> m2(12);
    m2.set(5, 1);
    EXPECT_EQ(m2.get(5), 1);

    m2 = m1;
    EXPECT_EQ(m2.get(47), 1);
    EXPECT_EQ(m2.get(48), -1);
    EXPECT_EQ(m2.get(5), -1);
    m1.reset_hard();
    EXPECT_EQ(m1.get(47), -1);
    EXPECT_EQ(m1.get(48), -1);
    EXPECT_EQ(m2.get(47), 1);
    EXPECT_EQ(m2.get(48), -1);
}

TEST(workset_test, inc_nr) {
    workset_t<int> m;
    m.initialize(4);
    EXPECT_EQ(m.get(1), -1);
    for(long i = 0; i < INT16_MAX; ++i) {
        m.inc_nr(1);
        EXPECT_EQ(m.get(1), i);
    }
}

TEST(workset_test, inc) {
    workset_t<int> m;
    m.initialize(4);
    EXPECT_EQ(m.get(1), -1);
    for(long i = 0; i < INT16_MAX; ++i) {
        m.inc(1);
        EXPECT_EQ(m.get(1), i);
    }
}

TEST(workset_test, inc_and_reset) {
    workset_t<int> m;
    m.initialize(4);
    EXPECT_EQ(m.get(1), -1);
    for(long i = 0; i < INT16_MAX; ++i) {
        m.inc(1);
        EXPECT_EQ(m.get(1), i);
    }
    m.reset();
    EXPECT_EQ(m.get(1), -1);
}


TEST(workset_test, copy_reset_queue) {
    workset_t<int> m1;
    m1.initialize(128);
    EXPECT_EQ(m1.get(47), -1);
    m1.inc(47);
    m1.inc(47);
    EXPECT_EQ(m1.get(47), 1);

    workset_t<int> m2 = m1;
    EXPECT_EQ(m2.get(47), 1);
    EXPECT_EQ(m2.get(48), -1);
    m2.reset();
    EXPECT_EQ(m2.get(47), -1);
    EXPECT_EQ(m2.get(48), -1);
    EXPECT_EQ(m1.get(47), 1);
    EXPECT_EQ(m1.get(48), -1);
    m1.reset();
    EXPECT_EQ(m1.get(47), -1);
    EXPECT_EQ(m1.get(48), -1);
}

