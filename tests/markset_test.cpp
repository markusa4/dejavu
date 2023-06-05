#include "gtest/gtest.h"
#include "../ds.h"

using dejavu::ds::mark_set;

TEST(markset_test, basic_construction) {
    mark_set m;
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

    mark_set m2(128);
    for(int i = 0; i < 128; ++i) {
        EXPECT_EQ(m2.get(i), false);
    }
}

TEST(markset_test, basic_reset) {
    mark_set m;
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
    mark_set m;
    m.initialize(128);
    EXPECT_EQ(m.get(47), false);
    m.set(47);
    EXPECT_EQ(m.get(47), true);
    m.unset(47);
    EXPECT_EQ(m.get(47), false);
}

TEST(markset_test, many_resets) {
    mark_set m;
    m.initialize(4);
    for(long i = 0; i < INT16_MAX; ++i) {
        EXPECT_EQ(m.get(1), false);
        m.set(1);
        EXPECT_EQ(m.get(1), true);
        m.reset();
    }
}