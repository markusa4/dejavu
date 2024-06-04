// Copyright 2024 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../dejavu.h"

thread_local bool bulk_domain_reset = false;

dejavu::ir::refinement dgtest_test_r;
dejavu::sgraph         dgtest_graph;
int*   dej_test_col;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}