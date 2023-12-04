// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_HELPER_FUNCTIONS_TEST_H
#define DEJAVU_HELPER_FUNCTIONS_TEST_H

#include "gtest/gtest.h"
#include "../dejavu.h"

extern      dejavu::ir::refinement dgtest_test_r;
extern      dejavu::sgraph         dgtest_graph;
extern int* dgtest_col;

static void gtest_certify_hook([[maybe_unused]] int n, [[maybe_unused]] const int *p, [[maybe_unused]] int nsupp,
                      [[maybe_unused]] const int *supp) {
    EXPECT_TRUE(dgtest_test_r.certify_automorphism_sparse(&dgtest_graph, p, nsupp, supp));
    //EXPECT_TRUE(dgtest_test_r->certify_automorphism(dgtest_graph, dgtest_col, p));
    //EXPECT_TRUE(dgtest_test_r->certify_automorphism_sparse(dgtest_graph, dgtest_col, p, nsupp, supp));
}

#endif //DEJAVU_HELPER_FUNCTIONS_TEST_H
