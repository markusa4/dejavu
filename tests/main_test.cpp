#include "gtest/gtest.h"
#include "../dejavu.h"

thread_local bool bulk_domain_reset = false;

dejavu::ir::refinement test_r;
dejavu::sgraph dej_test_graph;
int*   dej_test_col;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}