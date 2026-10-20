// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../dejavu.h"
#include "../parser.h"
#include "helper_functions_test.h"
#include <filesystem>

static void write_temp_file(std::string content, std::string fname) {
    std::ofstream file(fname);

    // Check if the file opened successfully
    EXPECT_TRUE(file);
    file << content;
    file.close();
}

// out of range vertex (1)
TEST(parser_test, invalid_input1) {
    write_temp_file("p edge 3 2\ne 1 2\ne 3 4","____test.dimacs");
    dejavu::sgraph g;
    int*   colmap;
    EXPECT_FALSE(dejavu::parse_dimacs("____test.dimacs", &g, &colmap));
    EXPECT_TRUE(std::remove("____test.dimacs") == 0);
}

// underclaim edges
TEST(parser_test, invalid_input2) {
    write_temp_file("p edge 3 1\ne 1 2\ne 3 2","____test.dimacs");
    dejavu::sgraph g;
    int*   colmap;
    EXPECT_FALSE(dejavu::parse_dimacs("____test.dimacs", &g, &colmap));
    EXPECT_TRUE(std::remove("____test.dimacs") == 0);
}


// double edge
TEST(parser_test, invalid_input3) {
    write_temp_file("p edge 3 2\ne 1 2\ne 2 1","____test.dimacs");
    dejavu::sgraph g;
    int*   colmap;
    EXPECT_FALSE(dejavu::parse_dimacs("____test.dimacs", &g, &colmap));
    EXPECT_TRUE(std::remove("____test.dimacs") == 0);
}


// invalid character
TEST(parser_test, invalid_input4) {
    write_temp_file("p edge 3 2\ne 1 2\nm 2 3","____test.dimacs");
    dejavu::sgraph g;
    int*   colmap;
    EXPECT_FALSE(dejavu::parse_dimacs("____test.dimacs", &g, &colmap));
    EXPECT_TRUE(std::remove("____test.dimacs") == 0);
}


// overclaim edges
TEST(parser_test, invalid_input5) {
    write_temp_file("p edge 7 9\ne 1 2\ne 4 3","____test.dimacs");
    dejavu::sgraph g;
    int*   colmap;
    EXPECT_FALSE(dejavu::parse_dimacs("____test.dimacs", &g, &colmap));
    EXPECT_TRUE(std::remove("____test.dimacs") == 0);
}


// out of range vertex (2)
TEST(parser_test, invalid_input6) {
    write_temp_file("p edge 3 2\ne 1 2\ne 2 3\n 4 0","____test.dimacs");
    dejavu::sgraph g;
    int*   colmap;
    EXPECT_FALSE(dejavu::parse_dimacs("____test.dimacs", &g, &colmap));
    EXPECT_TRUE(std::remove("____test.dimacs") == 0);
}


TEST(parser_test, valid_input1) {
    write_temp_file("p edge 3 2\ne 1 2\ne 2 3","____test.dimacs");
    dejavu::sgraph g;
    int*   colmap;
    EXPECT_TRUE(dejavu::parse_dimacs("____test.dimacs", &g, &colmap));
    EXPECT_EQ(g.degree_of(0), 1);
    EXPECT_EQ(g.degree_of(1), 2);
    EXPECT_EQ(g.degree_of(2), 1);
    EXPECT_TRUE(std::remove("____test.dimacs") == 0);
}


TEST(parser_test, valid_input2) {
    write_temp_file("p edge 3 2\nc this is a comment       \ne 1 2\ne 2 3\nc another comment over here","____test.dimacs");
    dejavu::sgraph g;
    int*   colmap;
    EXPECT_TRUE(dejavu::parse_dimacs("____test.dimacs", &g, &colmap));
    EXPECT_EQ(g.degree_of(0), 1);
    EXPECT_EQ(g.degree_of(1), 2);
    EXPECT_EQ(g.degree_of(2), 1);
    EXPECT_TRUE(std::remove("____test.dimacs") == 0);
}


TEST(parser_test, valid_input3) {
    write_temp_file("p edge 3 2\nc this is a comment       \ne 1 2\ne 2 3\nc another comment over here\nn 1 0\nn 2 1\nn 3 2","____test.dimacs");
    dejavu::sgraph g;
    int*   colmap = nullptr;
    EXPECT_TRUE(dejavu::parse_dimacs("____test.dimacs", &g, &colmap));
    EXPECT_EQ(g.degree_of(0), 1);
    EXPECT_EQ(g.degree_of(1), 2);
    EXPECT_EQ(g.degree_of(2), 1);
    EXPECT_TRUE(std::remove("____test.dimacs") == 0);
}

