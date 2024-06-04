// Copyright 2024 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../dejavu.h"
#include "helper_functions_test.h"
#include <filesystem>


dejavu::groups::orbit orbits1;
dejavu::groups::orbit orbits2;

static void orbit1_test_hook(int n, const int *p, int nsupp, const int *supp) {
    assert(dgtest_test_r.certify_automorphism_sparse(&dgtest_graph, p, nsupp, supp));
    orbits1.add_automorphism_to_orbit(p, nsupp, supp);
}

static void orbit2_test_hook(int n, const int *p, int nsupp, const int *supp) {
    assert(dgtest_test_r.certify_automorphism_sparse(&dgtest_graph, p, nsupp, supp));
    orbits2.add_automorphism_to_orbit(p, nsupp, supp);
}

void read_symmetry_file(const std::string& filename, int domain_size, dejavu_hook* hook) {
    dejavu::groups::automorphism_workspace auto_ws(domain_size);
    std::ifstream infile(filename);
    std::string line;
    while (std::getline(infile, line)){
        auto_ws.reset();
        std::string from_str;
        std::string to_str;
        int read_mode = 0;
        for(std::string::size_type i = 0; i < line.size(); ++i) {
            if (line[i] == '(' && read_mode == 0) {
                read_mode = 1;
                from_str = "";
                to_str   = "";
            } else if (read_mode == 1 && line[i] != '-') {
                from_str += line[i];
            } else if (read_mode == 1 && line[i] == '-') {
                read_mode = 2;
            } else if (read_mode == 2 && line[i] == '>') {
                read_mode = 3;
            } else if (read_mode == 3 && line[i] != ')') {
                to_str += line[i];
            } else if (read_mode == 3 && line[i] == ')') {
                read_mode = 0;
                const int from = std::stoi(from_str);
                const int to   = std::stoi(to_str);
                auto_ws.write_single_map(from, to);
            }
        }
        (*hook)(domain_size, auto_ws.p(), auto_ws.nsupp(), auto_ws.supp());
        auto_ws.reset();
    }
}

dejavu::big_number read_grp_sz_file(const std::string& filename) {
    std::ifstream infile(filename);
    std::string line;
    dejavu::big_number result;
    while (std::getline(infile, line)){
        std::string man_str;
        std::string exp_str;
        int read_mode = 0;
        for(std::string::size_type i = 0; i < line.size(); ++i) {
            if (read_mode == 0 && line[i] != '*') {
                man_str += line[i];
            } else if (read_mode == 0 && line[i] == '*') {
                read_mode = 1;
            } else if (read_mode == 1 && line[i] == '^') {
                read_mode = 2;
            } else if (read_mode == 2) {
                exp_str += line[i];
            }
        }
        result.mantissa = std::stod(man_str);
        result.exponent = std::stoi(exp_str);
        break;
    }

    return result;
}

void test_graph_orbit_check(std::string filename) {
    std::string sym_filename = filename + ".sym";
    std::string grp_sz_filename = filename + ".grp_sz";
    dejavu::sgraph *g = new dejavu::sgraph();
    std::cout << "Parsing " << filename << "..." << std::endl;
    int* colmap = nullptr;
    parse_dimacs(filename, g, &colmap);

    dgtest_graph.copy_graph(g);
    const int domain_size = g->v_size;
    orbits1.initialize(domain_size);
    orbits2.initialize(domain_size);

    auto test_hook_dejavu = dejavu_hook(orbit1_test_hook);

    std::cout << "Running dejavu " << filename << "..." << std::endl;
    dejavu::solver d;
    d.set_error_bound(20); // in theory, due to the probabilistic nature of the algorithm, some tests could fail without
                          // there being a bug -- but this "almost ensures" that it will be a bug...
    d.set_print(false);
    d.automorphisms(g, colmap, &test_hook_dejavu);

    auto test_hook_sym_file = dejavu_hook(orbit2_test_hook);
    std::cout << "Reading symmetry file " << sym_filename << "..." << std::endl;
    read_symmetry_file(sym_filename, domain_size, &test_hook_sym_file);

    std::cout << "Comparing orbits..." << std::endl;
    EXPECT_TRUE(orbits1 == orbits2);

    dejavu::big_number file_grp_sz = read_grp_sz_file(grp_sz_filename);
    file_grp_sz.multiply(1);
    dejavu::big_number dejavu_grp_sz = d.get_automorphism_group_size();
    std::cout << "Comparing group size " << file_grp_sz << ":" << dejavu_grp_sz << "..." << std::endl;
    EXPECT_EQ(dejavu_grp_sz.exponent, file_grp_sz.exponent);
    EXPECT_NEAR(dejavu_grp_sz.mantissa, file_grp_sz.mantissa, 0.1);
}

TEST(graphs_test, graph_suite_comb) {
    std::string directory = TEST_RESOURCE_DIR;
    test_graph_orbit_check(directory + "X2.bliss");
    //test_graph(directory + "Ranreg1024.bliss");
    test_graph_orbit_check(directory + "pp-16-10.dimacs");
    test_graph_orbit_check(directory + "pp-25-100");
    test_graph_orbit_check(directory + "ag2-49.dimacs");
    test_graph_orbit_check(directory + "pg2-47.dimacs");
    test_graph_orbit_check(directory + "NCKL900K.bliss");
    test_graph_orbit_check(directory + "kef14.dimacs");
    test_graph_orbit_check(directory + "had-188.dimacs");
    test_graph_orbit_check(directory + "had-sw-40-6.dimacs");
    test_graph_orbit_check(directory + "f-lex-reg-20-1.dimacs");
    test_graph_orbit_check(directory + "f-lex-srg-14-1.dimacs");
    test_graph_orbit_check(directory + "latin-20.dimacs");
    test_graph_orbit_check(directory + "latin-30.dimacs");
    test_graph_orbit_check(directory + "latin-sw-23-4.dimacs");
    test_graph_orbit_check(directory + "cfi-172.dimacs");
    //test_graph_orbit_check(directory + "usr4_116-1.dimacs");
    //test_graph_orbit_check(directory + "CHH_cc7-7_1078-1.dimacs");
    test_graph_orbit_check(directory + "sts-79.dimacs");
    test_graph_orbit_check(directory + "sts-sw-49-8.dimacs");
    test_graph_orbit_check(directory + "triang-20.dimacs");
    test_graph_orbit_check(directory + "cmz-50.dimacs");
    test_graph_orbit_check(directory + "mz-aug-50.dimacs");
    test_graph_orbit_check(directory + "mz-aug2-50.dimacs");
    test_graph_orbit_check(directory + "grid-3-19.dimacs");
    test_graph_orbit_check(directory + "17cube.bliss");
    test_graph_orbit_check(directory + "tran_2880_12960.bliss");
    test_graph_orbit_check(directory + "paley-461.dimacs");
    test_graph_orbit_check(directory + "cfi-rigid-t2-0384-01-2");
    test_graph_orbit_check(directory + "sat_cfi_mult_4000_a.dmc");
    test_graph_orbit_check(directory + "group_128_164.dimacs");
    test_graph_orbit_check(directory + "flag6.bliss");
}

TEST(graphs_test, graph_suite_prep) {
    std::string directory = TEST_RESOURCE_DIR;
    test_graph_orbit_check(directory + "k-100.dimacs");
    test_graph_orbit_check(directory + "e-100.dimacs");
    test_graph_orbit_check(directory + "rantree-100.bliss");
    test_graph_orbit_check(directory + "highschool1-aigio.dimacs");
    test_graph_orbit_check(directory + "AS.bliss");
    test_graph_orbit_check(directory + "s4-4-3-5.bliss");
    test_graph_orbit_check(directory + "academictimetablesmall.dimacs");
    test_graph_orbit_check(directory + "heuristic_199.dimacs");
    test_graph_orbit_check(directory + "vlsat2_11_26.cnf.dimacs");
    test_graph_orbit_check(directory + "vlsat2_44_545.cnf.dimacs");
}

TEST(graphs_test, graph_suite_ref) {
    std::string directory = TEST_RESOURCE_DIR;
    test_graph_orbit_check(directory + "ransq_1000_a.bliss");
    test_graph_orbit_check(directory + "ran10_500_a.bliss");
    test_graph_orbit_check(directory + "ran2_200_a.bliss");
}