// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#include "gtest/gtest.h"
#include "../dejavu.h"
#include "../parser.h"
#include "helper_functions_test.h"
#include <string>
#include <fstream>
#include <dirent.h>
#include <algorithm>
#include <vector>


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
    while (std::getline(infile, line)) {
        auto_ws.reset();
        std::string vert_str;
        int first_in_cycle = -1;
        int prev_in_cycle  = -1;
        for(size_t i = 0; i < line.size(); ++i) {
            if (line[i] == '(') {
                vert_str = "";
                first_in_cycle = -1;
                prev_in_cycle  = -1;
            } else if ((line[i] != ' ') && (line[i] != ')')) {
                vert_str += line[i];
            } else {
                if(vert_str != "") {
                    const int vert = std::stoi(vert_str);
                    if(first_in_cycle == -1) first_in_cycle = vert;
                    else {
                        const int from = prev_in_cycle;
                        const int to   = vert;
                        auto_ws.write_single_map(from-1, to-1);
                    }
                    prev_in_cycle = vert;
                    vert_str = "";
                }
                if (line[i] == ')') {
                    auto_ws.write_single_map(prev_in_cycle-1, first_in_cycle-1);
                }
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
    result.set(0, 0);
    double divide_by = 1;
    while (std::getline(infile, line)){
        for(size_t i = 0; i < line.size(); ++i) {
            if(line[i] >= '0' && line[i] <= '9') {
                if(i > 0) result.multiply(10);
                const int dec_place = line[i] - '0';
                if(divide_by < 1000000000) {
                    result.mantissa += 1.0*dec_place / divide_by;
                    divide_by *= 10;
                }
            }
        }
        break;
    }

    return result;
}

bool test_graph_orbit_check(std::string filename) {
    std::string sym_filename = filename + ".sym";
    std::string grp_sz_filename = filename + ".grpsz";
    dejavu::sgraph *g = new dejavu::sgraph();
    int* colmap = nullptr;
    std::cout << console_blue << filename << console_neutral << std::flush;
    parse_dimacs(filename, g, &colmap);
    EXPECT_TRUE(g->v_size > 0);

    dgtest_graph.copy_graph(g);
    const int domain_size = g->v_size;
    orbits1.initialize(domain_size);
    orbits2.initialize(domain_size);

    auto test_hook_dejavu = dejavu_hook(orbit1_test_hook);

    //std::cout << "Running dejavu " << filename << "..." << std::endl;
    dejavu::solver d;
    d.set_error_bound(20); // in theory, due to the probabilistic nature of the algorithm, some tests could fail without
                          // there being a bug -- but this "almost ensures" that it will be a bug...
    d.set_print(false);
    d.automorphisms(g, colmap, &test_hook_dejavu);

    auto test_hook_sym_file = dejavu_hook(orbit2_test_hook);
    //std::cout << "Reading symmetry file " << sym_filename << "..." << std::endl;
    read_symmetry_file(sym_filename, domain_size, &test_hook_sym_file);

    //std::cout << "Comparing orbits..." << std::endl;
    EXPECT_TRUE(orbits1 == orbits2);
    const bool check1 = orbits1 == orbits2;

    dejavu::big_number file_grp_sz = read_grp_sz_file(grp_sz_filename);
    file_grp_sz.multiply(1);
    dejavu::big_number dejavu_grp_sz = d.get_automorphism_group_size();
    //std::cout << "Comparing group size " << file_grp_sz << ":" << dejavu_grp_sz << "..." << std::endl;
    EXPECT_EQ(dejavu_grp_sz.exponent, file_grp_sz.exponent);
    const bool check2 = dejavu_grp_sz.exponent == file_grp_sz.exponent;
    EXPECT_NEAR(dejavu_grp_sz.mantissa, file_grp_sz.mantissa, 0.1);
    if(check1 && check2)
        std::cout << console_green << " OK " << console_neutral << std::endl;
    else
        std::cout << console_red << "FAIL" << console_neutral << std::endl;

    return (check1 && check2);
}

TEST(graphs_test, graph_suite) {
    std::string directory = TEST_RESOURCE_DIR;

    std::vector<std::string> graphs;

    DIR* dir = opendir(directory.c_str());
    ASSERT_NE(dir, nullptr) << "Could not open directory: " << directory;

    while (dirent* entry = readdir(dir)) {
        std::string name = entry->d_name;

        if (name.size() < 7 || name.substr(name.size() - 7) != ".dimacs") {
            continue;
        }

        std::string graph = directory + name;

        if (std::ifstream(graph + ".sym").good()
            && std::ifstream(graph + ".grpsz").good()) {
            graphs.push_back(graph);
        }
    }

    closedir(dir);

    std::sort(graphs.begin(), graphs.end());

    ASSERT_FALSE(graphs.empty()) << "No graph/reference triples found in " << directory;

    size_t c = 1;
    int fails = 0;
    for (const std::string& graph : graphs) {
        SCOPED_TRACE(graph);
        if(fails > 0) std::cout << console_red;
        std::cout << "[" << c << "/" << graphs.size() << "] "; 
        std::cout << console_neutral;
        fails += !test_graph_orbit_check(graph);
        ++c;
    }
}
