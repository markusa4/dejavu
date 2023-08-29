// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#include <iostream>
#include "parser.h"
#include <thread>
#include "dejavu.h"
#include <chrono>
#include <string>

typedef std::chrono::high_resolution_clock Clock;

dejavu::ir::refinement test_r;
dejavu::sgraph dej_test_graph;
int*   dej_test_col;

volatile int dejavu_kill_request = 0;

bool finished = false;

/*void kill_thread(volatile int* kill_switch, int timeout) {
    Clock::time_point start = Clock::now();
    while(!finished) {
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        if(((std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count()) / 1000000.0) > timeout * 1000.0) {
            std::cout << "Killing" << std::endl;
            *kill_switch = 1;
        }
    }
}*/

void empty_hook(int, const int*, int, const int *) {}

void bench_dejavu(dejavu::sgraph* g, int* colmap, double* dejavu_solve_time) {
    // touch the graph (mitigate cache variance)
    volatile int acc = 0;
    for(int i = 0; i < g->v_size; ++i) {
        acc += g->v[i] + g->d[i];
    }
    for(int i = 0; i < g->e_size; ++i) {
        acc += g->e[i];
    }

    Clock::time_point timer = Clock::now();
    auto empty_hook_func = dejavu_hook(empty_hook);
    //dejavu_automorphisms(g, colmap, &empty_hook_func);
    dejavu::dejavu2 d;
    bool del = false;
    if(colmap == nullptr) {
        colmap = (int*) calloc(g->v_size, sizeof(int));
        del = true;
    }

#ifndef NDEBUG
    auto test_hook_func = dejavu_hook(dejavu::test_hook);
    dej_test_graph.copy_graph(g);
    dej_test_col = (int*) calloc(g->v_size, sizeof(int));
    memcpy(dej_test_col, colmap, sizeof(int) * g->v_size);
    d.automorphisms(g, colmap, &test_hook_func);
    free(dej_test_col);
#else
    d.automorphisms(g, colmap, &empty_hook_func);
#endif
    *dejavu_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "deterministic: " << d.s_deterministic_termination << ", error: 1/2^"<< d.h_error_bound << std::endl;
    std::cout << std::setprecision(4) << "#symmetries: " << d.s_grp_sz << std::endl;
     finished = true;
    if(del) free(colmap);
}

int commandline_mode(int argc, char **argv) {
    std::string filename = "";
    bool entered_file = false;
    int  timeout = -1;
    bool permute_graph = false;
    bool permute_graph_have_seed  = false;
    int  permute_graph_given_seed = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = std::string(argv[i]);
        std::transform(arg.begin(), arg.end(), arg.begin(), ::toupper);
        std::replace(arg.begin(), arg.end(), '-', '_');

        if (arg == "__FILE") {
            if (i + 1 < argc) {
                i++;
                filename = argv[i];
                entered_file = true;
            } else {
                std::cerr << "--file option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "__TIMEOUT") {
            if (i + 1 < argc) {
                i++;
                timeout = atoi(argv[i]);
            } else {
                std::cerr << "--timeout option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "__ERR") {
            if (i + 1 < argc) {
                i++;
                //config.CONFIG_RAND_ABORT = atoi(argv[i]);
            } else {
                std::cerr << "--err option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "__LEAF_LIMIT") {
            if (i + 1 < argc) {
                i++;
                //config.CONFIG_IR_LEAF_STORE_LIMIT = atoi(argv[i]);
            } else {
                std::cerr << "--leaf-limit option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "__WRITE_AUTO") {
            //config.CONFIG_WRITE_AUTOMORPHISMS = true;
        } else if (arg == "__WRITE_AUTO_GAP") {
            //config.CONFIG_WRITE_AUTOMORPHISMS_GAP = true;
        }  else if (arg == "__THREADS") {
            if (i + 1 < argc) {
                i++;
                //config.CONFIG_THREADS_REFINEMENT_WORKERS = atoi(argv[i]);
            } else {
                std::cerr << "--threads option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "__PERMUTE") {
            permute_graph = true;
        }  else if (arg == "__PERMUTE_SEED") {
            if (i + 1 < argc) {
                i++;
                permute_graph_have_seed  = true;
                permute_graph_given_seed = atoi(argv[i]);
            } else {
                std::cerr << "--permute_seed option requires one argument." << std::endl;
                return 1;
            }
        } else if (argv[i][0] != '-') {
            if(!entered_file) {
                filename = argv[i];
                entered_file = true;
            } else {
                std::cout << "Extraneous file '" << argv[i] << "'. Only 1 file required." << std::endl;
                return 1;
            }
        } else {
            std::cout << "Invalid commandline option '" << argv[i] << "'." << std::endl;
            return 1;
        }
    }

    if (!entered_file) {
        std::cerr << "--file not specified" << std::endl;
        return 1;
    }

    if(!file_exists(filename)) {
        std::cerr << "File '" << filename << "' does not exist." << std::endl;
        return 1;
    }

    dejavu::sgraph *g = new dejavu::sgraph();
    std::cout << "Parsing " << filename << "..." << std::endl;
    int* colmap = nullptr;

    int permute_seed = 0;
    if(permute_graph) {
        permute_seed = permute_graph_given_seed;
        if(!permute_graph_have_seed) {
            std::random_device r;
            permute_seed = static_cast<int>(r());
        }
        std::cout << "permutation(" << permute_seed << ")" << std::endl;
    }

    parse_dimacs_file_fast(filename, g, &colmap, permute_seed);

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "dejavu 2.0beta" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    double dejavu_solve_time;

    finished = false;
    std::thread killer;
    //if(timeout > 0)
    //    killer = std::thread(kill_thread, &dejavu_kill_request, timeout);
    bench_dejavu(g, colmap, &dejavu_solve_time);
    if(timeout > 0)
        killer.join();

    std::cout << "Solve time: " << dejavu_solve_time / 1000000.0 << "ms" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    return 0;
}


int main(int argc, char *argv[]) {
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "parser" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    return commandline_mode(argc, argv);
}
