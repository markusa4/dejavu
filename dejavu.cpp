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

bool finished = false;

void empty_hook(int, const int*, int, const int *) {}

dejavu::big_number run_dejavu(dejavu::sgraph* g, int* colmap, double* dejavu_solve_time, bool print = true,
                              bool true_random = false, bool randomize_seed = false, dejavu_hook* hook = nullptr) {
    // touch the graph (mitigate cache variance)
    volatile int acc = 0;
    for(int i = 0; i < g->v_size; ++i) acc += g->v[i] + g->d[i];
    for(int i = 0; i < g->e_size; ++i) acc += g->e[i];

    bool del = false;
    if(colmap == nullptr) {
        colmap = (int*) calloc(g->v_size, sizeof(int));
        del = true;
    }

    Clock::time_point timer = Clock::now();
    dejavu::dejavu2 d;
    d.set_print(print);
    if(randomize_seed) d.randomize_seed();
    d.set_true_random(true_random);

    d.automorphisms(g, colmap, hook);

/*#ifndef NDEBUG
    auto test_hook_func = dejavu_hook(dejavu::test_hook);
    dej_test_graph.copy_graph(g);
    dej_test_col = (int*) calloc(g->v_size, sizeof(int));
    memcpy(dej_test_col, colmap, sizeof(int) * g->v_size);
    d.automorphisms(g, colmap, &test_hook_func);
    free(dej_test_col);
#else
    d.automorphisms(g, colmap, &empty_hook_func);
#endif*/
    *dejavu_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    dejavu::big_number grp_sz = d.get_automorphism_group_size();
    if(print) std::cout << "------------------------------------------------------------------" << std::endl;
    if(print) std::cout << std::setprecision(4) << "symmetries=" << grp_sz
                        << ", deterministic=" << (d.get_deterministic_termination()?"true":"false")
                        << ", error=1/2^"<< d.get_error_bound() << "," << std::endl;
     finished = true;
    if(del) free(colmap);
    return grp_sz;
}

int commandline_mode(int argc, char **argv) {
    std::string filename = "";
    bool entered_file = false;
    int  timeout = -1;
    bool permute_graph = false;
    bool permute_graph_have_seed  = false;
    int  permute_graph_given_seed = 0;
    bool print = true;

    bool true_random = false;
    bool true_random_seed = false;

    bool benchmark_mode = false;
    bool write_grp_sz = false;
    bool write_auto_stdout = false;
    bool        write_auto_file      = false;
    std::string write_auto_file_name;

    for (int i = 1; i < argc; ++i) {
        std::string arg = std::string(argv[i]);
        std::transform(arg.begin(), arg.end(), arg.begin(), ::toupper);
        std::replace(arg.begin(), arg.end(), '-', '_');

        if (arg == "__HELP" || arg == "_H") {
            std::cout << "Usage: dejavu [file] [options]" << std::endl;
            std::cout << "Computes the automorphism group of the graph described in the given FILE." << std::endl;
            std::cout << "The FILE must describe a graph using the DIMACS format." << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "    "  << std::left << std::setw(20) <<
            "--err [n]" << std::setw(16) <<
            "Sets the error to be bounded by 1/2^N" << std::endl;
            std::cout << "    " << std::left << std::setw(20) <<
            "--silent" << std::setw(16) <<
            "Does not print progress of the solver" << std::endl;
            std::cout << "    " << std::left << std::setw(20) <<
            "--gens" << std::setw(16) <<
            "Prints found generators line-by-line to console" << std::endl;
            std::cout << "    " << std::left << std::setw(20) <<
            "--gens-file [f]" << std::setw(16) <<
           "Writes found generators line-by-line to file F" << std::endl;
            std::cout << "    " << std::left << std::setw(20) <<
            "--grp-sz" << std::setw(16) <<
            "Prints group size to console (even if --silent)" << std::endl;
            std::cout << "    "  << std::left << std::setw(20) <<
            "--pseudo-random" << std::setw(16) <<
            "Uses pseudo random numbers (default)" << std::endl;
            std::cout << "    " << std::left << std::setw(20) <<
            "--true-random" << std::setw(16) <<
            "Uses random device of OS" << std::endl;
            std::cout << "    "  << std::left << std::setw(20) <<
            "--true-random-seed" << std::setw(16) <<
            "Seeds pseudo random with random device of OS" << std::endl;
            std::cout << "    "  << std::left << std::setw(20) <<
            "--permute" << std::setw(16) <<
            "Randomly permutes the given graph" << std::endl;
            std::cout << "    "  << std::left << std::setw(20) <<
            "--permute-seed [n]" << std::setw(16) <<
            "Seed for the previous option with N" << std::endl;
            return 0;
        } else if (arg == "__VERSION" || arg == "_V") {
            std::cout << DEJAVU_VERSION_MAJOR << "." << DEJAVU_VERSION_MINOR <<
                        (DEJAVU_VERSION_IS_BETA?"beta":"") << std::endl;
            return 0;
        } else if (arg == "__FILE") {
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
        } else if (arg == "__GRP_SZ") {
            write_grp_sz = true;
        } else if (arg == "__GENS") {
            write_auto_stdout = true;
        }  else if (arg == "__GENS_FILE") {
            if (i + 1 < argc) {
                i++;
                write_auto_file = true;
                write_auto_file_name = argv[i];
            } else {
                std::cerr << "--write-gens-file option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "__TRUE_RANDOM") {
            if(true_random_seed) {
                std::cerr << "--true-random and --true-random-seed can not be activated at the same time:" <<
                          "--true-random-seed seeds a pseudo random number generator with a random number" << std::endl;
                return 1;
            }
            true_random = true;
        } else if (arg == "__PSEUDO_RANDOM") {
            true_random = false;
        }  else if (arg == "__TRUE_RANDOM_SEED") {
            if(true_random) {
                std::cerr << "--true-random and --true-random-seed can not be activated at the same time:" <<
                "--true-random-seed seeds a pseudo random number generator with a random number" << std::endl;
                return 1;
            }
            true_random_seed = true;
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
        }  else if (arg == "__SILENT") {
            print = false;
        }  else if (argv[i][0] != '-') {
            if(!entered_file) {
                filename = argv[i];
                entered_file = true;
            } else {
                std::cerr << "Extraneous file '" << argv[i] << "'. Only 1 file required." << std::endl;
                return 1;
            }
        } else {
            std::cerr << "Invalid commandline option '" << argv[i] << "'." << std::endl;
            return 1;
        }
    }

    if (!entered_file) {
        std::cerr << "no file was specified, usage: dejavu [file] [options]" << std::endl;
        return 1;
    }

    if(!file_exists(filename)) {
        std::cerr << "File '" << filename << "' does not exist." << std::endl;
        return 1;
    }

    //if(print) std::cout << "------------------------------------------------------------------" << std::endl;
    if(print) std::cout << "dejavu version=" << DEJAVU_VERSION_MAJOR << "." << DEJAVU_VERSION_MINOR <<
                        (DEJAVU_VERSION_IS_BETA?"beta":"") << std::endl;
    if(print) std::cout << "------------------------------------------------------------------" << std::endl;

    dejavu::sgraph *g = new dejavu::sgraph();
    if(print) std::cout << "parsing '" << filename << "'" << std::endl;
    int* colmap = nullptr;

    int permute_seed = 0;
    if(permute_graph) {
        permute_seed = permute_graph_given_seed;
        if(!permute_graph_have_seed) {
            std::random_device r;
            permute_seed = static_cast<int>(r());
        }
        if(print) std::cout << (true_random?"true_random=true, ":"") << (true_random_seed?"true_random_seed=true":"");
        if(print) std::cout << "permutation_seed=" << permute_seed << ", ";
    }
    parse_dimacs_file_fast(filename, g, &colmap, !print, permute_seed);
    if(print) std::cout << ", n=" << g->v_size << ", " << "m=" << g->e_size/2 << std::endl << std::endl;

    // manage hooks
    auto empty_hook_func = dejavu_hook(empty_hook);
    dejavu::multi_hook hooks;
    std::ofstream output_file;
    dejavu::ostream_hook file_hook(output_file);
    dejavu::ostream_hook cout_hook(std::cout);
    dejavu_hook* hook;

    // write automorphism to file or cout
    if(write_auto_stdout) hooks.add_hook(cout_hook.get_hook());
    if(write_auto_file) {
        output_file.open(write_auto_file_name);
        hooks.add_hook(file_hook.get_hook());
    }

    // debug hook
#ifndef NDEBUG
    auto test_hook_func = dejavu_hook(dejavu::test_hook);
    dej_test_graph.copy_graph(g);
    hooks.add_hook(&test_hook_func);
#endif

    // use multi-hook or empty hook
    if(hooks.size() == 0) hook = &empty_hook_func; // using empty_hook_func for fair benchmarks, 'nullptr' is faster
    else hook = hooks.get_hook();

    // now run the solver with the given options...
    double dejavu_solve_time;
    finished = false;
    auto grp_sz = run_dejavu(g, colmap, &dejavu_solve_time, print, true_random, true_random_seed, hook);
    if(print) std::cout << "solve_time=" << dejavu_solve_time / 1000000.0 << "ms" << std::endl;
    if(!print && write_grp_sz) std::cout << grp_sz << std::endl;
    return 0;
}

int main(int argc, char *argv[]) {
    return commandline_mode(argc, argv);
}
