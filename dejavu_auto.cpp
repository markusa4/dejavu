#include <iostream>
#include "parser.h"
#include "dejavu_auto.h"
#include <assert.h>
#include "configuration.h"
#include <chrono>
#include <string>
#include <fstream>

typedef std::chrono::high_resolution_clock Clock;

configstruct config;
volatile int dejavu_kill_request = 0;
thread_local int numnodes;
thread_local int colorcost;

bool finished = false;

void kill_thread(volatile int* kill_switch, int timeout) {
    Clock::time_point start = Clock::now();
    while(!finished) {
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        if(((std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count()) / 1000000.0) > timeout * 1000.0) {
            std::cout << "Killing" << std::endl;
            *kill_switch = 1;
        }
    }
}

void empty_hook(int n, const int * p, int support, const int *) {
    /*bijection<int> test_p;
    test_p.read_from_array(p, n);
    const bool test_auto = test_R.certify_automorphism(&test_graph, &test_p); // TODO: sparse automorphism certification?
    assert(test_auto);*/
    //std::cout << "automorphism support " << support << std::endl;
    return;
}

void bench_dejavu(sgraph* g, int* colmap, double* dejavu_solve_time) {
    test_graph.copy_graph(g);

    // touch the graph (mitigate cache variance)
    Clock::time_point timer = Clock::now();
    dejavu_automorphisms(g, colmap, empty_hook);
    //dejavu d;
    //d.automorphisms(g, nullptr);
    *dejavu_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;
}

int commandline_mode(int argc, char **argv) {
    std::string filename = "";
    bool entered_file = false;
    int  timeout = -1;
    bool comp_dejavu = true;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    bool permute_graph = false;

    config.CONFIG_IR_WRITE_GROUPORDER = true;
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
                config.CONFIG_RAND_ABORT = atoi(argv[i]);
            } else {
                std::cerr << "--err option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "__LEAF_LIMIT") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_IR_LEAF_STORE_LIMIT = atoi(argv[i]);
            } else {
                std::cerr << "--leaf-limit option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "__WRITE_AUTO") {
            config.CONFIG_WRITE_AUTOMORPHISMS = true;
        } else if (arg == "__WRITE_AUTO_GAP") {
            config.CONFIG_WRITE_AUTOMORPHISMS_GAP = true;
        }  else if (arg == "__THREADS") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_THREADS_REFINEMENT_WORKERS = atoi(argv[i]);
            } else {
                std::cerr << "--threads option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "__PERMUTE") {
            permute_graph = true;
        } else if (arg == "__KDEVIATION") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_IR_EXPAND_DEVIATION = atoi(argv[i]);
            } else {
                std::cerr << "--kdeviation option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "__NO_IDLESKIP") {
            config.CONFIG_IR_IDLE_SKIP = false;
        }  else if (arg == "__EARLY_IR") {
            config.CONFIG_IR_INDIVIDUALIZE_EARLY = true;
        } else if (arg == "__COMPRESS") {
            config.CONFIG_PREPROCESS_COMPRESS      = true;
            config.CONFIG_PREPROCESS_EDGELIST_SORT = true;
        } else if (arg == "__PERMUTE_SEED") {
            if (i + 1 < argc) {
                i++;
                permute_graph = true;
                seed = atoi(argv[i]);
            } else {
                std::cerr << "--permute_seed option requires one argument." << std::endl;
                return 1;
            }
        }  else if (arg == "__ONLY_COLOR_REF_INVARIANT") {
                config.CONFIG_ONLY_COLOR_REF_INVARIANT = true;
        }  else if (arg == "__FORCE_SELECTOR") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_IR_CELL_SELECTOR  = atoi(argv[i]);
                config.CONFIG_IR_FORCE_SELECTOR = true;
            } else {
                std::cerr << "--force_selector option requires one argument." << std::endl;
                return 1;
            }
        }   else if (argv[i][0] != '-'){
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

    parser p;
    sgraph *g = new sgraph;
    std::cout << "Parsing " << filename << "..." << std::endl;
    int* colmap = nullptr;
    p.parse_dimacs_file(filename, g, &colmap);
    sgraph *_g = new sgraph;
    if(permute_graph) {
        std::cout << "Permuting graph..." << std::endl;
        bijection pr;
        bijection::random_bijection(&pr, g->v_size, seed);
        g->permute_graph(_g, &pr); // permute graph
        int* rmap = pr.extract_map();
        if(colmap != nullptr)
            permute_colmap(&colmap, g->v_size, rmap);
        delete g;
        delete[] rmap;
    } else {
        _g = g;
    }

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "dejavu" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    double dejavu_solve_time;

    finished = false;
    std::thread killer;
    if(timeout > 0)
        killer = std::thread(kill_thread, &dejavu_kill_request, timeout);
    bench_dejavu(_g, colmap, &dejavu_solve_time);
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
