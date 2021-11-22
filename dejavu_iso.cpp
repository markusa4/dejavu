#include "dejavu_iso.h"
#include <iostream>
#include "parser.h"
#include <assert.h>
#include "configuration.h"
#include <chrono>
#include <string>

typedef std::chrono::high_resolution_clock Clock;
bool finished = false;

configstruct config;
volatile int dejavu_kill_request = 0;
thread_local int numnodes;
thread_local int colorcost;

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

bool bench_dejavu_iso(sgraph *g1, sgraph *g2, double* dejavu_solve_time) {
    Clock::time_point timer = Clock::now();
    bool res = dejavu_isomorphic(g1, g2);
    *dejavu_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;
    return res;
}

int commandline_mode(int argc, char **argv) {
    std::string filename1 = "";
    std::string filename2 = "";
    bool entered_file1 = false;
    bool entered_file2 = false;

    int  timeout = -1;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    bool permute_graph = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = std::string(argv[i]);
        std::transform(arg.begin(), arg.end(), arg.begin(), ::toupper);
        std::replace(arg.begin(), arg.end(), '-', '_');

        if ((arg == "__FILE1") || (arg == "__F1")) {
            if (i + 1 < argc) {
                i++;
                filename1 = argv[i];
                entered_file1 = true;
            } else {
                std::cerr << "--file1 requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "__FILE2") || (arg == "__F2")) {
            if (i + 1 < argc) {
                i++;
                filename2 = argv[i];
                entered_file2 = true;
            } else {
                std::cerr << "--file2 requires one argument." << std::endl;
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
        } else if (arg == "__WRITE_ISO") {
            config.CONFIG_WRITE_AUTOMORPHISMS = true;
        } else if (arg == "__THREADS") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_THREADS_REFINEMENT_WORKERS = atoi(argv[i]);
            } else {
                std::cerr << "--threads option requires one argument." << std::endl;
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
        } else if (arg == "__FORCE_SELECTOR") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_IR_CELL_SELECTOR  = atoi(argv[i]);
                config.CONFIG_IR_FORCE_SELECTOR = true;
            } else {
                std::cerr << "--force_selector option requires one argument." << std::endl;
                return 1;
            }
        } else if (argv[i][0] != '-'){
            if(!entered_file1) {
                filename1 = argv[i];
                entered_file1 = true;
            } else if(!entered_file2) {
                filename2 = argv[i];
                entered_file2 = true;
            }  else {
                std::cout << "Extraneous file '" << argv[i] << "'. Only 2 files required." << std::endl;
                return 1;
            }
        } else {
            std::cout << "Invalid commandline option '" << argv[i] << "'." << std::endl;
            return 1;
        }
    }

    if (!entered_file1 || !entered_file2) {
        std::cerr << "--file1 or --file2 not specified" << std::endl;
        return 1;
    }

    if(!file_exists(filename1)) {
        std::cerr << "File '" << filename1 << "' does not exist." << std::endl;
        return 1;
    }

    if(!file_exists(filename2)) {
        std::cerr  << "File '" << filename2 << "' does not exist." << std::endl;
        return 1;
    }

    parser p;
    int* colmap1 = nullptr;
    int* colmap2 = nullptr;
    sgraph *g1 = new sgraph;
    std::cout << "Parsing " << filename1 << "..." << std::endl;
    p.parse_dimacs_file(filename1, g1, &colmap1);
    sgraph *g2 = new sgraph;
    std::cout << "Parsing " << filename2 << "..." << std::endl;
    p.parse_dimacs_file(filename2, g2, &colmap2);
    sgraph *_g1 = new sgraph;
    sgraph *_g2 = new sgraph;
    if(permute_graph) {
        std::cout << "Permuting graphs..." << std::endl;
        bijection<int> pr1;
        bijection<int>::random_bijection(&pr1, g1->v_size, seed);
        g1->permute_graph(_g1, &pr1); // permute graph
        int* rmap1 = pr1.extract_map();
        std::cout << "Colmap" << std::endl;
        if(colmap1 != nullptr)
            permute_colmap(&colmap1, g1->v_size, rmap1);
        delete g1;
        delete[] rmap1;

        bijection<int> pr2;
        bijection<int>::random_bijection(&pr2, g2->v_size, seed * 123);
        g2->permute_graph(_g2, &pr2); // permute graph
        int* rmap2 = pr2.extract_map();
        if(colmap2 != nullptr)
            permute_colmap(&colmap2, g2->v_size, rmap2);
        delete g2;
        delete[] rmap2;
    } else {
        _g1 = g1;
        _g2 = g2;
    }

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "dejavu-iso" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    double dejavu_solve_time;

    finished = false;
    std::thread killer;
    if(timeout > 0)
        killer = std::thread(kill_thread, &dejavu_kill_request, timeout);
    bool res = bench_dejavu_iso(_g1, _g2, &dejavu_solve_time);
    if(timeout > 0)
        killer.join();

    if(res) {
        std::cout << "ISOMORPHIC" << std::endl;
    } else {
        std::cout << "NON_ISOMORPHIC" << std::endl;
    }

    std::cout << "Solve time: " << dejavu_solve_time / 1000000.0 << "ms" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    return 0;
}

int main(int argc, char *argv[]) {
    return commandline_mode(argc, argv);
}
