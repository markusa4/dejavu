#include "dejavu_api.h"
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

bool bench_dejavu_api(sgraph *g1, int max_length, int num, double* dejavu_solve_time) {
    // touch the graph (mitigate cache variance)
    Clock::time_point timer = Clock::now();
    std::set<std::tuple<int*, int, int*, long>> paths;
    //random_paths(g1, nullptr, max_length, num, &paths);
    dejavu_automorphisms(g1, nullptr, nullptr);
    *dejavu_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;
    return true;
}

int commandline_mode(int argc, char **argv) {
    std::string filename1 = "";
    bool entered_file1 = false;

    int  timeout = -1;
    int max_length = 1;
    int num = 1;
    bool comp_dejavu = true;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    bool permute_graph = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = std::string(argv[i]);
        std::transform(arg.begin(), arg.end(), arg.begin(), ::toupper);
        std::replace(arg.begin(), arg.end(), '-', '_');

        if ((arg == "__FILE") || (arg == "__F")) {
            if (i + 1 < argc) {
                i++;
                filename1 = argv[i];
                entered_file1 = true;
            } else {
                std::cerr << "--file requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "__TIMEOUT") {
            if (i + 1 < argc) {
                i++;
                timeout = atoi(argv[i]);
            } else {
                std::cerr << "--timeout option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "__WRITE_AUTO") {
            config.CONFIG_WRITE_AUTOMORPHISMS = true;
        }

        if (arg == "__THREADS") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_THREADS_REFINEMENT_WORKERS = atoi(argv[i]);
            } else {
                std::cerr << "--threads option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "__MAX_LENGTH") {
            if (i + 1 < argc) {
                i++;
                max_length = atoi(argv[i]);
            } else {
                std::cerr << "--max-length option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "__NUM") {
            if (i + 1 < argc) {
                i++;
                num = atoi(argv[i]);
            } else {
                std::cerr << "--num option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "__ERR") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_RAND_ABORT = atoi(argv[i]);
            } else {
                std::cerr << "--err option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "__LEAF_LIMIT") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_IR_LEAF_STORE_LIMIT = atoi(argv[i]);
            } else {
                std::cerr << "--leaf-limit option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "__PERMUTE") {
            permute_graph = true;
        }

        if (arg == "__KDEVIATION") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_IR_EXPAND_DEVIATION = atoi(argv[i]);
            } else {
                std::cerr << "--kdeviation option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "__NO_IDLESKIP") {
            config.CONFIG_IR_IDLE_SKIP = false;
        }

        if (arg == "__COMPRESS") {
            config.CONFIG_PREPROCESS_COMPRESS      = true;
            config.CONFIG_PREPROCESS_EDGELIST_SORT = true;
        }

        if (arg == "__PERMUTE_SEED") {
            if (i + 1 < argc) {
                i++;
                permute_graph = true;
                seed = atoi(argv[i]);
            } else {
                std::cerr << "--permute_seed option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "__FORCE_SELECTOR") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_IR_CELL_SELECTOR  = atoi(argv[i]);
                config.CONFIG_IR_FORCE_SELECTOR = true;
            } else {
                std::cerr << "--force_selector option requires one argument." << std::endl;
                return 1;
            }
        }
    }

    if (!entered_file1) {
        std::cerr << "--file not specified" << std::endl;
        return 1;
    }

    if(!file_exists(filename1)) {
        std::cerr << "File '" << filename1 << "' does not exist." << std::endl;
        return 1;
    }

    parser p;
    int* colmap1 = nullptr;
    int* colmap2 = nullptr;
    sgraph *g1 = new sgraph;
    std::cout << "Parsing " << filename1 << "..." << std::endl;
    p.parse_dimacs_file(filename1, g1, &colmap1);
    sgraph *_g1;
    _g1 = g1;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "dejavu-api benchmark" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    double dejavu_solve_time;

    finished = false;
    std::thread killer;
    if(timeout > 0)
        killer = std::thread(kill_thread, &dejavu_kill_request, timeout);
    for(int i = 0; i < num; ++i ) {
        std::cout << "Run " << i << ":" << std::endl;
        bool res = bench_dejavu_api(_g1, max_length, num, &dejavu_solve_time);
    }
    if(timeout > 0)
        killer.join();

    std::cout << "Solve time: " << dejavu_solve_time / 1000000.0 << "ms" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    getchar();
    delete g1;
    if(colmap1)
        delete[] colmap1;
    return 0;
}


int main(int argc, char *argv[]) {
    return commandline_mode(argc, argv);
}
