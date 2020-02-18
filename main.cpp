#include <iostream>
#include "parser.h"
#include "dejavu.h"
#include <assert.h>
#include "configuration.h"
#include <chrono>
#include <string>
#include <fstream>

typedef std::chrono::high_resolution_clock Clock;

class time_point;

configstruct config;
volatile int dejavu_kill_request = 0;

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

void bench_dejavu(sgraph *g, double* dejavu_solve_time) {
    // touch the graph (mitigate cache variance)
    Clock::time_point timer = Clock::now();
    dejavu_automorphisms(g, nullptr);
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

    for (int i = 1; i < argc; ++i) {
        std::string arg = std::string(argv[i]);
        std::transform(arg.begin(), arg.end(), arg.begin(), ::toupper);

        if (arg == "--FILE") {
            if (i + 1 < argc) {
                i++;
                filename = argv[i];
                entered_file = true;
            } else {
                std::cerr << "--file option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "--TIMEOUT") {
            if (i + 1 < argc) {
                i++;
                timeout = atoi(argv[i]);
            } else {
                std::cerr << "--timeout option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "--WRITE_AUTO") {
            config.CONFIG_WRITE_AUTOMORPHISMS = true;
        }

        if (arg == "--THREADS") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_THREADS_REFINEMENT_WORKERS = atoi(argv[i]);
            } else {
                std::cerr << "--threads option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "--PERMUTE") {
            permute_graph = true;
        }

        if (arg == "--COMPRESS") {
            config.CONFIG_PREPROCESS_COMPRESS      = true;
            config.CONFIG_PREPROCESS_EDGELIST_SORT = true;
        }

        if (arg == "--PERMUTE_SEED") {
            if (i + 1 < argc) {
                i++;
                permute_graph = true;
                seed = atoi(argv[i]);
            } else {
                std::cerr << "--permute_seed option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "--FORCE_SELECTOR") {
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

    if (!entered_file) {
        std::cerr << "--file not specified" << std::endl;
        return 1;
    }
    parser p;
    sgraph *g = new sgraph;
    std::cout << "Parsing " << filename << "..." << std::endl;
    p.parse_dimacs_file(filename, g);
    sgraph *_g = new sgraph;
    if(permute_graph) {
        std::cout << "Permuting graph..." << std::endl;
        bijection_temp<int> pr;
        bijection_temp<int>::random_bijection(&pr, g->v_size, seed);
        g->permute_graph(_g, &pr); // permute graph
        delete g;
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
    bench_dejavu(_g, &dejavu_solve_time);
    if(timeout > 0)
        killer.join();

    std::cout << "Solve time: " << dejavu_solve_time / 1000000.0 << "ms" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    return 0;
}


int main(int argc, char *argv[]) {
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "dejavu" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    return commandline_mode(argc, argv);
}
