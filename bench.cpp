#include <iostream>
#include "parser.h"
#include "dejavu_auto.h"
#include <assert.h>

extern "C" {
#include "saucy/saucy.h"
#include "nauty/traces.h"
}

#include "nauty/naugroup.h"
#include "configuration.h"
#include <chrono>
#include <string>
#include <fstream>

typedef std::chrono::high_resolution_clock Clock;

class time_point;

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

struct colorComparator {
    colorComparator(const int* colmap) : colmap(colmap) {}
    const int* colmap;

    bool operator()(const int & v1, const int & v2) {
        return (colmap[v1] < colmap[v2]);
    }
};

void make_lab_ptn_from_colmap(int* lab, int* ptn, int* colmap, int colmap_sz) {
    for(int i = 0; i < colmap_sz; ++i) {
        lab[i] = i;
        ptn[i] = 1;
    }
    int last_new_cell = 0;
    std::sort(lab, lab + colmap_sz, colorComparator(colmap));
    for(int i = 0; i < colmap_sz; i++) {
        if (i + 1 == colmap_sz) {
            ptn[last_new_cell] = (i - last_new_cell) > 0;
            ptn[i] = 0;
            break;
        }
        if (colmap[lab[i]] < colmap[lab[i + 1]]) {
            ptn[i] = 0;
            ptn[last_new_cell] = (i - last_new_cell) > 0;
            last_new_cell = i + 1;
            continue;
        }
    }
}

void make_small_colmap(int* smallcolmap, int* colmap, int colmap_sz) {
    int* lab = new int[colmap_sz];
    for(int i = 0; i < colmap_sz; ++i) {
        lab[i] = i;
    }

    int last_new_cell = 0;
    std::sort(lab, lab + colmap_sz, colorComparator(colmap));

    int col_small = -1;
    int last_col_orig  = -1;

    for(int i = 0; i < colmap_sz; i++) {
        if (colmap[lab[i]] != last_col_orig) {
            last_col_orig = colmap[lab[i]];
            ++col_small;
        }
        smallcolmap[lab[i]] = col_small;
    }

    delete[] lab;
}

void bench_nauty(sgraph_t<int, int, int> *g, int* colmap, double* nauty_solve_time) {
    //TracesStats stats;
    statsblk stats;
    sparsegraph sg;

    DYNALLSTAT(int, lab, lab_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, orbits, orbits_sz);
    //static DEFAULTOPTIONS_SPARSEGRAPH(options);

    SG_INIT(sg);
    int m = SETWORDSNEEDED(g->v_size);
    nauty_check(WORDSIZE, m, g->v_size, NAUTYVERSIONID);
    SG_ALLOC(sg, g->v_size, g->e_size, "malloc");
    sg.nv = g->v_size;
    sg.nde = g->e_size;
    //static DEFAULTOPTIONS_TRACES(options);
    static DEFAULTOPTIONS_SPARSEGRAPH(options);
    options.schreier = true;
    //schreier_fails(10);
    options.defaultptn = false;
    //options.writeautoms = true;

    DYNALLOC1(int, lab, lab_sz, sg.nv, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, sg.nv, "malloc");
    DYNALLOC1(int, orbits, orbits_sz, sg.nv, "malloc");

    for (int i = 0; i < g->v_size; ++i) {
        lab[i] = i;
        ptn[i] = 1;
        sg.v[i] = g->v[i];
        sg.d[i] = g->d[i];
    }
    ptn[g->v_size - 1] = 0;
    if(colmap != nullptr) {
        make_lab_ptn_from_colmap(lab, ptn, colmap, g->v_size);
    }

    for (int i = 0; i < g->e_size; ++i) {
        sg.e[i] = g->e[i];
    }

    // touch the graph (mitigate cache variance
    int acc = 0;
    for(int i = 0; i < g->v_size; ++i) {
        acc += sg.v[i] + sg.d[i];
    }
    for(int i = 0; i < g->e_size; ++i) {
        acc += sg.e[i];
    }
    Clock::time_point timer = Clock::now();

    sparsenauty(&sg, lab, ptn, orbits, &options, &stats, NULL);
    std::cout << "Group size: ";
    writegroupsize(stdout, stats.grpsize1, stats.grpsize2);
    std::cout << std::endl;
    //Traces(&sg,lab,ptn,orbits,&options,&stats,NULL);
    *nauty_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;
}

void empty_consumer(int n, const int * p, int support, const int *) {
    /*bijection<int> test_p;
    test_p.read_from_array(p, n);
    const bool test_auto = test_R.certify_automorphism(&test_graph, &test_p); // TODO: sparse automorphism certification?
    assert(test_auto);*/
    //std::cout << "automorphism support " << support << std::endl;
    return;
}

int bench_dejavu(sgraph_t<int, int, int> *g, int* colmap, double* dejavu_solve_time) {
    // touch the graph (mitigate cache variance)
    int acc = 0;
    for(int i = 0; i < g->v_size; ++i) {
        acc += g->v[i] + g->d[i];
    }
    for(int i = 0; i < g->e_size; ++i) {
        acc += g->e[i];
    }

    if(colmap != nullptr) {
        for(int i = 0; i < g->v_size; ++i) {
            acc += colmap[i];
        }
    }

    Clock::time_point timer = Clock::now();
    config.CONFIG_IR_WRITE_GROUPORDER = true;
    dejavu_automorphisms(g, colmap, empty_consumer);
    *dejavu_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;

    return acc;
}

int test_cycle = 0;

static int saucyConsume(int n, const int *perm, int nsupp, int * support, void *) {
    int i, k;
    for (i = 0; i < nsupp; ++i) {
        k = support[i];
        test_cycle += perm[k];
    }
    return true;
}


void bench_saucy(sgraph_t<int, int, int> *g, int* colmap, double* saucy_solve_time) {
    saucy_graph sg;
    sg.n = g->v_size;
    sg.e = g->e_size;
    sg.edg = new int[g->e_size];
    sg.adj = new int[sg.n + 1];
    //sg.colors = new int[g->v_size];

    int epos = 0;

    for(int i = 0; i < g->v_size; ++i) {
        const int npt = g->v[i];
        const int nd  = g->d[i];
        sg.adj[i] = epos;
        for(int j = 0; j < nd; ++j) {
            const int adj_n = g->e[npt + j];
            sg.edg[epos] = adj_n;
            ++epos;
        }
    }

    sg.adj[g->v_size] = epos;
    assert(epos == g->e_size);

    if(colmap == nullptr) {
        colmap = new int[g->v_size];
        for(int i = 0; i < g->v_size; ++i) {
            colmap[i] = 0;
        }
    }

    int* small_colmap = new int[g->v_size];
    make_small_colmap(small_colmap, colmap, g->v_size);

    Clock::time_point timer = Clock::now();
    // automorphisms
    struct saucy *s = saucy_alloc(sg.n); // 100000
    struct saucy_stats stats;
    saucy_search(s, &sg, 0, colmap, saucyConsume, 0, &stats);
    saucy_free(s);
    *saucy_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Group size: ";
    std::cout << stats.grpsize_base << "*10^" << stats.grpsize_exp << std::endl;
    finished = true;
}

void test_auto(int i, int* x1, int x2) {
    return;
}

void bench_traces(sgraph_t<int, int, int> *g, int* colmap, double* traces_solve_time) {
    //sleep(1);
    TracesStats stats;
    //statsblk stats;
    sparsegraph sg;

    DYNALLSTAT(int, lab, lab_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, orbits, orbits_sz);
    //static DEFAULTOPTIONS_SPARSEGRAPH(options);

    SG_INIT(sg);
    int m = SETWORDSNEEDED(g->v_size);
    nauty_check(WORDSIZE, m, g->v_size, NAUTYVERSIONID);
    SG_ALLOC(sg, g->v_size, g->e_size, "malloc");
    sg.nv = g->v_size;
    sg.nde = g->e_size;
    static DEFAULTOPTIONS_TRACES(options);
    // static DEFAULTOPTIONS_SPARSEGRAPH(options);
    //options.schreier = true;
    // options.writeautoms = true;
    //schreier_fails(10);
    options.defaultptn = false;
    options.userautomproc = test_auto;

    DYNALLOC1(int, lab, lab_sz, sg.nv, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, sg.nv, "malloc");
    DYNALLOC1(int, orbits, orbits_sz, sg.nv, "malloc");

    for (int i = 0; i < g->v_size; ++i) {
        lab[i] = i;
        ptn[i] = 1;
        sg.v[i] = g->v[i];
        sg.d[i] = g->d[i];
    }
    ptn[g->v_size - 1] = 0;
    if(colmap != nullptr) {
        make_lab_ptn_from_colmap(lab, ptn, colmap, g->v_size);
    }

    for (int i = 0; i < g->e_size; ++i) {
        sg.e[i] = g->e[i];
    }

    // touch the graph (mitigate cache variance
    int acc = 0;
    for(int i = 0; i < g->v_size; ++i) {
        acc += sg.v[i] + sg.d[i];
    }
    for(int i = 0; i < g->e_size; ++i) {
        acc += sg.e[i];
    }
    Clock::time_point timer = Clock::now();

    //sparsenauty(&sg,lab,ptn,orbits,&options,&stats,NULL);
    Traces(&sg, lab, ptn, orbits, &options, &stats, NULL);
    std::cout << "Group size: ";
    writegroupsize(stdout, stats.grpsize1, stats.grpsize2);
    std::cout << "Gens: " << stats.numgenerators << std::endl;
    std::cout << std::endl;
    std::cout << "Num nodes: " << stats.numnodes << std::endl;
    *traces_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;
}

int commandline_mode(int argc, char **argv) {
    std::string filename = "";
    bool entered_file = false;
    int  timeout = -1;
    std::fstream stat_file;
    std::string stat_filename = "test.dat";
    bool entered_stat_file = true;
    bool comp_nauty = false;
    bool comp_traces = false;
    bool comp_saucy = false;
    bool comp_dejavu = false;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    bool permute_graph = false;

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
        }

        if (arg == "__STAT_FILE") {
            if (i + 1 < argc) {
                i++;
                stat_filename = argv[i];
                entered_stat_file = true;
            } else {
                std::cerr << "--stat_file option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "__TIMEOUT") {
            if (i + 1 < argc) {
                i++;
                timeout = atoi(argv[i]);
                entered_stat_file = true;
            } else {
                std::cerr << "--timeout option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "__COMPRESS") {
            config.CONFIG_PREPROCESS_COMPRESS      = true;
            config.CONFIG_PREPROCESS_EDGELIST_SORT = true;
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

        if (arg == "__NAUTY") {
            comp_nauty = true;
        }

        if (arg == "__TRACES") {
            comp_traces = true;
        }

        if (arg == "__DEJAVU") {
            comp_dejavu = true;
        }

        if (arg == "__SAUCY") {
            comp_saucy = true;
        }


        if (arg == "__NO_IDLESKIP") {
            config.CONFIG_IR_IDLE_SKIP = false;
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

        if (arg == "__PERMUTE") {
            permute_graph = true;
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

    if (!entered_file) {
        std::cerr << "--file not specified" << std::endl;
        return 1;
    }
    parser p;
    sgraph *g = new sgraph;
   // p.parse_dimacs_file_digraph(filename, g);
    std::cout << "Parsing " << filename << "..." << std::endl;
    int* colmap = nullptr;
    //p.parse_dimacs_file(filename, g, &colmap);
    p.parse_dimacs_file_fast(filename, g, &colmap);
    sgraph *_g = new sgraph;
    if(permute_graph) {
        bijection<int> pr;
        std::cout << "Generating random bijection (seed " << seed << ")..." << std::endl;
        bijection<int>::random_bijection(&pr, g->v_size, seed);
        std::cout << "Permuting graph..." << std::endl;
        g->permute_graph(_g, &pr); // permute graph
        if(colmap != nullptr)
            permute_colmap(&colmap, g->v_size, pr.map);
        _g = g;
    } else {
        _g = g;
    }

    if(comp_nauty) {
        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << "nauty" << std::endl;
        std::cout << "------------------------------------------------------------------" << std::endl;
    }
    double nauty_solve_time;
    if(comp_nauty) {
        finished = false;
        nauty_kill_request = 0;
        std::thread killer;
        if(timeout > 0)
            killer = std::thread(kill_thread, &nauty_kill_request, timeout);
        bench_nauty(_g, colmap, &nauty_solve_time);
        if(timeout > 0)
            killer.join();
    }
    if(comp_nauty) {
        std::cout << "Solve time: " << nauty_solve_time / 1000000.0 << "ms" << std::endl;
    }
    if(comp_traces) {
        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << "Traces" << std::endl;
        std::cout << "------------------------------------------------------------------" << std::endl;
    }
    double traces_solve_time;
    if(comp_traces) {
        finished = false;
        nauty_kill_request = 0;
        std::thread killer;
        if(timeout > 0)
            killer = std::thread(kill_thread, &nauty_kill_request, timeout);
        bench_traces(_g, colmap, &traces_solve_time);
        if(timeout > 0)
            killer.join();
    }
    if(comp_traces) {
        std::cout << "Solve time: " << traces_solve_time / 1000000.0 << "ms" << std::endl;
    }

    if(comp_saucy) {
        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << "saucy" << std::endl;
        std::cout << "------------------------------------------------------------------" << std::endl;
    }
    double saucy_solve_time;
    if(comp_saucy) {
        finished = false;
        nauty_kill_request = 0;
        std::thread killer;
        if(timeout > 0)
            killer = std::thread(kill_thread, &nauty_kill_request, timeout);
        bench_saucy(_g, colmap, &saucy_solve_time);
        if(timeout > 0)
            killer.join();
    }
    if(comp_saucy) {
        std::cout << "Solve time: " << saucy_solve_time / 1000000.0 << "ms" << std::endl;
    }

    if(comp_dejavu) {
        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << "dejavu" << std::endl;
        std::cout << "------------------------------------------------------------------" << std::endl;
    }

    double dejavu_solve_time;
    if(comp_dejavu) {
        finished = false;
        std::thread killer;
        if(timeout > 0)
            killer = std::thread(kill_thread, &dejavu_kill_request, timeout);
        bench_dejavu(_g, colmap, &dejavu_solve_time);
        if(timeout > 0)
            killer.join();
    }
    if(comp_dejavu) {
        std::cout << "Solve time: " << dejavu_solve_time / 1000000.0 << "ms" << std::endl;
    }

    // std::cout << "Compare (nauty): " << nauty_solve_time / dejavu_solve_time << std::endl;
    // std::cout << "Compare (Traces): " << traces_solve_time / dejavu_solve_time << std::endl;

    if (entered_stat_file) {
        stat_file.open(stat_filename, std::fstream::in | std::fstream::out | std::fstream::app);

        std::cout << "Appending to " << stat_filename << ".\n";
        stat_file << filename << " " << _g->v_size << " ";
        if(comp_dejavu)
            stat_file << dejavu_solve_time / 1000000.0 << " ";
        if(comp_nauty)
            stat_file << nauty_solve_time  / 1000000.0 << " ";
        if(comp_traces)
            stat_file << traces_solve_time / 1000000.0 << " ";
        if(comp_saucy)
            stat_file << saucy_solve_time / 1000000.0 << " ";
        stat_file << "\n";
        stat_file.close();
    }

    delete _g;

    return 0;
}


int main(int argc, char *argv[]) {
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "bench" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
     return commandline_mode(argc, argv);
}
