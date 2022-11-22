#include <iostream>
#include "parser.h"
#include "dejavu_auto.h"
#include <assert.h>

extern "C" {
#include "saucy/saucy.h"
#include "nauty/traces.h"
#include "bliss/bliss_C.h"
}

#include "nauty/naugroup.h"
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
        if (colmap[lab[i]] != colmap[lab[i + 1]]) {
            ptn[i] = 0;
            ptn[last_new_cell] = (i - last_new_cell) > 0;
            last_new_cell = i + 1;
            continue;
        }
    }
    std::cout << last_new_cell << std::endl;
}

void make_small_colmap(int* lab, int* smallcolmap, int* colmap, int colmap_sz) {
    for(int i = 0; i < colmap_sz; ++i) {
        lab[i] = i;
    }

    int last_new_cell = 0;
    std::sort(lab, lab + colmap_sz, colorComparator(colmap));

    int col_small = -1;
    int last_col_orig  = INT32_MAX;

    for(int i = 0; i < colmap_sz; i++) {
        if (colmap[lab[i]] != last_col_orig) {
            last_col_orig = colmap[lab[i]];
            ++col_small;
        }
        smallcolmap[lab[i]] = col_small;
    }
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
    if(g->v_size > 0)
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

    if(g->v_size > 0) {
        sparsenauty(&sg, lab, ptn, orbits, &options, &stats, NULL);
    }

    std::cout << "Group size: ";
    writegroupsize(stdout, stats.grpsize1, stats.grpsize2);
    std::cout << std::endl;
    //Traces(&sg,lab,ptn,orbits,&options,&stats,NULL);
    *nauty_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;
}

sgraph* test_g;
int*    test_colmap;
refinement test_R;
bool all_certified = true;

void certifying_consumer_dense(int n, const int * p, int support, const int *) {
    //assert(test_auto);
    return;
}

void certifying_hook_sparse(int n, const int * p, int support, const int* support_arr) {
    const bool test_auto = test_R.certify_automorphism_sparse(test_g, test_colmap, p, support, support_arr);
    all_certified = all_certified && test_auto;
    assert(test_auto);
    return;
}

void empty_hook(int n, const int * p, int support, const int *) {
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
    dejavu_automorphisms(g, colmap, empty_hook);
    *dejavu_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;

    return acc;
}

int bench_dejavu_certify(sgraph_t<int, int, int> *g, int* colmap, double* dejavu_solve_time) {
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
    sgraph _test_g;
    _test_g.copy_graph(g);
    test_g = &_test_g;
    test_colmap = new int[g->v_size];
    test_R = refinement();
    for(int j = 0; j < g->v_size; ++j)
        test_colmap[j] = colmap[j];

    dejavu_automorphisms(g, colmap, certifying_hook_sparse);
    *dejavu_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;

    return acc;
}

static void blissHook(void *user_param, unsigned int N, const unsigned int *aut) {
    return;
}


int bench_bliss(sgraph_t<int, int, int> *g, int* colmap, double* dejavu_solve_time) {
    // touch the graph (mitigate cache variance)

    BlissGraph* bgraph = bliss_new(0);
    for(int i = 0; i < g->v_size; ++i) {
        bliss_add_vertex(bgraph, colmap[i]);
    }

    for(int i = 0; i < g->v_size; ++i) {
        const int npt = g->v[i];
        const int nd  = g->d[i];
        for(int j = 0; j < nd; ++j) {
            const int adj_n = g->e[npt + j];
            if(i < adj_n) {
                bliss_add_edge(bgraph, i, adj_n);
            }
        }
    }

    Clock::time_point timer = Clock::now();
    BlissStats stats;
    bliss_find_automorphisms(bgraph, &blissHook, nullptr, &stats);
    *dejavu_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());

    std::cout << "Group size: ";
    std::cout << stats.group_size_approx << std::endl;

    finished = true;
    return true;
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
    int* small_colmap         = new int[g->v_size];
    //int* small_colmap_ordered = new int[g->v_size];
    int* lab                  = new int[g->v_size];
    if(colmap == nullptr) {
        colmap = new int[g->v_size];
        for(int i = 0; i < g->v_size; ++i) {
            colmap[i] = 0;
        }
    }
    make_small_colmap(lab, small_colmap, colmap, g->v_size);

    saucy_graph sg;
    sg.n = g->v_size;
    sg.e = g->e_size;
    sg.edg = new int[g->e_size];
    sg.adj = new int[sg.n + 1];
    //sg.colors = new int[g->v_size];

    int epos = 0;

    for(int i = 0; i < g->v_size; ++i) {
        //const int i = lab[ii];
        const int npt = g->v[i];
        const int nd  = g->d[i];
        sg.adj[i] = epos;
        for(int j = 0; j < nd; ++j) {
            const int adj_n = g->e[npt + j];
            sg.edg[epos] = adj_n;
            ++epos;
        }
    }

    /*for(int ii = 0; ii < g->v_size; ++ii) {
        const int i = lab[ii];
        small_colmap_ordered[]
    }*/

    sg.adj[g->v_size] = epos;
    assert(epos == g->e_size);


    Clock::time_point timer = Clock::now();
    // automorphisms
    struct saucy_stats stats;
    if(g->v_size > 0) {
        struct saucy *s = saucy_alloc(sg.n); // 100000
        saucy_search(s, &sg, 0, small_colmap, saucyConsume, 0, &stats);
        saucy_free(s);
    }
    *saucy_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Group size: ";
    std::cout << stats.grpsize_base << "*10^" << stats.grpsize_exp << std::endl;
    finished = true;

    delete[] small_colmap;
    //delete[] small_colmap_ordered;
    delete[] lab;
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
    if(g->v_size > 0)
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

    if(g->v_size > 0) {
        //sparsenauty(&sg,lab,ptn,orbits,&options,&stats,NULL);
        Traces(&sg, lab, ptn, orbits, &options, &stats, NULL);
    }
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
    bool preprocess = false;
    bool comp_nauty = false;
    bool comp_traces = false;
    bool comp_saucy = false;
    bool comp_dejavu = false;
    bool comp_bliss  = false;
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

        if (arg == "__PREP_DEACT_DEG01") {
            config.CONFIG_PREP_DEACT_DEG01 = true;
        }
        if (arg == "__PREP_DEACT_DEG2") {
            config.CONFIG_PREP_DEACT_DEG2 = true;
        }
        if (arg == "__PREP_DEACT_PROBE") {
            config.CONFIG_PREP_DEACT_PROBE = true;
        }

        if (arg == "__PREP_ALT_SCHEDULE") {
            config.CONFIG_PREP_ALT_SCHEDULE = true;
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

        if (arg == "__PREP") {
            preprocess = true;
        }

        if (arg == "__DEJAVU_NOPREP") {
            config.CONFIG_PREPROCESS = false;
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
        if (arg == "__BLISS") {
            comp_bliss = true;
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
        bijection pr;
        std::cout << "Generating random bijection (seed " << seed << ")..." << std::endl;
        bijection::random_bijection(&pr, g->v_size, seed);
        std::cout << "Permuting graph..." << std::endl;
        g->permute_graph(_g, &pr); // permute graph
        if(colmap != nullptr)
            permute_colmap(&colmap, g->v_size, pr.map);
        _g = g;
    } else {
        _g = g;
    }

    if(colmap == nullptr) {
        colmap = new int[g->v_size];
        for(int i = 0; i < g->v_size; ++i) {
            colmap[i] = 0;
        }
    }

    double prep_time;
    if(preprocess) {
        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << "sassy" << std::endl;
        std::cout << "------------------------------------------------------------------" << std::endl;

        Clock::time_point timer = Clock::now();
        sassy preprocessor;
        preprocessor.reduce(_g, colmap, empty_hook);
        prep_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());

        std::cout << "Prep time: " << prep_time / 1000000.0 << "ms" << std::endl;
        std::cout << "graph sz: " << g->v_size << ", " << g->e_size << std::endl;
        std::cout << "(partial) group sz: " << preprocessor.base << "*10^" << preprocessor.exp << std::endl;
        std::cout << "(reported group sizes below must be multiplied with group size above)" << std::endl;
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
        std::cout << "Solve time: " << (nauty_solve_time + prep_time) / 1000000.0 << "ms" << std::endl;
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
        std::cout << "Solve time: " << (traces_solve_time + prep_time) / 1000000.0 << "ms" << std::endl;
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
        std::cout << "Solve time: " << (saucy_solve_time + prep_time) / 1000000.0 << "ms" << std::endl;
    }

    if(comp_bliss) {
        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << "bliss" << std::endl;
        std::cout << "------------------------------------------------------------------" << std::endl;
    }
    double bliss_solve_time;
    if(comp_bliss) {
        finished = false;
        nauty_kill_request = 0;
        std::thread killer;
        if(timeout > 0)
            killer = std::thread(kill_thread, &nauty_kill_request, timeout);
        bench_bliss(_g, colmap, &bliss_solve_time);
        if(timeout > 0)
            killer.join();
    }
    if(comp_bliss) {
        std::cout << "Solve time: " << (bliss_solve_time + prep_time) / 1000000.0 << "ms" << std::endl;
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
        std::cout << "Solve time: " << (dejavu_solve_time + prep_time) / 1000000.0 << "ms" << std::endl;
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
