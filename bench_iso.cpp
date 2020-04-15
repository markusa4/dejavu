#include <iostream>
#include "parser.h"
#include "vujade.h"
#include <assert.h>

extern "C" {
#include "nauty/traces.h"
}

extern "C" {
#include "conauto/conauto.h"
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

void sgraph_to_adjgraph(sgraph* g, adjgraph* ag) {
    ag->num_vert = g->v_size;
    ag->num_arc  = g->e_size;
    ag->deg = new uint64_t[g->v_size];
    for(int i = 0; i < g->v_size; ++i) {
        ag->deg[i] = (uint64_t)0;
        /*for (int j = 0; j < g->d[i]; ++j)
            ag->deg[i] += (uint64_t)UINT64_C(0x1000000000000);*/
    }

    ag->adj = new int8_t*[g->v_size];

    for (int from = 0; from < g->v_size; from++ ) {
        ag->adj[from] = new int8_t[g->v_size];
        for (int to = 0; to < g->v_size; to++ )
            ag->adj[from][to] = NOT_ADJ;
    }

    for(int i = 0; i < g->v_size; ++i) {
        for (int j = g->v[i]; j < g->v[i] + g->d[i]; ++j) {
            uint16_t from = i;
            uint16_t to   = g->e[j];
            if ( ag->adj[from][to] == NOT_ADJ )
            {
                ag->deg[from] = ag->deg[from] + UINT64_C(0x100000000) + (uint64_t)(from!=to);
                ag->deg[to] = ag->deg[to] + UINT64_C(0x10000) + (uint64_t)(from!=to);
            }
            else
            {
                ag->deg[from] = ag->deg[from] - UINT64_C(0x10000) + UINT64_C(0x1000000000000);
                ag->deg[to]   = ag->deg[to] - UINT64_C(0x100000000) + UINT64_C(0x1000000000000);
            }
            ag->adj[from][to] = ag->adj[from][to] | ARC_OUT;
            ag->adj[to][from] = ag->adj[to][from] | ARC_IN;
        }
    }
    std::cout << "Converted" << std::endl;
}

void bench_conauto(sgraph *g1, sgraph *g2, double* dejavu_solve_time) {
    // touch the graph (mitigate cache variance)
    adjgraph adj_g1;
    adjgraph adj_g2;
    sgraph_to_adjgraph(g1, &adj_g1);
    sgraph_to_adjgraph(g2, &adj_g2);

    Clock::time_point timer = Clock::now();
    int iso = are_isomorphic(&adj_g1, &adj_g2);
    *dejavu_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    if(iso) {
        std::cout << "ISOMORPHIC" << std::endl;
    } else {
        std::cout << "NON_ISOMORPHIC" << std::endl;
    }
    finished = true;
}

void bench_nauty(sgraph* g1, sgraph* g2, double* nauty_solve_time) {
    statsblk stats;
    sparsegraph sg2;

    DYNALLSTAT(int, lab2, lab_sz2);
    DYNALLSTAT(int, ptn2, ptn_sz2);
    DYNALLSTAT(int, orbits2, orbits_sz2);

    SG_INIT(sg2);
    int m1 = SETWORDSNEEDED(g2->v_size);
    nauty_check(WORDSIZE, m1, g2->v_size, NAUTYVERSIONID);
    SG_ALLOC(sg2, g2->v_size, g2->e_size, "malloc");
    sg2.nv = g2->v_size;
    sg2.nde = g2->e_size;
    static DEFAULTOPTIONS_SPARSEGRAPH(options1);
    options1.schreier = true;
    options1.defaultptn = false;
    options1.getcanon = true;


    DYNALLOC1(int, lab2, lab_sz2, sg2.nv, "malloc");
    DYNALLOC1(int, ptn2, ptn_sz2, sg2.nv, "malloc");
    DYNALLOC1(int, orbits2, orbits_sz2, sg2.nv, "malloc");

    for (int i = 0; i < g2->v_size; ++i) {
        lab2[i] = i;
        ptn2[i] = 1;
        sg2.v[i] = g2->v[i];
        sg2.d[i] = g2->d[i];
    }

    ptn2[g2->v_size - 1] = 0;

    for (int i = 0; i < g2->e_size; ++i) {
        sg2.e[i] = g2->e[i];
    }

    sparsegraph sg1;

    DYNALLSTAT(int, lab1, lab_sz1);
    DYNALLSTAT(int, ptn1, ptn_sz1);
    DYNALLSTAT(int, orbits1, orbits_sz1);
    //static DEFAULTOPTIONS_SPARSEGRAPH(options);

    SG_INIT(sg1);
    int m2 = SETWORDSNEEDED(g1->v_size);
    nauty_check(WORDSIZE, m2, g1->v_size, NAUTYVERSIONID);
    SG_ALLOC(sg1, g1->v_size, g1->e_size, "malloc");
    sg1.nv = g1->v_size;
    sg1.nde = g1->e_size;
    static DEFAULTOPTIONS_SPARSEGRAPH(options2);
    options2.schreier = true;
    options2.defaultptn = false;
    options2.getcanon = true;


    DYNALLOC1(int, lab1, lab_sz1, sg1.nv, "malloc");
    DYNALLOC1(int, ptn1, ptn_sz1, sg1.nv, "malloc");
    DYNALLOC1(int, orbits1, orbits_sz1, sg1.nv, "malloc");

    for (int i = 0; i < g1->v_size; ++i) {
        lab1[i] = i;
        ptn1[i] = 1;
        sg1.v[i] = g1->v[i];
        sg1.d[i] = g1->d[i];
    }

    ptn1[g1->v_size - 1] = 0;

    for (int i = 0; i < g1->e_size; ++i) {
        sg1.e[i] = g1->e[i];
    }

    Clock::time_point timer = Clock::now();

    SG_DECL(cg1);
    SG_DECL(cg2);

    sparsenauty(&sg1, lab1, ptn1, orbits1, &options1, &stats, &cg1);
    sparsenauty(&sg2, lab2, ptn2, orbits2, &options2, &stats, &cg2);
    if (aresame_sg(&cg1,&cg2)) {
        std::cout << "ISOMORPHIC" << std::endl;
    } else {
        std::cout << "NON_ISOMORPHIC" << std::endl;
    }
    //Traces(&sg,lab,ptn,orbits,&options,&stats,NULL);
    *nauty_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;
}


void bench_vujade(sgraph *g1, sgraph *g2, double* dejavu_solve_time) {
    // touch the graph (mitigate cache variance)
    Clock::time_point timer = Clock::now();
    vujade_iso(g1, g2);
    *dejavu_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;
}

void bench_traces(sgraph *g1, sgraph *g2, double* traces_solve_time) {
    TracesStats stats;
    sparsegraph sg1;
    sparsegraph sg2;
    DYNALLSTAT(int, lab2, lab_sz2);
    DYNALLSTAT(int, ptn2, ptn_sz2);
    DYNALLSTAT(int, orbits2, orbits_sz2);
    SG_INIT(sg2);
    int m2 = SETWORDSNEEDED(g2->v_size);
    nauty_check(WORDSIZE, m2, g2->v_size, NAUTYVERSIONID);
    SG_ALLOC(sg2, g2->v_size, g2->e_size, "malloc");
    sg2.nv = g2->v_size;
    sg2.nde = g2->e_size;
    static DEFAULTOPTIONS_TRACES(options);
    options.defaultptn = false;
    options.getcanon   = true;

    DYNALLOC1(int, lab2, lab_sz2, sg2.nv, "malloc");
    DYNALLOC1(int, ptn2, ptn_sz2, sg2.nv, "malloc");
    DYNALLOC1(int, orbits2, orbits_sz2, sg2.nv, "malloc");

    for (int i = 0; i < g2->v_size; ++i) {
        lab2[i] = i;
        ptn2[i] = 1;
        sg2.v[i] = g2->v[i];
        sg2.d[i] = g2->d[i];
    }

    ptn2[g2->v_size - 1] = 0;

    for (int i = 0; i < g2->e_size; ++i) {
        sg2.e[i] = g2->e[i];
    }

    DYNALLSTAT(int, lab1, lab_sz1);
    DYNALLSTAT(int, ptn1, ptn_sz1);
    DYNALLSTAT(int, orbits1, orbits_sz1);
    SG_INIT(sg1);
    int m = SETWORDSNEEDED(g1->v_size);
    nauty_check(WORDSIZE, m, g1->v_size, NAUTYVERSIONID);
    SG_ALLOC(sg1, g1->v_size, g1->e_size, "malloc");
    sg1.nv = g1->v_size;
    sg1.nde = g1->e_size;
    static DEFAULTOPTIONS_TRACES(options1);
    options1.defaultptn = false;
    options1.getcanon   = true;

    DYNALLOC1(int, lab1, lab_sz1, sg1.nv, "malloc");
    DYNALLOC1(int, ptn1, ptn_sz1, sg1.nv, "malloc");
    DYNALLOC1(int, orbits1, orbits_sz1, sg1.nv, "malloc");

    for (int i = 0; i < g1->v_size; ++i) {
        lab1[i] = i;
        ptn1[i] = 1;
        sg1.v[i] = g1->v[i];
        sg1.d[i] = g1->d[i];
    }

    ptn1[g1->v_size - 1] = 0;

    for (int i = 0; i < g1->e_size; ++i) {
        sg1.e[i] = g1->e[i];
    }
    Clock::time_point timer = Clock::now();

    SG_DECL(cg1);
    SG_DECL(cg2);
    Traces(&sg1, lab1, ptn1, orbits1, &options1, &stats, &cg1);
    Traces(&sg2, lab2, ptn2, orbits2, &options, &stats, &cg2);
    if (aresame_sg(&cg1,&cg2)) {
        std::cout << "ISOMORPHIC" << std::endl;
    } else {
        std::cout << "NON_ISOMORPHIC" << std::endl;
    }
    *traces_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;
}

int commandline_mode(int argc, char **argv) {
    std::string filename1 = "";
    std::string filename2 = "";
    bool entered_file1 = false;
    bool entered_file2 = false;
    int  timeout = -1;
    std::fstream stat_file;
    std::string stat_filename = "test.dat";
    bool entered_stat_file = true;
    bool comp_nauty = true;
    bool comp_traces = true;
    bool comp_dejavu = true;
    bool comp_conauto = true;
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
        }

        if ((arg == "__FILE2") || (arg == "__F2")) {
            if (i + 1 < argc) {
                i++;
                filename2 = argv[i];
                entered_file2 = true;
            } else {
                std::cerr << "--file2 requires one argument." << std::endl;
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

        if (arg == "__ERR") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_RAND_ABORT = atoi(argv[i]);
            } else {
                std::cerr << "--err option requires one argument." << std::endl;
                return 1;
            }
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

        if (arg == "__NO_NAUTY") {
            comp_nauty = false;
        }

        if (arg == "__NO_TRACES") {
            comp_traces = false;
        }

        if (arg == "__NO_DEJAVU") {
            comp_dejavu = false;
        }

        if (arg == "__NO_CONAUTO") {
            comp_conauto = false;
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

    if (!entered_file1 || !entered_file2) {
        std::cerr << "--file1 or --file2 not specified" << std::endl;
        return 1;
    }
    parser p;
    sgraph *g1 = new sgraph;
    std::cout << "Parsing " << filename1 << "..." << std::endl;
    p.parse_dimacs_file(filename1, g1);
    sgraph *g2 = new sgraph;
    std::cout << "Parsing " << filename2 << "..." << std::endl;
    p.parse_dimacs_file(filename2, g2);
    sgraph *_g1 = new sgraph;
    sgraph *_g2 = new sgraph;
    if(permute_graph) {
        std::cout << "Permuting graphs..." << std::endl;
        bijection<int> pr1;
        bijection<int>::random_bijection(&pr1, g1->v_size, seed);
        g1->permute_graph(_g1, &pr1); // permute graph
        delete g1;

        bijection<int> pr2;
        bijection<int>::random_bijection(&pr2, g2->v_size, seed * 123);
        g2->permute_graph(_g2, &pr2); // permute graph
        delete g2;
    } else {
        _g1 = g1;
        _g2 = g2;
    }

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "vujade" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    double dejavu_solve_time;

    if(comp_dejavu) {
        finished = false;
        std::thread killer;
        if(timeout > 0)
            killer = std::thread(kill_thread, &dejavu_kill_request, timeout);
        bench_vujade(_g1, _g2, &dejavu_solve_time);
        if(timeout > 0)
            killer.join();
    }

    std::cout << "Solve time: " << dejavu_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "nauty" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    double nauty_solve_time;
    if(comp_nauty) {
        finished = false;
        nauty_kill_request = 0;
        std::thread killer;
        if(timeout > 0)
            killer = std::thread(kill_thread, &nauty_kill_request, timeout);
        bench_nauty(_g1, _g2, &nauty_solve_time);
        if(timeout > 0)
            killer.join();
    }
    std::cout << "Solve time: " << nauty_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Traces" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    double traces_solve_time;
    if(comp_traces) {
        finished = false;
        nauty_kill_request = 0;
        std::thread killer;
        if(timeout > 0)
            killer = std::thread(kill_thread, &nauty_kill_request, timeout);
        bench_traces(_g1, _g2, &traces_solve_time);
        if(timeout > 0)
            killer.join();
    }
    std::cout << "Solve time: " << traces_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "conauto" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    double conauto_solve_time;
    if(comp_conauto) {
        finished = false;
        nauty_kill_request = 0;
        std::thread killer;
        if(timeout > 0)
            killer = std::thread(kill_thread, &nauty_kill_request, timeout);
        bench_conauto(_g1, _g2, &conauto_solve_time);
        if(timeout > 0)
            killer.join();
    }
    std::cout << "Solve time: " << conauto_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;

    if (entered_stat_file) {
        stat_file.open(stat_filename, std::fstream::in | std::fstream::out | std::fstream::app);

        std::cout << "Appending to " << stat_filename << ".\n";
        stat_file << filename1 << " " << _g1->v_size << " ";
        if(comp_dejavu)
            stat_file << dejavu_solve_time / 1000000.0 << " ";
        if(comp_nauty)
            stat_file << nauty_solve_time  / 1000000.0 << " ";
        if(comp_traces)
            stat_file << traces_solve_time / 1000000.0 << " ";
        if(comp_conauto)
            stat_file << conauto_solve_time / 1000000.0 << " ";
        stat_file << "\n";
        stat_file.close();
    }

    delete _g1;
    delete _g2;
    return 0;
}


int main(int argc, char *argv[]) {
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "dejavu" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
     return commandline_mode(argc, argv);
}
