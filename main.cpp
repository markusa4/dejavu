#include <iostream>
#include "parser.h"
#include "dejavu.h"
#include <assert.h>

extern "C" {
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
int dejavu_kill_request = 0;

bool finished = false;

void bench_nauty(sgraph *g, double* nauty_solve_time) {
    Clock::time_point timer = Clock::now();
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

    for (int i = 0; i < g->e_size; ++i) {
        sg.e[i] = g->e[i];
    }
    sparsenauty(&sg, lab, ptn, orbits, &options, &stats, NULL);
    std::cout << "Group size: ";
    writegroupsize(stdout, stats.grpsize1, stats.grpsize2);
    std::cout << std::endl;
    //Traces(&sg,lab,ptn,orbits,&options,&stats,NULL);
    *nauty_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;
}


void bench_dejavu(sgraph *g, double* dejavu_solve_time) {
    Clock::time_point timer = Clock::now();
    dejavu d;
    d.automorphisms(g);
    *dejavu_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    finished = true;
}

void bench_traces(sgraph *g, double* traces_solve_time) {
    //sleep(1);
    Clock::time_point timer = Clock::now();
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

    for (int i = 0; i < g->e_size; ++i) {
        sg.e[i] = g->e[i];
    }
    //sparsenauty(&sg,lab,ptn,orbits,&options,&stats,NULL);
    Traces(&sg, lab, ptn, orbits, &options, &stats, NULL);
    std::cout << "Group size: ";
    writegroupsize(stdout, stats.grpsize1, stats.grpsize2);
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
    bool comp_nauty = true;
    bool comp_traces = true;
    bool comp_dejavu = true;

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
        if (arg == "--STAT_FILE") {
            if (i + 1 < argc) {
                i++;
                stat_filename = argv[i];
                entered_stat_file = true;
            } else {
                std::cerr << "--stat_file option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "--TIMEOUT") {
            if (i + 1 < argc) {
                i++;
                timeout = atoi(argv[i]);
                entered_stat_file = true;
            } else {
                std::cerr << "--timeout option requires one argument." << std::endl;
                return 1;
            }
        }

        if (arg == "--NO_NAUTY") {
            comp_nauty = false;
        }

        if (arg == "--NO_TRACES") {
            comp_traces = false;
        }

        if (arg == "--NO_DEJAVU") {
            comp_dejavu = false;
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
        if (arg == "--IR_CELL_SELECTOR") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_IR_CELL_SELECTOR = atoi(argv[i]);
            } else {
                std::cerr << "--ir_cell_selector option requires one argument." << std::endl;
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
   // p.parse_dimacs_file_digraph(filename, g);
    std::cout << "Parsing " << filename << "..." << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    p.parse_dimacs_file(filename, g);
    //sleep(1);
    std::cout << "Permuting graph..." << std::endl;
    sgraph _g;
    bijection pr;
    bijection::random_bijection(&pr, g->v_size);
    g->permute_graph(&_g, &pr); // permute graph
    delete g;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "dejavu" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    double dejavu_solve_time;

    if(comp_dejavu) {
        finished = false;
        std::thread bench(bench_dejavu, &_g, &dejavu_solve_time);
        Clock::time_point start = Clock::now();
        while(!finished) {
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            //std::cout << (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count()) / 1000000.0 << ", " << timeout * 1000.0 << std::endl;
            if(((std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count()) / 1000000.0) > timeout * 1000.0) {
                std::cout << "Killing" << std::endl;
                dejavu_kill_request = 1;
            }
        }
        bench.join();
    }

    std::cout << "Solve time: " << dejavu_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "nauty" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    double nauty_solve_time;
    if(comp_nauty) {
        finished = false;
        nauty_kill_request = 0;
        std::thread bench(bench_nauty, &_g, &nauty_solve_time);
        //bench_nauty(&_g, &nauty_solve_time);
        Clock::time_point start = Clock::now();
        while(!finished) {
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            //std::cout << (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count()) / 1000000.0 << ", " << timeout * 1000.0 << std::endl;
            if(((std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count()) / 1000000.0) > timeout * 1000.0) {
                std::cout << "Killing" << std::endl;
                nauty_kill_request = 1;
            }
        }
        bench.join();
    }
    std::cout << "Solve time: " << nauty_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Traces" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    double traces_solve_time;
    if(comp_traces) {
        finished = false;
        nauty_kill_request = 0;
        //bench_traces(&_g, &traces_solve_time);
        std::thread bench(bench_traces, &_g, &traces_solve_time);
        Clock::time_point start = Clock::now();
        while(!finished) {
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            //std::cout << (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count()) / 1000000.0 << ", " << timeout * 1000.0 << std::endl;
            if(((std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count()) / 1000000.0) > timeout * 1000.0) {
                std::cout << "Killing" << std::endl;
                nauty_kill_request = 1;
            }
        }
        bench.join();
    }
    std::cout << "Solve time: " << traces_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    // std::cout << "Compare (nauty): " << nauty_solve_time / dejavu_solve_time << std::endl;
    // std::cout << "Compare (Traces): " << traces_solve_time / dejavu_solve_time << std::endl;

    if (entered_stat_file) {
        stat_file.open(stat_filename, std::fstream::in | std::fstream::out | std::fstream::app);

        std::cout << "Appending to " << stat_filename << ".\n";
        stat_file << filename << " " << _g.v_size << " ";
        if(comp_dejavu)
            stat_file << dejavu_solve_time / 1000000.0 << " ";
        if(comp_nauty)
            stat_file << nauty_solve_time  / 1000000.0 << " ";
        if(comp_traces)
            stat_file << traces_solve_time / 1000000.0 << " ";
        stat_file << "\n";
        stat_file.close();
    }

    return 0;
}


int main(int argc, char *argv[]) {
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "dejavu" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    // return commandline_mode(argc, argv);

    // parse a sgraph
    parser p;
    sgraph g;
    //p.parse_dimacs_file(argv[1], &g);
     //  p.parse_dimacs_file("/home/markus/Downloads/graphs/rantree/rantree/rantree-50000.bliss", &g);
     // p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/lattice/lattice/lattice-30", &g);
    // p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/grid/grid/grid-3-20", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/k/k/k2500.dimacs", &g);
     //  p.parse_dimacs_file("/home/markus/Downloads/graphs/combinatorial/combinatorial/FTWKB11.bliss", &g);
      // p.parse_dimacs_file("/home/markus/Downloads/graphs/tran/tran/tran_56832_255744.bliss", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/cmz/cmz-50", &g);
      p.parse_dimacs_file("/home/markus/Downloads/graphs/hypercubes/18cube.bliss", &g);
      // p.parse_dimacs_file("/home/markus/CLionProjects/dejavu/graph_tools/k2000.dimacs", &g);
    // p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/mz-aug2/mz-aug2/mz-aug2-50", &g);
    // p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/ag/ag/ag2-49", &g);
      //p.parse_dimacs_file("/home/markus/Downloads/graphs/ranreg/ranreg/Ranreg131072.bliss", &g);
    // p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/sat_cfi/sat_cfi_dim/sat_cfi_mult_5000_d.dmc", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/cfi/cfi/cfi-200", &g);
  // .parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/paley/paley/paley-461", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/sts-sw/sts-sw/sts-sw-79-11", &g);
      //p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/sts/sts/sts-67", &g);
      // p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/had-sw/had-sw/had-sw-32-1", &g);
      // p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/had-sw/had-sw/had-sw-44-11", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/had/had/had-156", &g);
    //p.parse_dimacs_file_digraph("/home/markus/Downloads/graphs/rnd-3-reg_cfi/rnd-3-reg-4000-2", &g);
    // p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/latin/latin/latin-30", &g); // skiplevels / base_size thing
     // p.parse_dimacs_file_g("/home/markus/Downloads/graphssaucy/states/AS.g", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/ranreg/ranreg/32768.bliss", &g);
      //p.parse_dimacs_file("/home/markus/Downloads/graphs/ranreg/ranreg/Ranreg32768.bliss", &g);
       //  p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/pp/pp/pp-25-170", &g);
     // p.parse_dimacs_file("/home/markus/Downloads/graphs/ranreg/ranreg/Ranreg131072.bliss", &g);
        // p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/latin-sw/latin-sw/latin-sw-30-5", &g);
   // p.parse_dimacs_file("/home/markus/Downloads/graphs/cfi-rigid-t2-tar/cfi-rigid-t2/cfi-rigid-t2-0504-01-1", &g); // <- significantly faster here!
      // p.parse_dimacs_file("/home/markus/Downloads/graphs/ran2/ran2/ran2_5000_a.bliss", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/ransq/ransq/ransq_10000_a.bliss", &g);
    // p.parse_dimacs_file("/home/markus/Downloads/hypercubes/18cube.bliss", &g);
    // p.parse_dimacs_file("/home/markus/Downloads/graphs/dac/dac/pipe/7pipe.bliss", &g);
    // p.parse_dimacs_file("/home/markus/Downloads/graphs/dac/dac/other/hole12.bliss", &g);
     // p.parse_dimacs_file("/home/markus/Downloads/graphs/dac/dac/other/pret150_25_ms.bliss", &g); // smallest is crazy here...
      // p.parse_dimacs_file("/home/markus/Downloads/graphs/dac/dac/other/s4-4-3-9.bliss", &g);
    // p.parse_dimacs_file("/home/markus/Downloads/graphs/dac/dac/other/fpga11_20.bliss", &g);

    std::cout << "Permuting graph---------------------------------------------------" << std::endl;
    sgraph _g;
    bijection pr;
    bijection::random_bijection(&pr, g.v_size);
    g.permute_graph(&_g, &pr); // permute graph
    //sgraph backup_graph;
    //backup_graph.copy_graph(&_g);
     //_g = g;
    std::cout << "Path Sampling-----------------------------------------------------" << std::endl;
    int repeat = 1;
    double avg = 0;
    Clock::time_point timer;
    for (int i = 0; i < repeat; ++i) {
        timer = Clock::now();
        dejavu d;
        d.automorphisms(&_g);
        double solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
        avg += solve_time;
        std::cout << "Solve time: " << solve_time / 1000000.0 << "ms" << std::endl;
    }
    std::cout << "Avg solve time: " << avg / 1000000.0 << "ms" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "nauty" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    timer = Clock::now();
    double nauty_solve_time;
    //bench_nauty(&_g, &nauty_solve_time);
    std::cout << "Solve time: " << nauty_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Traces" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    double traces_solve_time;
    finished = false;
    nauty_kill_request = 0;
    std::thread bench(bench_traces, &_g, &traces_solve_time);
    Clock::time_point start = Clock::now();
    while(!finished) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        //std::cout << (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count()) / 1000000.0 << ", " << timeout * 1000.0 << std::endl;
        if(((std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count()) / 1000000.0) > 100 * 1000.0) {
            std::cout << "Killing" << std::endl;
            nauty_kill_request = 1;
        }
    }
    bench.join();

    std::cout << "Solve time: " << traces_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Compare (nauty): " << nauty_solve_time / avg << std::endl;
    std::cout << "Compare (Traces): " << traces_solve_time / avg << std::endl;
    return 0;
}
