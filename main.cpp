#include <iostream>
#include "parser.h"
#include "refinement.h"
#include "selector.h"
#include "ir_tools.h"
#include "auto_blaster.h"
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

void label_graph(sgraph* g, bijection* canon_p) {
    coloring c;
    g->initialize_coloring(&c);
    std::cout << "------------------" << std::endl;

    ir_tools IR;
    IR.label_graph(g, canon_p);
}

void bench_nauty(sgraph* g) {
    //TracesStats stats;
    statsblk stats;
    sparsegraph sg;

    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    //static DEFAULTOPTIONS_SPARSEGRAPH(options);

    SG_INIT(sg);
    int m = SETWORDSNEEDED(g->v_size);
    nauty_check(WORDSIZE,m,g->v_size,NAUTYVERSIONID);
    SG_ALLOC(sg,g->v_size,g->e_size,"malloc");
    sg.nv  = g->v_size;
    sg.nde = g->e_size;
    //static DEFAULTOPTIONS_TRACES(options);
    static DEFAULTOPTIONS_SPARSEGRAPH(options);
    options.schreier = true;
    schreier_fails(10);
    options.defaultptn = false;
    //options.writeautoms = true;

    DYNALLOC1(int,lab,lab_sz,sg.nv,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,sg.nv,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,sg.nv,"malloc");

    for(int i = 0; i < g->v_size; ++i) {
        lab[i] = i;
        ptn[i] = 1;
        sg.v[i] = g->v[i];
        sg.d[i] = g->d[i];
    }

    ptn[g->v_size - 1] = 0;

    for(int i = 0; i < g->e_size; ++i) {
        sg.e[i] = g->e[i];
    }
    sparsenauty(&sg,lab,ptn,orbits,&options,&stats,NULL);
    std::cout << "Group size: ";
    writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
    std::cout << std::endl;
    //Traces(&sg,lab,ptn,orbits,&options,&stats,NULL);
}

void bench_traces(sgraph* g) {
    TracesStats stats;
    //statsblk stats;
    sparsegraph sg;

    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    //static DEFAULTOPTIONS_SPARSEGRAPH(options);

    SG_INIT(sg);
    int m = SETWORDSNEEDED(g->v_size);
    nauty_check(WORDSIZE,m,g->v_size,NAUTYVERSIONID);
    SG_ALLOC(sg,g->v_size,g->e_size,"malloc");
    sg.nv  = g->v_size;
    sg.nde = g->e_size;
    static DEFAULTOPTIONS_TRACES(options);
   // static DEFAULTOPTIONS_SPARSEGRAPH(options);
    //options.schreier = true;
    //options.writeautoms = true;
    schreier_fails(10);
    options.defaultptn = false;

    DYNALLOC1(int,lab,lab_sz,sg.nv,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,sg.nv,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,sg.nv,"malloc");

    for(int i = 0; i < g->v_size; ++i) {
        lab[i] = i;
        ptn[i] = 1;
        sg.v[i] = g->v[i];
        sg.d[i] = g->d[i];
    }

    ptn[g->v_size - 1] = 0;

    for(int i = 0; i < g->e_size; ++i) {
        sg.e[i] = g->e[i];
    }
    //sparsenauty(&sg,lab,ptn,orbits,&options,&stats,NULL);
    Traces(&sg,lab,ptn,orbits,&options,&stats,NULL);
    std::cout << "Group size: ";
    writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
    std::cout << std::endl;
}

int commandline_mode(int argc, char** argv) {
    std::string filename = "";
    bool entered_file = false;
    std::fstream stat_file;
    std::string stat_filename = "test.dat";
    bool entered_stat_file = true;
    bool comp_nauty = true;
    bool comp_traces = true;

    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--file") {
            if (i + 1 < argc) {
                i++;
                filename = argv[i];
                entered_file = true;
            } else {
                std::cerr << "--file option requires one argument." << std::endl;
                return 1;
            }
        }
        if (std::string(argv[i]) == "--stat_file") {
            if (i + 1 < argc) {
                i++;
                stat_filename = argv[i];
                entered_stat_file = true;
            } else {
                std::cerr << "--stat_file option requires one argument." << std::endl;
                return 1;
            }
        }
        if (std::string(argv[i]) == "--no_nauty") {
            comp_nauty = false;
        }
        if (std::string(argv[i]) == "--no_traces") {
            comp_traces = false;
        }
        if (std::string(argv[i]) == "--THREADS_REFINEMENT_WORKERS") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_THREADS_REFINEMENT_WORKERS = atoi(argv[i]);
            } else {
                std::cerr << "--THREADS_REFINEMENT_WORKERS option requires one argument." << std::endl;
                return 1;
            }
        }
        if (std::string(argv[i]) == "--THREADS_PIPELINE_DEPTH") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_THREADS_PIPELINE_DEPTH = atoi(argv[i]);
            } else {
                std::cerr << "--THREADS_PIPELINE_DEPTH option requires one argument." << std::endl;
                return 1;
            }
        }
        if (std::string(argv[i]) == "--IR_CELL_SELECTOR") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_IR_CELL_SELECTOR = atoi(argv[i]);
            } else {
                std::cerr << "--IR_CELL_SELECTOR option requires one argument." << std::endl;
                return 1;
            }
        }
        if (std::string(argv[i]) == "--IR_BACKTRACK") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_IR_BACKTRACK = (atoi(argv[i]) == 1);
            } else {
                std::cerr << "--IR_BACKTRACK option requires one argument." << std::endl;
                return 1;
            }
        }
        if (std::string(argv[i]) == "--IR_REFINEMENT") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_IR_REFINEMENT = atoi(argv[i]);
            } else {
                std::cerr << "--IR_REFINEMENT option requires one argument." << std::endl;
                return 1;
            }
        }
        if (std::string(argv[i]) == "--THREADS_COPYG") {
            if (i + 1 < argc) {
                i++;
                config.CONFIG_THREADS_COPYG = (atoi(argv[i]) == 1);
            } else {
                std::cerr << "--THREADS_COPYG option requires one argument." << std::endl;
                return 1;
            }
        }
    }

    if(!entered_file) {
        std::cerr << "--file not specified" << std::endl;
        return 1;
    }
    parser p;
    sgraph* g = new sgraph;
    p.parse_dimacs_file_digraph(filename, g);

    //sleep(1);
    std::cout << "Permuting graph---------------------------------------------------" << std::endl;
    sgraph _g;
    bijection pr;
    bijection::random_bijection(&pr, g->v_size);
    g->permute_graph(&_g, &pr); // permute graph
    delete g;
    std::cout << "Path Sampling-----------------------------------------------------" << std::endl;
    int repeat = 1;
    double avg = 0;
    Clock::time_point timer;
    for(int i = 0; i < repeat; ++i) {
        timer = Clock::now();
        auto_blaster A;
        shared_switches switches;
        if (config.CONFIG_THREADS_PIPELINE_DEPTH <= 0) {
            //A.sample(&_g, true, &done);
        } else {
            if (config.CONFIG_IR_REFINEMENT == 0) {
                //A.sample_pipelined(&_g, true, &switches, nullptr, nullptr, nullptr, nullptr, nullptr, -1, nullptr,  nullptr, nullptr,  nullptr);
                A.sample_shared(&_g, true, &switches, nullptr, nullptr, nullptr, nullptr, nullptr, -1, nullptr,  nullptr, nullptr,  nullptr);
            } else if (config.CONFIG_IR_REFINEMENT == 1) {
               // A.sample_pipelined_bucket(&_g, true, &done, nullptr, nullptr, nullptr);
            } else {
                std::cout << "Unknown IR_REFINEMENT." << std::endl;
            }
        }
        double solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
        avg += solve_time / repeat;
        std::cout << "Solve time: " << solve_time / 1000000.0 << "ms" << std::endl;
    }
    std::cout << "Avg solve time: " << avg / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "nauty" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    timer = Clock::now();
    bench_nauty(&_g);
    double nauty_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Solve time: " << nauty_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Traces" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    timer = Clock::now();
    bench_traces(&_g);
    double traces_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Solve time: " << traces_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Compare (nauty): "  << nauty_solve_time / avg << std::endl;
    std::cout << "Compare (Traces): " << traces_solve_time / avg << std::endl;

    if(entered_stat_file) {
        stat_file.open(stat_filename, std::fstream::in | std::fstream::out | std::fstream::app);

        // If file does not exist, Create new file
        if (!stat_file )
        {
            std::cout << "Creating new stat_file..";

            stat_file.open(stat_filename,  std::fstream::in | std::fstream::out | std::fstream::trunc);
            stat_file <<"V dejavu nauty Traces\n";

        }// use existing file
        std::cout<<"Appending to " << stat_filename << ".\n";
        stat_file << g->v_size << " " << avg / 1000000.0 << " " << nauty_solve_time / 1000000.0 << " " << traces_solve_time / 1000000.0 << "\n";
        stat_file.close();
    }

    return 0;
}



int main(int argc, char *argv[]) {
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "dejavu" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    return commandline_mode(argc, argv);

    // parse a sgraph
    parser p;
    sgraph g;
    //p.parse_dimacs_file(argv[1], &g);
     //p.parse_dimacs_file("/home/markus/Downloads/rantree/rantree/rantree-5000.bliss", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/lattice/lattice/lattice-30", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/k/k/k-100", &g);
      //p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/mz/mz/mz-50", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/ag/ag/ag2-49", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/ranreg/ranreg/Ranreg65536.bliss", &g);
      //p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/sat_cfi/sat_cfi_dim/sat_cfi_mult_5000_d.dmc", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/cfi/cfi/cfi-200", &g);
     p.parse_dimacs_file_digraph("/home/markus/Downloads/graphs/rnd-3-reg_cfi/rnd-3-reg-3000-3", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/latin/latin/latin-20", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/rantree/rantree/rantree-10000.bliss", &g);
   // p.parse_dimacs_file("/home/markus/Downloads/graphs/ranreg/ranreg/Ranreg65536.bliss", &g);
      //p.parse_dimacs_file("/home/markus/Downloads/cfi/cfi/cfi-200", &g);
   // p.parse_dimacs_file("/home/markus/Downloads/graphs/cfi-rigid-t2-tar/cfi-rigid-t2/cfi-rigid-t2-0408-03-2", &g); // <- significantly faster here!
     //p.parse_dimacs_file("C:\\Users\\Markus\\Downloads\\undirected_dim\\undirected_dim\\cfi\\cfi-200", &g);
    //p.parse_dimacs_file("C:\\Users\\Markus\\Downloads\\undirected_dim\\undirected_dim\\mz-aug2\\mz-aug2\\mz-aug2-22", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/ran2/ran2/ran2_3000_a.bliss", &g);
      //p.parse_dimacs_file("/home/markus/Downloads/ransq/ransq/ransq_2000_a.bliss", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/hypercubes/15cube.bliss", &g);
    // p.parse_dimacs_file("/home/markus/Downloads/graphs/dac/dac/4pipe.bliss", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/dac/dac/fpga11_20.bliss", &g);

    std::cout << "Permuting graph---------------------------------------------------" << std::endl;
    sgraph _g;
    bijection pr;
    bijection::random_bijection(&pr, g.v_size);
    g.permute_graph(&_g, &pr); // permute graph
   // _g = g;
    std::cout << "Path Sampling-----------------------------------------------------" << std::endl;
    int repeat = 1;
    double avg = 0;
    Clock::time_point timer;
    for(int i = 0; i < repeat; ++i) {
        timer = Clock::now();
        auto_blaster A;
        shared_switches switches;
        if (config.CONFIG_THREADS_PIPELINE_DEPTH <= 0) {
         //   A.sample(&_g, true, &done);
        } else {
            if (config.CONFIG_IR_REFINEMENT == 0) {
                //A.sample_pipelined(&_g, true, &switches, nullptr, nullptr, nullptr, nullptr, nullptr, -1, nullptr,  nullptr, nullptr,  nullptr);
                A.sample_shared(&_g, true, &switches, nullptr, nullptr, nullptr, nullptr, nullptr, -1, nullptr,  nullptr, nullptr,  nullptr);
            } else if (config.CONFIG_IR_REFINEMENT == 1) {
               // A.sample_pipelined_bucket(&_g, true, &done, nullptr);
            } else {
                std::cout << "Unknown IR_REFINEMENT." << std::endl;
            }
        }
        double solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
        avg += solve_time / repeat;
        std::cout << "Solve time: " << solve_time / 1000000.0 << "ms" << std::endl;
    }
    std::cout << "Avg solve time: " << avg / 1000000.0 << "ms" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "nauty" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    timer = Clock::now();
    //bench_nauty(&_g);
    double nauty_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Solve time: " << nauty_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Traces" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    timer = Clock::now();
    bench_traces(&_g);
    double traces_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Solve time: " << traces_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Compare (nauty): "  << nauty_solve_time / avg << std::endl;
    std::cout << "Compare (Traces): " << traces_solve_time / avg << std::endl;
    return 0;
}
