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
    int m = SETWORDSNEEDED(g->v.size());
    nauty_check(WORDSIZE,m,g->v.size(),NAUTYVERSIONID);
    SG_ALLOC(sg,g->v.size(),g->e.size(),"malloc");
    sg.nv  = g->v.size();
    sg.nde = g->e.size();
    //static DEFAULTOPTIONS_TRACES(options);
    static DEFAULTOPTIONS_SPARSEGRAPH(options);
    options.schreier = true;
    schreier_fails(10);
    options.defaultptn = false;
    //options.writeautoms = true;

    DYNALLOC1(int,lab,lab_sz,sg.nv,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,sg.nv,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,sg.nv,"malloc");

    for(int i = 0; i < g->v.size(); ++i) {
        lab[i] = i;
        ptn[i] = 1;
        sg.v[i] = g->v[i];
        sg.d[i] = g->d[i];
    }

    ptn[g->v.size() - 1] = 0;

    for(int i = 0; i < g->e.size(); ++i) {
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
    int m = SETWORDSNEEDED(g->v.size());
    nauty_check(WORDSIZE,m,g->v.size(),NAUTYVERSIONID);
    SG_ALLOC(sg,g->v.size(),g->e.size(),"malloc");
    sg.nv  = g->v.size();
    sg.nde = g->e.size();
    static DEFAULTOPTIONS_TRACES(options);
   // static DEFAULTOPTIONS_SPARSEGRAPH(options);
    //options.schreier = true;
    //options.writeautoms = true;
    schreier_fails(10);
    options.defaultptn = false;

    DYNALLOC1(int,lab,lab_sz,sg.nv,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,sg.nv,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,sg.nv,"malloc");

    for(int i = 0; i < g->v.size(); ++i) {
        lab[i] = i;
        ptn[i] = 1;
        sg.v[i] = g->v[i];
        sg.d[i] = g->d[i];
    }

    ptn[g->v.size() - 1] = 0;

    for(int i = 0; i < g->e.size(); ++i) {
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
    }

    if(!entered_file) {
        std::cerr << "--file not specified" << std::endl;
        return 1;
    }
    parser p;
    sgraph g;
    p.parse_dimacs_file(filename, &g);

    std::cout << "Path Sampling-----------------------------------------------------" << std::endl;
    Clock::time_point timer = Clock::now();
    auto_blaster A;
    bool done = false;
    if(config.CONFIG_THREADS_PIPELINE_DEPTH <= 0) {
        A.sample(&g, true, &done);
    } else {
        A.sample_pipelined(&g, true, &done, nullptr);
    }
    double solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Solve time: " << solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "nauty" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    timer = Clock::now();
    bench_nauty(&g);
    double nauty_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Solve time: " << nauty_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Traces" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    timer = Clock::now();
    bench_traces(&g);
    double traces_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Solve time: " << traces_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Compare (nauty): "  << nauty_solve_time / solve_time << std::endl;
    std::cout << "Compare (Traces): " << traces_solve_time / solve_time << std::endl;

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
        stat_file << g.v.size() << " " << solve_time / 1000000.0 << " " << nauty_solve_time / 1000000.0 << " " << traces_solve_time / 1000000.0 << "\n";
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
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/rantree/rantree/rantree-5000.bliss", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/lattice/lattice/lattice-30", &g);
     p.parse_dimacs_file("/home/markus/Downloads/graphs/k/k/k150.dimacs", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/mz/mz/mz-50", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/ag/ag/ag2-47", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/cfi/cfi/cfi-200", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/dac/dac/4pipe.bliss", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/dac/dac/fpga11_20.bliss", &g);
     //g = g.permute_graph(bijection::random_bijection(g.v.size())); // permute graph
    // canonically label the sgraph

    //bijection canon_p;
    //label_graph(&g, &canon_p);

    std::cout << "Path Sampling-----------------------------------------------------" << std::endl;
    Clock::time_point timer = Clock::now();
    auto_blaster A;
    bool done = false;
    if(config.CONFIG_THREADS_PIPELINE_DEPTH <= 0) {
        A.sample(&g, true, &done);
    } else {
        A.sample_pipelined(&g, true, &done, nullptr);
    }
    double solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Solve time: " << solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "nauty" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    timer = Clock::now();
    bench_nauty(&g);
    double nauty_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Solve time: " << nauty_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Traces" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    timer = Clock::now();
    bench_traces(&g);
    double traces_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Solve time: " << traces_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Compare (nauty): "  << nauty_solve_time / solve_time << std::endl;
    std::cout << "Compare (Traces): " << traces_solve_time / solve_time << std::endl;
    return 0;
}