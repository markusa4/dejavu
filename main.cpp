#include <iostream>
#include "parser.h"
#include "refinement.h"
#include "selector.h"
#include "ir_tools.h"
#include "auto_blaster.h"
#include <assert.h>
#include "nauty/traces.h"
#include "nauty/naugroup.h"
#include "configuration.h"
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

class time_point;

void label_graph(sgraph* g, bijection* canon_p) {
    coloring c;
    g->initialize_coloring(&c);
    std::cout << "------------------" << std::endl;

    ir_tools IR;
    IR.label_graph(g, canon_p);
}

void bench_traces(sgraph* g) {
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


int main(int argc, char *argv[]) {
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "dejavu" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    // parse a sgraph
    parser p;
    sgraph g;
    p.parse_dimacs_file(argv[1], &g);
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/rantree/rantree/rantree-10000.bliss", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/lattice/lattice/lattice-30", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/k/k/k-100", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/mz/mz/mz-50", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/ag/ag/ag2-47", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/cfi/cfi/cfi-50", &g);
     //p.parse_dimacs_file("/home/markus/Downloads/graphs/dac/dac/5pipe.bliss", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/dac/dac/fpga11_20.bliss", &g);
     //g = g.permute_graph(bijection::random_bijection(g.v.size())); // permute graph
    // canonically label the sgraph

    //bijection canon_p;
    //label_graph(&g, &canon_p);

    std::cout << "Path Sampling-----------------------------------------------------" << std::endl;
    Clock::time_point timer = Clock::now();
    auto_blaster A;
    bool done = false;
    if(CONFIG_SIFT_PIPELINE_DEPTH <= 0) {
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
    bench_traces(&g);
    double nauty_solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Solve time: " << nauty_solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Compare: " << nauty_solve_time / solve_time << std::endl;
    return 0;
}