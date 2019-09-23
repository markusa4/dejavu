#include <iostream>
#include "parser.h"
#include "refinement.h"
#include "selector.h"
#include "ir_tools.h"
#include "auto_blaster.h"
#include <assert.h>
#include "nauty/traces.h"
#include "nauty/naugroup.h"
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


int main() {
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "dejavu" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    // parse a sgraph
    parser p;
    sgraph g;
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/rantree/rantree/rantree-500.bliss", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/k/k/k-100", &g);
     p.parse_dimacs_file("/home/markus/Downloads/graphs/ag/ag/ag2-5", &g);
    //p.parse_dimacs_file("/home/markus/Downloads/graphs/dac/dac/hole11.bliss", &g);
    // canonically label the sgraph


    //bijection test;
    //test.map = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    //std::cout << g.certify_automorphism(test) << std::endl;

    //bijection canon_p;
    //label_graph(&g, &canon_p);

    std::cout << "Sampling..." << std::endl;
    Clock::time_point timer = Clock::now();
    auto_blaster A;
    bool done;
    A.sample(&g, true, &done);
    double solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Solve time: " << solve_time / 1000000.0 << "ms" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "nauty" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;

    timer = Clock::now();
    bench_traces(&g);
    solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer).count());
    std::cout << "Solve time: " << solve_time / 1000000.0 << "ms" << std::endl;
    return 0;
}