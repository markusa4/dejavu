//
// Created by markus on 11/11/2019.
//

#include <iostream>
#include "lowdeg.h"
#include "refinement.h"

std::pair<sgraph*, coloring*> lowdeg::preprocess(sgraph* g) {
    sgraph* reduced_graph = new sgraph;
    reduced_graph->v = new int[g->v_size];
    reduced_graph->e = new int[g->e_size];
    reduced_graph->d = new int[g->d_size];
    reduced_graph->v_size = g->v_size;
    reduced_graph->e_size = g->e_size;
    reduced_graph->d_size = g->d_size;
    reduced_graph->max_degree = g->max_degree;
    //std::cout << reduced_graph->max_degree << std::endl;

    std::cout << "dejavu_lowdeg_pre-------------------------------------------------" << std::endl;
    work_set_int deg1_counter;
    deg1_counter.initialize(g->v_size);
    work_set_int g_to_reduced_g;
    g_to_reduced_g.initialize(g->v_size);

    int deg1_vert = 0;
    int deg2_vert = 0;
    int pairs = 0;

    coloring* c = new coloring;
    g->initialize_coloring(c);

    int v, vpos, epos, reduce_vpos;

    reduce_vpos = 0;
    for(int i = 0; i < c->lab_sz; ++i) {
        v = c->lab[i];
        if(g->d[v] != 1) {
            if(g->d[v] == 2)
                deg2_vert += 1;

            g_to_reduced_g.set(v, reduce_vpos++);
            continue;
        }

        if(g->d[v] == 1) {
            deg1_counter.inc_nr(g->e[g->v[v]]); // ToDo: collect more information here, worklists?
                                                // ToDo: also for enriching generators in the end
            deg1_vert += 1;

            if(deg1_counter.get(v) == 1) { // detect pairs
                pairs += 1;
            }
        }
    }

    std::cout << "deg1_vert: " << deg1_vert << std::endl;
    std::cout << "deg2_vert: " << deg2_vert << std::endl;
    std::cout << "pairs: "     << pairs << std::endl;

    std::cout << "Writing graph..." << std::endl;
    vpos = 0;
    epos = 0;

    int con_v, map_v;

    // ToDo: also write new coloring
    for(int i = 0; i < c->lab_sz; ++i) {
        v = c->lab[i]; // ToDo: need to be aware of "current cell" and split it according to deg1_counter (this is basically color refinement)
        if(g->d[v] == 1)
            continue;

        map_v = g_to_reduced_g.get(v);
        assert(map_v >= 0);
        assert(map_v < reduce_vpos);

        reduced_graph->v[map_v] = epos;
        int d = 0;
        for(int j = g->v[v]; j < g->v[v] + g->d[v]; ++j) {
            con_v = g->e[j];
            if(g->d[con_v] == 1)
                continue;
            reduced_graph->e[epos++] = g_to_reduced_g.get(con_v);
            ++d;
        }

        reduced_graph->d[map_v] = d;
    }

    reduced_graph->v_size = reduce_vpos;
    reduced_graph->e_size = epos;
    reduced_graph->d_size = reduce_vpos;
    std::cout << vpos << ", " << reduce_vpos << std::endl;

    std::cout << "Reduced (" << g->v_size << ", " << g->e_size << ") -> (" << reduced_graph->v_size << ", " << reduced_graph->e_size << ")" << std::endl;


    std::cout << "------------------------------------------------------------------" << std::endl;
    return std::pair<sgraph*, coloring*>(reduced_graph, c);
}

void lowdeg::postprocess(mschreier *gp, mpermnode *gens) {
    std::cout << "dejavu_lowdeg_post------------------------------------------------" << std::endl;

    std::cout << "------------------------------------------------------------------" << std::endl;
}
