// currently not used

#ifndef DEJAVU_LOWDEG_H
#define DEJAVU_LOWDEG_H


#include "sgraph.h"
#include "schreier_shared.h"
#include "refinement.h"
#include "group_shared.h"

class lowdeg {
    lowdeg* nest = nullptr;
    sgraph* g;
    coloring* c;
    int domain_size;
    bool abort = false;
    int reduced_domain_size;
    work_set_int g_to_reduced_g, reduced_g_to_g, deg1_counter;

    work_set_int deg1_to_counter;
    work_set_int* counters;
    int counters_sz;
public:
    std::pair<sgraph*, coloring*> preprocess(coloring* c, sgraph* g, refinement* R);
    long double postprocess(group_shared* G);

    std::pair<sgraph *, coloring *> preprocess2(coloring *c, sgraph *g, refinement *R);
};


#endif //DEJAVU_LOWDEG_H
