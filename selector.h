//
// Created by markus on 19.09.19.
//

#ifndef BRUTUS_SELECTOR_H
#define BRUTUS_SELECTOR_H


#include "coloring.h"
#include "sgraph.h"
#include "configuration.h"
#include "refinement.h"
#include <list>

class selector {
    int skipstart = 0;
    //std::list<std::pair<int, int>> largest_cache;
    ring_pair largest_cache;
    int init = false;

public:
    int select_color(sgraph *g, coloring *c, int seed);

    int seeded_select_color(sgraph *g, coloring *c, int seed);

    int select_color_largest(coloring *c);

    int select_color_smallest(sgraph *g, coloring *c);

    int select_color_first(sgraph *g, coloring *c);

    int select_color_largest_degseq2(sgraph *g, coloring *c);

    void empty_cache();

    void pop_cache();
};


#endif //BRUTUS_SELECTOR_H
