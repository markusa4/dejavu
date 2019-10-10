//
// Created by markus on 19.09.19.
//

#ifndef BRUTUS_SELECTOR_H
#define BRUTUS_SELECTOR_H


#include "coloring.h"
#include "sgraph.h"
#include "configuration.h"
#include <list>

class selector {
    int skipstart = 0;
    std::list<std::pair<int, int>> largest_cache;

public:
    std::pair<int, int> select_color_bucket(sgraph *g, coloring_bucket *c, int seed, int level);

    int select_color(sgraph *g, coloring *c, int seed);

    int seeded_select_color(sgraph *g, coloring *c, int seed);

    int select_color_largest(sgraph *g, coloring *c);

    int select_color_smallest(sgraph *g, coloring *c);

    int select_color_first(sgraph *g, coloring *c);

    int select_color_largest_degseq2(sgraph *g, coloring *c);

    void empty_cache();
};


#endif //BRUTUS_SELECTOR_H
