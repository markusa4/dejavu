#ifndef DEJAVU_SELECTOR_H
#define DEJAVU_SELECTOR_H


#include "coloring.h"
#include "sgraph.h"
#include "configuration.h"
#include "refinement.h"
#include <list>

enum selector_type {SELECTOR_FIRST, SELECTOR_LARGEST, SELECTOR_SMALLEST, SELECTOR_TRACES, SELECTOR_RANDOM};

struct strategy {
    bijection* leaf = nullptr;
    invariant* I    = nullptr;
    selector_type cell_selector_type = SELECTOR_FIRST;
    int           cell_selector_seed = 0;

    strategy() = default;
    strategy(bijection* leaf, invariant* I, selector_type cell_selector_type, int seed) {
        this->leaf = leaf;
        this->I    = I;
        this->cell_selector_type = cell_selector_type;
        this->cell_selector_seed = seed;
    }

    void replace(strategy* s) {
        leaf = s->leaf;
        I    = s->I;
        cell_selector_type = s->cell_selector_type;
        cell_selector_seed = s->cell_selector_seed;
    }
};

class selector {
    int skipstart = 0;
    int hint      = -1;
    int hint_sz   = -1;

    //std::list<std::pair<int, int>> largest_cache;
    ring_pair largest_cache;
    work_list non_trivial_list;
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

    int select_color_dynamic(sgraph *g, coloring *c, strategy *s);

    int select_color_traces(coloring *c);
};

#endif // DEJAVU_SELECTOR_H
