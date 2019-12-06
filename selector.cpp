#include <iostream>
#include <assert.h>
#include "selector.h"

int selector::select_color_first(sgraph *g, coloring *c) {
    int first_cell  = -1;

    for(int i = skipstart; i < c->ptn_sz;){
        if(c->ptn[i] > 0) {
            skipstart = i;
            first_cell = i;
            break;
        }
        i += c->ptn[i] + 1;
    }
    return first_cell;
}


// "first largest", -1 if coloring is discrete
int selector::select_color_smallest(sgraph *g, coloring *c) {
    int smallest_cell  = -1;
    int smallest_cell_sz = c->lab_sz + 1;
    bool only_trivial = true;
    for(int i = skipstart; i < c->ptn_sz;){
        if (c->ptn[i] != 0 && only_trivial) {
            skipstart = i;
            only_trivial = false;
        }
        if(c->ptn[i] < smallest_cell_sz && c->ptn[i] > 0) {
            smallest_cell = i;
            smallest_cell_sz = c->ptn[i];
            if(smallest_cell_sz == 2)
                break;
        }
        i += c->ptn[i] + 1;
    }
    return smallest_cell;
}

void selector::empty_cache() {
    skipstart = 0;
    hint      = -1;
    hint_sz   = -1;
    largest_cache.reset();
}

void selector::pop_cache() {
    largest_cache.pop();
}

int selector::select_color_largest(coloring *c) {
    if(!init) {
        largest_cache.initialize(c->lab_sz);
        non_trivial_list.initialize(c->lab_sz);
        init = true;
    }

    //if(hint != -1 && c->ptn[hint] == hint_sz)
    //    return hint;

    while(!largest_cache.empty()) {
        std::pair<int, int>* elem = largest_cache.front();
        if(c->ptn[elem->first] == elem->second)
            return elem->first;
        largest_cache.pop();
    }

    int largest_cell  = -1;
    int largest_cell_sz = -1;
    bool only_trivial = true;

    assert(skipstart < c->ptn_sz);
    for(int i = skipstart; i < c->ptn_sz;) { // c->ptn[i] > largest_cell_sz &&
        assert(i < c->ptn_sz);
        if (c->ptn[i] != 0 && only_trivial) {
            skipstart = i;
            only_trivial = false;
        }
        if (c->ptn[i] > largest_cell_sz && c->ptn[i] > 0) {
            largest_cell = i;
            largest_cell_sz = c->ptn[i];
            largest_cache.reset();
            largest_cache.push_back(std::pair<int, int>(i, c->ptn[i]));
            //hint      = -1;
            //hint_sz   = -1;
        } else if(c->ptn[i] == largest_cell_sz) {
            largest_cache.push_back(std::pair<int, int>(i, c->ptn[i]));
            //hint    = i;
            //hint_sz = c->ptn[i];
        }

        i += c->ptn[i] + 1;
    }
    return largest_cell;
}

int selector::select_color_traces(coloring *c) {
    //if(!init) {
    //    largest_cache.initialize(c->lab_sz);
    //    non_trivial_list.initialize(c->lab_sz);
    //    init = true;
    //}

    // int stack_sz = c->color_choices.size();

    // if(stack_sz == 0) return select_color_largest(c);

    // int col, col_sz;
    // int largest_cell    = -1;
    // int largest_cell_sz = -1;

    //for(int j = 0; j < stack_sz; ++j) {
        // col    = c->color_choices[stack_sz - j - 1].first;
        // col_sz = c->color_choices[stack_sz - j - 1].second;

        // largest_cell    = -1;
        // largest_cell_sz = -1;

        // for (int i = col; i < col + col_sz;) { // c->ptn[i] > largest_cell_sz &&
        //     if(c->ptn[i] >= largest_cell_sz && c->ptn[i] > 0) {
        //         largest_cell_sz = c->ptn[i];
        //         largest_cell    = i;
        //     }
        //     i += c->ptn[i] + 1;
        // }

        // if(largest_cell_sz > 0) break;
    // }

    // if(largest_cell != -1)
    //     return largest_cell;
    // else
    //     return select_color_largest(c);

    return -1;
}

int selector::seeded_select_color(sgraph *g, coloring *c, int seed) {
    std::vector<int> cells;
    for(int i = 0; i < c->ptn_sz;){
        if(c->ptn[i] > 0) {
            cells.push_back(i);
        }
        i += c->ptn[i] + 1;
    }
    if(cells.size() == 0) {
        return -1;
    } else {
        int target_cell = seed % cells.size();
        return cells[target_cell];
    }
}

int selector::select_color(sgraph *g, coloring *c, int seed) {
   // assert(config.CONFIG_IR_CELL_SELECTOR == 1);
   if(c->cells == g->v_size)
       return -1;

    switch(config.CONFIG_IR_CELL_SELECTOR) {
        case 0:
            return seeded_select_color(g, c, seed);
        case 1:
            return select_color_largest(c);
        case 2:
            return select_color_smallest(g, c);
        case 3:
        default:
            return select_color_first(g, c);
    }
}

int selector::select_color_dynamic(sgraph *g, coloring *c, strategy* s) {
    switch(s->cell_selector_type) {
        case SELECTOR_RANDOM:
            return seeded_select_color(g, c, s->cell_selector_seed);
        case SELECTOR_LARGEST:
            //return select_color_traces(c);
            return select_color_largest(c);
        case SELECTOR_SMALLEST:
            return select_color_smallest(g, c);
        case SELECTOR_TRACES:
            return select_color_largest(c);
            //return select_color_traces(c);
        case SELECTOR_FIRST:
        default:
            //return select_color_traces(c);
            return select_color_first(g, c);
    }
}