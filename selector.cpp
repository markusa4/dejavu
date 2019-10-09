//
// Created by markus on 19.09.19.
//

#include <iostream>
#include "selector.h"

int selector::select_color_first(sgraph *g, coloring *c) {
    int first_cell  = -1;
    for(int i = 0; i < c->ptn_sz;){
        if(c->ptn[i] > 0) {
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
    for(int i = 0; i < c->ptn_sz;){
        if(c->ptn[i] < smallest_cell_sz && c->ptn[i] > 0) {
            smallest_cell = i;
            smallest_cell_sz = c->ptn[i];
        }
        i += c->ptn[i] + 1;
    }
    return smallest_cell;
}

int selector::select_color_largest(sgraph *g, coloring *c) {
    int largest_cell  = -1;
    int largest_cell_sz = -1;
    for(int i = 0; i < c->ptn_sz;){ // c->ptn[i] > largest_cell_sz &&
        if(c->ptn[i] > largest_cell_sz && c->ptn[i] > 0) {
            largest_cell = i;
            largest_cell_sz = c->ptn[i];
        }
        i += c->ptn[i] + 1;
    }
    return largest_cell;
}

// ToDo: special code to make this happen?
// rantree, pipe, etc. only have few basepoints left
// "special code" simple version: just randomly permute leftover color class (since they are K)
// base has to be properly fixed though!
int selector::select_color_largest_degseq2(sgraph *g, coloring *c) {
    int largest_cell  = -1;
    int largest_cell_sz = -1;
    int largest_cell_deg = -1;
    for(int i = 0; i < c->ptn_sz;) { // c->ptn[i] > largest_cell_sz &&
        if((c->ptn[i] > largest_cell_sz && (g->d[c->lab[i]] > 1) &&  c->ptn[i] > 0)) {
            largest_cell = i;
            largest_cell_sz = c->ptn[i];
            largest_cell_deg = g->d[c->lab[i]];
        }
        i += c->ptn[i] + 1;
    }
   // if(largest_cell_deg == 1)
     //   std::cout << "picked deg 1" << std::endl;
    return largest_cell;
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
    switch(config.CONFIG_IR_CELL_SELECTOR) {
        case 0:
            return seeded_select_color(g, c, seed);
        case 1:
            return select_color_largest(g, c);
        case 2:
            return select_color_smallest(g, c);
        case 3:
        default:
            return select_color_first(g, c);
    }
}

std::pair<int, int> selector::select_color_bucket(sgraph *g, coloring_bucket *c, int seed, int level) {
    int last_start = -1;

    int largest_cell_sz  = g->v.size()+ 1;
    int largest_cell_pos = -1;

    for(int i = 0; i < c->lab_sz; ++i) {
        if(c->ptn[i] <= level) {
            int cell_sz = i - last_start;
            if(cell_sz > 1 && cell_sz < largest_cell_sz) {
                largest_cell_sz = cell_sz;
                largest_cell_pos = last_start + 1;
            }
            last_start = i;
        }
    }

    return std::pair<int, int>(largest_cell_pos, largest_cell_sz);
}

