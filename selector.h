#ifndef DEJAVU_SELECTOR_H
#define DEJAVU_SELECTOR_H


#include "coloring.h"
#include "sgraph.h"
#include "configuration.h"
#include "refinement.h"
#include <list>

enum selector_type {SELECTOR_FIRST, SELECTOR_LARGEST, SELECTOR_SMALLEST, SELECTOR_TRACES, SELECTOR_RANDOM};

template<class vertex_t>
struct strategy {
    bijection<vertex_t>* leaf = nullptr;
    invariant* I    = nullptr;
    selector_type cell_selector_type = SELECTOR_FIRST;
    int           cell_selector_seed = 0;

    strategy() = default;
    strategy(bijection<vertex_t>* leaf, invariant* I, selector_type cell_selector_type, int seed) {
        this->leaf = leaf;
        this->I    = I;
        this->cell_selector_type = cell_selector_type;
        this->cell_selector_seed = seed;
    }

    void replace(strategy<vertex_t>* s) {
        leaf = s->leaf;
        I    = s->I;
        cell_selector_type = s->cell_selector_type;
        cell_selector_seed = s->cell_selector_seed;
    }
};

template<class vertex_t, class degree_t, class edge_t>
class selector {
    int skipstart = 0;
    int hint      = -1;
    int hint_sz   = -1;

    ring_pair largest_cache;
    work_list non_trivial_list;
    int init = false;

public:
    int seeded_select_color(coloring<vertex_t> *c, int seed) {
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

    int select_color_largest(coloring<vertex_t> *c) {
        if(!init) {
            largest_cache.initialize(c->lab_sz);
            non_trivial_list.initialize(c->lab_sz);
            init = true;
        }

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
        for(int i = skipstart; i < c->ptn_sz;) {
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
            } else if(c->ptn[i] == largest_cell_sz) {
                largest_cache.push_back(std::pair<int, int>(i, c->ptn[i]));
            }

            i += c->ptn[i] + 1;
        }
        return largest_cell;
    }

    int select_color_smallest(coloring<vertex_t> *c) {
        int smallest_cell    = -1;
        int smallest_cell_sz = c->lab_sz + 1;
        bool only_trivial    = true;
        for(int i = skipstart; i < c->ptn_sz; i += c->ptn[i] + 1){
            if(c->ptn[i] != 0) {
                if(only_trivial) {
                    skipstart = i;
                    only_trivial = false;
                }
                if ((c->ptn[i] + 1) < smallest_cell_sz) {
                    smallest_cell = i;
                    smallest_cell_sz = (c->ptn[i] + 1);
                    if (smallest_cell_sz == 2) {
                        break;
                    }
                    if (smallest_cell_sz == c->smallest_cell_lower_bound + 1) {
                        break;
                    }
                }
            }
        }
        return smallest_cell;
    }

    int select_color_first(coloring<vertex_t> *c) {
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

    void empty_cache() {
        skipstart = 0;
        hint      = -1;
        hint_sz   = -1;
        largest_cache.reset();
    }

    int select_color(sgraph_t<vertex_t, degree_t, edge_t> *g, coloring<vertex_t> *c, int seed) {
        if(c->cells == g->v_size)
            return -1;
        switch(config.CONFIG_IR_CELL_SELECTOR) {
            case 0:
                return seeded_select_color(c, seed);
            case 1:
                return select_color_largest(c);
            case 2:
                return select_color_smallest(c);
            case 3:
            default:
                return select_color_first(c);
        }
    }

    int select_color_dynamic(sgraph_t<vertex_t, degree_t, edge_t> *g, coloring<vertex_t>  *c, strategy<vertex_t> *s) {
        if(c->cells == g->v_size)
            return -1;
        switch(s->cell_selector_type) {
            case SELECTOR_RANDOM:
                return seeded_select_color(c, s->cell_selector_seed);
            case SELECTOR_LARGEST:
                return select_color_largest(c);
            case SELECTOR_SMALLEST:
                return select_color_smallest(c);
            case SELECTOR_TRACES:
                return select_color_largest(c);
            case SELECTOR_FIRST:
            default:
                return select_color_first(c);
        }
    }
};

#endif // DEJAVU_SELECTOR_H
