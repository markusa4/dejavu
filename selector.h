#ifndef SASSY_SELECTOR_H
#define SASSY_SELECTOR_H

#include "coloring.h"
#include "sgraph.h"

namespace sassy {
    // ring queue for pairs of integers
    class ring_pair {
        std::pair<int, int> *arr = 0;
        bool init = false;
        int arr_sz = -1;
        int front_pos = -1;
        int back_pos = -1;

    public:
        void initialize(int size) {
            if(init)
                delete[] arr;
            arr = new std::pair<int, int>[size];
            arr_sz = size;
            back_pos = 0;
            front_pos = 0;
            init = true;
        }

        void push_back(std::pair<int, int> value) {
            arr[back_pos] = value;
            back_pos = (back_pos + 1) % arr_sz;
        }

        std::pair<int, int> *front() {
            return &arr[front_pos];
        }

        void pop() {
            front_pos = (front_pos + 1) % arr_sz;
        }

        bool empty() {
            return (front_pos == back_pos);
        }

        ~ring_pair() {
            if (init)
                delete[] arr;
        }

        void reset() {
            front_pos = back_pos;
        }
    };

    enum selector_type {
        SELECTOR_FIRST, SELECTOR_LARGEST, SELECTOR_SMALLEST, SELECTOR_TRACES, SELECTOR_RANDOM
    };

    struct strategy {
        selector_type cell_selector_type = SELECTOR_FIRST;
        int cell_selector_seed = 0;
        bool init = false;

        strategy() = default;
    };

    class selector {
        int skipstart = 0;
        int hint = -1;
        int hint_sz = -1;

        ring_pair largest_cache;
        dejavu::work_list non_trivial_list;
        int init = false;

    public:
        int seeded_select_color(coloring *c, int seed) {
            std::vector<int> cells;
            for (int i = 0; i < c->ptn_sz;) {
                if (c->ptn[i] > 0) {
                    cells.push_back(i);
                }
                i += c->ptn[i] + 1;
            }
            if (cells.size() == 0) {
                return -1;
            } else {
                int target_cell = seed % cells.size();
                return cells[target_cell];
            }
        }

        int select_color_largest(coloring *c) {
            if (!init) {
                largest_cache.initialize(c->lab_sz);
                non_trivial_list.allocate(c->lab_sz);
                init = true;
            }

            while (!largest_cache.empty()) {
                std::pair<int, int> *elem = largest_cache.front();
                if (c->ptn[elem->first] == elem->second) {
                    assert(c->ptn[elem->first] > 0);
                    return elem->first;
                }
                largest_cache.pop();
            }

            int largest_cell = -1;
            int largest_cell_sz = -1;
            bool only_trivial = true;

            assert(skipstart < c->ptn_sz);
            for (int i = skipstart; i < c->ptn_sz;) {
                assert(c->vertex_to_col[c->lab[i]] == i);
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
                } else if (c->ptn[i] == largest_cell_sz) {
                    largest_cache.push_back(std::pair<int, int>(i, c->ptn[i]));
                }

                i += c->ptn[i] + 1;
            }
            assert(c->ptn[largest_cell] > 0);
            return largest_cell;
        }

        int select_color_smallest(coloring *c) {
            int smallest_cell = -1;
            int smallest_cell_sz = c->lab_sz + 1;
            bool only_trivial = true;
            for (int i = skipstart; i < c->ptn_sz; i += c->ptn[i] + 1) {
                if (c->ptn[i] != 0) {
                    if (only_trivial) {
                        skipstart = i;
                        only_trivial = false;
                    }
                    if ((c->ptn[i] + 1) < smallest_cell_sz) {
                        smallest_cell = i;
                        smallest_cell_sz = (c->ptn[i] + 1);
                        if (smallest_cell_sz == 2) {
                            break;
                        }
                    }
                }
            }
            return smallest_cell;
        }

        int select_color_first(coloring *c) {
            int first_cell = -1;

            for (int i = skipstart; i < c->ptn_sz;) {
                if (c->ptn[i] > 0) {
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
            hint = -1;
            hint_sz = -1;
            largest_cache.reset();
        }

        int
        select_color_dynamic(dejavu::sgraph *g, coloring *c, strategy *s) {
            if (c->cells == g->v_size)
                return -1;
            switch (s->cell_selector_type) {
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
}
#endif // SASSY_SELECTOR_H
