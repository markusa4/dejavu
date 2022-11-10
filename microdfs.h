#ifndef DEJAVU_MICRODFS_H
#define DEJAVU_MICRODFS_H

#include <random>
#include <chrono>
#include "sgraph.h"
#include "invariant.h"
#include "sassy.h"
#include "selector.h"
#include "bfs.h"
#include "schreier_sequential.h"

namespace dejavu {

    struct splitmap {

    };

    // heuristics to attempt to find a good selector, as well as creating a split-map for microdfs
    class cellwizard {
        selector find_selector(sgraph* g, coloring* c) {

        }

        splitmap find_splitmap(sgraph* g, coloring* c) {

        }
    };

    struct ir_state {
        coloring*    c;
        invariant*   I;
        mark_set*    touched_color;
        work_list*   touched_color_list;
        work_list*   prev_color_list;
        std::vector<int>  base_vertex;
        std::vector<int>  base_color;
        std::vector<int>  base_touched_color_list_pt;
        int base_pos = 0;
    };

    // minmal DFS IR-solver with limited failure, parallelizes according to given split-map
    class microdfs {
        int         fail_cnt = 0;
        int         fail_lim = 0;
        int         threads  = 1;
        refinement* R        = nullptr;
        selector*   S        = nullptr;
        splitmap*   SM       = nullptr;
    public:
        void setup(int fail_lim, int threads, refinement* R, selector* S, splitmap* SM) {
            this->R        = R;
            this->fail_lim = fail_lim;
            this->threads  = threads;
        }

        void move_to_leaf(sgraph* g, ir_state* state) {
            while(state->c->cells != g->v_size) {
                move_to_child(g, state);
            }
        }

        void move_to_child(sgraph* g, ir_state* state) {
            if (state->c->cells == g->v_size) {
                return;
            }
            ++state->base_pos;
            const int color_class = S->select_color_largest(state->c);
            const int ind_v = state->c->lab[color_class];
            const int init_color_class = R->individualize_vertex(state->c, ind_v);
            state->base_vertex.push_back(ind_v);
            state->base_color.push_back(color_class);
            state->base_touched_color_list_pt.push_back(state->touched_color_list->cur_pos);
            R->refine_coloring(g, state->c, state->I, init_color_class, nullptr, -1, -1,
                               nullptr, state->touched_color, state->touched_color_list, state->prev_color_list);
        }

        void move_to_parent(sgraph* g, ir_state* state) {
            // ToDo
            --state->base_pos;
            while(state->prev_color_list->cur_pos > state->base_touched_color_list_pt[state->base_pos]) {
                const int old_color = state->prev_color_list->pop_back();
                const int new_color = state->touched_color_list->pop_back();

                state->touched_color->unset(new_color);

                const int new_color_sz = state->c->ptn[new_color] + 1;
                state->c->ptn[old_color] += new_color_sz;
                state->c->ptn[new_color]  = 0;

                for(int j = 0; j < new_color_sz; ++j) {
                    const int v = state->c->lab[new_color + j];
                    state->c->vertex_to_col[v] = old_color;
                }

                --state->c->cells;
            }
        }

        void dfs(sgraph* g, coloring* c) {
            tiny_orbit orbs;
            orbs.initialize(g->v_size);

            selector Se;
            Se.empty_cache();
            S = &Se;

            invariant I;
            mark_set  touched_color(g->v_size);
            work_list touched_color_list(g->v_size);
            work_list prev_color_list(g->v_size);

            ir_state local_state;
            local_state.c = c;
            local_state.I = &I;
            local_state.touched_color = &touched_color;
            local_state.touched_color_list = &touched_color_list;
            local_state.prev_color_list = &prev_color_list;

            move_to_leaf(g, &local_state);

            for(int i = 0; i < local_state.base_color.size(); ++i) {
                std::cout << local_state.base_color[i] << " ";
            }
            std::cout << std::endl;
        }
    };
}

#endif //DEJAVU_MICRODFS_H
