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

    enum ir_mode { IR_MODE_COMPARE, IR_MODE_RECORD};

    struct ir_state {
        coloring*    c;
        trace*       T;
        mark_set*    touched_color;
        work_list*   touched_color_list;
        work_list*   prev_color_list;
        std::vector<int> singletons;

        std::vector<int>  base_vertex;
        std::vector<int>  base_color;
        std::vector<int>  base_touched_color_list_pt;
        std::vector<int>  base_cells;
        std::vector<int>  base_singleton_pt;

        std::vector<int>  compare_base_color;
        std::vector<int>  compare_base_cells;
        std::vector<int>  compare_singletons;

        ir_mode mode = IR_MODE_RECORD;

        int base_pos = 0;

    public:
        coloring* get_coloring() {
            return c;
        }

        // move this to static functions
        void  move_to_child(refinement* R, sgraph* g, int v) {
            if (c->cells == g->v_size) {
                return;
            }

            ++base_pos;
            base_singleton_pt.push_back(singletons.size());
            base_vertex.push_back(v);
            base_color.push_back(c->vertex_to_col[v]);
            base_touched_color_list_pt.push_back(touched_color_list->cur_pos);
            const int init_color_class = R->individualize_vertex(c, v, T, touched_color,
                                                                 touched_color_list, prev_color_list, &singletons);


            if(mode == IR_MODE_RECORD) {
                R->refine_coloring(g, c, init_color_class, -1, touched_color,
                                   touched_color_list, prev_color_list, T, &singletons, nullptr, nullptr);
            } else {
                R->refine_coloring(g, c, init_color_class, compare_base_cells[base_pos-1],
                                   touched_color, touched_color_list, prev_color_list, T, &singletons, nullptr, nullptr);
            }

            base_cells.push_back(c->cells);
        }
        void  move_to_parent() {
            // unwind invariant
            T->rewind_to_individualization();

            --base_pos;
            while(prev_color_list->cur_pos > base_touched_color_list_pt[base_pos]) {
                const int old_color = prev_color_list->pop_back();
                const int new_color = touched_color_list->pop_back();

                touched_color->unset(new_color);

                const int new_color_sz = c->ptn[new_color] + 1;
                c->ptn[old_color] += new_color_sz;
                c->ptn[new_color]  = 1;

                for(int j = 0; j < new_color_sz; ++j) {
                    const int v = c->lab[new_color + j];
                    c->vertex_to_col[v] = old_color;
                    assert(state->c->vertex_to_lab[v] == new_color + j);
                }

                --c->cells;
            }

            int const new_singleton_pos = base_singleton_pt.back();
            singletons.resize(new_singleton_pos);

            base_vertex.pop_back();
            base_color.pop_back();
            base_touched_color_list_pt.pop_back();
            base_cells.pop_back();
            base_singleton_pt.pop_back();
        }
    };

    struct splitmap {

    };

    // heuristics to attempt to find a good selector, as well as creating a split-map for microdfs
    class cellmagic {
        int locked_lim = 512;

    public:
        mark_set test_set;

        void  __attribute__ ((noinline)) find_base(refinement* R, sgraph* g, ir_state* state) {
            std::vector<int> candidates;
            test_set.initialize(g->v_size);
            candidates.reserve(locked_lim);
            int prev_color = -1;

            mark_set neighbour_color;
            neighbour_color.initialize(g->v_size);

            while(state->get_coloring()->cells != g->v_size) {
                int best_color = -1;

                // pick previous color if possible
                if(prev_color >= 0 && state->get_coloring()->ptn[prev_color] > 0) {
                    //std::cout << "previous color" << std::endl;
                    best_color = prev_color;
                } else if(prev_color >= 0) { // pick neighbour of previous color if possible
                    const int test_vertex = state->get_coloring()->lab[prev_color];
                    for(int i = 0; i < g->d[test_vertex]; ++i) {
                        const int other_vertex = g->e[g->v[test_vertex] + i];
                        const int other_color = state->get_coloring()->vertex_to_col[other_vertex];
                        if(state->get_coloring()->ptn[other_color] > 0) {
                            //std::cout << "neighbour color" << std::endl;
                            neighbour_color.set(other_color);
                            best_color = other_color;
                            //break;
                        }
                    }
                }

                // heuristic, try to pick "good" color
                if(best_color == -1) {
                    //std::cout << "heuristic color" << std::endl;
                    candidates.clear();
                    int best_score = -1;

                    int alt_best_color = -1;
                    int alt_best_score = -1;

                    for (int i = 0; i < state->get_coloring()->ptn_sz;) {
                        if (state->get_coloring()->ptn[i] > 0) {
                            candidates.push_back(i);
                        }

                        //if(candidates.size() >= (size_t) locked_lim)
                        //    break;

                        i += state->get_coloring()->ptn[i] + 1;
                    }

                    const int num_tested = candidates.size();

                    while (!candidates.empty()) {
                        const int test_color = candidates.back();
                        candidates.pop_back();

                        int test_score = color_score(g, state, test_color);
                        if(neighbour_color.get(test_color)) {
                            test_score *= 10;
                        }
                        if (test_score > best_score) {
                            best_color = test_color;
                            best_score = test_score;
                        }

                        int alt_test_score = alt_score(g, state, test_color);
                        if(neighbour_color.get(test_color)) {
                            alt_test_score *= 10;
                        }
                        if (alt_test_score > alt_best_score) {
                            alt_best_color = test_color;
                            alt_best_score = test_score;
                        }
                    }

                    if (best_color != alt_best_color) {
                        //std::cout << best_color << " vs. " << alt_best_color << std::endl;
                        state->move_to_child(R, g, state->get_coloring()->lab[best_color]);
                        const int score1 = state_score(g, state);
                        state->move_to_parent();
                        state->move_to_child(R, g, state->get_coloring()->lab[alt_best_color]);
                        const int score2 = state_score(g, state);
                        state->move_to_parent();

                        //std::cout << score1 << ":" << score2 << std::endl;
                        if (score2 > score1) {
                            best_color = alt_best_color;
                        }
                    }
                }

                prev_color = best_color;
                state->move_to_child(R, g, state->get_coloring()->lab[best_color]);
                //std::cout << "picked color " << best_color << " score " << best_score << " among " << num_tested << " colors" << ", cells: " << state->c->cells << "/" << g->v_size << std::endl;
            }
            return;
        }

        int color_score(sgraph* g, ir_state* state, int color) {
            test_set.reset();
            const int v     = state->get_coloring()->lab[color];
            const int d     = g->d[v];
            const int ept   = g->v[v];
            int non_triv_col_d = 1;
            for(int i = 0; i < d; ++i) {
                const int test_col = state->get_coloring()->vertex_to_col[g->e[ept + i]];
                if(!test_set.get(test_col)) {
                    non_triv_col_d += 1;
                    test_set.set(test_col);
                }
            }
            //return state->c->ptn[color] + d + non_triv_col_d;
            return state->get_coloring()->ptn[color] * non_triv_col_d;
        }

        int alt_score(sgraph* g, ir_state* state, int color) {
            test_set.reset();
            const int v     = state->get_coloring()->lab[color];
            const int d     = g->d[v];
            const int ept   = g->v[v];
            int non_triv_col_d = 1;
            for(int i = 0; i < d; ++i) {
                const int test_col = state->get_coloring()->vertex_to_col[g->e[ept + i]];
                if(!test_set.get(test_col)) {
                    non_triv_col_d += 1;
                    test_set.set(test_col);
                }
            }
            //return state->c->ptn[color] + d + non_triv_col_d;
            return non_triv_col_d;
        }

        int state_score(sgraph* g, ir_state* state) {
            return state->get_coloring()->cells;
        }

        splitmap find_splitmap(sgraph* g, coloring* c) {
            return splitmap();
        }
    };

    // small DFS IR-solver with limited failure, parallelizes according to given split-map
    // good at solving Tinhofer graphs, not good at pruning difficult combinatorial graphs
    class microdfs {
        int         fail_cnt = 0;
        int         fail_lim = 0;
        int         threads  = 1;
        refinement* R        = nullptr;
        selector*   S        = nullptr;
        splitmap*   SM       = nullptr;
    public:
        long double grp_sz   = 1.0;

        void setup(int fail_lim, int threads, refinement* R, selector* S, splitmap* SM) {
            this->R        = R;
            this->fail_lim = fail_lim;
            this->threads  = threads;
        }

        void  __attribute__ ((noinline)) touch_initial_colors(coloring* c, mark_set* touched_color) {
            int i = 0;
            while(i < c->lab_sz) {
                touched_color->set(i);
                i += c->ptn[i] + 1;
            }
        }

        std::pair<bool, bool>  __attribute__ ((noinline)) recurse_to_equal_leaf(sgraph* g, work_list* initial_colors, ir_state* state, int fails, work_list* automorphism, work_list* automorphism_supp) {
            if(fails == 0) {
                int depth = 0;
                bool prev_fail     = false;
                int  prev_fail_pos = -1;
                int cert_pos   = 0;

                while((size_t) state->base_pos < state->compare_base_color.size()) {
                    ++depth;
                    const int col = state->compare_base_color[state->base_pos];
                    const int col_sz = state->c->ptn[col] + 1;
                    if(col_sz < 2)
                        return {false, false};
                    const int ind_v = state->c->lab[col];
                    state->move_to_child(R, g, ind_v);
                    write_singleton_automorphism(&state->compare_singletons, &state->singletons, state->base_singleton_pt[state->base_singleton_pt.size()-1], state->singletons.size(),
                                                 automorphism->get_array(), automorphism_supp);
                    //bool found_auto = R->certify_automorphism_sparse(g, initial_colors->get_array(), automorphism->get_array(),
                    //                                                 automorphism_supp->cur_pos, automorphism_supp->get_array());

                    bool prev_cert = true;

                    if(prev_fail) {
                        prev_cert = R->check_single_failure(g, initial_colors->get_array(), automorphism->get_array(), prev_fail_pos);
                    }

                    //if(state->c->cells == g->v_size) {
                    if(prev_cert) {
                        auto cert_res = R->certify_automorphism_sparse_report_fail_resume(g, initial_colors->get_array(),
                                                                              automorphism->get_array(),
                                                                              automorphism_supp->cur_pos,
                                                                              automorphism_supp->get_array(), cert_pos);
                        cert_pos = std::get<2>(cert_res);
                        if (std::get<0>(cert_res)) {
                            cert_res = R->certify_automorphism_sparse_report_fail_resume(g, initial_colors->get_array(),
                                                                                         automorphism->get_array(),
                                                                                         automorphism_supp->cur_pos,
                                                                                         automorphism_supp->get_array(),
                                                                                         0);
                        }

                        if (std::get<0>(cert_res)) {
                            //std::cout << "deep_automorphism" << "@" << state->base_pos << "/" << state->compare_base_color.size() << "(" << state->c->cells << "/" << g->v_size << ")" << std::endl;
                            return {true, true};
                        } else {
                            prev_fail = true;
                            prev_fail_pos = std::get<1>(cert_res);
                        }
                    }
                }
                if(state->c->cells == g->v_size) {
                    //std::cout << "fail" << std::endl;
                    return {true, false};
                } else {
                    return {false, false};
                }
            }

            if((size_t) state->base_pos == state->compare_base_color.size()) {
                // at leaf?
                if(state->c->cells == g->v_size) {
                    return {true, false}; // check leaf first?
                } else {
                    return {false, false};
                }
            }

            const int col    = state->compare_base_color[state->base_pos];
            const int col_sz = state->c->ptn[col] + 1;
            if(col_sz < 2) {
                return {false, false};
            }

            for(int i = 0; i < col_sz; ++i) {
                const int ind_v = state->c->lab[col + i];
                state->move_to_child(R, g, ind_v);
                const auto res = recurse_to_equal_leaf(g, initial_colors, state, fails, automorphism, automorphism_supp);
                if(res.first)
                    return res;
                state->move_to_parent();
            }

            return {true, false};
        }

        void  color_diff_automorphism(int n, int* vertex_to_col, int* col_to_vertex, int* automorphism_map, work_list* automorphism_supp) {
            for(int v1 = 0; v1 < n; ++v1) {
                const int col = vertex_to_col[v1];
                const int v2  = col_to_vertex[col];
                if(v1 != v2) {
                    automorphism_map[v1] = v2;
                    automorphism_supp->push_back(v1);
                }
            }
        }

        // reset internal automorphism structure to the identity
        static void  reset_automorphism(int* rautomorphism, work_list* automorphism_supp) {
            for(int i = 0; i < automorphism_supp->cur_pos; ++i) {
                rautomorphism[(*automorphism_supp)[i]] = (*automorphism_supp)[i];
            }
            automorphism_supp->reset();
        };

        // adds a given automorphism to the tiny_orbit structure
        void  __attribute__ ((noinline)) add_automorphism_to_orbit(tiny_orbit* orbit, int* automorphism, int nsupp, int* supp) {
            for(int i = 0; i < nsupp; ++i) {
                orbit->combine_orbits(automorphism[supp[i]], supp[i]);
            }
        }

        void  __attribute__ ((noinline)) write_singleton_automorphism(std::vector<int>* singletons1, std::vector<int>* singletons2, int pos_start, int pos_end, int* rautomorphism, work_list* automorphism_supp) {
            for(int i = pos_start; i < pos_end; ++i) {
                const int from = (*singletons1)[i];
                const int to   = (*singletons2)[i];
                assert(rautomorphism[from] == from);
                if(from != to) {
                    automorphism_supp->push_back(from);
                    rautomorphism[from] = to;
                }
            }
        }

        trace compare_T;

        void  make_compare_snapshot(ir_state* state) {
            compare_T.set_compare(true);
            compare_T.set_record(false);
            compare_T.set_compare_trace(state->T);
            compare_T.set_position(state->T->get_position());
            state->T = &compare_T;
            state->compare_base_color.resize(state->base_color.size());
            std::copy(state->base_color.begin(), state->base_color.end(), state->compare_base_color.begin());
            state->compare_base_cells.resize(state->base_cells.size());
            std::copy(state->base_cells.begin(), state->base_cells.end(), state->compare_base_cells.begin());
            state->compare_singletons.resize(state->singletons.size());
            std::copy(state->singletons.begin(), state->singletons.end(), state->compare_singletons.begin());


            state->mode = IR_MODE_COMPARE;
        }

        // returns base-level reached (from leaf)
        int do_dfs(sgraph* g, coloring* c) {
            int gens = 0;
            work_list initial_colors(g->v_size);
            for(int i = 0; i < g->v_size; ++i) {
                initial_colors[i] = c->vertex_to_col[i];
            }

            tiny_orbit orbs;
            orbs.initialize(g->v_size);

            work_list automorphism(g->v_size);
            work_list automorphism_supp(g->v_size);

            for(int i = 0; i < g->v_size; ++i)
                automorphism[i] = i;

            selector Se;

            mark_set  touched_color(g->v_size);
            work_list touched_color_list(g->v_size);
            Se.empty_cache();
            S = &Se;
            work_list prev_color_list(g->v_size);

            trace T;
            T.set_compare(false);
            T.set_record(true);

            mark_set color_class(g->v_size);

            ir_state local_state;
            local_state.c = c;
            local_state.T = &T;
            local_state.touched_color = &touched_color;
            local_state.touched_color_list = &touched_color_list;
            local_state.prev_color_list = &prev_color_list;
            touch_initial_colors(c, &touched_color);

            // start DFS from a leaf!
            //move_to_leaf(g, &local_state);
            cellmagic wizard;
            wizard.find_base(R, g, &local_state);

            std::cout << "trace length: " << T.get_position() << std::endl;

            // make a snapshot of the leaf to compare to!
            make_compare_snapshot(&local_state);

            coloring leaf_color;
            leaf_color.copy(c);

            std::vector<int> individualize;

            bool fail = false;

            // loop that serves to optimize Tinhofer graphs
            while(local_state.base_pos > 0 && !fail) {
                local_state.move_to_parent();
                const int col    = local_state.base_color[local_state.base_pos]; // TODO: detect stack of "same color"?
                const int col_sz = local_state.c->ptn[col] + 1;
                const int vert   = local_state.base_vertex[local_state.base_pos];

                /*individualize.clear();
                for(int i = 0; i < col_sz; ++i) { // TODO: this is quadratic when it shouldn't be!
                    // TODO: do left-most, right-most, check sizes, stuff?
                    // TODO: I can save col and col_sz, and use compare snapshot leaf_color
                    const int v = local_state.c->lab[col + i];
                    if(v != vert && orbs.represents_orbit(v))
                        individualize.push_back(v);
                }*/

                //std::cout << "do_dfs@" << local_state.base_pos << ", col " <<  col << " sz " << col_sz << " con " << individualize.size() << std::endl;

                int count_leaf = 0;
                int count_orb  = 0;

                for(int i = col_sz - 1; i >= 0; --i) {
                    const int ind_v = leaf_color.lab[col + i];
                    if(ind_v == vert || !orbs.represents_orbit(ind_v))
                        continue;
                    if(orbs.are_in_same_orbit(ind_v, vert)) { // TODO somehow skip to ones not in same orbit?
                        ++count_orb;
                        continue;
                    }

                    ++count_leaf;

                    // call actual DFS with limited fails
                    // Tinhofer graphs will never fail, and have recursion width = 1
                    const int prev_base_pos = local_state.base_pos;
                    local_state.T->reset_trace_equal(); // probably need to re-adjust position?
                    local_state.move_to_child(R, g, ind_v);
                    write_singleton_automorphism(&local_state.compare_singletons, &local_state.singletons, local_state.base_singleton_pt[local_state.base_singleton_pt.size()-1], local_state.singletons.size(),
                                                 automorphism.get_array(), &automorphism_supp);
                    bool found_auto = R->certify_automorphism_sparse(g, initial_colors.get_array(), automorphism.get_array(),
                                                                  automorphism_supp.cur_pos, automorphism_supp.get_array());
                    assert(automorphism[vert] == ind_v);
                    //std::cout << found_auto << ", " << automorphism_supp.cur_pos << std::endl;
                    //reset_automorphism(automorphism.get_array(), &automorphism_supp);

                    // try proper recursion
                    if(!found_auto) {
                        auto rec_succeeded = recurse_to_equal_leaf(g, &initial_colors, &local_state, 0, &automorphism, &automorphism_supp);
                        found_auto = (rec_succeeded.first && rec_succeeded.second);
                        if (rec_succeeded.first && !rec_succeeded.second) {
                            reset_automorphism(automorphism.get_array(), &automorphism_supp);
                            color_diff_automorphism(g->v_size, local_state.c->vertex_to_col, leaf_color.lab,
                                                    automorphism.get_array(), &automorphism_supp);
                            found_auto = R->certify_automorphism_sparse(g, initial_colors.get_array(),
                                                                          automorphism.get_array(),
                                                                          automorphism_supp.cur_pos,
                                                                          automorphism_supp.get_array());
                        }
                    }

                    if (found_auto) {
                        assert(automorphism[vert] == ind_v);
                        ++gens;
                        add_automorphism_to_orbit(&orbs, automorphism.get_array(),
                                                  automorphism_supp.cur_pos, automorphism_supp.get_array());
                    }
                    reset_automorphism(automorphism.get_array(), &automorphism_supp);

                    while(prev_base_pos < local_state.base_pos) {
                        local_state.move_to_parent();
                    }

                    if(!found_auto) {
                        std::cout << ind_v << "(F) " << std::endl;
                        fail = true;
                        break;
                    }

                    if(orbs.orbit_size(vert) == col_sz) {
                        break;
                    }
                }

                if(!fail) {
                    grp_sz *= col_sz;
                }

                //std::cout << "leaf/orb: " << count_leaf << "/" << count_orb << std::endl;
            }

            //grp_sz = 10;
            std::cout << "dfs: gens " << gens << ", levels covered " << local_state.base_pos << "-" << local_state.compare_base_color.size() << ", grp_sz covered: " << grp_sz << std::endl;
            return local_state.base_pos;
        }
    };
}

#endif //DEJAVU_MICRODFS_H
