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
    };

    struct splitmap {

    };

    // heuristics to attempt to find a good selector, as well as creating a split-map for microdfs
    class cellmagic {

        int locked_lim = 512;

        // move this to static functions
        void  move_to_child(refinement* R, sgraph* g, ir_state* state, int v) {
            if (state->c->cells == g->v_size) {
                return;
            }

            ++state->base_pos;
            state->base_vertex.push_back(v);
            state->base_color.push_back(state->c->vertex_to_col[v]);
            const int touched_col_prev = state->touched_color_list->cur_pos;
            state->base_touched_color_list_pt.push_back(touched_col_prev);
            state->base_singleton_pt.push_back(state->singletons.size());

            const int init_color_class = R->individualize_vertex(state->c, v, state->T, state->touched_color,
                                                                 state->touched_color_list, state->prev_color_list, &state->singletons);

            if(state->mode == IR_MODE_RECORD) {
                R->refine_coloring(g, state->c, init_color_class, nullptr, -1, -1, nullptr,
                                   state->touched_color, state->touched_color_list, state->prev_color_list, state->T, &state->singletons);
            } else {
                R->refine_coloring(g, state->c, init_color_class, nullptr, state->compare_base_cells[state->base_pos-1], -1, nullptr,
                                   state->touched_color, state->touched_color_list, state->prev_color_list, state->T, &state->singletons);
            }

            state->base_cells.push_back(state->c->cells);

        }
        void  move_to_parent(ir_state* state) {
            // unwind invariant
            state->T->rewind_to_individualization();

            --state->base_pos;
            while(state->prev_color_list->cur_pos > state->base_touched_color_list_pt[state->base_pos]) {
                const int old_color = state->prev_color_list->pop_back();
                const int new_color = state->touched_color_list->pop_back();

                state->touched_color->unset(new_color);

                const int new_color_sz = state->c->ptn[new_color] + 1;
                state->c->ptn[old_color] += new_color_sz;
                state->c->ptn[new_color]  = 1;

                for(int j = 0; j < new_color_sz; ++j) {
                    const int v = state->c->lab[new_color + j];
                    state->c->vertex_to_col[v] = old_color;
                    assert(state->c->vertex_to_lab[v] == new_color + j);
                }

                --state->c->cells;
            }

            int const new_singleton_pos = state->base_singleton_pt.back();
            state->singletons.resize(new_singleton_pos);

            state->base_vertex.pop_back();
            state->base_color.pop_back();
            state->base_touched_color_list_pt.pop_back();
            state->base_cells.pop_back();
            state->base_singleton_pt.pop_back();
        }
    public:
        mark_set test_set;

        void  __attribute__ ((noinline)) find_base(refinement* R, sgraph* g, ir_state* state) {
            std::vector<int> candidates;
            test_set.initialize(g->v_size);
            candidates.reserve(locked_lim);
            while(state->c->cells != g->v_size) {
                candidates.clear();
                int best_color = -1;
                int best_score = -1;

                int alt_best_color = -1;
                int alt_best_score = -1;

                for(int i = 0; i < state->c->ptn_sz;) {
                    if(state->c->ptn[i] > 0) {
                        candidates.push_back(i);
                    }

                    //if(candidates.size() >= (size_t) locked_lim)
                    //    break;

                    i += state->c->ptn[i] + 1;
                }

                const int num_tested = candidates.size();

                while(!candidates.empty()) {
                    const int test_color = candidates.back();
                    candidates.pop_back();

                    const int test_score = color_score(g, state, test_color);
                    if(test_score > best_score) {
                        best_color = test_color;
                        best_score = test_score;
                    }

                    const int alt_test_score = alt_score(g, state, test_color);
                    if(alt_test_score > alt_best_score) {
                        alt_best_color = test_color;
                        alt_best_score = test_score;
                    }
                }

                if(best_color != alt_best_color) {
                    //std::cout << best_color << " vs. " << alt_best_color << std::endl;
                    move_to_child(R, g, state, state->c->lab[best_color]);
                    const int score1 = state_score(g, state);
                    move_to_parent(state);
                    move_to_child(R, g, state, state->c->lab[alt_best_color]);
                    const int score2 = state_score(g, state);
                    move_to_parent(state);

                    //std::cout << score1 << ":" << score2 << std::endl;
                    if(score2 > score1) {
                        best_color = alt_best_color;
                    }
                }


                move_to_child(R, g, state, state->c->lab[best_color]);
                //std::cout << "picked color " << best_color << " score " << best_score << " among " << num_tested << " colors" << ", cells: " << state->c->cells << "/" << g->v_size << std::endl;
            }
            return;
        }

        int color_score(sgraph* g, ir_state* state, int color) {
            test_set.reset();
            const int v     = state->c->lab[color];
            const int d     = g->d[v];
            const int ept   = g->v[v];
            int non_triv_col_d = 1;
            for(int i = 0; i < d; ++i) {
                const int test_col = state->c->vertex_to_col[g->e[ept + i]];
                if(!test_set.get(test_col)) {
                    non_triv_col_d += 1;
                    test_set.set(test_col);
                }
            }
            //return state->c->ptn[color] + d + non_triv_col_d;
            return state->c->ptn[color] * non_triv_col_d;
        }

        int alt_score(sgraph* g, ir_state* state, int color) {
            test_set.reset();
            const int v     = state->c->lab[color];
            const int d     = g->d[v];
            const int ept   = g->v[v];
            int non_triv_col_d = 1;
            for(int i = 0; i < d; ++i) {
                const int test_col = state->c->vertex_to_col[g->e[ept + i]];
                if(!test_set.get(test_col)) {
                    non_triv_col_d += 1;
                    test_set.set(test_col);
                }
            }
            //return state->c->ptn[color] + d + non_triv_col_d;
            return non_triv_col_d;
        }

        int state_score(sgraph* g, ir_state* state) {
            return state->c->cells;
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
        void setup(int fail_lim, int threads, refinement* R, selector* S, splitmap* SM) {
            this->R        = R;
            this->fail_lim = fail_lim;
            this->threads  = threads;
        }

        void  move_to_leaf(sgraph* g, ir_state* state) {
            while(state->c->cells != g->v_size) {
                move_to_child(g, state, -1);
            }
        }

        void   __attribute__ ((noinline)) move_to_child(sgraph* g, ir_state* state, int v) {
            if (state->c->cells == g->v_size) {
                return;
            }
            ++state->base_pos;
            if(v == -1) {
                const int color_class = S->select_color_largest(state->c);
                v = state->c->lab[color_class];
            }
            state->base_singleton_pt.push_back(state->singletons.size());
            state->base_vertex.push_back(v);
            state->base_color.push_back(state->c->vertex_to_col[v]);
            state->base_touched_color_list_pt.push_back(state->touched_color_list->cur_pos);
            const int init_color_class = R->individualize_vertex(state->c, v, state->T, state->touched_color,
                                                                 state->touched_color_list, state->prev_color_list, &state->singletons);

            if(state->mode == IR_MODE_RECORD) {
                R->refine_coloring(g, state->c, init_color_class, nullptr, -1, -1, nullptr,
                                   state->touched_color, state->touched_color_list, state->prev_color_list, state->T, &state->singletons);
            } else {
                R->refine_coloring(g, state->c, init_color_class, nullptr, state->compare_base_cells[state->base_pos-1], -1, nullptr,
                                   state->touched_color, state->touched_color_list, state->prev_color_list, state->T, &state->singletons);
            }

            state->base_cells.push_back(state->c->cells);
        }

        void   __attribute__ ((noinline)) move_to_parent(ir_state* state) {
            // unwind invariant
            state->T->rewind_to_individualization();

            --state->base_pos;
            while(state->prev_color_list->cur_pos > state->base_touched_color_list_pt[state->base_pos]) {
                const int old_color = state->prev_color_list->pop_back();
                const int new_color = state->touched_color_list->pop_back();

                state->touched_color->unset(new_color);

                const int new_color_sz = state->c->ptn[new_color] + 1;
                state->c->ptn[old_color] += new_color_sz;
                state->c->ptn[new_color]  = 1;

                for(int j = 0; j < new_color_sz; ++j) {
                    const int v = state->c->lab[new_color + j];
                    state->c->vertex_to_col[v] = old_color;
                    assert(state->c->vertex_to_lab[v] == new_color + j);
                }

                --state->c->cells;
            }
            int const new_singleton_pos = state->base_singleton_pt.back();
            state->singletons.resize(new_singleton_pos);

            state->base_singleton_pt.pop_back();
            state->base_vertex.pop_back();
            state->base_color.pop_back();
            state->base_touched_color_list_pt.pop_back();
            state->base_cells.pop_back();
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
                    move_to_child(g, state, ind_v);
                    write_singleton_automorphism(&state->compare_singletons, &state->singletons, state->base_singleton_pt[state->base_singleton_pt.size()-1], state->singletons.size(),
                                                 automorphism->get_array(), automorphism_supp);
                    //bool found_auto = R->certify_automorphism_sparse(g, initial_colors->get_array(), automorphism->get_array(),
                    //                                                 automorphism_supp->cur_pos, automorphism_supp->get_array());

                    bool prev_cert = true;

                    if(prev_fail) {
                        prev_cert = R->check_single_failure(g, initial_colors->get_array(), automorphism->get_array(), prev_fail_pos);
                    }

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
                            //std::cout << "success" << "@" << state->base_pos << "/" << state->compare_base_color.size() << std::endl;
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
                move_to_child(g, state, ind_v);
                const auto res = recurse_to_equal_leaf(g, initial_colors, state, fails, automorphism, automorphism_supp);
                if(res.first)
                    return res;
                move_to_parent(state);
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
        void  __attribute__ ((noinline)) add_automorphism_to_orbit(tiny_orbit* orbit, mark_set* interest, int* automorphism, int nsupp, int* supp) {
            for(int i = 0; i < nsupp; ++i) {
                //if(interest->get(supp[i]))
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
                move_to_parent(&local_state);
                //orbs.reset(); // TODO no need to reset when going up in tree!
                const int col    = local_state.base_color[local_state.base_pos];
                const int col_sz = local_state.c->ptn[col] + 1;
                const int vert   = local_state.base_vertex[local_state.base_pos];
                //std::cout << "do_dfs@" << local_state.base_pos << ", col " <<  col << " sz " << col_sz << std::endl;

                color_class.reset();
                individualize.clear();
                for(int i = 0; i < col_sz; ++i) {
                    const int v = local_state.c->lab[col + i];
                    color_class.set(v);
                    individualize.push_back(v);
                }

                int count_leaf = 0;
                int count_orb  = 0;

                for(auto ind_v : individualize) {
                    if(orbs.are_in_same_orbit(ind_v, vert)) {
                        ++count_orb;
                        continue;
                    }

                    ++count_leaf;

                    // call actual DFS with limited fails
                    // Tinhofer graphs will never fail, and have recursion width = 1
                    const int prev_base_pos = local_state.base_pos;
                    local_state.T->reset_trace_equal(); // probably need to re-adjust position?
                    move_to_child(g, &local_state, ind_v);
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
                        add_automorphism_to_orbit(&orbs, &color_class, automorphism.get_array(),
                                                  automorphism_supp.cur_pos, automorphism_supp.get_array());
                    }

                    reset_automorphism(automorphism.get_array(), &automorphism_supp);

                    while(prev_base_pos < local_state.base_pos) {
                        move_to_parent(&local_state);
                    }

                    if(!found_auto) {
                        std::cout << ind_v << "(F) " << std::endl;
                        fail = true;
                        break;
                    }
                }
                //std::cout << "leaf/orb: " << count_leaf << "/" << count_orb << std::endl;
            }

            std::cout << "dfs: gens " << gens << ", levels covered " << local_state.base_pos << "-" << local_state.compare_base_color.size() << std::endl;
            return local_state.base_pos;
        }
    };
}

#endif //DEJAVU_MICRODFS_H
