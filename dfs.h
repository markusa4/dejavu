#ifndef DEJAVU_DFS_H
#define DEJAVU_DFS_H

#include <random>
#include <chrono>
#include "sgraph.h"
#include "invariant.h"
#include "selector.h"
#include "bfs.h"
#include "schreier_sequential.h"

namespace dejavu {

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// UTILITY                                                                                                     ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * The ir_mode determines in which mode the trace is used: whether a new trace is recorded, or whether the current
     * computation is compared to a stored trace.
     */
    enum ir_mode { IR_MODE_COMPARE_TRACE, IR_MODE_RECORD_TRACE, IR_MODE_RECORD_HASH};

    static void color_diff_automorphism(int n, int* vertex_to_col, int* col_to_vertex,
                                        int* automorphism_map, work_list* automorphism_supp) {
        for(int v1 = 0; v1 < n; ++v1) {
            const int col = vertex_to_col[v1];
            const int v2  = col_to_vertex[col];
            if(v1 != v2) {
                automorphism_map[v1] = v2;
                automorphism_supp->push_back(v1);
            }
        }
    }


    // TODO enable a "compressed state" only consisting of base (but code outside should be oblivious to it)
    // TODO this is only supposed to be an "incomplete" state -- should there be complete states?
    /**
     * Using this class, partial information of a state of an IR computation can be stored. Using this information,
     * IR computations can be resumed from this state either using BFS or random walks. The state in particular does not
     * keep enough information to resume using DFS.
     */
    class ir_reduced_save {
        std::vector<int> base_vertex;
        coloring c;
        long     invariant = 0;
    public:
        void set_state(std::vector<int>& base_vertex, coloring& c, long invariant) {
            this->base_vertex = base_vertex;
            this->c.copy(&c);
            this->invariant = invariant;
        }

        coloring* get_coloring() {
            return &c;
        }

        long get_invariant_hash() {
            return invariant;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// IR CONTROL                                                                                                  ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    struct ir_controller {
        coloring*    c = nullptr;
        trace*       T = nullptr;
        mark_set     touched_color;
        work_list    touched_color_list;
        work_list    prev_color_list;
        std::vector<int>  singletons;
        std::vector<int>  base_vertex;
        std::vector<int>  base_color;
        std::vector<int>  base_touched_color_list_pt;
        std::vector<int>  base_cells;
        std::vector<int>  base_singleton_pt;

        std::vector<int>  compare_base_color;
        std::vector<int>  compare_base_cells;
        std::vector<int>  compare_singletons;

        std::function<type_split_color_hook>    my_split_hook;
        std::function<type_worklist_color_hook> my_worklist_hook;

        // settings
        ir_mode mode = IR_MODE_RECORD_TRACE;

        // heuristics
        bool h_last_refinement_singleton_only = true;
        int  h_hint_color                     = -1;
        bool h_hint_color_is_singleton_now    = true;
        bool h_cell_active                    = false;
        bool h_individualize                  = false;

        int base_pos = 0;
    private:
        void touch_initial_colors() {
            int i = 0;
            while(i < c->lab_sz) {
                touched_color.set(i);
                i += c->ptn[i] + 1;
            }
        }

    public:
        ir_controller(coloring* c, trace* T) {
            this->c = c;
            this->T = T;

            touched_color.initialize(c->lab_sz);
            touched_color_list.initialize(c->lab_sz);
            prev_color_list.initialize(c->lab_sz);

            touch_initial_colors();

            my_split_hook    = self_split_hook();
            my_worklist_hook = self_worklist_hook();
        }

        // TODO incomplete state save for BFS & random reset
        void save_state(ir_reduced_save* state) {
            state->set_state(base_vertex, *c, T->get_hash());
        }

        // TODO incomplete state load for BFS & random reset
        void load_state(ir_reduced_save* state) {
            c->copy(state->get_coloring());
            T->set_hash(state->get_invariant_hash());
        }

        coloring* get_coloring() {
            return c;
        }

        bool split_hook(const int old_color, const int new_color, const int new_color_sz) {
            // update some heuristic values
            if(new_color_sz > 1) {
                h_last_refinement_singleton_only = false;
                h_hint_color = new_color;
                h_hint_color_is_singleton_now = false;
            }
            if(new_color == h_hint_color && new_color_sz == 1) {
                h_hint_color_is_singleton_now = true;
            }

            // write singletons to singleton list
            if(new_color_sz == 1) {
                singletons.push_back(c->lab[new_color]);
            }

            // record colors that were changed
            if(!touched_color.get(new_color)) {
                touched_color.set(new_color);
                prev_color_list.push_back(old_color);
                touched_color_list.push_back(new_color);
            }

            // record split into trace invariant, unless we are individualizing
            if(T && !h_individualize) T->op_refine_cell_record(new_color, new_color_sz, 1);

            return true;
        }

        bool worklist_hook(const int color, const int color_sz) {
            if(h_cell_active) {
                if(T) T->op_refine_cell_end();
                h_cell_active = false;
            }

            // update some heuristic values
            if(T) {
                if(!T->blueprint_is_next_cell_active()) {
                    if(config.CONFIG_IR_IDLE_SKIP) {
                        T->blueprint_skip_to_next_cell();
                        return false;
                    }
                }
            }

            if(T) T->op_refine_cell_start(color);

            h_cell_active = true;

            return true;
        }

        std::function<type_split_color_hook> self_split_hook() {
            return std::bind(&ir_controller::split_hook, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        }

        std::function<type_worklist_color_hook> self_worklist_hook() {
            return std::bind(&ir_controller::worklist_hook, this, std::placeholders::_1, std::placeholders::_2);
        }

        // move this to static functions
        void move_to_child(refinement* R, sgraph* g, int v) {
            ++base_pos;
            base_singleton_pt.push_back(singletons.size());
            base_vertex.push_back(v);
            base_color.push_back(c->vertex_to_col[v]);
            base_touched_color_list_pt.push_back(touched_color_list.cur_pos);

            assert(!h_cell_active);

            h_individualize = true;
            const int init_color_class = R->individualize_vertex(c, v, my_split_hook);
            T->op_individualize(c->vertex_to_col[v]);
            h_individualize = false;

            h_last_refinement_singleton_only = true;

            if(T) T->op_refine_start();

            if(mode == IR_MODE_RECORD_TRACE) {
                R->refine_coloring(g, c, init_color_class, -1, my_split_hook, my_worklist_hook);
                if(T && h_cell_active) T->op_refine_cell_end();
                if(T) T->op_refine_end();
            } else {
                R->refine_coloring(g, c, init_color_class, compare_base_cells[base_pos-1], my_split_hook, my_worklist_hook);
                if(T) T->skip_to_individualization();
            }

            h_cell_active = false;

            base_cells.push_back(c->cells);
        }

        void   __attribute__ ((noinline)) move_to_parent() {
            // unwind invariant
            if(T) T->rewind_to_individualization();

            --base_pos;
            while(prev_color_list.cur_pos > base_touched_color_list_pt[base_pos]) {
                const int old_color = prev_color_list.pop_back();
                const int new_color = touched_color_list.pop_back();

                touched_color.unset(new_color);

                const int new_color_sz = c->ptn[new_color] + 1;
                c->ptn[old_color] += new_color_sz;
                c->ptn[new_color]  = 1;

                if(c->ptn[old_color] > 0) {
                    h_hint_color = old_color;
                    h_hint_color_is_singleton_now = false;
                }

                for(int j = 0; j < new_color_sz; ++j) {
                    const int v = c->lab[new_color + j];
                    c->vertex_to_col[v] = old_color;
                    assert(c->vertex_to_lab[v] == new_color + j);
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

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// SELECTOR                                                                                                    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // heuristics to attempt to find a good selector, as well as creating a split-map for dfs_ir
    // TODO: enable different selector strategies for restarts
    // TODO: make a grab-able hook for dynamic selector which can also deviate from base color (blueprint selector)
    class selector_factory {
        int locked_lim = 512;

    public:
        mark_set test_set;

        void  __attribute__ ((noinline)) find_base(refinement* R, sgraph* g, ir_controller* state) {
            state->T->set_compare(false);
            state->T->set_record(true);

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

                        if(candidates.size() >= (size_t) locked_lim)
                            break;

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

                        /*int alt_test_score = alt_score(g, state, test_color);
                        if(neighbour_color.get(test_color)) {
                            alt_test_score *= 10;
                        }
                        if (alt_test_score > alt_best_score) {
                            alt_best_color = test_color;
                            alt_best_score = test_score;
                        }*/
                    }

                    /*if (best_color != alt_best_color) {
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
                    }*/
                }

                if(best_color == -1) {
                    std::cout << state->get_coloring()->cells << "/" << g->v_size << std::endl;
                 }

                assert(best_color >= 0);
                assert(best_color < g->v_size);
                prev_color = best_color;
                state->move_to_child(R, g, state->get_coloring()->lab[best_color]);
                //std::cout << "picked color " << best_color << " score " << best_score << " among " << num_tested << " colors" << ", cells: " << state->c->cells << "/" << g->v_size << std::endl;
            }
            return;
        }

        int color_score(sgraph* g, ir_controller* state, int color) {
            //return (state->c->ptn[color] +1)*100; //+ d + non_triv_col_d;

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

        int alt_score(sgraph* g, ir_controller* state, int color) {
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

        int state_score(sgraph* g, ir_controller* state) {
            return state->get_coloring()->cells;
        }

        splitmap find_splitmap(sgraph* g, coloring* c) {
            return splitmap();
        }
    };

    // TODO tree structure for BFS + random walk
    class ir_tree {

    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// SCHREIER-SIMS                                                                                               ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // TODO implement sparse schreier-sims for fixed base
    // TODO support for sparse tables, sparse automorphism, maybe mix sparse & dense
    class schreier {

    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// BFS                                                                                                         ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // TODO implement new, simpler bfs
    // TODO depends on ir_tree, and selector
    class bfs_ir {

    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// RANDOM IR                                                                                                   ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // TODO implement dejavu strategy, more simple
    // TODO depends on ir_tree, selector, and given base (no need to sift beyond base!)
    class random_ir {

    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// DFS                                                                                                         ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // DFS IR which does not backtrack, parallelizes according to given split-map
    // not able to solve difficult combinatorial graphs
    // TODO should depend on given selector
    class dfs_ir {
        int         fail_cnt = 0;
        int         threads  = 1;
        refinement* R        = nullptr;
        splitmap*   SM       = nullptr;
    public:
        long double grp_sz_man   = 1.0;
        int grp_sz_exp = 0;

        void setup(int threads, refinement* R, selector* S) {
            this->R        = R;
            this->threads  = threads;
        }

        std::pair<bool, bool> recurse_to_equal_leaf(sgraph* g, int* initial_colors, ir_controller* state, work_list* automorphism, work_list* automorphism_supp) {
            bool prev_fail     = false;
            int  prev_fail_pos = -1;
            int cert_pos   = 0;

            while((size_t) state->base_pos < state->compare_base_color.size()) {
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

                assert(state->h_hint_color_is_singleton_now?state->h_last_refinement_singleton_only:true);

                if(prev_fail && state->h_last_refinement_singleton_only && state->h_hint_color_is_singleton_now) {
                    prev_cert = R->check_single_failure(g, initial_colors, automorphism->get_array(), prev_fail_pos);
                }

                //if(state->c->cells == g->v_size) {
                if(prev_cert && state->h_last_refinement_singleton_only && state->h_hint_color_is_singleton_now) {
                    // TODO: add better heuristic to not always do this check, too expensive!
                    auto cert_res = R->certify_automorphism_sparse_report_fail_resume(g, initial_colors,
                                                                          automorphism->get_array(),
                                                                          automorphism_supp->cur_pos,
                                                                          automorphism_supp->get_array(), cert_pos);
                    cert_pos = std::get<2>(cert_res);
                    if (std::get<0>(cert_res)) {
                        cert_res = R->certify_automorphism_sparse_report_fail_resume(g, initial_colors,
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

        // reset internal automorphism structure to the identity
        static void  reset_automorphism(int* rautomorphism, work_list* automorphism_supp) {
            for(int i = 0; i < automorphism_supp->cur_pos; ++i) {
                rautomorphism[(*automorphism_supp)[i]] = (*automorphism_supp)[i];
            }
            automorphism_supp->reset();
        };

        // adds a given automorphism to the tiny_orbit structure
        void  add_automorphism_to_orbit(tiny_orbit* orbit, int* automorphism, int nsupp, int* supp) {
            for(int i = 0; i < nsupp; ++i) {
                orbit->combine_orbits(automorphism[supp[i]], supp[i]);
            }
        }

        void  write_singleton_automorphism(std::vector<int>* singletons1, std::vector<int>* singletons2, int pos_start, int pos_end, int* rautomorphism, work_list* automorphism_supp) {
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

        void  make_compare_snapshot(ir_controller* state) {
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


            state->mode = IR_MODE_COMPARE_TRACE;
        }

        // returns base-level reached (from leaf)
        int do_dfs(sgraph* g, int* initial_colors, ir_controller& local_state) {
            int gens = 0;

            tiny_orbit orbs;
            orbs.initialize(g->v_size);

            work_list automorphism(g->v_size);
            work_list automorphism_supp(g->v_size);

            for(int i = 0; i < g->v_size; ++i)
                automorphism[i] = i;

            mark_set color_class(g->v_size);

            //ir_controller local_state(c, &T);

            // start DFS from a leaf!
            //move_to_leaf(g, &local_state);
            // TODO move this outside DFS, only give state in leaf node (we dont even need a selector here!)

            // make a snapshot of the leaf to compare to!
            make_compare_snapshot(&local_state);

            coloring leaf_color;
            leaf_color.copy(local_state.get_coloring());

            std::vector<int> individualize;

            bool fail = false;

            // loop that serves to optimize Tinhofer graphs
            while(local_state.base_pos > 0 && !fail) {
                local_state.move_to_parent();
                const int col    = local_state.base_color[local_state.base_pos]; // TODO: detect stack of "same color"?
                const int col_sz = local_state.c->ptn[col] + 1;
                const int vert   = local_state.base_vertex[local_state.base_pos];

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
                    bool found_auto = false;
                    //if(local_state.h_last_refinement_singleton_only) {
                    write_singleton_automorphism(&local_state.compare_singletons, &local_state.singletons,
                                                 local_state.base_singleton_pt[
                                                         local_state.base_singleton_pt.size() - 1],
                                                 local_state.singletons.size(),
                                                 automorphism.get_array(), &automorphism_supp);
                    found_auto = R->certify_automorphism_sparse(g, initial_colors,
                                                                     automorphism.get_array(),
                                                                     automorphism_supp.cur_pos,
                                                                     automorphism_supp.get_array());
                    assert(automorphism[vert] == ind_v);
                    //}
                    //std::cout << found_auto << ", " << automorphism_supp.cur_pos << std::endl;
                    //reset_automorphism(automorphism.get_array(), &automorphism_supp);

                    // try proper recursion
                    if(!found_auto) {
                        auto rec_succeeded = recurse_to_equal_leaf(g, initial_colors, &local_state, &automorphism, &automorphism_supp);
                        found_auto = (rec_succeeded.first && rec_succeeded.second);
                        if (rec_succeeded.first && !rec_succeeded.second) {
                            reset_automorphism(automorphism.get_array(), &automorphism_supp);
                            color_diff_automorphism(g->v_size, local_state.c->vertex_to_col, leaf_color.lab,
                                                    automorphism.get_array(), &automorphism_supp);
                            found_auto = R->certify_automorphism_sparse(g, initial_colors,
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
                    grp_sz_man *= col_sz;
                    while(grp_sz_man > 10) {
                        grp_sz_man /= 10;
                        grp_sz_exp += 1;
                    }
                }
            }

            std::cout << "dfs_ir: gens " << gens << ", levels covered " << local_state.base_pos << "-" << local_state.compare_base_color.size() << ", grp_sz covered: " << grp_sz_man << "*10^" << grp_sz_exp << std::endl;
            return local_state.base_pos;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// MAIN ALGORITHM                                                                                              ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // TODO high-level strategy
    //  - full restarts, change selector, but make use of previously computed automorphisms (maybe inprocess) -- goal:
    //    stable, single-threaded performance
    //  - exploit strengths of DFS, random walks, BFS, and their synergies
    //  - try to use no-pruning DFS to fill up available leafs more efficiently

    /**
     * The dejavu solver.
     */
    class dejavu2 {
    private:
        // high-level modules of the algorithm
        sassy::preprocessor m_prep;   /**< preprocessor */
        dfs_ir    m_dfs;              /**< depth-first search */
        bfs_ir    m_bfs;              /**< breadth-first search */
        random_ir m_rand;             /**< randomized search */

        // utility tools used by other modules
        refinement m_refinement;      /**< workspace for color refinement and other utilities */
        selector_factory m_selectors; /**< cell selector creation */
        schreier  m_schreier;         /**< Schreier-Sims algorithm */
        ir_tree   m_tree;             /**< IR-tree */

        // TODO: should not be necessary in the end!
        void transfer_sgraph_to_sassy_sgraph(sgraph* g, sassy::sgraph* gg) {
            gg->v = g->v;
            gg->d = g->d;
            gg->e = g->e;
            gg->v_size = g->v_size;
            gg->d_size = g->d_size;
            gg->e_size = g->e_size;
        }
        // TODO: should not be necessary in the end!
        void transfer_sassy_sgraph_to_sgraph(sgraph* g, sassy::sgraph* gg) {
            g->v = gg->v;
            g->d = gg->d;
            g->e = gg->e;
            g->v_size = gg->v_size;
            g->d_size = gg->d_size;
            g->e_size = gg->e_size;
        }

    public:
        /**
         * Compute the automorphisms of the graph \p g colored with vertex colors \p colmap. Automorphisms are returned
         * using the function pointer \p hook.
         * @param g The graph.
         * @param colmap The vertex coloring of \p g. A null pointer is admissible as the trivial coloring.
         * @param hook The hook used for returning automorphisms. A null pointer is admissible if this is not needed.
         *
         * @todo return automorphism group size
         */
        void automorphisms(sgraph* g, int* colmap = nullptr, dejavu_hook* hook = nullptr) {
            std::cout << "dejavu2" << std::endl;

            // TODO facilities for restarts
            
            // preprocess the graph using sassy
            sassy::sgraph gg;
            transfer_sgraph_to_sassy_sgraph(g, &gg);
            m_prep.reduce(&gg, colmap, hook);
            transfer_sassy_sgraph_to_sgraph(g, &gg);

            // initialize a coloring using colors of preprocessed graph
            coloring local_coloring;
            g->initialize_coloring(&local_coloring, colmap);

            // setup a local state for IR computations
            trace local_trace;
            ir_controller local_state(&local_coloring, &local_trace);

            // find a selector, moves local_state to a leaf
            m_selectors.find_base(&m_refinement, g, &local_state);
            std::cout << "trace length: " << local_trace.get_position() << std::endl;

            // depth-first search starting from the computed leaf in local_state
            m_dfs.setup(0, &m_refinement, nullptr);
            const int dfs_reached_level = m_dfs.do_dfs(g, colmap, local_state);
            if(dfs_reached_level == 0) {
                /*long double add_grp_sz = D.grp_sz_man;
                while(add_grp_sz > 10) {
                    a.grp_sz_exp += 1;
                    add_grp_sz = add_grp_sz / 10;
                }
                a.grp_sz_exp += D.grp_sz_exp;
                a.grp_sz_man *= add_grp_sz;*/
                std::cout << "DFS finished graph" << std::endl;
                return;
            }

            // intertwined random automorphisms and breadth-first search

            // random automorphisms

            // breadth-first search

        }
    };
}

#endif //DEJAVU_DFS_H
