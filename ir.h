// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_IR_H
#define DEJAVU_IR_H

#include "refinement.h"
#include "coloring.h"
#include "sgraph.h"
#include "trace.h"

namespace dejavu {
    /**
     * \brief IR fundamentals.
     *
     * Contains fundamental algorithms and data structures to implement individualization-refinement algorithms. This
     * includes graphs, colorings, color refinement, cell selectors as well as higher level control mechanisms.
     */
    namespace ir {
        /**
         * \brief Mode of trace for IR search
         *
         * The ir_mode determines in which mode the trace is used: whether a new trace is recorded, or whether the current
         * computation is compared to a stored trace.
         *
         */
        enum ir_mode {
            IR_MODE_COMPARE_TRACE, IR_MODE_RECORD_TRACE, IR_MODE_RECORD_HASH_IRREVERSIBLE
        };

        /**
         * int type_selector_hook(coloring* c, const int base_pos);
         */
        typedef int type_selector_hook(const coloring *, const int);

        /**
         * \brief Reduced IR save state
         *
         * Using this class, partial information of a state of an IR computation can be stored. Using this information,
         * IR computations can be resumed from this state either using BFS or random walks. The state in particular does not
         * keep enough information to resume using DFS.
         */
        class reduced_save {
            // TODO enable a "compressed state" only consisting of base (but code outside should be oblivious to it)
            // TODO this is only supposed to be an "incomplete" state -- should there be complete states?

            std::vector<int> base_vertex; /**< base of vertices of this IR node (optional) */
            coloring c; /**< vertex coloring of this IR node */
            long invariant = 0; /**< hash of invariant of this IR node */
            int trace_position = 0; /**< position of trace of this IR node */
            int base_position = 0; /**< length of base of this IR node */
        public:
            void set_state(std::vector<int> &base_vertex, coloring &c, long invariant, int trace_position,
                           int base_position) {
                this->base_vertex = base_vertex;
                this->c.copy_force(&c);
                this->invariant = invariant;
                this->trace_position = trace_position;
                this->base_position = base_position;
            }

            coloring *get_coloring() {
                return &c;
            }

            /**
             * @return hash of invariant of this IR node
             */
            long get_invariant_hash() {
                return invariant;
            }

            /**
             * @return position of trace of this IR node
             */
            int get_trace_position() {
                return trace_position;
            }

            /**
             * @return length of base of this IR node
             */
            int get_base_position() {
                return base_position;
            }

            /**
             * @return base of this IR node
             */
            std::vector<int>& get_base() {
                return base_vertex;
            }
        };

        /**
         * \brief Controls movement in IR tree
         *
         * Keeps a state of an IR node. Enables the movement to a child of the node, or to the parent of the node.
         * The controller manages data structures and functions, which facilitate the trace \a T as well as the reversal of
         * color refinement.
         *
         * Has different modes (managed by \a mode) depending on whether color refinement should be reversible or not.
         */
        struct controller {
            coloring      *c  = nullptr;
            trace         *T  = nullptr;
            trace         *cT = nullptr;

            trace _T1;
            trace _T2;

            refinement* R;


            mark_set touched_color;
            work_list touched_color_list;
            work_list prev_color_list;
            std::vector<int> singletons;
            std::vector<int> base_vertex;
            std::vector<int> base_color;
            std::vector<int> base_color_sz;
            std::vector<int> base_touched_color_list_pt;
            std::vector<int> base_cells;
            std::vector<int> base_singleton_pt;

            std::vector<int> compare_base_color;
            std::vector<int> compare_base_cells;
            std::vector<int> compare_singletons;
            std::vector<int> compare_base;

            coloring leaf_color;

            std::function<type_split_color_hook> my_split_hook;
            std::function<type_worklist_color_hook> my_worklist_hook;
            std::function<type_additional_info_hook> my_add_hook;

            work_list workspace;

            // settings
            ir_mode mode = IR_MODE_RECORD_TRACE;

            // heuristics
            bool h_last_refinement_singleton_only = true;
            int h_hint_color = -1;
            bool h_hint_color_is_singleton_now = true;
            bool h_cell_active = false;
            bool h_individualize = false;
            bool h_trace_early_out = false;

            bool  h_deviation_inc_active  = false;
            int   h_deviation_inc         = 48;
            int   h_deviation_inc_current = 0;

            int base_pos = 0;
        private:
            void touch_initial_colors() {
                int i = 0;
                while (i < c->lab_sz) {
                    touched_color.set(i);
                    i += c->ptn[i] + 1;
                }
            }

            void reset_touched() {
                touched_color.reset();
                touched_color_list.reset();
                prev_color_list.reset();
                touch_initial_colors();
            }

        public:

            void mode_search_for_base() {
                T->reset();
                cT->reset();

                reset_touched();
                mode = IR_MODE_RECORD_TRACE;
                T->set_compare(false);
                T->set_record(true);
            }

            void mode_dfs() {
                mode = ir::IR_MODE_COMPARE_TRACE;
            }

            void mode_bfs() {

            }

            void mode_random_walk() {

            }

            void flip_trace() {
                trace* t_flip = T;
                T  = cT;
                cT = t_flip;
            }

            void compare_to_this() {
                cT->set_compare(true);
                cT->set_record(false);
                cT->set_compare_trace(T);
                cT->set_position(T->get_position());
                flip_trace();

                compare_base_color.clear();
                compare_base_color.resize(base_color.size());
                std::copy(base_color.begin(), base_color.end(), compare_base_color.begin());
                compare_base_cells.clear();
                compare_base_cells.resize(base_cells.size());
                std::copy(base_cells.begin(), base_cells.end(), compare_base_cells.begin());
                compare_singletons.clear();
                compare_singletons.resize(singletons.size());
                std::copy(singletons.begin(), singletons.end(), compare_singletons.begin());
                compare_base.clear();
                compare_base.resize(base_vertex.size());
                std::copy(base_vertex.begin(), base_vertex.end(), compare_base.begin());

                leaf_color.copy_force(c);

                mode = ir::IR_MODE_COMPARE_TRACE;
            }

            controller(coloring *c) {
                this->c = c;

                T  = &_T1;
                cT = &_T2;

                workspace.allocate(c->lab_sz);
                touched_color.initialize(c->lab_sz);
                touched_color_list.allocate(c->lab_sz);
                prev_color_list.allocate(c->lab_sz);

                touch_initial_colors();

                my_split_hook = self_split_hook();
                my_worklist_hook = self_worklist_hook();
                my_add_hook   = self_add_hook();
            }

            void use_limited_reversible_for_next() {
                reset_touched();
                mode = IR_MODE_COMPARE_TRACE;
            }

            void reset_trace_equal() {
                T->reset_trace_equal();
                h_deviation_inc_current = 0;
            }

            void use_trace_early_out(bool trace_early_out) {
                this->h_trace_early_out = trace_early_out;
                if(!trace_early_out) use_increase_deviation_hash(false);
            }

            void use_increase_deviation_hash(bool deviation_inc_active) {
                h_deviation_inc_active = deviation_inc_active;
            }

            /**
             * Save a partial state of this controller.
             *
             * @param state A reference to the reduced_save in which the state will be stored.
             */
            void save_reduced_state(reduced_save &state) {
                state.set_state(base_vertex, *c, T->get_hash(), T->get_position(), base_pos);
            }

            /**
             * Load a partial state into this controller.
             *
             * @param state A reference to the reduced_save from which the state will be loaded.
             */
            void __attribute__ ((noinline)) load_reduced_state(reduced_save &state) {
                bool partial_base = false;
                /*if(state.get_base().size() <= base_vertex.size()) {
                    int i;
                    for (i = 0; i < state.get_base().size(); ++i) if (state.get_base()[i] != base_vertex[i]) break;
                    partial_base = (i == state.get_base().size());
                }

                if(partial_base) {
                    c->copy(state.get_coloring());
                } else {*/
                    c->copy_force(state.get_coloring());
                //}

                T->set_hash(state.get_invariant_hash());
                T->set_position(state.get_trace_position());
                T->reset_trace_equal();
                T->set_compare(true);
                base_pos = state.get_base_position();
                base_vertex = state.get_base();

                // these become meaningless
                base_singleton_pt.clear();
                base_color.clear();
                base_color_sz.clear();
                base_touched_color_list_pt.clear();
                base_cells.clear();

                // deactivate reversability
                mode = IR_MODE_RECORD_HASH_IRREVERSIBLE;
            }

            void load_reduced_state_without_coloring(reduced_save &state) {
                T->set_hash(state.get_invariant_hash());
                T->set_position(state.get_trace_position());
                T->reset_trace_equal();
                T->set_compare(true);
                base_pos = state.get_base_position();
                base_vertex = state.get_base();

                // these become meaningless
                base_singleton_pt.clear();
                base_color.clear();
                base_color_sz.clear();
                base_touched_color_list_pt.clear();
                base_cells.clear();

                // deactivate reversability
                mode = IR_MODE_RECORD_HASH_IRREVERSIBLE;
            }

            coloring *get_coloring() {
                return c;
            }

            int get_base_pos() {
                return base_pos;
            }

            void set_base_pos(int pos) {
                base_pos = pos;
            }

            void add_hook(const long d) {
                if (T) T->op_additional_info(d);
            }

            bool split_hook(const int old_color, const int new_color, const int new_color_sz) {
                        if (mode != IR_MODE_RECORD_HASH_IRREVERSIBLE) {
                            // update some heuristic values
                            if (new_color_sz > 1) {
                                h_last_refinement_singleton_only = false;
                                h_hint_color = new_color;
                                h_hint_color_is_singleton_now = false;
                            }
                            if (new_color == h_hint_color && new_color_sz == 1) {
                                h_hint_color_is_singleton_now = true;
                            }
                            // write singletons to singleton list
                            if (new_color_sz == 1) {
                                singletons.push_back(c->lab[new_color]);
                            }

                            // record colors that were changed
                            if (!touched_color.get(new_color)) {
                                touched_color.set(new_color);
                                prev_color_list.push_back(old_color);
                                touched_color_list.push_back(new_color);
                            }
                        }

                        // record split into trace invariant, unless we are individualizing
                        if (!h_individualize && old_color != new_color) T->op_refine_cell_record(new_color, new_color_sz, 1);

                        const bool cont = !h_trace_early_out || T->trace_equal();
                        h_deviation_inc_current += (!cont);
                        const bool deviation_override = h_deviation_inc_active && (h_deviation_inc_current <= h_deviation_inc);

                        //const bool ndone = !((mode != IR_MODE_RECORD_TRACE) && T->trace_equal() && compare_base_cells[base_pos - 1] == c->cells);
                        const bool ndone = true;
                        return ndone && (cont || deviation_override);
            }

            void __attribute__ ((noinline)) write_strong_invariant(sgraph* g) {
                for(int l = 0; l < g->v_size/2; ++l) { // "half of them should be enough"...
                    const int v = c->lab[l];
                    unsigned int inv1 = 0;
                    const int start_pt = g->v[v];
                    const int end_pt   = start_pt + g->d[v];
                    for(int pt = start_pt; pt < end_pt; ++pt) {
                        const int other_v = g->e[pt];
                        inv1 += hash((unsigned int) c->vertex_to_col[other_v]);
                    }
                    T->op_additional_info((int) inv1);
                }
            }

            bool worklist_hook(const int color, const int color_sz) {
                if (h_cell_active) {
                    T->op_refine_cell_end();
                    h_cell_active = false;
                }

                // update some heuristic values
                // TODO: only activate blueprints on first few restarts!
                if (T->trace_equal() && !T->blueprint_is_next_cell_active()) {
                    T->blueprint_skip_to_next_cell();
                    return false;
                }

                T->op_refine_cell_start(color);
                if (!h_individualize) T->op_additional_info(color_sz);

                h_cell_active = true;
                return true;
            }

            std::function<type_split_color_hook> self_split_hook() {
                return std::bind(&controller::split_hook, this, std::placeholders::_1, std::placeholders::_2,
                                 std::placeholders::_3);
            }

            std::function<type_worklist_color_hook> self_worklist_hook() {
                return std::bind(&controller::worklist_hook, this, std::placeholders::_1, std::placeholders::_2);
            }

            std::function<type_additional_info_hook> self_add_hook() {
                return std::bind(&controller::add_hook, this, std::placeholders::_1);
            }

            /**
             * Move IR node kept in this controller to a child, specified by a vertex to be individualized.
             *
             * @param R a refinement workspace
             * @param g the graph
             * @param v the vertex to be individualized
             */
            void move_to_child(refinement *R, sgraph *g, int v) {
                ++base_pos;
                base_vertex.push_back(v); // always keep track of base (needed for BFS)

                if (mode != IR_MODE_RECORD_HASH_IRREVERSIBLE) {
                    base_singleton_pt.push_back(singletons.size());
                    base_color.push_back(c->vertex_to_col[v]);
                    base_color_sz.push_back(c->ptn[c->vertex_to_col[v]] + 1);
                    base_touched_color_list_pt.push_back(touched_color_list.cur_pos);
                }

                assert(!h_cell_active);

                h_individualize = true;
                const int prev_col = c->vertex_to_col[v];
                const int init_color_class = R->individualize_vertex(c, v, my_split_hook);
                if (T) T->op_individualize(prev_col, c->vertex_to_col[v]);
                h_individualize = false;
                h_last_refinement_singleton_only = true;

                if (T) T->op_refine_start();

                if (mode == IR_MODE_RECORD_TRACE) {
                    R->refine_coloring(g, c, init_color_class, -1, my_split_hook,
                                       my_worklist_hook);
                    if (T && h_cell_active) T->op_refine_cell_end();
                    if (T) T->op_refine_end();

                    //std::cout << T->get_position() <<  ", " << c->cells << std::endl;
                } else {
                    // TODO compare_base_cells not necessarily applicable
                    // compare_base_cells[base_pos - 1]
                    // T->trace_equal()?compare_base_cells[base_pos - 1]:-1

                    R->refine_coloring(g, c, init_color_class, T->trace_equal()?compare_base_cells[base_pos - 1]:-1, my_split_hook,
                                       my_worklist_hook);
                    assert(T->trace_equal()?c->cells==compare_base_cells[base_pos-1]:true);
                    if (T && T->trace_equal()) {T->skip_to_individualization();}
                }

                h_cell_active = false;

                if (mode != IR_MODE_RECORD_HASH_IRREVERSIBLE) base_cells.push_back(c->cells);
            }

            /**
             * Move IR node kept in this controller to a child, specified by a vertex to be individualized.
             *
             * @param R a refinement workspace
             * @param g the graph
             * @param v the vertex to be individualized
             */
            void move_to_child_no_trace(refinement *R, sgraph *g, int v) {
                const int init_color_class = R->individualize_vertex(c, v);
                R->refine_coloring_first(g, c, init_color_class);
            }

            void refine(refinement *R, sgraph *g) {
                R->refine_coloring_first(g, c);
            }

            /**
             * Move IR node kept in this controller back to its parent.
             */
            void move_to_parent() {
                assert(mode != IR_MODE_RECORD_HASH_IRREVERSIBLE);

                // unwind invariant
                if (T) T->rewind_to_individualization();

                --base_pos;
                while (prev_color_list.cur_pos > base_touched_color_list_pt[base_touched_color_list_pt.size()-1]) {
                    const int old_color = prev_color_list.pop_back();
                    const int new_color = touched_color_list.pop_back();

                    touched_color.unset(new_color);

                    const int new_color_sz = c->ptn[new_color] + 1;
                    c->ptn[old_color] += new_color_sz;
                    c->ptn[new_color] = 1;

                    if (c->ptn[old_color] > 0) {
                        h_hint_color = old_color;
                        h_hint_color_is_singleton_now = false;
                    }

                    for (int j = 0; j < new_color_sz; ++j) {
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
                base_color_sz.pop_back();
                base_touched_color_list_pt.pop_back();
                base_cells.pop_back();
                base_singleton_pt.pop_back();

                T->reset_trace_equal();
            }
        };

        /**
         * \brief Creates cell selectors.
         *
         * Heuristics which enable the creation of different cell selectors, as well as moving an \ref controller to a
         * leaf of the IR tree.
         */
        class selector_factory {
            // heuristics to attempt to find a good selector, as well as creating a split-map for dfs_ir
            // TODO: enable different selector strategies for restarts
            // TODO: make a grab-able hook for dynamic selector which can also deviate from base color (blueprint selector)
            int locked_lim = 512;

            std::vector<int> saved_color_base;
            std::function<type_selector_hook> dynamic_seletor;

        public:
            mark_set test_set;

            // TODO: save the base of state for blueprint selection
            void save_base(controller *state) {

            }

            // TODO configuration of dynamic selector when deviating from base

            /**
             * Dynamic cell selector, chooses first non-trivial color class, unless color class from stored base is
             * applicable.
             *
             * @param c Coloring from which a color class shall be selected.
             * @param base_pos Current position in base.
             */
            int dynamic_selector_first(const coloring *c, const int base_pos) {
                if (base_pos >= 0 && base_pos < saved_color_base.size() && c->ptn[saved_color_base[base_pos]] > 0) {
                    return saved_color_base[base_pos];
                }
                for (int i = 0; i < c->ptn_sz;) {
                    if (c->ptn[i] > 0) {
                        return i;
                    }
                    i += c->ptn[i] + 1;
                }
                return -1;
            }

            /**
             * @return A selector hook based on the saved base and configured dynamic selector.
             */
            std::function<type_selector_hook> *get_selector_hook() {
                dynamic_seletor = std::bind(&selector_factory::dynamic_selector_first, this, std::placeholders::_1,
                                            std::placeholders::_2);
                return &dynamic_seletor;
            }


            void find_base(const int m_choose, refinement *R, sgraph *g, controller *state) {
                switch(m_choose % 3) {
                    case 0:
                        find_sparse_optimized_base(R, g, state);
                        break;
                    case 1:
                        find_combinatorial_optimized_base(R, g, state);
                        break;
                    case 2:
                        find_small_optimized_base(R, g, state);
                        break;
                }
            }
            /**
             * Find a base/selector for a given graph \p g from state \p state. Attemtps to find a good base for simple,
             * sparse graphs.
             *
             * @param R A refinement workspace.
             * @param g The graph.
             * @param state The IR state from which a base is created.
             */
            void find_sparse_optimized_base(refinement *R, sgraph *g, controller *state) {
                state->mode_search_for_base();

                std::vector<int> candidates;
                test_set.initialize(g->v_size);
                candidates.reserve(locked_lim);
                int prev_color = -1;

                mark_set neighbour_color;
                neighbour_color.initialize(g->v_size);

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;

                    // pick previous color if possible
                    if (prev_color >= 0 && state->get_coloring()->ptn[prev_color] > 0) {
                        //std::cout << "previous color" << std::endl;
                        best_color = prev_color;
                    } else if (prev_color >= 0) { // pick neighbour of previous color if possible
                        const int test_vertex = state->get_coloring()->lab[prev_color];
                        for (int i = 0; i < g->d[test_vertex]; ++i) {
                            const int other_vertex = g->e[g->v[test_vertex] + i];
                            const int other_color = state->get_coloring()->vertex_to_col[other_vertex];
                            if (state->get_coloring()->ptn[other_color] > 0) {
                                //std::cout << "neighbour color" << std::endl;
                                neighbour_color.set(other_color);
                                best_color = other_color;
                                //break;
                            }
                        }
                    }

                    // heuristic, try to pick "good" color
                    if (best_color == -1) {
                        candidates.clear();
                        int best_score = -1;

                        for (int i = 0; i < state->get_coloring()->ptn_sz;) {
                            if (state->get_coloring()->ptn[i] > 0) {
                                candidates.push_back(i);
                            }

                            if (candidates.size() >= (size_t) locked_lim)
                                break;

                            i += state->get_coloring()->ptn[i] + 1;
                        }

                        while (!candidates.empty()) {
                            const int test_color = candidates.back();
                            candidates.pop_back();

                            int test_score = color_score(g, state, test_color);
                            if (neighbour_color.get(test_color)) {
                                test_score *= 10;
                            }
                            if (test_score > best_score) {
                                best_color = test_color;
                                best_score = test_score;
                            }
                        }
                    }

                    assert(best_color >= 0);
                    assert(best_color < g->v_size);
                    prev_color = best_color;
                    state->move_to_child(R, g, state->get_coloring()->lab[best_color]);
                }

                saved_color_base = state->base_color;
            }

            void find_test_base(refinement *R, sgraph *g, controller *state) {
                state->mode_search_for_base();

                std::vector<int> candidates;
                candidates.reserve(locked_lim);
                int prev_color = -1;

                mark_set neighbour_color;
                neighbour_color.initialize(g->v_size);

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;

                    // pick previous color if possible
                    if (prev_color >= 0 && state->get_coloring()->ptn[prev_color] > 0) {
                        //std::cout << "previous color" << std::endl;
                        best_color = prev_color;
                    }

                    // heuristic, try to pick "good" color
                    if (best_color == -1) {
                        candidates.clear();
                        int best_score = -1;

                        for (int i = 0; i < state->get_coloring()->ptn_sz;) {
                            if (state->get_coloring()->ptn[i] > 0) {
                                candidates.push_back(i);
                            }

                            i += state->get_coloring()->ptn[i] + 1;
                        }

                        while (!candidates.empty()) {
                            const int test_color = candidates.back();
                            candidates.pop_back();

                            int test_score = color_score_size(g, state, test_color);
                            if (neighbour_color.get(test_color)) {
                                test_score *= 10;
                            }
                            if (test_score > best_score) {
                                best_color = test_color;
                                best_score = test_score;
                            }
                        }
                    }

                    assert(best_color >= 0);
                    assert(best_color < g->v_size);
                    prev_color = best_color;
                    state->move_to_child(R, g, state->get_coloring()->lab[best_color]);
                }

                saved_color_base = state->base_color;
            }

            /**
             * Find a base/selector for a given graph \p g from state \p state. Attempts to find a with small color
             * classes.
             *
             * @param R A refinement workspace.
             * @param g The graph.
             * @param state The IR state from which a base is created.
             */
            void find_small_optimized_base(refinement *R, sgraph *g, controller *state) {
                state->mode_search_for_base();

                std::vector<int> candidates;
                test_set.initialize(g->v_size);
                candidates.reserve(locked_lim);
                int prev_color = -1;

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;

                    // heuristic, try to pick "good" color
                    if (best_color == -1) {
                        candidates.clear();
                        int best_score = INT32_MIN;

                        for (int i = 0; i < state->get_coloring()->ptn_sz;) {
                            if (state->get_coloring()->ptn[i] > 0) {
                                candidates.push_back(i);
                            }

                            if (candidates.size() >= (size_t) locked_lim)
                                break;

                            i += state->get_coloring()->ptn[i] + 1;
                        }

                        while (!candidates.empty()) {
                            const int test_color = candidates.back();
                            candidates.pop_back();

                            int test_score = color_score_anti_size(g, state, test_color);
                            if (test_score > best_score || best_color == -1) {
                                best_color = test_color;
                                best_score = test_score;
                            }
                        }
                    }

                    assert(best_color >= 0);
                    assert(best_color < g->v_size);
                    prev_color = best_color;
                    state->move_to_child(R, g, state->get_coloring()->lab[best_color]);
                }

                saved_color_base = state->base_color;
            }

            /**
             * Find a base/selector for a given graph \p g from state \p state. Attemtps to find a good base for
             * combinatorial graphs solved by bfs/random walks (i.e., by choosing large colors).
             *
             * @param R A refinement workspace.
             * @param g The graph.
             * @param state The IR state from which a base is created.
             */
            void find_combinatorial_optimized_base(refinement *R, sgraph *g, controller *state) {
                state->mode_search_for_base();

                std::vector<int> candidates;
                test_set.initialize(g->v_size);
                candidates.reserve(locked_lim);

                mark_set neighbour_color;
                neighbour_color.initialize(g->v_size);

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;

                    // heuristic, try to pick "good" color
                    if (best_color == -1) {
                        candidates.clear();
                        int best_score = -1;

                        for (int i = 0; i < state->get_coloring()->ptn_sz;) {
                            if (state->get_coloring()->ptn[i] > 0) {
                                candidates.push_back(i);
                            }

                            //if (candidates.size() >= (size_t) locked_lim)
                            //    break;

                            i += state->get_coloring()->ptn[i] + 1;
                        }

                        while (!candidates.empty()) {
                            const int test_color = candidates.back();
                            candidates.pop_back();

                            int test_score = color_score_size(g, state, test_color);
                            if (neighbour_color.get(test_color)) {
                                test_score *= 10;
                            }
                            if (test_score > best_score) {
                                best_color = test_color;
                                best_score = test_score;
                            }
                        }
                    }

                    assert(best_color >= 0);
                    assert(best_color < g->v_size);
                    state->move_to_child(R, g, state->get_coloring()->lab[best_color]);
                }

                saved_color_base = state->base_color;
            }

            int color_score(sgraph *g, controller *state, int color) {
                test_set.reset();
                const int v = state->get_coloring()->lab[color];
                const int d = g->d[v];
                const int ept = g->v[v];
                int non_triv_col_d = 1;
                for (int i = 0; i < d; ++i) {
                    const int test_col = state->get_coloring()->vertex_to_col[g->e[ept + i]];
                    if (!test_set.get(test_col)) {
                        non_triv_col_d += 1;
                        test_set.set(test_col);
                    }
                }
                return state->get_coloring()->ptn[color] * non_triv_col_d;
            }

            int color_score_size(sgraph *g, controller *state, int color) {
                return state->get_coloring()->ptn[color];
            }

            int color_score_anti_size(sgraph *g, controller *state, int color) {
                return INT32_MAX - state->get_coloring()->ptn[color];
            }


            int state_score(sgraph *g, controller *state) {
                return state->get_coloring()->cells;
            }
        };

        class tree_node {
            std::mutex    lock;
            reduced_save* data;
            tree_node*    next;
            bool          is_base = false;
            bool          is_pruned = false;
        public:
            tree_node(reduced_save* data, tree_node* next) {
                this->data = data;
                this->next = next;
                if(next == nullptr) {
                    next = this;
                }
            }
            tree_node* get_next() {
                return next;
            }
            void set_next(tree_node* next) {
                this->next = next;
            }
            reduced_save* get_save() {
                return data;
            }

            void prune() {
                is_pruned = true;
            }
            bool get_prune() {
                return is_pruned;
            }

            void base() {
                is_base = true;
            }
            bool get_base() {
                return is_base;
            }
        };

        typedef std::pair<ir::tree_node*, int> missing_node;

        // TODO shared_tree structure for BFS + random walk
        class shared_tree {
            shared_queue_t<missing_node>         missing_nodes;
            std::vector<tree_node*>              tree_data;
            std::vector<std::vector<tree_node*>> tree_data_jump_map;
            std::vector<int>        tree_level_size;
            int                     finished_up_to = 0;

            std::vector<long>       node_invariant;
        public:

            std::vector<long>* get_node_invariant() {
                return &node_invariant;
            }

            void initialize(int base_size, ir::reduced_save* root) {
                tree_data.resize(base_size + 1);
                tree_level_size.resize(base_size + 1);
                tree_data_jump_map.resize(base_size + 1);
                add_node(0, root, true);
                node_invariant.resize(root->get_coloring()->lab_sz);
            }

            void reset(int new_base_size, ir::reduced_save* root, int keep_level) {

            }


            void queue_reserve(const int n) {
                missing_nodes.reserve(n);
            }

            void queue_missing_node(missing_node node) {
                missing_nodes.add(node);
            }

            bool queue_missing_node_empty() {
                return missing_nodes.empty();
            }

            missing_node queue_missing_node_pop() {
                return missing_nodes.pop();
            }

            void mark_first_level(mark_set& marks) {
                if(tree_data[1] == nullptr) return;

                tree_node * first = tree_data[1];
                tree_node * next = first;
                do {
                    marks.set(next->get_save()->get_base()[0]);
                    next = next->get_next();
                } while (next != first);
            }

            void record_invariant(int v, long inv) {
                node_invariant[v] = inv;
            }

            void add_node(int level, reduced_save* data, bool is_base = false) {
                // TODO use locks
                if(tree_data[level] == nullptr) {
                    tree_level_size[level] = 0;
                    tree_data[level] = new tree_node(data, nullptr);
                    tree_data[level]->set_next(tree_data[level]);
                    if(is_base) tree_data[level]->base();
                } else {
                    tree_node* a_node    = tree_data[level];
                    tree_node* next_node = a_node->get_next();
                    auto       new_node  = new tree_node(data, next_node);
                    if(is_base) new_node->base();
                    a_node->set_next(new_node);
                    tree_data[level] = new_node;
                }
                ++tree_level_size[level];
            }

            void finish_level(int level) {
                if(tree_data_jump_map[level].size() == 0) {
                    tree_node * first = tree_data[level];
                    tree_data_jump_map[level].reserve(tree_level_size[level]);
                    tree_node * next = first;
                    do {
                        tree_data_jump_map[level].push_back(next);
                        next = next->get_next();
                    } while (next != first);
                    assert(tree_data_jump_map[level].size() == tree_level_size[level]);
                }
            }

            ir::tree_node* pick_node_from_level(const int level, int num) {
                finish_level(level);
                // TODO very sure this is too slow!
                num = num % tree_level_size[level];
                return tree_data_jump_map[level][num];
            }

            int get_finished_up_to() {
                return finished_up_to;
            }

            void set_finished_up_to(const int finished_up_to) {
                this->finished_up_to = finished_up_to;
            }

            tree_node* get_level(int level) {
                return tree_data[level];
            }

            int get_level_size(int level) {
                return tree_level_size[level];
            }

            ~shared_tree() {
                //TODO: write this
            };
        };
    }
}

#endif //DEJAVU_IR_H
