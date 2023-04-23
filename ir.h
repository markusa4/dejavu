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
            IR_MODE_RECORD_TRACE, IR_MODE_COMPARE_TRACE_REVERSIBLE, IR_MODE_COMPARE_TRACE_IRREVERSIBLE
        };

        /**
         * int type_selector_hook(coloring* c, const int s_base_pos);
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
        public:
            coloring      *c  = nullptr;
            trace         *T  = nullptr;

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

            work_list workspace;

            // statistics
            int   s_base_pos = 0;                          /**< how large the base of the current IR node is*/
            bool  s_last_refinement_singleton_only = true; /**< whether in the last refinement, only singleton cells
                                                            * appeared*/
            int   s_hint_color                     = -1;   /**< a color that was not singleton in last refinement*/
            bool  s_hint_color_is_singleton_now    = true; /**< whether the \a s_hint_color is singleton now */

            // settings for heuristics
            int   h_deviation_inc = 48; /**< how many additional splits to wait before terminating whenever
                                         * use_increase_deviation is used  */
        private:
            refinement* R;
            trace *cT = nullptr;
            trace _T1;
            trace _T2;

            mark_set  touched_color;
            work_list touched_color_list;
            work_list prev_color_list;

            std::function<type_split_color_hook>     my_split_hook;
            std::function<type_worklist_color_hook>  my_worklist_hook;
            std::function<type_additional_info_hook> my_add_hook;

            // settings
            ir_mode mode = IR_MODE_RECORD_TRACE;

            // internal flags for heuristics
            bool h_cell_active     = false;
            bool h_individualize   = false;
            bool h_trace_early_out = false;

            bool  h_deviation_inc_active  = false;
            int   h_deviation_inc_current = 0;

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

            void flip_trace() {
                trace* t_flip = T;
                T  = cT;
                cT = t_flip;
            }

            bool split_hook(const int old_color, const int new_color, const int new_color_sz) {
                if (mode != IR_MODE_COMPARE_TRACE_IRREVERSIBLE) {
                    // update some heuristic values
                    if (new_color_sz > 1) {
                        s_last_refinement_singleton_only = false;
                        s_hint_color = new_color;
                        s_hint_color_is_singleton_now = false;
                    } else if (new_color == s_hint_color && new_color_sz == 1) {
                        s_hint_color_is_singleton_now = true;
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

                // TODO: should this be in split hook, or in color refinement?
                const bool ndone = !((mode != IR_MODE_RECORD_TRACE) && T->trace_equal() && compare_base_cells[s_base_pos - 1] == c->cells);
                //const bool ndone = true;
                return ndone && (cont || deviation_override);
            }

            std::function<type_split_color_hook> self_split_hook() {
                return std::bind(&controller::split_hook, this, std::placeholders::_1, std::placeholders::_2,
                                 std::placeholders::_3);
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

            std::function<type_worklist_color_hook> self_worklist_hook() {
                return std::bind(&controller::worklist_hook, this, std::placeholders::_1, std::placeholders::_2);
            }

            std::function<type_additional_info_hook> self_add_hook() {
                return std::bind(&controller::write_to_trace, this, std::placeholders::_1);
            }

        public:
            /**
             * Initialize this controller using a refinement and a graph coloring for the initial state.
             *
             * @param R The refinement workspace to use.
             * @param c The initial coloring.
             */
            controller(refinement* R, coloring *c) {
                this->c = c;
                this->R = R;

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

            /**
             * Sets internal trace into recording mode. We write a trace which we might want to compare to later.
             *
             * Always reversible (see also \a use_reversible).
             */
            void mode_write_base() {
                T->reset();
                cT->reset();

                reset_touched();
                mode = IR_MODE_RECORD_TRACE;
                T->set_compare(false);
                T->set_record(true);
            }

            /**
             * We compare to a pre-existing trace, recorded earlier using mode_write_base. Must use \a compare_to_this or
             * provide a comparison trace in another manner before being able to use this mode in the intended manner.
             */
            void mode_compare_base() {
                mode = ir::IR_MODE_COMPARE_TRACE_REVERSIBLE;
            }

            /**
             * Compare all following computations to this IR node. The \a mode must be set to `IR_MODE_RECORD_TRACE`
             * (using \a mode_write_base) to call this function. Changes the \a mode to `IR_MODE_COMPARE_TRACE_REVERSIBLE`
             * (i.e., such as calling \a mode_compare_base).
             */
            void compare_to_this() {
                assert(mode == ir::IR_MODE_RECORD_TRACE);

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

                mode = ir::IR_MODE_COMPARE_TRACE_REVERSIBLE;
            }

            /**
             * Enables or disables whether following \a move_to_child calls will be reversible using \a move_to_parent,
             * or not. Using non-reversible calls to \a move_to_child is faster.
             *
             * When recording a trace using \a mode_write_base, computations are always reversible and calling this
             * function will have no effect.
             * @param reversible
             */
            void use_reversible(const bool reversible) {
                if(mode == IR_MODE_RECORD_TRACE) return;
                if(reversible) {
                    reset_touched();
                    mode = IR_MODE_COMPARE_TRACE_REVERSIBLE;
                } else {
                    mode = IR_MODE_COMPARE_TRACE_IRREVERSIBLE;
                }
            }

            /**
             * Whether to terminate color refinement whenever a deviation to its comparison trace is found.
             *
             * @param trace_early_out Flag that determines whether the early out is used.
             */
            void use_trace_early_out(bool trace_early_out) {
                this->h_trace_early_out = trace_early_out;
                if(!trace_early_out) use_increase_deviation(false);
            }

            /**
             * Whether to record additional deviation information once a deviation from the comparison trace is found.
             * Only applicable when \a use_trace_early_out is set. Essentially delays the termination of color refinement
             * to record more information into the trace invariant.
             *
             * @param deviation_inc_active Whether increased deviation is recorded or not.
             */
            void use_increase_deviation(bool deviation_inc_active) {
                h_deviation_inc_active = deviation_inc_active;
            }

            /**
             * Resets whether the trace is deemed equal to its comparison trace.
             */
            void reset_trace_equal() {
                T->reset_trace_equal();
                h_deviation_inc_current = 0;
            }

            /**
             * Save a partial state of this controller.
             *
             * @param state A reference to the reduced_save in which the state will be stored.
             */
            void save_reduced_state(reduced_save &state) {
                state.set_state(base_vertex, *c, T->get_hash(), T->get_position(), s_base_pos);
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
                s_base_pos = state.get_base_position();
                base_vertex = state.get_base();

                // these become meaningless
                base_singleton_pt.clear();
                base_color.clear();
                base_color_sz.clear();
                base_touched_color_list_pt.clear();
                base_cells.clear();

                // if reversible, need to reset touched colors
                if(mode == IR_MODE_COMPARE_TRACE_REVERSIBLE) reset_touched();
            }

            // TODO: hopefully can be deprecated
            void load_reduced_state_without_coloring(reduced_save &state) {
                T->set_hash(state.get_invariant_hash());
                T->set_position(state.get_trace_position());
                T->reset_trace_equal();
                T->set_compare(true);
                s_base_pos = state.get_base_position();
                base_vertex = state.get_base();

                // these become meaningless
                base_singleton_pt.clear();
                base_color.clear();
                base_color_sz.clear();
                base_touched_color_list_pt.clear();
                base_cells.clear();

                // deactivate reversability
                mode = IR_MODE_COMPARE_TRACE_IRREVERSIBLE;
            }

            coloring *get_coloring() {
                return c;
            }

            int get_base_pos() {
                return s_base_pos;
            }

            /**
             * Write additional information into the internal trace. Care must be taken that the written data is
             * isomorphism-invariant.
             *
             * @param d Data to be written to the trace.
             */
            void write_to_trace(const int d) {
                T->op_additional_info(d);
            }

            /**
             * Writes a stronger invariant using the internal coloring and graph to the trace.
             *
             * @param g The graph.
             */
            void write_strong_invariant(const sgraph* g) const {
                for(int l = 0; l < g->v_size; ++l) { // "half of them should be enough"...
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

            /**
             * Move IR node kept in this controller to a child, specified by a vertex to be individualized.
             *
             * @param R a refinement workspace
             * @param g the graph
             * @param v the vertex to be individualized
             */
            void move_to_child(sgraph *g, int v) {
                ++s_base_pos;
                base_vertex.push_back(v); // always keep track of base (needed for BFS)

                if (mode != IR_MODE_COMPARE_TRACE_IRREVERSIBLE) {
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
                s_last_refinement_singleton_only = true;

                if (T) T->op_refine_start();

                if (mode == IR_MODE_RECORD_TRACE) {
                    R->refine_coloring(g, c, init_color_class, -1, my_split_hook,
                                       my_worklist_hook);
                    if (T && h_cell_active) T->op_refine_cell_end();
                    if (T) T->op_refine_end();
                } else {
                    // TODO compare_base_cells currently implemented in split hook, could reconsider this...
                    // T->trace_equal()?compare_base_cells[s_base_pos - 1]:-1
                    R->refine_coloring(g, c, init_color_class, -1, my_split_hook,
                                       my_worklist_hook);
                    assert(T->trace_equal() ?c->cells==compare_base_cells[s_base_pos - 1] : true);
                    if (T && T->trace_equal()) {T->skip_to_individualization();}
                }

                h_cell_active = false;

                if (mode != IR_MODE_COMPARE_TRACE_IRREVERSIBLE) base_cells.push_back(c->cells);
            }

            /**
             * Move IR node kept in this controller to a child, specified by a vertex to be individualized.
             *
             * @param R a refinement workspace
             * @param g the graph
             * @param v the vertex to be individualized
             */
            void move_to_child_no_trace(sgraph *g, int v) {
                const int init_color_class = R->individualize_vertex(c, v);
                R->refine_coloring_first(g, c, init_color_class);
            }

            /**
             * Perform color refinement on the internal coloring based on the given graph.
             *
             * @param g The graph.
             */
            void refine(sgraph *g) {
                R->refine_coloring_first(g, c);
            }

            /**
             * Move IR node kept in this controller back to its parent.
             */
            void move_to_parent() {
                assert(mode != IR_MODE_COMPARE_TRACE_IRREVERSIBLE);
                assert(base_touched_color_list_pt.size() > 0);

                // unwind invariant
                if (T) T->rewind_to_individualization();

                --s_base_pos;
                while (prev_color_list.cur_pos > base_touched_color_list_pt[base_touched_color_list_pt.size()-1]) {
                    const int old_color = prev_color_list.pop_back();
                    const int new_color = touched_color_list.pop_back();

                    touched_color.unset(new_color);

                    const int new_color_sz = c->ptn[new_color] + 1;
                    c->ptn[old_color] += new_color_sz;
                    c->ptn[new_color] = 1;

                    if (c->ptn[old_color] > 0) {
                        s_hint_color = old_color;
                        s_hint_color_is_singleton_now = false;
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

            void walk(sgraph *g, ir::reduced_save &start_from, std::vector<int>& vertices) {
                load_reduced_state(start_from);

                while (s_base_pos < (int) vertices.size()) {
                    int v = vertices[s_base_pos];
                    move_to_child(g, v);
                }
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
                state->mode_write_base();

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
                    state->move_to_child(g, state->get_coloring()->lab[best_color]);
                }

                saved_color_base = state->base_color;
            }

            void find_test_base(refinement *R, sgraph *g, controller *state) {
                state->mode_write_base();

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
                    state->move_to_child(g, state->get_coloring()->lab[best_color]);
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
                state->mode_write_base();

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
                    state->move_to_child(g, state->get_coloring()->lab[best_color]);
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
                state->mode_write_base();

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
                    state->move_to_child(g, state->get_coloring()->lab[best_color]);
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

        /**
         * \brief An IR leaf
         *
         * A stored leaf of an IR tree. The leaf can be stored in a dense manner (coloring of the leaf), or a sparse
         * manner (base of the walk that leads to this leaf).
         */
        class stored_leaf {
        public:
            enum stored_leaf_type { STORE_LAB, ///< stores coloring of the leaf
                STORE_BASE ///< stores only base of the leaf
            };

            stored_leaf(int* arr, int arr_sz, stored_leaf_type store_type) : store_type(store_type) {
                lab_or_base.allocate(arr_sz);
                memcpy(lab_or_base.get_array(), arr, arr_sz * sizeof(int));
                lab_or_base.set_size(arr_sz);
            }

            stored_leaf(std::vector<int>& arr, stored_leaf_type store_type) : store_type(store_type) {
                lab_or_base.allocate(arr.size());
                std::copy(arr.begin(), arr.end(), lab_or_base.get_array());
                lab_or_base.set_size(arr.size());
            }

            const int* get_lab_or_base() {
                return lab_or_base.get_array();
            }

            int get_lab_or_base_size() {
                return lab_or_base.size();
            }

            [[nodiscard]] stored_leaf_type get_store_type() const {
                return store_type;
            }

        private:
            work_list lab_or_base;
            stored_leaf_type store_type;
        };

        /**
         * \brief Collection of leaves
         *
         * Can be used across multiple threads.
         *
         */
        class shared_leaves {
            std::mutex lock;
            std::unordered_multimap<long, stored_leaf*> leaf_store;
            std::vector<stored_leaf*> garbage_collector;

        public:
            int s_leaves = 0; /**< number of leaves stored */

            /**
             * Free up all memory.
             */
            ~shared_leaves() {
                for(int i = 0; i < garbage_collector.size(); ++i) {
                    delete garbage_collector[i];
                }
            }

            /**
             * Lookup whether a leaf with the given hash already exists.
             *
             * @param hash
             * @return
             */
            stored_leaf* lookup_leaf(long hash) {
                lock.lock();
                auto find = leaf_store.find(hash);
                if(find != leaf_store.end()) {
                    auto result = find->second;
                    lock.unlock();
                    return result;
                } else {
                    lock.unlock();
                    return nullptr;
                }
            }

            /**
             * Add leaf with the given hash. Does not add the leaf, if a leaf with the given hash already exists.
             *
             * @param hash
             * @param ptr
             */
            void add_leaf(long hash, coloring& c, std::vector<int>& base) {
                lock.lock();

                // check whether hash already exists
                if(leaf_store.contains(hash)) {
                    lock.unlock();
                    return;
                }

                // if not, add the leaf
                if(s_leaves < 100) {

                } else {

                }
                const bool full_save = s_leaves < 100;
                auto type = full_save?stored_leaf::stored_leaf_type::STORE_LAB:stored_leaf::stored_leaf_type::STORE_BASE;
                auto new_leaf = full_save?new stored_leaf(c.lab,c.lab_sz, type):new stored_leaf(base, type);
                leaf_store.insert(std::pair<long, stored_leaf*>(hash, new_leaf));
                garbage_collector.push_back(new_leaf);
                ++s_leaves;
                lock.unlock();
            }

            /**
             * Empty this leaf container.
             */
            void clear() {
                s_leaves = 0;
                leaf_store.clear();
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

        /**
         * \brief IR tree structure
         *
         * Datastructure to explicitly store parts of an IR tree, such as a level-wise store, leaf store, as well as
         * further information used for pruning in BFS.
         *
         * Can be used across multiple threads.
         */
        class shared_tree {
            shared_queue_t<missing_node>         missing_nodes;
            std::vector<tree_node*>              tree_data;
            std::vector<std::vector<tree_node*>> tree_data_jump_map;
            std::vector<int>        tree_level_size;
            std::vector<tree_node*> garbage_collector;
            int                     finished_up_to = 0;

            std::vector<int> current_base;

            std::vector<long>       node_invariant;
        public:
            shared_leaves stored_leaves;

            std::vector<long>* get_node_invariant() {
                return &node_invariant;
            }

            void initialize(std::vector<int> &base, ir::reduced_save* root) {
                tree_data.resize(base.size() + 1);
                tree_level_size.resize(base.size() + 1);
                tree_data_jump_map.resize(base.size() + 1);
                add_node(0, root, true);
                node_invariant.resize(root->get_coloring()->lab_sz);

                current_base = base;
            }

            void clear_leaves() {
                stored_leaves.clear();
            }

            /**
             * @return How many leaves were stored during random search.
             */
            int stat_leaves() {
                return stored_leaves.s_leaves;
            }

            bool reset(std::vector<int> &new_base, ir::reduced_save* root, bool keep_old) {
                const int old_size = current_base.size();
                const int new_size = new_base.size();

                // compare with stored base, keep whatever is possible
                int keep_until = 0;
                if(keep_old) {
                    for (; keep_until < old_size && keep_until < new_size; ++keep_until) {
                        if (current_base[keep_until] != new_base[keep_until]) break;
                    }
                } else {
                }

                if(keep_until == new_size && new_size == old_size) return false;
                finished_up_to = std::min(keep_until, finished_up_to);

                tree_data.resize(new_size + 1);
                tree_level_size.resize(new_size + 1);
                tree_data_jump_map.resize(new_size + 1);

                for (int i = keep_until; i < old_size; ++i) {
                    tree_level_size[i] = 0;
                    tree_data_jump_map[i].clear();
                    tree_data[i] = nullptr;
                }

                if(keep_until == 0) add_node(0, root, true);
                assert(missing_nodes.empty());

                return true;
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

                    garbage_collector.push_back( tree_data[level]);
                    if(is_base) tree_data[level]->base();
                } else {
                    tree_node* a_node    = tree_data[level];
                    tree_node* next_node = a_node->get_next();
                    auto       new_node  = new tree_node(data, next_node);
                    garbage_collector.push_back( new_node);
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
                for(int i = 0; i < garbage_collector.size(); ++i) delete garbage_collector[i];
            };
        };
    }
}

#endif //DEJAVU_IR_H
