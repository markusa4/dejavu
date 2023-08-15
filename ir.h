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
        class limited_save {
            // TODO enable a "compressed state" only consisting of base (but code outside should be oblivious to it)
            // TODO this is only supposed to be an "incomplete" state -- should there be complete states?

            std::vector<int> base_vertex; /**< base of vertices of this IR node (optional) */
            coloring c;                   /**< vertex coloring of this IR node   */
            long invariant     = 0;       /**< hash of invariant of this IR node */
            int trace_position = 0;       /**< position of trace of this IR node */
            int base_position  = 0;       /**< length of base of this IR node    */
        public:
            void save(std::vector<int> &s_base_vertex, coloring &s_c, long s_invariant, int s_trace_position,
                      int s_base_position) {
                this->base_vertex = s_base_vertex;
                this->c.copy_any(&s_c);
                this->invariant = s_invariant;
                this->trace_position = s_trace_position;
                this->base_position = s_base_position;
            }

            /**
             * @return coloring of this IR node
             */
            coloring *get_coloring() {
                return &c;
            }

            /**
             * @return hash of invariant of this IR node
             */
            [[nodiscard]] long get_invariant_hash() const {
                return invariant;
            }

            /**
             * @return position of trace of this IR node
             */
            [[nodiscard]] int get_trace_position() const {
                return trace_position;
            }

            /**
             * @return length of base of this IR node
             */
            [[nodiscard]] int get_base_position() const {
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

            // TODO: should these be public? or rather accesible through const references or something?

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

            // statistics
            int   s_base_pos = 0;                          /**< how large the base of the current IR node is*/
            bool  s_last_refinement_singleton_only = true; /**< whether in the last refinement, only singleton cells
                                                            * appeared*/
            int   s_hint_color                     = -1;   /**< a color that was not singleton in last refinement*/
            bool  s_hint_color_is_singleton_now    = true; /**< whether the \a s_hint_color is singleton now */

            // settings
            int   h_deviation_inc = 48; /**< how many additional splits to wait before terminating whenever
                                         * use_increase_deviation is used  */
        private:
            refinement* R;
            trace *cT = nullptr;
            trace internal_T1;
            trace internal_T2;

            mark_set  touched_color;
            work_list touched_color_list;
            work_list prev_color_list;

            mark_set  diff_tester;

            mark_set  diff_vertices;
            mark_set  diff_is_singleton;
            work_list diff_vertices_list;
            work_list diff_vertices_list_pt;

            bool diff_diverge = false;
            int singleton_pt_start = 0;

            std::function<type_split_color_hook>     my_split_hook;
            std::function<type_worklist_color_hook>  my_worklist_hook;

            // settings
            ir_mode mode = IR_MODE_RECORD_TRACE;

            // internal flags for heuristics
            bool s_cell_active     = false;
            bool s_individualize   = false;

            bool h_trace_early_out = false;       /*< use trace early out                */
            bool h_deviation_inc_active  = false; /*< use increased trace deviation      */
            int  s_deviation_inc_current = 0;     /*< number of current trace deviations */



            void touch_initial_colors() {
                int i = 0;
                while (i < c->domain_size) {
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
                    if (new_color_sz == 1) singletons.push_back(c->lab[new_color]);

                    // record colors that were changed
                    if (!touched_color.get(new_color)) {
                        touched_color.set(new_color);
                        prev_color_list.push_back(old_color);
                        touched_color_list.push_back(new_color);
                    }
                }

                // record split into trace invariant, unless we are individualizing
                if (!s_individualize && old_color != new_color) T->op_refine_cell_record(new_color);

                const bool cont = !h_trace_early_out || T->trace_equal();
                s_deviation_inc_current += (!cont);
                const bool deviation_override = h_deviation_inc_active && (s_deviation_inc_current <= h_deviation_inc);

                const bool ndone = !((mode != IR_MODE_RECORD_TRACE) && T->trace_equal()
                                    && compare_base_cells[s_base_pos - 1] == c->cells); // compare_base_cells[s_base_pos - 1] == c->cells
                return ndone && (cont || deviation_override);
            }

            std::function<type_split_color_hook> self_split_hook() {
                return [this](auto && PH1, auto && PH2, auto && PH3) { return
                        split_hook(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2),
                                std::forward<decltype(PH3)>(PH3)); };
            }

            bool worklist_hook(const int color, const int color_sz) {
                if (s_cell_active) {
                    T->op_refine_cell_end();
                    s_cell_active = false;
                }

                // blueprints
                if (T->trace_equal() && !T->blueprint_is_next_cell_active()) {
                    T->blueprint_skip_to_next_cell();
                    return false;
                }

                T->op_refine_cell_start(color);
                if (!s_individualize) T->op_additional_info(color_sz);

                s_cell_active = true;
                return true;
            }

            std::function<type_worklist_color_hook> self_worklist_hook() {
                return [this](auto && PH1, auto && PH2) {
                    return worklist_hook(std::forward<decltype(PH1)>(PH1),
                                       std::forward<decltype(PH2)>(PH2));
                };
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

                T  = &internal_T1;
                cT = &internal_T2;

                touched_color.initialize(c->domain_size);
                touched_color_list.allocate(c->domain_size);
                prev_color_list.allocate(c->domain_size);

                touch_initial_colors();

                my_split_hook = self_split_hook();
                my_worklist_hook = self_worklist_hook();

                diff_tester.initialize(c->domain_size);

                diff_vertices.initialize(c->domain_size);
                diff_vertices_list.resize(c->domain_size);
                diff_vertices_list_pt.resize(c->domain_size);
                diff_is_singleton.initialize(c->domain_size);

                singleton_pt_start = 0;
            }

            work_list* get_touched_colors() {
                return &touched_color_list;
            }

            mark_set* get_touched_colors_test() {
                return &touched_color;
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

            void add_diff_vertex(int v) {
                if(!diff_vertices.get(v) && !diff_is_singleton.get(v)) {
                    diff_vertices.set(v);
                    diff_vertices_list.push_back(v);
                    diff_vertices_list_pt[v] = diff_vertices_list.cur_pos - 1;
                }
            }

            void remove_diff_vertex(int v) {
                if(diff_vertices.get(v)) {
                    const int pt      = diff_vertices_list_pt[v];
                    assert(diff_vertices_list[pt] == v);

                    // which element is in the back?
                    const int back_pt  = diff_vertices_list.cur_pos-1;
                    const int back_v = diff_vertices_list[back_pt];

                    // now swap back_col to col
                    diff_vertices_list[pt] = back_v;
                    assert(diff_vertices.get(back_v));
                    assert(diff_vertices_list_pt[back_v] == back_pt);
                    diff_vertices_list_pt[back_v] = pt;
                    diff_vertices_list_pt[v] = -1; // I don't trust myself

                    diff_vertices.unset(v);
                    --diff_vertices_list.cur_pos;
                }
            }

            void __attribute__((noinline)) color_fix_difference(controller* state, int col) {
                const int col_sz = c->ptn[col] + 1;
                if(col_sz != state->c->ptn[col] + 1) {
                    assert(false);
                    return;
                }

                if(col_sz == 1) {
                    add_diff_vertex(state->c->lab[col]);
                    remove_diff_vertex(c->lab[col]);
                    if(!diff_is_singleton.get(c->lab[col])) {
                        diff_is_singleton.set(c->lab[col]);
                    }
                    return;
                }

                diff_tester.reset();
                for(int i = 0; i < col_sz; ++i) diff_tester.set(state->c->lab[col + i]);
                for(int i = 0; i < col_sz; ++i) if(!diff_tester.get(c->lab[col + i]))        add_diff_vertex(c->lab[col + i]);

                diff_tester.reset();
                for(int i = 0; i < col_sz; ++i) diff_tester.set(c->lab[col + i]);
                for(int i = 0; i < col_sz; ++i) if(!diff_tester.get(state->c->lab[col + i])) add_diff_vertex(state->c->lab[col + i]);
            }

            std::pair<int, int> diff_pair(controller* state) {
                diff_tester.reset();
                assert(diff_vertices_list.cur_pos > 0);

                const int right = diff_vertices_list[0];
                const int right_col = c->vertex_to_col[right];
                const int left = state->c->lab[right_col];

                return {left, right};
            }

            void reset_diff() {
                diff_vertices.reset();
                diff_vertices_list.reset();
                diff_is_singleton.reset();

                diff_tester.reset();
                singleton_pt_start = singletons.size();

                diff_diverge = false;
            }

            void __attribute__((noinline))  singleton_automorphism(controller* state, groups::automorphism_workspace*  automorphism) {
                automorphism->reset();
                for(int i = singleton_pt_start; i < singletons.size(); ++i) {
                    automorphism->write_single_map(singletons[i], state->singletons[i]);
                }
            }


            int get_diff_diverge() {
                return diff_diverge;
            }


            // returns whether there is diff cells
            bool __attribute__((noinline))  update_diff_vertices_last_individualization(controller* other_state) {
                int i = base_touched_color_list_pt[base_touched_color_list_pt.size()-1];
                assert(touched_color_list.cur_pos == other_state->touched_color_list.cur_pos);
                for(;i < touched_color_list.cur_pos; ++i) {
                    const int old_color = prev_color_list[i];
                    const int new_color = touched_color_list[i];

                    assert(old_color != new_color);

                    const int old_color_sz = c->ptn[old_color] + 1;
                    const int new_color_sz = c->ptn[new_color] + 1;

                    const int old_color_sz_other = other_state->c->ptn[old_color] + 1;
                    const int new_color_sz_other = other_state->c->ptn[new_color] + 1;

                    diff_diverge = diff_diverge || (old_color_sz != old_color_sz_other)
                                   || (new_color_sz != new_color_sz_other);

                    color_fix_difference(other_state, new_color);
                    //color_fix_difference(other_state, old_color);
                    if(old_color_sz == 1) color_fix_difference(other_state, old_color);
                }

                return (diff_vertices_list.cur_pos != 0) || diff_diverge; // there is a diff!
            }

            void copy(controller* state) {
                this->c->copy_any(state->c);

                T->set_compare(true);
                T->set_record(false);
                T->set_compare_trace(state->T->get_compare_trace());
                T->set_position(state->T->get_position());
                T->set_hash(state->T->get_hash());

                compare_base_color = state->compare_base_color;
                compare_base_cells = state->compare_base_cells;
                compare_singletons = state->compare_singletons;
                compare_base       = state->compare_base;

                leaf_color.copy_any(&state->leaf_color);
                mode = state->mode;

                base_color  = state->base_color;
                base_vertex = state->base_vertex;
                base_cells  = state->base_cells;
                base_color_sz = state->base_color_sz;
                base_touched_color_list_pt = state->base_touched_color_list_pt;
                base_singleton_pt = state->base_singleton_pt;
                touched_color.copy(&state->touched_color);
                prev_color_list.copy(&state->prev_color_list);
                touched_color_list.copy(&state->touched_color_list);
                base_touched_color_list_pt = state->base_touched_color_list_pt;
                singletons = state->singletons;

                diff_tester.copy(&state->diff_tester);
                s_base_pos = state->s_base_pos;

                this->R = state->R;
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

                leaf_color.copy_any(c);
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
                s_deviation_inc_current = 0;
            }

            /**
             * Save a partial state of this controller.
             *
             * @param state A reference to the limited_save in which the state will be stored.
             */
            void save_reduced_state(limited_save &state) {
                state.save(base_vertex, *c, T->get_hash(), T->get_position(),
                           s_base_pos);
            }

            /**
             * Load a partial state into this controller.
             *
             * @param state A reference to the limited_save from which the state will be loaded.
             */
            void load_reduced_state(limited_save &state) {
                c->copy_any(state.get_coloring());

                T->set_hash(state.get_invariant_hash());
                T->set_position(state.get_trace_position());
                T->reset_trace_equal();
                T->set_compare(true);
                s_base_pos = state.get_base_position();
                base_vertex = state.get_base();

                // these become meaningless, so clear them out
                base_singleton_pt.clear();
                base_color.clear();
                base_color_sz.clear();
                base_touched_color_list_pt.clear();
                base_cells.clear();

                // if reversible, need to reset touched colors
                if(mode == IR_MODE_COMPARE_TRACE_REVERSIBLE) reset_touched();
            }

            // TODO: hopefully can be deprecated
            void load_reduced_state_without_coloring(limited_save &state) {
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

            [[nodiscard]] coloring *get_coloring() const {
                return c;
            }

            [[nodiscard]] int get_base_pos() const {
                return s_base_pos;
            }

            /**
             * Write additional information into the internal trace. Care must be taken that the written data is
             * isomorphism-invariant.
             *
             * @param d Data to be written to the trace.
             */
            void write_to_trace(const int d) const {
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

            void write_strong_invariant_singleton(const sgraph* g) const {
                for(int l = 0; l < g->v_size; ++l) { // "half of them should be enough"...
                    if(c->ptn[l] != 0) continue;
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
            void __attribute__((noinline))  move_to_child(sgraph *g, int v) {
                ++s_base_pos;
                base_vertex.push_back(v); // always keep track of base (needed for BFS)

                if (mode != IR_MODE_COMPARE_TRACE_IRREVERSIBLE) {
                    base_singleton_pt.push_back((int) singletons.size());
                    base_color.push_back(c->vertex_to_col[v]);
                    base_color_sz.push_back(c->ptn[c->vertex_to_col[v]] + 1);
                    base_touched_color_list_pt.push_back(touched_color_list.cur_pos);
                }

                assert(!s_cell_active);

                s_individualize = true;
                const int prev_col = c->vertex_to_col[v];
                const int init_color_class = R->individualize_vertex(c, v, my_split_hook);
                T->op_individualize(prev_col, c->vertex_to_col[v]);
                s_individualize = false;
                s_last_refinement_singleton_only = true;
                T->op_refine_start();

                if (mode == IR_MODE_RECORD_TRACE) {
                    R->refine_coloring(g, c, init_color_class, -1, &my_split_hook,
                                       my_worklist_hook);
                    if (s_cell_active) T->op_refine_cell_end();
                    T->op_refine_end();
                } else {
                    R->refine_coloring(g, c, init_color_class, -1, &my_split_hook,
                                       my_worklist_hook);
                    if (T->trace_equal() && c->cells==compare_base_cells[s_base_pos - 1]) {
                        T->skip_to_individualization();
                    }
                }

                s_cell_active = false;

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
            void __attribute__((noinline)) move_to_parent() {
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

            /**
             * Perform the given walk from the given IR node.
             *
             * @param g The graph.
             * @param start_from A limited_save describing the IR node from which the walk should be performed.
             * @param vertices A vector of vertices that describes the walk, i.e., these vertices will be individualized
             *                 in the given order.
             */
            void walk(sgraph *g, ir::limited_save &start_from, std::vector<int>& vertices) {
                load_reduced_state(start_from);

                while (s_base_pos < (int) vertices.size()) {
                    int v = vertices[s_base_pos];
                    move_to_child(g, v);
                }
            }

            /**
             * Certifies whether the provided \p automorphism is indeed an automorphism of the graph \p g. Runs in
             * time roughly in the support of \p automorphism.
             *
             * Uses internal workspace to perform this operation efficiently.
             *
             * @param g The graph.
             * @param automorphism The automorphism.
             * @return Whether \p automorphism is an automorphism of \p g.
             */
            bool __attribute__((noinline))  certify(sgraph* g, groups::automorphism_workspace& automorphism) {
                if(automorphism.nsupport() > g->v_size/4) {
                    return R->certify_automorphism(g, automorphism.perm());
                } else {
                    return R->certify_automorphism_sparse(g, automorphism.perm(), automorphism.nsupport(),
                                                          automorphism.support());
                }
            }

            std::tuple<bool, int, int> certify_sparse_report_fail_resume(const sgraph *g, const int *colmap,
                                              groups::automorphism_workspace& automorphism, int pos_start) {
                return R->certify_automorphism_sparse_report_fail_resume(g, colmap,
                                                                         automorphism.perm(),
                                                                         automorphism.nsupport(),
                                                                         automorphism.support(),
                                                                         pos_start);
            }

            bool check_single_failure(const sgraph *g, const int *colmap, const groups::automorphism_workspace& automorphism,
                                      int failure) {
                return R->check_single_failure(g, colmap, automorphism.perm(), failure);
            }
        };

        /**
         * \brief Creates cell selectors.
         *
         * Heuristics which enable the creation of different cell selectors, as well as moving an \ref controller to a
         * leaf of the IR tree.
         */
        class base_selector {
            const int locked_lim = 512;

            std::vector<int> saved_color_base;
            std::vector<int> saved_color_sizes;
            std::function<type_selector_hook> dynamic_seletor;
            mark_set test_set;
            std::vector<int> candidates;

            big_number ir_tree_size_estimate;;

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

            static int color_score_size(controller *state, int color) {
                return state->get_coloring()->ptn[color];
            }

            static int color_score_anti_size(controller *state, int color) {
                return INT32_MAX - state->get_coloring()->ptn[color];
            }

        public:

            /**
             * Dynamic cell selector, chooses first non-trivial color class, unless color class from stored base is
             * applicable.
             *
             * @param c Coloring from which a color class shall be selected.
             * @param base_pos Current position in base.
             */
            int dynamic_selector_first(const coloring *c, const int base_pos) {
                if (base_pos >= 0 && base_pos < (int) saved_color_base.size() && c->ptn[saved_color_base[base_pos]] > 0 &&
                    c->vertex_to_col[c->lab[saved_color_base[base_pos]]] == saved_color_base[base_pos] &&
                        c->ptn[saved_color_base[base_pos]] + 1 == saved_color_sizes[base_pos]) {
                    return saved_color_base[base_pos];
                }
                int best_score = -1;
                int best_color = -1;
                for (int i = 0; i < c->domain_size;) {
                    if (c->ptn[i] > best_score && c->ptn[i] > 0) {
                        best_score = c->ptn[i];
                        best_color = i;
                    }
                    if(best_color >= 0) break;
                    i += c->ptn[i] + 1;
                }
                return best_color;
            }

            /**
             * @return A selector hook based on the saved base and configured dynamic selector.
             */
            std::function<type_selector_hook> *get_selector_hook() {
                dynamic_seletor = std::bind(&base_selector::dynamic_selector_first, this, std::placeholders::_1,
                                            std::placeholders::_2);
                return &dynamic_seletor;
            }

            big_number get_ir_size_estimate() {
                return ir_tree_size_estimate;
            }

            void find_base(sgraph *g, controller *state, const int h_choose) {
                int perturbe = 0;
                if(h_choose > 6) {
                    perturbe = (int) (hash(h_choose) % 256);
                }
                switch(h_choose % 4) {
                    case 0:
                        find_sparse_optimized_base(g, state, 0);
                        //find_first_base(g, state);
                        break;
                    case 1:
                        find_combinatorial_optimized_base_recurse(g, state, perturbe);
                        break;
                    case 2:
                        find_small_optimized_base(g, state, perturbe);
                        break;
                    case 3:
                        find_combinatorial_optimized_base(g, state);
                        break;
                }


                saved_color_sizes = state->base_color_sz;
                ir_tree_size_estimate.mantissa = 1.0;
                ir_tree_size_estimate.exponent = 0;
                for(auto col_sz : state->base_color_sz) {
                    ir_tree_size_estimate.multiply(col_sz);
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
            void find_sparse_optimized_base(sgraph *g, controller *state, int perturbe) {
                state->mode_write_base();

                int prev_color = -1;
                assert(state->s_base_pos == 0);
                assert(state->base_vertex.size() == 0);
                assert(state->base_color.size() == 0);
                assert(state->base_cells.size() == 0);
                assert(state->base_touched_color_list_pt.size() == 0);
                assert(state->base_singleton_pt.size() == 0);

                test_set.initialize(g->v_size);

                int start_test_from = 0;

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;
                    int test_lim   = 512;

                    // pick previous color if possible
                    if (prev_color >= 0 && state->get_coloring()->ptn[prev_color] > 0) {
                        best_color = prev_color;
                    } else if (prev_color >= 0) { // pick neighbour of previous color if possible
                        const int test_vertex = state->get_coloring()->lab[prev_color];
                        for (int i = 0; i < g->d[test_vertex]; ++i) {
                            const int other_vertex = g->e[g->v[test_vertex] + i];
                            const int other_color = state->get_coloring()->vertex_to_col[other_vertex];
                            if (state->get_coloring()->ptn[other_color] > 0) {
                                best_color = other_color;
                            }
                        }
                    }

                    // heuristic, try to pick "good" color
                    if (best_color == -1) {
                        int best_score = -1;

                        for (int i = start_test_from; i < state->get_coloring()->domain_size && test_lim >= 0;) {
                            const int col_sz = state->get_coloring()->ptn[i];

                            if (col_sz > 0) {
                                --test_lim;
                                const int test_score = color_score(g, state, i);
                                if (test_score > best_score) {
                                    best_color = i;
                                    best_score = test_score;
                                }
                            }

                            if(col_sz == 0 && start_test_from == i) start_test_from += col_sz + 1;
                            i += col_sz + 1;
                        }
                    }

                    assert(best_color >= 0);
                    assert(best_color < g->v_size);
                    prev_color = best_color;
                    state->move_to_child(g, state->get_coloring()->lab[best_color]);
                }

                saved_color_base  = state->base_color;
            }

            void find_split_driven_base(sgraph *g, controller *state, int perturbe) {
                state->mode_write_base();

                int prev_color = -1;
                assert(state->s_base_pos == 0);
                assert(state->base_vertex.size() == 0);
                assert(state->base_color.size() == 0);
                assert(state->base_cells.size() == 0);
                assert(state->base_touched_color_list_pt.size() == 0);
                assert(state->base_singleton_pt.size() == 0);

                test_set.initialize(g->v_size);

                int start_test_from    = state->get_touched_colors()->cur_pos;
                int start_test_from_cl = 0;

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;
                    int best_score = INT32_MIN;
                    int test_lim   = 64;

                    // pick previous color if possible
                    if (prev_color >= 0 && state->get_coloring()->ptn[prev_color] > 0) {
                        best_color = prev_color;
                    } else if (prev_color >= 0) { // pick neighbour of previous color if possible
                        const int test_vertex = state->get_coloring()->lab[prev_color];
                        for (int i = 0; i < g->d[test_vertex]; ++i) {
                            const int other_vertex = g->e[g->v[test_vertex] + i];
                            const int other_color = state->get_coloring()->vertex_to_col[other_vertex];
                            if (state->get_coloring()->ptn[other_color] > 0 &&
                                state->get_touched_colors_test()->get(other_color)) {
                                best_color = other_color;
                            }
                        }
                    }

                    if(best_color == -1) {
                        for(int i = start_test_from; i < state->get_touched_colors()->cur_pos && test_lim >= 0;) {
                            const int col    = (*state->get_touched_colors())[i];
                            const int col_sz = state->get_coloring()->ptn[col];

                            if (col_sz > 0) {
                                --test_lim;
                                const int test_score = color_score_size(state, col);
                                if(test_score > best_score) {
                                    best_color = col;
                                    best_score = test_score;
                                }
                                break;
                            }

                            if(col_sz == 0 && start_test_from == i) start_test_from += 1;
                            i += 1;
                        }
                    }

                    // heuristic, try to pick "good" color
                    if (best_color == -1) {
                        for (int i = start_test_from_cl; i < state->get_coloring()->domain_size && test_lim >= 0;) {
                            const int col_sz = state->get_coloring()->ptn[i];

                            if (col_sz > 0) {
                                --test_lim;
                                const int test_score = color_score_size(state, i);
                                if (test_score > best_score) {
                                    best_color = i;
                                    best_score = test_score;
                                }
                            }

                            if(col_sz == 0 && start_test_from_cl == i) start_test_from_cl += col_sz + 1;
                            i += col_sz + 1;
                        }
                    }

                    assert(best_color >= 0);
                    assert(best_color < g->v_size);
                    prev_color = best_color;
                    state->move_to_child(g, state->get_coloring()->lab[best_color]);
                }

                saved_color_base  = state->base_color;
            }

            void find_first_base(sgraph *g, controller *state) {
                state->mode_write_base();

                int prev_color = -1;
                assert(state->s_base_pos == 0);
                assert(state->base_vertex.size() == 0);
                assert(state->base_color.size() == 0);
                assert(state->base_cells.size() == 0);
                assert(state->base_touched_color_list_pt.size() == 0);
                assert(state->base_singleton_pt.size() == 0);

                test_set.initialize(g->v_size);

                int start_test_from = 0;

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;
                    int test_lim   = 512;


                    for (int i = start_test_from; i < state->get_coloring()->domain_size && test_lim >= 0;) {
                        const int col_sz = state->get_coloring()->ptn[i];

                        if (col_sz > 0) {
                                best_color = i;
                                break;
                        }

                        if(col_sz == 0 && start_test_from == i) start_test_from += col_sz + 1;
                        i += col_sz + 1;
                    }

                    assert(best_color >= 0);
                    assert(best_color < g->v_size);
                    prev_color = best_color;
                    state->move_to_child(g, state->get_coloring()->lab[best_color]);
                }

                saved_color_base  = state->base_color;
            }

            /**
             * Find a base/selector for a given graph \p g from state \p state. Attemtps to find a good base for simple,
             * sparse graphs.
             *
             * @param R A refinement workspace.
             * @param g The graph.
             * @param state The IR state from which a base is created.
             */
            void find_sparse_optimized_base_recurse(sgraph *g, controller *state, int perturbe) {
                state->mode_write_base();

                test_set.initialize(g->v_size);
                int prev_color    = -1;
                int prev_color_sz = 0;

                assert(state->s_base_pos == 0);
                assert(state->base_vertex.size() == 0);
                assert(state->base_color.size() == 0);
                assert(state->base_cells.size() == 0);
                assert(state->base_touched_color_list_pt.size() == 0);
                assert(state->base_singleton_pt.size() == 0);

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;
                    int test_lim   = 512;

                    // recurse into previous color if possible
                    if(prev_color >= 0) {
                        for (int i = 0; i < prev_color + prev_color_sz;) {
                            if (state->get_coloring()->ptn[i] > 0) {
                                best_color = i;
                                break;
                            }
                            i += state->get_coloring()->ptn[i] + 1;
                        }
                    }

                    if (prev_color >= 0 && best_color == -1) { // pick neighbour of previous color if possible
                        const int test_vertex = state->get_coloring()->lab[prev_color];
                        for (int i = 0; i < g->d[test_vertex]; ++i) {
                            const int other_vertex = g->e[g->v[test_vertex] + i];
                            const int other_color = state->get_coloring()->vertex_to_col[other_vertex];
                            if (state->get_coloring()->ptn[other_color] > 0) {
                                best_color = other_color;
                            }
                        }
                    }

                    // heuristic, try to pick "good" color
                    if (best_color == -1) {
                        int best_score = -1;

                        for (int i = 0; i < state->get_coloring()->domain_size && test_lim >= 0;) {
                            if (state->get_coloring()->ptn[i] > 0) {
                                int test_score = color_score(g, state, i);
                                if (test_score > best_score) {
                                    best_color = i;
                                    best_score = test_score;
                                }
                                --test_lim;
                            }

                            i += state->get_coloring()->ptn[i] + 1;
                        }
                    }

                    assert(best_color >= 0);
                    assert(best_color < g->v_size);
                    prev_color    = best_color;
                    prev_color_sz = state->get_coloring()->ptn[best_color] + 1;
                    state->move_to_child(g, state->get_coloring()->lab[(best_color +
                                                                        (perturbe%state->get_coloring()->ptn[best_color]))]);
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
            void find_small_optimized_base(sgraph *g, controller *state, int perturbe) {
                state->mode_write_base();

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;
                    int test_lim = 256;

                    // heuristic, try to pick "good" color
                    int best_score = INT32_MIN;
                    for (int i = 0; i < state->get_coloring()->domain_size && test_lim >= 0;) {
                        if (state->get_coloring()->ptn[i] > 0) {
                            int test_score = color_score_anti_size(state, i);
                            if (test_score > best_score || best_color == -1) {
                                best_color = i;
                                best_score = test_score;
                                --test_lim;
                                if(best_score == INT32_MAX - 1) break;
                            }
                        }

                        i += state->get_coloring()->ptn[i] + 1;
                    }

                    assert(best_color >= 0);
                    assert(best_color < g->v_size);
                    state->move_to_child(g, state->get_coloring()->lab[(best_color +
                                                                    (perturbe%state->get_coloring()->ptn[best_color]))]);
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
            void find_combinatorial_optimized_base(sgraph *g, controller *state, int perturbe) {
                state->mode_write_base();

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;

                    // heuristic, try to pick "good" color
                    int best_score = -1;
                    int skip = perturbe % 3;
                    for (int i = 0; i < state->get_coloring()->domain_size;) {
                        if (state->get_coloring()->ptn[i] > 0) {
                            int test_score = color_score_size(state, i);
                            if (test_score > best_score) {
                                best_color = i;
                                best_score = test_score;
                            }
                        }
                        i += state->get_coloring()->ptn[i] + 1;
                    }

                    assert(best_color >= 0);
                    assert(best_color < g->v_size);
                    state->move_to_child(g, state->get_coloring()->lab[best_color]);
                }

                saved_color_base = state->base_color;
            }

            void find_combinatorial_optimized_base(sgraph *g, controller *state) {
                state->mode_write_base();

                test_set.initialize(g->v_size);
                candidates.clear();
                candidates.reserve(locked_lim);

                mark_set neighbour_color;
                neighbour_color.initialize(g->v_size);

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;

                    // heuristic, try to pick "good" color
                    candidates.clear();
                    int best_score = -1;
                    for (int i = 0; i < state->get_coloring()->domain_size;) {
                        if (state->get_coloring()->ptn[i] > 0) {
                            candidates.push_back(i);
                        }
                        i += state->get_coloring()->ptn[i] + 1;
                    }
                    while (!candidates.empty()) {
                        const int test_color = candidates.back();
                        candidates.pop_back();

                        int test_score = color_score_size(state, test_color);
                        if (neighbour_color.get(test_color)) {
                            test_score *= 10;
                        }
                        if (test_score >= best_score) {
                            best_color = test_color;
                            best_score = test_score;
                        }
                    }

                    assert(best_color >= 0);
                    assert(best_color < g->v_size);
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
            void find_combinatorial_optimized_base_recurse(sgraph *g, controller *state, int perturbe) {
                state->mode_write_base();

                int prev_color    = -1;
                int prev_color_sz = 0;

                while (state->get_coloring()->cells != g->v_size) {
                    int best_color = -1;
                    int best_score = -1;

                    // recurse into previous color if possible
                    if(prev_color >= 0) {
                        for (int i = 0; i < prev_color + prev_color_sz;) {
                            const int col_sz = state->get_coloring()->ptn[i];
                            int test_score = color_score_size(state, i);
                            if (test_score >= best_score && col_sz > 0) {
                                best_color = i;
                                best_score = test_score;
                            }
                            i += state->get_coloring()->ptn[i] + 1;
                        }
                    }

                    if(best_color == -1) {
                        // heuristic, try to pick "good" color
                        for (int i = 0; i < state->get_coloring()->domain_size;) {
                            const int col_sz = state->get_coloring()->ptn[i];
                            int test_score = color_score_size(state, i);
                            if (test_score > best_score && col_sz > 0) {
                                best_color = i;
                                best_score = test_score;
                            }

                            i += col_sz + 1;
                        }
                    }

                    assert(best_color >= 0);
                    assert(best_color < g->v_size);
                    prev_color    = best_color;
                    prev_color_sz = state->get_coloring()->ptn[best_color] + 1;
                    state->move_to_child(g, state->get_coloring()->lab[(best_color + (perturbe% prev_color_sz))]);
                }

                saved_color_base = state->base_color;
            }
        };

        /**
         *  \brief Store deviations for a BFS level
         */
        class deviation_map {
        private:
            std::unordered_set<long> deviation_map;
            int computed_for_base = 0;
            int expected_for_base = 0;
            bool deviation_done = false;

            void check_finished() {
                if(computed_for_base == expected_for_base) deviation_done = true;
                assert(computed_for_base <= expected_for_base);
            }

        public:
            void start(const int h_expected_for_base) {
                computed_for_base = 0;
                expected_for_base = h_expected_for_base;

                deviation_map.clear();
                deviation_done = false;
            }

            void record_deviation(long deviation) {
                deviation_map.insert(deviation);
                ++computed_for_base;
                check_finished();
            }

            void record_no_deviation() {
                ++computed_for_base;
                check_finished();
            }

            bool check_deviation(long deviation) {
                return !deviation_done || deviation_map.contains(deviation);
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
                lab_or_base.allocate((int) arr.size());
                std::copy(arr.begin(), arr.end(), lab_or_base.get_array());
                lab_or_base.set_size((int) arr.size());
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
            int s_leaves          = 0;   /**< number of leaves stored */
            int h_full_save_limit = 100; /**< number of leaves which will be stored fully */

            shared_leaves() {
                leaf_store.reserve(20);
            }

            /**
             * Free up all memory.
             */
            ~shared_leaves() {
                for(auto & l : garbage_collector) {
                    delete l;
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
                const bool full_save = s_leaves < h_full_save_limit;
                auto type
                       = full_save?stored_leaf::stored_leaf_type::STORE_LAB:stored_leaf::stored_leaf_type::STORE_BASE;
                auto new_leaf
                       = full_save?new stored_leaf(c.lab,c.domain_size, type):new stored_leaf(base, type);
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
                // TODO clear garbage collector?
            }
        };

        class tree_node {
            std::mutex    lock;
            limited_save* data;
            tree_node*    next;
            tree_node*    parent;
            bool          is_base = false;
            bool          is_pruned = false;
            long          hash = 0;
        public:
            tree_node(limited_save* data, tree_node* next, tree_node* parent) {
                this->data = data;
                this->next = next;
                if(next == nullptr) {
                    next = this; // TODO supposed to be this->next?
                }
                this->parent = parent;
            }
            tree_node* get_next() {
                return next;
            }
            [[maybe_unused]] tree_node* get_parent() {
                return parent;
            }
            void set_next(tree_node* new_next) {
                this->next = new_next;
            }
            limited_save* get_save() {
                return data;
            }

            void prune() {
                is_pruned = true;
            }
            [[nodiscard]] bool get_prune() const {
                return is_pruned;
            }
            void add_hash(long add) {
                this->hash += add;
            }

            [[nodiscard]] long get_hash() const {
                return hash;
            }

            void base() {
                is_base = true;
            }
            [[nodiscard]] bool get_base() const {
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

            std::vector<unsigned long> node_invariant; // TODO: move this to inprocessor

            bool init = false;
        public:
            groups::orbit h_bfs_top_level_orbit;
            int h_bfs_automorphism_pw = 0;
            shared_leaves stored_leaves;    /**< stores leaves of the IR tree */
            deviation_map stored_deviation; /**< stores trace deviations of a BFS level*/

            explicit shared_tree(int domain_size) {
                h_bfs_top_level_orbit.initialize(domain_size);
            };

            // TODO: move this to inprocessor
            void make_node_invariant(int domain_size) {
                if(finished_up_to > 1) {
                    for(int j = finished_up_to; j >= 1; --j) {
                        finish_level(j);
                        for(auto node : tree_data_jump_map[j]) {
                            //node->get_parent()->add_hash(node->get_hash() + 1);
                            const int v = node->get_save()->get_base()[0];
                            //node_invariant[v] += 1;//+node->get_hash();
                            node_invariant[v] += hash(j);
                        }
                    }
                    /*for(int i = 0; i < domain_size; ++i) {
                        auto node = tree_data_jump_map[1][i];
                        const int v = node->get_save()->get_base()[0];
                        node_invariant[v] += node->get_hash();
                    }*/
                }
            }

            std::vector<unsigned long>* get_node_invariant() {
                return &node_invariant;
            }

            void initialize(std::vector<int> &base, ir::limited_save* root) {
                tree_data.resize(base.size() + 1);
                tree_level_size.resize(base.size() + 1);
                tree_data_jump_map.resize(base.size() + 1);
                add_node(0, root, nullptr, true);
                node_invariant.resize(root->get_coloring()->domain_size);
                current_base = base;
                init = true;
            }

            void clear_leaves() {
                stored_leaves.clear();
            }

            /**
             * @return How many leaves were stored during random search.
             */
            int stat_leaves() const {
                return stored_leaves.s_leaves;
            }

            bool reset(std::vector<int> &new_base, ir::limited_save* root, bool keep_old) {
                if(!init) {
                    initialize(new_base, root);
                    return false;
                }

                const int old_size = (int) current_base.size();
                const int new_size = (int) new_base.size();

                // compare with stored base, keep whatever is possible
                int keep_until = 0;
                if(keep_old) {
                    for (; keep_until < old_size && keep_until < new_size; ++keep_until) {
                        if (current_base[keep_until] != new_base[keep_until]) break;
                    }
                }

                if(keep_until == 0) {
                    if(garbage_collector.size() > 1024) {
                        for(auto & i : garbage_collector) delete i;
                        garbage_collector.clear();
                    }
                }


                if(keep_until == new_size && new_size == old_size) return false;

                finished_up_to = std::min(keep_until, finished_up_to);

                tree_data.resize(new_size + 1);
                tree_level_size.resize(new_size + 1);
                tree_data_jump_map.resize(new_size + 1);

                for (int i = keep_until+1; i < new_size+1; ++i) {
                    tree_data[i] = nullptr;
                    tree_level_size[i] = 0;
                    tree_data_jump_map[i].clear();
                }

                if(keep_until == 0) {
                    tree_level_size[0] = 0;
                    tree_data_jump_map[0].clear();
                    tree_data[0] = nullptr;
                    add_node(0, root, nullptr, true);
                }
                assert(missing_nodes.empty());

                current_base = new_base;

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

            void record_add_invariant(int v, long inv) {
                node_invariant[v] += inv;
            }

            void add_node(int level, limited_save* data, tree_node* parent, bool is_base = false) {
                // TODO use locks
                if(tree_data[level] == nullptr) {
                    tree_level_size[level] = 0;
                    tree_data[level] = new tree_node(data, nullptr, parent);
                    tree_data[level]->set_next(tree_data[level]);

                    garbage_collector.push_back( tree_data[level]);
                    if(is_base) tree_data[level]->base();
                } else {
                    tree_node* a_node    = tree_data[level];
                    tree_node* next_node = a_node->get_next();
                    auto       new_node  = new tree_node(data, next_node, parent);
                    garbage_collector.push_back( new_node);
                    if(is_base) new_node->base();
                    a_node->set_next(new_node);
                    tree_data[level] = new_node;
                }
                ++tree_level_size[level];
            }

            void finish_level(int level) {
                if(tree_data_jump_map[level].empty()) {
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
                num = num % tree_level_size[level];
                return tree_data_jump_map[level][num];
            }

            [[nodiscard]] int get_finished_up_to() const {
                return finished_up_to;
            }

            void set_finished_up_to(const int new_finished_up_to) {
                this->finished_up_to = new_finished_up_to;
            }

            tree_node* get_level(int level) {
                return tree_data[level];
            }

            int get_current_level_size() {
                return tree_level_size[finished_up_to];
            }

            int get_current_level_tracepos() {
                return tree_data[finished_up_to]->get_save()->get_trace_position();
            }

            int get_level_size(int level) {
                return tree_level_size[level];
            }

            ~shared_tree() {
                for(auto & i : garbage_collector) delete i;
            };
        };
    }
}

#endif //DEJAVU_IR_H
