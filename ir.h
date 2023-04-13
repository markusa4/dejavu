#ifndef DEJAVU_IR_H
#define DEJAVU_IR_H

#include "refinement.h"
#include "bijection.h"
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
                this->c.copy(&c);
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
            coloring *c = nullptr;
            trace *T = nullptr;
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

            std::function<type_split_color_hook> my_split_hook;
            std::function<type_worklist_color_hook> my_worklist_hook;

            // settings
            ir_mode mode = IR_MODE_RECORD_TRACE;

            // heuristics
            bool h_last_refinement_singleton_only = true;
            int h_hint_color = -1;
            bool h_hint_color_is_singleton_now = true;
            bool h_cell_active = false;
            bool h_individualize = false;

            int base_pos = 0;
        private:
            void touch_initial_colors() {
                int i = 0;
                while (i < c->lab_sz) {
                    touched_color.set(i);
                    i += c->ptn[i] + 1;
                }
            }

        public:
            controller(coloring *c, trace *T) {
                this->c = c;
                this->T = T;

                touched_color.initialize(c->lab_sz);
                touched_color_list.initialize(c->lab_sz);
                prev_color_list.initialize(c->lab_sz);

                touch_initial_colors();

                my_split_hook = self_split_hook();
                my_worklist_hook = self_worklist_hook();
            }

            // TODO incomplete state save for BFS & random reset
            void save_reduced_state(reduced_save &state) {
                state.set_state(base_vertex, *c, T->get_hash(), T->get_position(), base_pos);
            }

            // TODO incomplete state load for BFS & random reset
            void load_reduced_state(reduced_save &state) {
                c->copy(state.get_coloring());
                T->set_hash(state.get_invariant_hash());
                T->set_position(state.get_trace_position());
                base_pos = state.get_base_position();
                base_vertex = state.get_base();

                // deactivate reversability
                mode = IR_MODE_RECORD_HASH_IRREVERSIBLE;
            }

            coloring *get_coloring() {
                return c;
            }

            int get_base_pos() {
                return base_pos;
            }

            int set_base_pos(int pos) {
                base_pos = pos;
            }

            bool split_hook(const int old_color, const int new_color, const int new_color_sz) {
                // update some heuristic values
                if (new_color_sz > 1) {
                    h_last_refinement_singleton_only = false;
                    h_hint_color = new_color;
                    h_hint_color_is_singleton_now = false;
                }
                if (new_color == h_hint_color && new_color_sz == 1) {
                    h_hint_color_is_singleton_now = true;
                }

                if (mode != IR_MODE_RECORD_HASH_IRREVERSIBLE) {
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
                if (T && !h_individualize) T->op_refine_cell_record(new_color, new_color_sz, 1);

                return true;
            }

            bool worklist_hook(const int color, const int color_sz) {
                if (h_cell_active) {
                    if (T) T->op_refine_cell_end();
                    h_cell_active = false;
                }

                // update some heuristic values
                if (T) {
                    if (T->trace_equal() && !T->blueprint_is_next_cell_active()) {
                        if (config.CONFIG_IR_IDLE_SKIP) {
                            T->blueprint_skip_to_next_cell();
                            return false;
                        }
                    }
                }

                if (T) T->op_refine_cell_start(color);

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
                const int init_color_class = R->individualize_vertex(c, v, my_split_hook);
                T->op_individualize(c->vertex_to_col[v]);
                h_individualize = false;

                h_last_refinement_singleton_only = true;

                if (T) T->op_refine_start();

                if (mode == IR_MODE_RECORD_TRACE) {
                    R->refine_coloring(g, c, init_color_class, -1, my_split_hook, my_worklist_hook);
                    if (T && h_cell_active) T->op_refine_cell_end();
                    if (T) T->op_refine_end();
                } else {
                    // TODO compare_base_cells not necessarily applicable
                    R->refine_coloring(g, c, init_color_class, compare_base_cells[base_pos - 1], my_split_hook,
                                       my_worklist_hook);
                    if (T) T->skip_to_individualization();
                }

                h_cell_active = false;

                if (mode != IR_MODE_RECORD_HASH_IRREVERSIBLE) base_cells.push_back(c->cells);
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
            }
        };

        struct splitmap {

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

            /**
             * Find a base/selector for a given graph \p g from state \p state. Attemtps to find a good base for simple,
             * sparse graphs.
             *
             * @param R A refinement workspace.
             * @param g The graph.
             * @param state The IR state from which a base is created.
             */
            void find_sparse_optimized_base(refinement *R, sgraph *g, controller *state) {
                state->T->set_compare(false);
                state->T->set_record(true);

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

            /**
             * Find a base/selector for a given graph \p g from state \p state. Attemtps to find a good base for
             * combinatorial graphs solved by bfs/random walks (i.e., by choosing large colors).
             *
             * @param R A refinement workspace.
             * @param g The graph.
             * @param state The IR state from which a base is created.
             */
            void find_combinatorial_optimized_base(refinement *R, sgraph *g, controller *state) {
                state->T->set_compare(false);
                state->T->set_record(true);

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


            int state_score(sgraph *g, controller *state) {
                return state->get_coloring()->cells;
            }

            splitmap find_splitmap(sgraph *g, coloring *c) {
                return splitmap();
            }
        };

        class tree_node {
            std::mutex    lock;
            reduced_save* data;
            tree_node*    next;
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
        };

        // TODO tree structure for BFS + random walk
        class tree {
            std::vector<tree_node*> tree_data;
            std::vector<int>        tree_level_size;
        public:
            void initialize(int base_size) {
                tree_data.resize(base_size);
                tree_level_size.resize(base_size);
            }

            void add_node(int level, reduced_save* data) {
                // TODO use locks
                ++tree_level_size[level];
                if(tree_data[level] == nullptr) {
                    tree_data[level] = new tree_node(data, nullptr);
                } else {
                    tree_node* a_node    = tree_data[level];
                    tree_node* next_node = a_node->get_next();
                    auto       new_node  = new tree_node(data, next_node);
                    a_node->set_next(new_node);
                    tree_data[level] = new_node;
                }
            }

            ~tree() {
                //TODO: write this
            };
        };
    }
}

#endif //DEJAVU_IR_H
