#ifndef DEJAVU_DFS_H
#define DEJAVU_DFS_H

#include <random>
#include <chrono>
#include "refinement.h"
#include "bijection.h"
#include "coloring.h"
#include "sgraph.h"
#include "trace.h"

namespace dejavu {

    /**
     * \brief Mode of trace for IR search
     *
     * The ir_mode determines in which mode the trace is used: whether a new trace is recorded, or whether the current
     * computation is compared to a stored trace.
     *
     */
    enum ir_mode { IR_MODE_COMPARE_TRACE, IR_MODE_RECORD_TRACE, IR_MODE_RECORD_HASH_IRREVERSIBLE};

    /**
     * int type_selector_hook(coloring* c, const int base_pos);
     */
    typedef int type_selector_hook(const coloring*, const int);


    static void progress_print_header() {

        PRINT("________________________________________________________________");
        PRINT(std::setw(16) << std::left <<"T (ms)"                                  << std::setw(16) << "proc"  << std::setw(16) << "P1"        << std::setw(16)        << "P2");
        PRINT("________________________________________________________________");
        PRINT(std::setw(16) << std::left << 0 << std::setw(16) << "start" << std::setw(16) << "_" << std::setw(16) << "_" );
    }

    static void progress_print(const std::string proc, const std::string p1, const std::string p2) {
        static std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
        PRINT(std::setw(16) << std::left << (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count()) / 1000000.0  << std::setw(16) << proc << std::setw(16) << p1 << std::setw(16) << p2);
    }

    // reset internal automorphism structure to the identity
    static void  reset_automorphism(int* rautomorphism, work_list* automorphism_supp) {
        for(int i = 0; i < automorphism_supp->cur_pos; ++i) {
            rautomorphism[(*automorphism_supp)[i]] = (*automorphism_supp)[i];
        }
        automorphism_supp->reset();
    };

    // make automorphism from colorings
    static void color_diff_automorphism(int n, const int* vertex_to_col, const int* col_to_vertex,
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

    /**
     * \brief Workspace for sparse automorphisms
     *
     * Enables O(1) lookup on a sparse automorphism by using an O(n) workspace.
     */
    class sparse_automorphism_workspace {
        work_list automorphism;
        work_list automorphism_supp;
        int domain_size;

        bool support01 = false;
    public:
        /**
         * Initializes the stored automorphism to the identity.
         *
         * @param domain_size Size of the domain on which automorphisms operate
         */
        sparse_automorphism_workspace(int domain_size) {
            automorphism.initialize(domain_size);
            for(int i = 0; i < domain_size; ++i)
                automorphism[i] = i;
            automorphism_supp.initialize(domain_size);
            this->domain_size = domain_size;
        }

        void set_support01(bool support01) {
            this->support01 = support01;
        }

        /**
         * Create automorphism from two given vertex colorings.
         *
         * @param vertex_to_col a vertex-to-color mapping which describes the first coloring
         * @param col_to_vertex a color-to-vertex mapping which describes the second coloring
         */
        void write_color_diff(const int* vertex_to_col, const int* col_to_vertex) {
            color_diff_automorphism(domain_size, vertex_to_col, col_to_vertex, automorphism.get_array(),
                                    &automorphism_supp);
        }

        /**
         * Apply another automorphism to the stored automorphism. Closely follows the implementation in nauty / Traces.
         *
         * @param other
         * @param pwr
         */
        void apply(work_list& scratch_apply1, work_list& scratch_apply2, mark_set& scratch_apply3, sparse_automorphism_workspace* other, int pwr = 1) {
            apply(scratch_apply1, scratch_apply2, scratch_apply3, other->perm(), pwr);
        }

        void __attribute__ ((noinline)) update_support() {
            // rewrite support
            if(!support01) {
                automorphism_supp.reset();
                for (int i = 0; i < domain_size; ++i) {
                    if (i != automorphism[i])
                        automorphism_supp.push_back(i);
                }
            } else {
                automorphism_supp.reset();
                int i;
                for (i = 0; i < domain_size; ++i) {
                    if (i != automorphism[i]) break;
                }
                automorphism_supp.cur_pos = (i != domain_size);
            }
        }

        /**
         * Apply another automorphism to the stored automorphism. Closely follows the implementation in nauty / Traces.
         *
         * @param other
         * @param pwr
         */
        void __attribute__ ((noinline)) apply(work_list& scratch_apply1, work_list& scratch_apply2, mark_set& scratch_apply3,
                                              const int* p, int pwr = 1) {
            if(pwr == 0)
                return;
            if(pwr <= 5) {
                if (pwr == 1)
                    for (int i = 0; i < domain_size; ++i) automorphism[i] = p[automorphism[i]];
                else if (pwr == 2)
                    for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[automorphism[i]]];
                else if (pwr == 3)
                    for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[p[automorphism[i]]]];
                else if (pwr == 4)
                    for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[p[p[automorphism[i]]]]];
                else if (pwr == 5)
                    for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[p[p[p[automorphism[i]]]]]];
            } else if(pwr <= 19) {
                // apply other automorphism
                if (pwr >= 6) {
                    for (int j = 0; j < domain_size; ++j) {
                        scratch_apply1[j] = p[p[p[j]]];
                    }
                    for (; pwr >= 6; pwr -= 6)
                        for (int j = 0; j < domain_size; ++j)
                            automorphism[j] = scratch_apply1[scratch_apply1[automorphism[j]]];
                }

                if (pwr == 1)
                    for (int i = 0; i < domain_size; ++i) automorphism[i] = p[automorphism[i]];
                else if (pwr == 2)
                    for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[automorphism[i]]];
                else if (pwr == 3)
                    for (int i = 0; i < domain_size; ++i) automorphism[i] = scratch_apply1[automorphism[i]];
                else if (pwr == 4)
                    for (int i = 0; i < domain_size; ++i) automorphism[i] = p[scratch_apply1[automorphism[i]]];
                else if (pwr == 5)
                    for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[scratch_apply1[automorphism[i]]]];
            } else {
                // 1 cycle at a time

                scratch_apply3.reset();
                for (int i = 0; i < domain_size; ++i) {
                    if (scratch_apply3.get(i)) continue;
                    if (p[i] == i)
                        scratch_apply2[i] = i;
                    else {
                        int cyclen = 1;
                        scratch_apply1[0] = i;
                        for (int j = p[i]; j != i; j = p[j]) {
                            scratch_apply1[cyclen++] = j;
                            scratch_apply3.set(j);
                        }
                        int kk = pwr % cyclen;
                        for (int j = 0; j < cyclen; ++j) {
                            scratch_apply2[scratch_apply1[j]] = scratch_apply1[kk];
                            if (++kk == cyclen) kk = 0;
                        }
                    }
                }
                for (int i = 0; i < domain_size; ++i) automorphism[i] = scratch_apply2[automorphism[i]];
            }

            // rewrite support
            update_support();
        }


        /**
         * Create automorphism from two canonically-ordered vectors of singletons. The resulting automorphism maps
         * \p singletons1[i] to \p singletons2[i] for i in \p pos_start, ..., \p pos_end.
         *
         * @param singletons1 first vector of singletons
         * @param singletons2 second vector of singletons
         * @param pos_start start reading the vectors at this position
         * @param pos_end stop reading the vecvtors at this position.
         */
        void  write_singleton(std::vector<int>* singletons1, std::vector<int>* singletons2, int pos_start, int pos_end) {
            for(int i = pos_start; i < pos_end; ++i) {
                const int from = (*singletons1)[i];
                const int to   = (*singletons2)[i];
                assert(automorphism[from] == from);
                if(from != to) {
                    automorphism_supp.push_back(from);
                    automorphism[from] = to;
                }
            }
        }

        void write_single_map(const int from, const int to) {
            assert(automorphism[from] == from);
            if(from != to) {
                automorphism_supp.push_back(from);
                automorphism[from] = to;
            }
        }

        /**
         * Reset the contained automorphism back to the identity.
         */
        void reset() {
            reset_automorphism(automorphism.get_array(), &automorphism_supp);
        }

        /**
         * @return Integer array \p p describing the stored automorphism, where point v is mapped to \p p[v].
         */
        int* perm() {
            return automorphism.get_array();
        }

        /**
         * @return Integer array which contains all vertices in the support of the contained automorphism.
         */
        int* support() {
            return automorphism_supp.get_array();
        }

        /**
         * @return Size of the support.
         */
        int nsupport() {
            return automorphism_supp.cur_pos;
        }
    };


    /**
     * \brief Reduced IR save state
     *
     * Using this class, partial information of a state of an IR computation can be stored. Using this information,
     * IR computations can be resumed from this state either using BFS or random walks. The state in particular does not
     * keep enough information to resume using DFS.
     */
    class ir_reduced_save {
        // TODO enable a "compressed state" only consisting of base (but code outside should be oblivious to it)
        // TODO this is only supposed to be an "incomplete" state -- should there be complete states?

        std::vector<int> base_vertex; /**< base of vertices of this IR node (optional) */
        coloring c; /**< vertex coloring of this IR node */
        long     invariant = 0; /**< hash of invariant of this IR node */
        int      trace_position = 0; /**< position of trace of this IR node */
        int      base_position = 0; /**< length of base of this IR node */
    public:
        void set_state(std::vector<int>& base_vertex, coloring& c, long invariant, int trace_position, int base_position) {
            this->base_vertex = base_vertex;
            this->c.copy(&c);
            this->invariant = invariant;
            this->trace_position = trace_position;
            this->base_position  = base_position;
        }

        coloring* get_coloring() {
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
    struct ir_controller {
        coloring*    c = nullptr;
        trace*       T = nullptr;
        mark_set     touched_color;
        work_list    touched_color_list;
        work_list    prev_color_list;
        std::vector<int>  singletons;
        std::vector<int>  base_vertex;
        std::vector<int>  base_color;
        std::vector<int>  base_color_sz;
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
        void save_reduced_state(ir_reduced_save& state) {
            state.set_state(base_vertex, *c, T->get_hash(), T->get_position(), base_pos);
        }

        // TODO incomplete state load for BFS & random reset
        void load_reduced_state(ir_reduced_save& state) {
            c->copy(state.get_coloring());
            T->set_hash(state.get_invariant_hash());
            T->set_position(state.get_trace_position());
            base_pos = state.get_base_position();

            // deactivate reversability
            mode = IR_MODE_RECORD_HASH_IRREVERSIBLE;
        }

        coloring* get_coloring() {
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
            if(new_color_sz > 1) {
                h_last_refinement_singleton_only = false;
                h_hint_color = new_color;
                h_hint_color_is_singleton_now = false;
            }
            if(new_color == h_hint_color && new_color_sz == 1) {
                h_hint_color_is_singleton_now = true;
            }

            if(mode != IR_MODE_RECORD_HASH_IRREVERSIBLE) {
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

        /**
         * Move IR node kept in this controller to a child, specified by a vertex to be individualized.
         *
         * @param R a refinement workspace
         * @param g the graph
         * @param v the vertex to be individualized
         */
        void move_to_child(refinement* R, sgraph* g, int v) {
            ++base_pos;

            if(mode != IR_MODE_RECORD_HASH_IRREVERSIBLE) {
                base_singleton_pt.push_back(singletons.size());
                base_vertex.push_back(v);
                base_color.push_back(c->vertex_to_col[v]);
                base_color_sz.push_back(c->ptn[c->vertex_to_col[v]]+1);
                base_touched_color_list_pt.push_back(touched_color_list.cur_pos);
            }

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
                // TODO compare_base_cells not necessarily applicable
                R->refine_coloring(g, c, init_color_class, compare_base_cells[base_pos-1], my_split_hook, my_worklist_hook);
                if(T) T->skip_to_individualization();
            }

            h_cell_active = false;

            if(mode != IR_MODE_RECORD_HASH_IRREVERSIBLE) base_cells.push_back(c->cells);
        }

        /**
         * Move IR node kept in this controller back to its parent.
         */
        void move_to_parent() {
            assert(mode != IR_MODE_RECORD_HASH_IRREVERSIBLE);

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
     * Heuristics which enable the creation of different cell selectors, as well as moving an \ref ir_controller to a
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
        void save_base(ir_controller* state) {

        }

        // TODO configuration of dynamic selector when deviating from base

        /**
         *
         */
         int dynamic_selector(const coloring * c, const int base_pos) {
             if(base_pos >= 0 && base_pos < saved_color_base.size() && c->ptn[saved_color_base[base_pos]] > 0) {
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
        std::function<type_selector_hook>* get_selector_hook() {
            dynamic_seletor = std::bind(&selector_factory::dynamic_selector, this, std::placeholders::_1, std::placeholders::_2);
            return &dynamic_seletor;
        }

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

            saved_color_base = state->base_color;
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

    namespace group_structure {

        typedef void type_unload_hook();

        class dynamic_auto_loader {
            int*                        loaded_automorphism;
            std::function<type_unload_hook>* unload_hook;
        public:
            void load(int* automorphism, std::function<type_unload_hook>* unload_hook) {
                loaded_automorphism = automorphism;
                this->unload_hook = unload_hook;
            }

            int* p() {
                return loaded_automorphism;
            }

            void unload() {
                if(unload_hook) (*unload_hook)();
            }
        };

        /**
         * \brief Stores an automorphism in a dense or sparse manner, dynamically.
         */
        class stored_automorphism {
        public:
            enum stored_automorphism_type {
                STORE_DENSE, STORE_SPARSE, STORE_BOTH
            };

        private:
            work_list data;
            int domain_size;
            stored_automorphism_type store_type = STORE_SPARSE;

        public:

            stored_automorphism_type get_store_type() {
                return store_type;
            }

            /**
             * Load this stored automorphism into workspace.
             *
             * @param automorphism
             */
            void __attribute__ ((noinline)) load(dynamic_auto_loader& loader, sparse_automorphism_workspace &automorphism) {
                if(store_type == STORE_SPARSE) {
                    automorphism.reset();
                    int first_of_cycle = 0;
                    for (int i = 0; i < data.size(); ++i) {
                        const int j = abs(data[i]) - 1;
                        const bool is_last = data[i] < 0;
                        assert(i == data.size() - 1 ? is_last : true);
                        if (is_last) {
                            automorphism.write_single_map(j, abs(data[first_of_cycle]) - 1);
                            first_of_cycle = i + 1;
                        } else {
                            assert(i + 1 < data.size());
                            automorphism.write_single_map(j, abs(data[i + 1]) - 1);
                        }
                    }
                    loader.load(automorphism.perm(), nullptr);
                } else {
                    loader.load(data.get_array(), nullptr);
                }
            }

            /**
             * Load inverse of stored automorphism into workspace.
             *
             * @param automorphism
             */
            void __attribute__ ((noinline)) load_inverse(dynamic_auto_loader& loader, sparse_automorphism_workspace &automorphism) {
                if(store_type == STORE_SPARSE) {
                    automorphism.reset();
                    int first_of_cycle = data.size() - 1;
                    for (int i = data.size() - 1; i >= 0; --i) {
                        if (data[i] < 0)
                            first_of_cycle = i;

                        const int j = abs(data[i]) - 1;
                        const bool is_last = i == 0 || (data[i - 1] < 0);
                        if (is_last) {
                            automorphism.write_single_map(j, abs(data[first_of_cycle]) - 1);
                        } else {
                            automorphism.write_single_map(j, abs(data[i - 1]) - 1);
                        }
                    }
                    loader.load(automorphism.perm(), nullptr);
                } else {
                    automorphism.reset();
                    for(int i = 0; i < domain_size; ++i) {
                        if(i != data[i]) {
                            automorphism.write_single_map(data[i], i);
                        }
                    }
                    loader.load(automorphism.perm(), nullptr);
                }
            }

            /**
             * Store the given automorphism workspace.
             *
             * @param automorphism
             */
            void __attribute__ ((noinline)) store(int domain_size, sparse_automorphism_workspace& automorphism, mark_set& helper) {
                this->domain_size = domain_size;
                assert(data.empty());

                int support = 0;
                for (int i = 0; i < domain_size; ++i) support += (automorphism.perm()[i] != i);

                if(support < domain_size / 4) { // domain_size / 4
                    store_type = STORE_SPARSE;
                    helper.reset();

                    data.initialize(support);
                    for (int i = 0; i < domain_size; ++i) {
                        if(automorphism.perm()[i] == i) continue;
                        const int j = i;
                        if (helper.get(j)) continue;
                        helper.set(j);
                        int map_j = automorphism.perm()[j];
                        assert(map_j != j);
                        while (!helper.get(map_j)) {
                            data.push_back(map_j + 1);
                            helper.set(map_j);
                            map_j = automorphism.perm()[map_j];
                        }
                        assert(map_j == j);
                        data.push_back(-(j + 1));
                    }
                    assert(data.size() == support);
                } else {
                    store_type = STORE_DENSE;
                    data.initialize(domain_size);
                    data.set_size(domain_size);
                    memcpy(data.get_array(), automorphism.perm(), domain_size * sizeof(int));
                    assert(data.size() == domain_size);
                }
            }
        };

        /**
         * A global (thread local) state used for computations in Schreier structures.
         */
        class schreier_workspace {
        public:
            schreier_workspace(int domain_size, refinement * R, sgraph* g) : scratch_auto(domain_size) {
                scratch1.initialize(domain_size);
                scratch2.initialize(domain_size);
                scratch_apply1.initialize(domain_size);
                scratch_apply2.initialize(domain_size);
                scratch_apply3.initialize(domain_size);
                this->R = R;
                this-> g = g;
            }
            dynamic_auto_loader loader;

            mark_set scratch1;
            mark_set scratch2;
            work_list scratch_apply1;
            work_list scratch_apply2;
            mark_set  scratch_apply3;
            sparse_automorphism_workspace scratch_auto;

            // TODO debug
            refinement* R;
            sgraph* g;
        };


        class stored_generators {
            std::mutex lock_generators;          /**< locks the generators */
            std::vector<stored_automorphism*> generators; /** list of generators */
            int domain_size;
        public:
            int s_stored_sparse = 0;
            int s_stored_dense  = 0;

            void setup(int domain_size) {
                this->domain_size = domain_size;
            }

            int add_generator(schreier_workspace& w, sparse_automorphism_workspace& automorphism) {
                lock_generators.lock();
                generators.emplace_back(new stored_automorphism);
                const int num = generators.size()-1;
                assert(w.R->certify_automorphism_sparse(w.g, automorphism.perm(), automorphism.nsupport(), automorphism.support()));
                generators[num]->store(domain_size, automorphism, w.scratch2);
                lock_generators.unlock();

                s_stored_sparse += (generators[num]->get_store_type() == stored_automorphism::stored_automorphism_type::STORE_SPARSE);
                s_stored_dense  += (generators[num]->get_store_type() == stored_automorphism::stored_automorphism_type::STORE_DENSE);

                return num;
            }

            stored_automorphism* get_generator(const int num) {
                // TODO locks...
                return generators[num];
            }

            int size() {
                return generators.size();
            }
        };

        class stored_transversal {
            enum stored_transversal_type {
                STORE_DENSE, STORE_SPARSE, STORE_BOTH
            };

            std::mutex lock_transversal;          /**< locks this transversal */

            int fixed;                            /**< vertex fixed by this transversal */
            int sz_upb = INT32_MAX;               /**< upper bound for size of the transversal (e.g. color class size) */
            int level;
            bool finished = false;

            std::vector<int> fixed_orbit;         /**< contains vertices of orbit at this schreier level */
            std::vector<int> fixed_orbit_to_perm; /**< maps fixed_orbit[i] to generators[i] in class \ref schreier. */
            std::vector<int> fixed_orbit_to_pwr;  /**< power that should to be applied to generators[i] in class \ref schreier. */

            stored_transversal_type store_type = STORE_SPARSE; /**< whether above structures are stored dense or sparse */

            void __attribute__ ((noinline)) load_orbit_to_scratch(schreier_workspace& w) {
                for(int i = 0; i < fixed_orbit.size(); ++i) {
                    w.scratch1.set(fixed_orbit[i]);
                }
            }

            void add_to_fixed_orbit(const int vertex, const int perm, const int pwr) {
                //std::cout << "add " << vertex << ", " << perm << ", " << pwr << std::endl;
                fixed_orbit.push_back(vertex);
                fixed_orbit_to_perm.push_back(perm);
                fixed_orbit_to_pwr.push_back(pwr);
            }


            void __attribute__ ((noinline)) apply_perm(schreier_workspace& w, sparse_automorphism_workspace &automorphism, stored_generators& generators, const int perm, const int pwr) {
                // load perm into workspace
                auto generator = generators.get_generator(perm);

                //if(level == 0) {
                //    std::cout << perm << "^" << pwr << std::endl;
                //}

                if(pwr < 0) {
                    generator->load_inverse(w.loader, w.scratch_auto);

                    // multiply
                    automorphism.apply(w.scratch_apply1, w.scratch_apply2, w.scratch_apply3, w.loader.p(), abs(pwr));
                } else if(pwr > 0) {
                    generator->load(w.loader, w.scratch_auto);

                    // multiply
                    automorphism.apply(w.scratch_apply1, w.scratch_apply2, w.scratch_apply3, w.loader.p(), pwr);
                }
            }

        public:
            int find_point(const int p) {
                for(int i = 0; i < fixed_orbit.size(); ++i) {
                    if(p == fixed_orbit[i]) {
                        return i;
                    }
                }
                return -1;
            }

            void reduce_to_unfinished(schreier_workspace& w, std::vector<int>& selection, int base_pos) {
                load_orbit_to_scratch(w);
                int back_swap = selection.size() - 1;
                int front_pt;
                for(front_pt = 0; front_pt <= back_swap;) {
                    if(!w.scratch1.get(selection[front_pt])) {
                        ++front_pt;
                    } else {
                        selection[front_pt] = selection[back_swap];
                        --back_swap;
                    }
                }
                selection.resize(front_pt);
                w.scratch1.reset();
            }

            bool is_finished() {
                return finished;
            }

            int fixed_point() {
                return fixed;
            }

            void setup(const int fixed_vertex, const int level, const int sz_upb) {
                fixed = fixed_vertex;
                this->level = level;
                this->sz_upb= sz_upb;
                add_to_fixed_orbit(fixed_vertex, -1, 0);
            }
            /**
             * @param generators
             * @param automorphism
             * @return whether transversal was extended
             */
            bool extend_with_automorphism(schreier_workspace& w, stored_generators& generators, sparse_automorphism_workspace &automorphism) {
                if(finished)
                    return false;

                load_orbit_to_scratch(w);
                bool changed = false;
                // TODO probably aquire lock already here
                // TODO how this is performed precisely should depend on fixed_orbit size versus support of generator
                int gen_num = -1;

                int self_apply_after = 0;

                for(int i = 0; i < fixed_orbit.size(); ++i) {
                    int j = automorphism.perm()[fixed_orbit[i]];
                    if(w.scratch1.get(j))
                        continue;

                    int ipwr = 0; // inverted power
                    int npwr = 1; //non-inverted power
                    int jj;
                    for (jj = j; !w.scratch1.get(jj); jj = automorphism.perm()[jj]) ++ipwr;
                    //bool reversible = (jj == fixed_orbit[i]) && (ipwr > 100);

                    /*if(!w.scratch1.get(j) && reversible && i == 0) {
                        assert(fixed_orbit[i] == fixed);
                        self_apply_after = ipwr;
                    }*/

                    //int pwr = 1;
                    while(!w.scratch1.get(j)) {
                        // we change this traversal
                        changed = true;

                        // add generator to generating set (once)
                        if(gen_num == -1) gen_num = generators.add_generator(w, automorphism);

                        // add entry to traversal
                        //if(ipwr <= npwr || !reversible) {
                            add_to_fixed_orbit(j, gen_num, ipwr);
                       // } else {
                        //    add_to_fixed_orbit(j, gen_num, -npwr);
                       // }
                        w.scratch1.set(j);

                        // we check out the entire cycle of j now
                        j = automorphism.perm()[j];
                        --ipwr;
                        ++npwr;
                    }
                }

                /*if(self_apply_after > 0) {
                    assert(changed);
                    assert(gen_num >= 0);
                    apply_perm(w, automorphism, generators, gen_num, self_apply_after);
                    assert(automorphism.perm()[fixed] == fixed);
                }*/

                if(sz_upb == fixed_orbit.size() && !finished) {
                    finished = true;
                }

                w.scratch1.reset();
                return changed;
            }

            /**
             * @param generators
             * @param automorphism
             * @return whether transversal was extended
             */
            bool fix_automorphism(schreier_workspace& w, stored_generators& generators, sparse_automorphism_workspace &automorphism) {
                int fixed_map = automorphism.perm()[fixed];
                while(fixed != fixed_map) {
                    const int pos = find_point(fixed_map);
                    const int perm = fixed_orbit_to_perm[pos];
                    const int pwr  = fixed_orbit_to_pwr[pos];
                    //std::cout << "fixing " << fixed_map << "->" << fixed << " using " << perm << "^" << pwr << std::endl;

                    assert(w.R->certify_automorphism_sparse(w.g, automorphism.perm(), automorphism.nsupport(), automorphism.support()));
                    apply_perm(w, automorphism, generators, perm, pwr);
                    assert(w.R->certify_automorphism_sparse(w.g, automorphism.perm(), automorphism.nsupport(), automorphism.support()));

                    fixed_map = automorphism.perm()[fixed];
                }

                assert(automorphism.perm()[fixed] == fixed);
                return automorphism.nsupport() == 0;
            }


            // TODO: flip to dense over certain threshold -- could also include "semi-dense" with sorted list or hash
        };

        /**
         * \brief Schreier structure with fixed base.
         *
         * Enables sifting of automorphisms into a Schreier structure with fixed base. Can be used across multiple threads
         * in a safe manner, i.e., the structure can lock appropriate parts of itself.
         *
         */
        class schreier {
            // TODO implement sparse schreier-sims for fixed base
            // TODO support for sparse tables, sparse automorphism, maybe mix sparse & dense
            // TODO could use color sizes for memory allocation, predicting dense/sparse?
            // TODO using inverses to reduce pwr seems like low-hanging fruit

            int domain_size;
            int finished_up_to = -1;

            stored_generators generators;
            work_list_t<stored_transversal*> transversals;

        public:

            int stat_sparsegen() {
                return generators.s_stored_sparse;
            }
            int stat_densegen() {
                return generators.s_stored_dense;
            }

            /**
             * Set up this Schreier structure using the given base. The base is then fixed and can not be adjusted
             * later on.
             *
             * @param base the base
             * @param stop integer which indicates to stop reading the base at this position
             */
            void setup(const int domain_size, std::vector<int>& base, std::vector<int>& base_sizes, const int stop) {
                assert(base.size() >= stop);
                this->domain_size = domain_size;
                generators.setup(domain_size);
                transversals.initialize(stop);
                transversals.set_size(stop);
                for(int i = 0; i < stop; ++i) {
                    transversals[i] = new stored_transversal();
                    transversals[i]->setup(base[i], i, base_sizes[i]); // TODO add proper upper bound
                }
            }

            int base_point(int pos) {
                return transversals[pos]->fixed_point();
            }

            int base_size() {
                return transversals.size();
            }

            bool is_in_base_orbit(const int base_pos, const int v) {
                if(base_pos >= transversals.size()) return false;
                assert(base_pos >= 0);
                assert(base_pos < transversals.size());
                const int search = transversals[base_pos]->find_point(v);
                return search != -1;
            }

            void reduce_to_unfinished(schreier_workspace& w, std::vector<int>& selection, int base_pos) {
                transversals[base_pos]->reduce_to_unfinished(w, selection, base_pos);
            }

            bool is_finished(const int base_pos) {
                return transversals[base_pos]->is_finished();
            }

            /**
             * Sift automorphism into the Schreier structure.
             *
             * @param automorphism Automorphism to be sifted. Will be manipulated by the method.
             * @return Whether automorphism was added to the Schreier structure or not.
             */
            bool sift(schreier_workspace& w, sgraph* g, refinement* R, sparse_automorphism_workspace &automorphism) {
                bool changed = false;

                automorphism.set_support01(true);
                for (int level = 0; level < transversals.size(); ++level) {
                    //std::cout << "extending transversal..." << std::endl;
                    changed = transversals[level]->extend_with_automorphism(w, generators, automorphism) || changed;

                    if(finished_up_to == level - 1 && transversals[level]->is_finished()) {
                        ++finished_up_to;
                    }

                    //std::cout << "fixing automorphism..." << std::endl;
                    const bool is_identity = transversals[level]->fix_automorphism(w, generators, automorphism);
                    if (is_identity) break;
                }
                automorphism.set_support01(false);
                automorphism.update_support();
                automorphism.reset();

                //std::cout << generators.size() << ", " << finished_up_to << "/" << transversals.size() << std::endl;
                return changed;
            }

            int finished_up_to_level() {
                return finished_up_to;
            }
        };
    }


    namespace search_strategy {
        /**
         * \brief Breadth-first search.
         */
        class bfs_ir {
            // TODO implement new, simpler bfs
            // TODO depends on ir_tree, and selector
        };

        class stored_leaf {
            coloring store_c;
        public:
            stored_leaf(coloring& c, std::vector<int>& base) {
                // TODO make the stored leaf more properly
                store_c.copy(&c);
            }

            coloring* get_coloring() {
                return &store_c;
            }
        };

        class stored_leafs {
            std::mutex lock;
            std::unordered_multimap<long, stored_leaf*> leaf_store;
            std::vector<stored_leaf*> garbage_collector;

        public:
            int s_leaves = 0;

            ~stored_leafs() {
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
                if(!leaf_store.contains(hash)) {
                    auto new_leaf = new stored_leaf(c, base);
                    leaf_store.insert(std::pair<long, stored_leaf*>(hash, new_leaf));
                    ++s_leaves;
                    lock.unlock();
                } else {
                    lock.unlock();
                }
            }
        };

        /**
         * \brief IR search using random walks.
         *
         * Performs random walks of the IR tree, sifting resulting automorphisms into the given Schreier structure. If the
         * Schreier structure is complete with respect to the base, or the probabilistic abort criterion satisfied, the
         * process terminates. The algorithm guarantees to find all automorphisms up to the specified error bound.
         *
         * Alternatively, a limit for the amount of discovered differing leafs can be set.
         */
        class random_ir {
            std::default_random_engine generator;
            stored_leafs leaf_storage;

            int consecutive_success = 0;
            int error_bound = 5;

        public:
            void setup(int error, int leaf_store_limit) {
            }

            void record_sift_result(const bool changed) {
                if(!changed) {
                    ++consecutive_success;
                } else {
                    if(consecutive_success > 0) {
                        ++error_bound;
                        consecutive_success = 0;
                    }
                }
            }

            bool probabilistic_abort_criterion() {
                return (consecutive_success > error_bound);
            }

            bool deterministic_abort_criterion(group_structure::schreier &group) {
                return (group.finished_up_to_level() + 1 == group.base_size());
            }

            int stat_leaves() {
                return leaf_storage.s_leaves;
            }

            void specific_walk(refinement &R, std::vector<int>& base_vertex, sgraph *g,
                               group_structure::schreier &group, ir_controller &local_state, ir_reduced_save &start_from) {
                local_state.load_reduced_state(start_from);
                int base_pos = local_state.base_pos;

                while (g->v_size != local_state.c->cells) {
                    int v = base_vertex[local_state.base_pos];
                    local_state.move_to_child(&R, g, v);
                    ++base_pos;
                }

                auto other_leaf = leaf_storage.lookup_leaf(local_state.T->get_hash());
                if(other_leaf == nullptr) {
                    std::cout << "adding leaf " << local_state.T->get_hash() << std::endl;
                    leaf_storage.add_leaf(local_state.T->get_hash(), *local_state.c, local_state.base_vertex);
                }
            }

            // TODO implement dejavu strategy, more simple
            // TODO depends on ir_tree, selector, and given base (no need to sift beyond base!)
            // TODO: swap out ir_reduced to weighted IR tree later? or just don't use automorphism pruning on BFS...?
            void random_walks(refinement &R, std::function<type_selector_hook> *selector, sgraph *g,
                              group_structure::schreier &group, ir_controller &local_state, ir_reduced_save &start_from) {
                sparse_automorphism_workspace automorphism(g->v_size);
                group_structure::schreier_workspace w(g->v_size, &R, g);
                std::vector<int> heuristic_reroll;

                while(!probabilistic_abort_criterion() && !deterministic_abort_criterion(group)) {
                    local_state.load_reduced_state(start_from);
                    int could_start_from = group.finished_up_to_level();
                    if(local_state.base_pos < could_start_from) {
                        while (local_state.base_pos <= could_start_from) {
                            //std::cout << local_state.base_pos << ", " << group.base_point(local_state.base_pos) << std::endl;
                            local_state.move_to_child(&R, g, group.base_point(local_state.base_pos));
                        }
                        local_state.save_reduced_state(start_from);
                    }

                    int base_pos = local_state.base_pos;

                    bool base_aligned = true;
                    bool uniform = true;

                    while (g->v_size != local_state.c->cells) {
                        const int col = (*selector)(local_state.c, base_pos);
                        const int col_sz = local_state.c->ptn[col] + 1;
                        const int rand = ((int) generator()) % col_sz;
                        int v = local_state.c->lab[col + rand];

                        if(group.finished_up_to_level() + 1 == base_pos && group.is_in_base_orbit(base_pos, v)) {
                            heuristic_reroll.clear();
                            for(int i = 0; i < col_sz; ++i) {
                                heuristic_reroll.push_back(local_state.c->lab[col + i]);
                            }
                            group.reduce_to_unfinished(w, heuristic_reroll, base_pos);
                            if(heuristic_reroll.size() > 0) {
                                const int rand = ((int) generator()) % heuristic_reroll.size();
                                v = heuristic_reroll[rand];
                                //std::cout << group.is_finished(base_pos) << "re-roll!"
                                //          << group.is_in_base_orbit(base_pos, v) << std::endl;
                            }
                        }

                        if(group.is_in_base_orbit(base_pos, v) && base_aligned) {
                            v = group.base_point(local_state.base_pos);
                            assert(local_state.c->vertex_to_col[v] == col);
                            uniform = false;
                        } else {
                            base_aligned = false;
                        }

                        local_state.move_to_child(&R, g, v);
                        ++base_pos;
                    }

                    auto other_leaf = leaf_storage.lookup_leaf(local_state.T->get_hash());
                    if(other_leaf == nullptr) {
                        std::cout << "adding leaf " << local_state.T->get_hash() << std::endl;
                        leaf_storage.add_leaf(local_state.T->get_hash(), *local_state.c, local_state.base_vertex);
                    } else {
                        //std::cout << "reading leaf " << local_state.T->get_hash() << std::endl;
                        automorphism.write_color_diff(local_state.c->vertex_to_col, other_leaf->get_coloring()->lab);
                        const bool cert = R.certify_automorphism_sparse(g, automorphism.perm(), automorphism.nsupport(), automorphism.support());
                        if(cert) {
                            //std::cout << "found automorphism, hash " << local_state.T->get_hash() << " support " << automorphism.nsupport() << std::endl;
                            const bool sift = group.sift(w, g, &R, automorphism);
                            if(uniform) record_sift_result(sift);
                        }

                        automorphism.reset();
                    }
                }
            }
        };


        /**
         * \brief Depth-first search without backtracking.
         *
         * Depth-first IR search which does not backtrack. Can parallelize along the base.
         * Due to the search not back-tracking, this module can not deal with difficult parts of combinatorial graphs.
         */
        class dfs_ir {
            int fail_cnt = 0;
            int threads = 1;
            refinement *R = nullptr;
            splitmap *SM = nullptr;
            trace compare_T;

        public:
            long double grp_sz_man = 1.0; /**< group size mantissa */
            int         grp_sz_exp = 0;   /**< group size exponent */

            int         cost_snapshot = 0;

            /**
             * Setup the DFS module.
             *
             * @param threads number of threads we are allowed to dispatch
             * @param R refinement workspace
             */
            void setup(int threads, refinement *R) {
                this->R = R;
                this->threads = threads;
            }

            std::pair<bool, bool> recurse_to_equal_leaf(sgraph *g, int *initial_colors, ir_controller *state,
                                                        sparse_automorphism_workspace &automorphism) {
                bool prev_fail = false;
                int prev_fail_pos = -1;
                int cert_pos = 0;

                while ((size_t) state->base_pos < state->compare_base_color.size()) {
                    const int col = state->compare_base_color[state->base_pos];
                    const int col_sz = state->c->ptn[col] + 1;
                    if (col_sz < 2)
                        return {false, false};
                    const int ind_v = state->c->lab[col];

                    state->move_to_child(R, g, ind_v);
                    automorphism.write_singleton(&state->compare_singletons, &state->singletons,
                                                 state->base_singleton_pt[state->base_singleton_pt.size() - 1],
                                                 state->singletons.size());
                    //bool found_auto = R->certify_automorphism_sparse(g, initial_colors->get_array(), automorphism->get_array(),
                    //                                                 automorphism_supp->cur_pos, automorphism_supp->get_array());

                    bool prev_cert = true;

                    assert(state->h_hint_color_is_singleton_now ? state->h_last_refinement_singleton_only : true);

                    if (prev_fail && state->h_last_refinement_singleton_only && state->h_hint_color_is_singleton_now) {
                        prev_cert = R->check_single_failure(g, initial_colors, automorphism.perm(), prev_fail_pos);
                    }

                    //if(state->c->cells == g->v_size) {
                    if (prev_cert && state->h_last_refinement_singleton_only && state->h_hint_color_is_singleton_now) {
                        // TODO: add better heuristic to not always do this check, too expensive!
                        auto cert_res = R->certify_automorphism_sparse_report_fail_resume(g, initial_colors,
                                                                                          automorphism.perm(),
                                                                                          automorphism.nsupport(),
                                                                                          automorphism.support(),
                                                                                          cert_pos);
                        cert_pos = std::get<2>(cert_res);
                        if (std::get<0>(cert_res)) {
                            cert_res = R->certify_automorphism_sparse_report_fail_resume(g, initial_colors,
                                                                                         automorphism.perm(),
                                                                                         automorphism.nsupport(),
                                                                                         automorphism.support(),
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
                if (state->c->cells == g->v_size) {
                    //std::cout << "fail" << std::endl;
                    return {true, false};
                } else {
                    return {false, false};
                }
            }

            // adds a given automorphism to the tiny_orbit structure
            void add_automorphism_to_orbit(tiny_orbit *orbit, int *automorphism, int nsupp, int *supp) {
                for (int i = 0; i < nsupp; ++i) {
                    orbit->combine_orbits(automorphism[supp[i]], supp[i]);
                }
            }

            void make_leaf_snapshot(ir_controller *state) {
                cost_snapshot = state->T->get_position();

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
            int do_dfs(sgraph *g, int *initial_colors, ir_controller &local_state) {
                // TODO should depend on given selector
                int gens = 0;

                tiny_orbit orbs;
                orbs.initialize(g->v_size);

                //work_list automorphism(g->v_size);
                //work_list automorphism_supp(g->v_size);
                sparse_automorphism_workspace pautomorphism(g->v_size);

                mark_set color_class(g->v_size);

                //ir_controller local_state(c, &T);

                // start DFS from a leaf!
                //move_to_leaf(g, &local_state);
                // TODO move this outside DFS, only give state in leaf node (we dont even need a selector here!)

                // make a snapshot of the leaf to compare to!
                make_leaf_snapshot(&local_state);

                coloring leaf_color;
                leaf_color.copy(local_state.get_coloring());

                std::vector<int> individualize;

                double recent_cost_snapshot = 0;

                bool fail = false;

                // loop that serves to optimize Tinhofer graphs
                while (recent_cost_snapshot < 0.25 && local_state.base_pos > 0 && !fail) {
                    local_state.move_to_parent();
                    const int col = local_state.base_color[local_state.base_pos]; // TODO: detect stack of "same color"?
                    const int col_sz = local_state.c->ptn[col] + 1;
                    const int vert = local_state.base_vertex[local_state.base_pos];

                    int count_leaf = 0;
                    int count_orb = 0;

                    for (int i = col_sz - 1; i >= 0; --i) {
                        int cost_start = local_state.T->get_position();

                        const int ind_v = leaf_color.lab[col + i];
                        if (ind_v == vert || !orbs.represents_orbit(ind_v))
                            continue;
                        if (orbs.are_in_same_orbit(ind_v, vert)) { // TODO somehow skip to ones not in same orbit?
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
                        pautomorphism.write_singleton(&local_state.compare_singletons, &local_state.singletons,
                                                      local_state.base_singleton_pt[
                                                              local_state.base_singleton_pt.size() - 1],
                                                      local_state.singletons.size());
                        found_auto = R->certify_automorphism_sparse(g, initial_colors, pautomorphism.perm(),
                                                                    pautomorphism.nsupport(),
                                                                    pautomorphism.support());
                        assert(pautomorphism.perm()[vert] == ind_v);
                        //}
                        //std::cout << found_auto << ", " << automorphism_supp.cur_pos << std::endl;
                        //reset_automorphism(automorphism.get_array(), &automorphism_supp);

                        // try proper recursion
                        if (!found_auto) {
                            // TODO: heuristic that once this becomes common, we start to copy away the state and recover it, and don't track
                            // TODO: touched stuff
                            auto rec_succeeded = recurse_to_equal_leaf(g, initial_colors, &local_state, pautomorphism);
                            found_auto = (rec_succeeded.first && rec_succeeded.second);
                            if (rec_succeeded.first && !rec_succeeded.second) {
                                pautomorphism.reset();
                                pautomorphism.write_color_diff(local_state.c->vertex_to_col, leaf_color.lab);
                                found_auto = R->certify_automorphism_sparse(g, initial_colors,
                                                                            pautomorphism.perm(),
                                                                            pautomorphism.nsupport(),
                                                                            pautomorphism.support());
                            }
                        }

                        int cost_end = local_state.T->get_position();
                        double cost_partial = (cost_end - cost_start) / (cost_snapshot*1.0);
                        recent_cost_snapshot = (cost_partial + recent_cost_snapshot * 3) / 4;

                        if (found_auto) {
                            assert(pautomorphism.perm()[vert] == ind_v);
                            ++gens;
                            add_automorphism_to_orbit(&orbs, pautomorphism.perm(),
                                                      pautomorphism.nsupport(), pautomorphism.support());
                        }
                        pautomorphism.reset();

                        while (prev_base_pos < local_state.base_pos) {
                            local_state.move_to_parent();
                        }

                        if (!found_auto) {
                            std::cout << ind_v << "(F) " << std::endl;
                            fail = true;
                            break;
                        }

                        if (orbs.orbit_size(vert) == col_sz) {
                            break;
                        }
                    }

                    if (!fail) {
                        grp_sz_man *= col_sz;
                        while (grp_sz_man > 10) {
                            grp_sz_man /= 10;
                            grp_sz_exp += 1;
                        }
                    }
                }

                return local_state.base_pos;
            }
        };
    }

    /**
     * \brief The dejavu solver.
     *
     * Contains the high-level strategy of the dejavu solver, controlling the interactions between different modules
     * of the solver.
     */
    class dejavu2 {
    private:
        // high-level modules of the algorithm
        sassy::preprocessor m_prep;        /**< preprocessor */
        search_strategy::dfs_ir    m_dfs;  /**< depth-first search */
        search_strategy::bfs_ir    m_bfs;  /**< breadth-first search */
        search_strategy::random_ir m_rand; /**< randomized search */

        // utility tools used by other modules
        refinement       m_refinement; /**< workspace for color refinement and other utilities */
        selector_factory m_selectors;  /**< cell selector creation */
        group_structure::schreier m_schreier; /**< Schreier-Sims algorithm */
        ir_tree   m_tree;              /**< IR-tree */

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

        void add_to_group_size(long double add_grp_sz_man, int add_grp_sz_exp) {
            while (add_grp_sz_man > 10) {
                grp_sz_exp += 1;
                add_grp_sz_man = add_grp_sz_man / 10;
            }
            grp_sz_exp += add_grp_sz_exp;
            grp_sz_man *= add_grp_sz_man;
            while (add_grp_sz_man > 10) {
                grp_sz_exp += 1;
                add_grp_sz_man = add_grp_sz_man / 10;
            }
        }

    public:
        long double grp_sz_man = 1.0; /**< group size mantissa, see also \a grp_sz_exp */
        int         grp_sz_exp = 0;   /**< group size exponent, see also \a grp_sz_man  */

        /**
         * Compute the automorphisms of the graph \p g colored with vertex colors \p colmap. Automorphisms are returned
         * using the function pointer \p hook.
         *
         * @param g The graph.
         * @param colmap The vertex coloring of \p g. A null pointer is admissible as the trivial coloring.
         * @param hook The hook used for returning automorphisms. A null pointer is admissible if this is not needed.
         *
         */
        void automorphisms(sgraph* g, int* colmap = nullptr, dejavu_hook* hook = nullptr) {
            // TODO high-level strategy
            //  - full restarts, change selector, but make use of previously computed automorphisms (maybe inprocess) -- goal:
            //    stable, single-threaded performance
            //  - exploit strengths of DFS, random walks, BFS, and their synergies
            //  - try to use no-pruning DFS to fill up available leafs more efficiently

            // TODO miscellaneous SAT-inspired stuff
            //      - consider calloc when appropriate
            //      - weigh variables for earlier individualization in restarts (maybe just non-sense, have to see how well
            //        restarts work out first anyway)

            // control values
            int h_limit_leaf = 0;
            int h_limit_fail = 0;
            int h_error_prob = 10;
            int h_restarts   = 0;

            // TODO facilities for restarts

            // preprocess the graph using sassy
            std::cout << "preprocessing..." << std::endl;
            sassy::sgraph gg;
            transfer_sgraph_to_sassy_sgraph(g, &gg);
            m_prep.reduce(&gg, colmap, hook);
            add_to_group_size(m_prep.base, m_prep.exp);
            transfer_sassy_sgraph_to_sgraph(g, &gg);

            std::cout << std::endl << "solving..." << std::endl;
            progress_print_header();

            // initialize a coloring using colors of preprocessed graph
            coloring local_coloring;
            g->initialize_coloring(&local_coloring, colmap);

            // set up a local state for IR computations
            trace local_trace;
            ir_controller local_state(&local_coloring, &local_trace);

            // save root state for random and BFS search
            ir_reduced_save root_save;
            local_state.save_reduced_state(root_save);

            // find a selector, moves local_state to a leaf of IR tree
            m_selectors.find_base(&m_refinement, g, &local_state);
            std::function<type_selector_hook>* current_selector = m_selectors.get_selector_hook();
            progress_print("selector", std::to_string(local_state.base_pos), std::to_string(local_trace.get_position()));
            int base_size = local_state.base_pos;

            std::vector<int> base       = local_state.base_vertex;
            std::vector<int> base_sizes = local_state.base_color_sz;

            // TODO fix this, should work...
            //local_state.T->update_blueprint_hash();
            //std::cout << local_state.T->get_hash() << std::endl;

            // depth-first search starting from the computed leaf in local_state
            m_dfs.setup(0, &m_refinement);
            //m_dfs.make_leaf_snapshot(&local_state);

           const int dfs_reached_level = m_dfs.do_dfs(g, colmap, local_state);
            progress_print("dfs", std::to_string(base_size) + "-" + std::to_string(dfs_reached_level), std::to_string(m_dfs.grp_sz_man)+"10^"+std::to_string(m_dfs.grp_sz_exp));
            if(dfs_reached_level == 0) {
                add_to_group_size(m_dfs.grp_sz_man, m_dfs.grp_sz_exp);
                //return;
            }

            // set up schreier structure
            m_schreier.setup(g->v_size, base, base_sizes, dfs_reached_level);

            // intertwined random automorphisms and breadth-first search

            // random automorphisms
            m_rand.setup(h_error_prob, h_limit_leaf);
            m_rand.specific_walk(m_refinement, base, g, m_schreier, local_state, root_save);
            m_rand.random_walks(m_refinement, current_selector, g, m_schreier, local_state, root_save);
            progress_print("urandom", std::to_string(m_rand.stat_leaves()), "_");
            progress_print("schreier", "s"+std::to_string(m_schreier.stat_sparsegen())+"/d"+std::to_string(m_schreier.stat_densegen()), "_");

            // breadth-first search

            std::cout << "#symmetries: " << grp_sz_man << "*10^" << grp_sz_exp << std::endl;
        }
    };
}

#endif //DEJAVU_DFS_H
