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
     * \brief Workspace to sparse automorphisms
     *
     * Enables O(1) lookup on a sparse automorphism by using an O(n) workspace.
     */
    class sparse_automorphism_workspace {
        work_list automorphism;
        work_list automorphism_supp;
        int domain_size;
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

    class stored_automorphism {
        enum stored_automorphism_type {STORE_DENSE, STORE_SPARSE, STORE_BOTH};

        std::vector<int> sparse_cycles;
        std::vector<int> dense_map;

        stored_automorphism_type store_type;

    };

    class stored_transversal {
        enum stored_transversal_type {STORE_DENSE, STORE_SPARSE, STORE_BOTH};

        std::mutex lock_transversal;          /**< locks this transversal */

        int fixed;

        std::vector<int> fixed_orbit;         /**< contains vertices of orbit at this schreier level */
        std::vector<int> fixed_orbit_to_perm; /**< maps fixed_orbit[i] to generators[i] in class \ref schreier. */
        std::vector<int> fixed_orbit_to_pwr;  /**< power that needs to be applied to generators[i] in class \ref schreier. */

        stored_transversal_type store_type = STORE_SPARSE; /**< whether above structures are stored dense or sparse */

        std::pair<int, int>& get_generator_pwr(int v) {

        }

    public:
        /**
         * @param automorphism
         * @return whether transversal was extended
         */
        bool extend_with_automorphism(sparse_automorphism_workspace& automorphism) {

        }

        /**
         * @param automorphism
         * @return whether transversal was extended
         */
        bool fix_automorphism(sparse_automorphism_workspace& automorphism) {

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

        std::mutex lock_generators;          /**< locks the generators */

        std::vector<stored_automorphism> generators; /** list of generators */
        std::vector<stored_transversal>  transversals;

    public:
        /**
         * Sift automorphism into the schreier structure.
         *
         * @param automorphism Automorphism to be sifted. Will be manipulated by the method.
         * @return Whether automorphism was added to the schreier structure or not.
         */
        bool sift(sparse_automorphism_workspace& automorphism) {
            bool changed = false;

            // TODO:
            for(int level = 0; level < transversals.size(); ++level) {
                if(transversals[level].extend_with_automorphism(automorphism)) {
                    changed = true;
                    // TODO: copy add to generators
                    // TODO: maybe re-arrange, dependencies are weird
                }
                const bool is_identity = transversals[level].fix_automorphism(automorphism);
                if(is_identity)
                    break;
            }
            return changed;
        }
    };



    /**
     * \brief Breadth-first search.
     */
    class bfs_ir {
        // TODO implement new, simpler bfs
        // TODO depends on ir_tree, and selector
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
        // TODO leaf storage, maybe do extra class that can be shared across threads

    public:
        void setup(int error, int leaf_store_limit) {

        }

        // TODO implement dejavu strategy, more simple
        // TODO depends on ir_tree, selector, and given base (no need to sift beyond base!)
        // TODO: swap out ir_reduced to weighted IR tree later? or just don't use automorphism pruning on BFS...?
        void random_walks(refinement& R, std::function<type_selector_hook>* selector, sgraph* g, schreier& group,
                          ir_controller& local_state, ir_reduced_save& start_from) {
            std::cout << "random walks" << std::endl;

            for(int i = 0; i < 50; ++i ) {
                local_state.load_reduced_state(start_from);
                int base_pos = local_state.base_pos;
                while (g->v_size != local_state.c->cells) {
                    const int col = (*selector)(local_state.c, base_pos);
                    const int col_sz = local_state.c->ptn[col] + 1;
                    const int rand = ((int) generator()) % col_sz;
                    const int v = local_state.c->lab[col + rand];
                    local_state.move_to_child(&R, g, v);
                    ++base_pos;
                }

                std::cout << "leaf " << local_state.T->get_hash() << std::endl;
            }
            // TODO continue...
        }
    };


    /**
     * \brief Depth-first search without backtracking.
     *
     * Depth-first IR search which does not backtrack. Can parallelize along the base.
     * Due to the search not back-tracking, this module can not deal with difficult parts of combinatorial graphs.
     */
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

        std::pair<bool, bool> recurse_to_equal_leaf(sgraph* g, int* initial_colors, ir_controller* state, sparse_automorphism_workspace& automorphism) {
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
                automorphism.write_singleton(&state->compare_singletons, &state->singletons,
                                             state->base_singleton_pt[state->base_singleton_pt.size() - 1],
                                             state->singletons.size());
                //bool found_auto = R->certify_automorphism_sparse(g, initial_colors->get_array(), automorphism->get_array(),
                //                                                 automorphism_supp->cur_pos, automorphism_supp->get_array());

                bool prev_cert = true;

                assert(state->h_hint_color_is_singleton_now?state->h_last_refinement_singleton_only:true);

                if(prev_fail && state->h_last_refinement_singleton_only && state->h_hint_color_is_singleton_now) {
                    prev_cert = R->check_single_failure(g, initial_colors, automorphism.perm(), prev_fail_pos);
                }

                //if(state->c->cells == g->v_size) {
                if(prev_cert && state->h_last_refinement_singleton_only && state->h_hint_color_is_singleton_now) {
                    // TODO: add better heuristic to not always do this check, too expensive!
                    auto cert_res = R->certify_automorphism_sparse_report_fail_resume(g, initial_colors,
                                                                          automorphism.perm(),
                                                                          automorphism.nsupport(),
                                                                          automorphism.support(), cert_pos);
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
            if(state->c->cells == g->v_size) {
                //std::cout << "fail" << std::endl;
                return {true, false};
            } else {
                return {false, false};
            }
        }

        // adds a given automorphism to the tiny_orbit structure
        void  add_automorphism_to_orbit(tiny_orbit* orbit, int* automorphism, int nsupp, int* supp) {
            for(int i = 0; i < nsupp; ++i) {
                orbit->combine_orbits(automorphism[supp[i]], supp[i]);
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
                    if(!found_auto) {
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

                    if (found_auto) {
                        assert(pautomorphism.perm()[vert] == ind_v);
                        ++gens;
                        add_automorphism_to_orbit(&orbs, pautomorphism.perm(),
                                                  pautomorphism.nsupport(), pautomorphism.support());
                    }
                    pautomorphism.reset();

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

    /**
     * \brief The dejavu solver.
     *
     * Contains the high-level strategy of the dejavu solver, controlling the interactions between the different modules
     * of the solver.
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
            // TODO high-level strategy
            //  - full restarts, change selector, but make use of previously computed automorphisms (maybe inprocess) -- goal:
            //    stable, single-threaded performance
            //  - exploit strengths of DFS, random walks, BFS, and their synergies
            //  - try to use no-pruning DFS to fill up available leafs more efficiently

            // control values
            int h_limit_leaf = 0;
            int h_limit_fail = 0;
            int h_error_prob = 10;
            int h_restarts   = 0;

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

            // save root state for random and BFS search
            ir_reduced_save root_save;
            local_state.save_reduced_state(root_save);

            // find a selector, moves local_state to a leaf of IR tree
            m_selectors.find_base(&m_refinement, g, &local_state);
            std::function<type_selector_hook>* current_selector = m_selectors.get_selector_hook();
            std::cout << "trace length: " << local_trace.get_position() << std::endl;
            std::cout << "leaf " << local_state.T->get_hash() << std::endl;

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
                //return;
            }

            // intertwined random automorphisms and breadth-first search

            // random automorphisms
            m_rand.setup(h_error_prob, h_limit_leaf);
            m_rand.random_walks(m_refinement, current_selector, g, m_schreier, local_state, root_save);

            // breadth-first search

        }
    };
}

#endif //DEJAVU_DFS_H
