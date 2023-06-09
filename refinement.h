// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_REFINEMENT_H
#define DEJAVU_REFINEMENT_H

#include "ds.h"
#include "coloring.h"
#include "sgraph.h"
#include "utility.h"

namespace dejavu {

    using namespace ds;

    namespace ir {

        // return whether to continue color refinement
        // bool split_color_hook(int color_initial, int new_color, int new_color_sz);
        typedef bool type_split_color_hook(const int, const int, const int);

        // return whether to continue splitting the respective cell, or skip it
        // bool worklist_color_hook(int color, int color_sz);
        typedef bool type_worklist_color_hook(const int, const int);

        // worklist implementation for color refinement
        class cell_worklist {
        public:
            void initialize(int domain_size) {
                arr = new int[domain_size];
                arr_sz = domain_size;
                init = true;
                cur_pos = 0;
            }

            ~cell_worklist() {
                delete[] arr;
            }

            int add_cell(work_set_int *queue_pointer, int col) {
                assert(init);
                assert(cur_pos >= 0 && cur_pos < arr_sz - 1);
                queue_pointer->set(col, cur_pos);
                arr[cur_pos] = col;
                cur_pos++;
                return 0;
            }

            int next_cell(work_set_int *queue_pointer, coloring *c) {
                // look at first 12 positions and pick the (first) smallest cell within these entries
                int sm_j = cur_pos - 1;
                for (int j = cur_pos - 1; j >= 0 && ((cur_pos - j) <= 12); --j) {
                    if (c->ptn[arr[j]] < c->ptn[arr[sm_j]]) {
                        sm_j = j;
                        if (c->ptn[arr[sm_j]] == 0)
                            break;
                    }
                }

                // swap sm_j and j
                const int sm_col = arr[sm_j];
                arr[sm_j] = arr[cur_pos - 1];
                queue_pointer->set(arr[sm_j], sm_j);

                cur_pos--;
                queue_pointer->set(sm_col, -1);
                return sm_col;
            }

            int next_cell(work_set_int *queue_pointer, coloring *c, work_list_t<int> *singleton_hint) {
                // use singleton_hint
                int sm_j = -1;
                while (!singleton_hint->empty()) {
                    const int next_hint = singleton_hint->pop_back();
                    sm_j = queue_pointer->get(next_hint);
                    if (sm_j == -1)
                        continue;
                    else
                        break;
                }

                // look at first 12 positions and pick the (first) smallest cell within these entries
                if (sm_j == -1) {
                    sm_j = cur_pos - 1;
                    for (int j = cur_pos - 1; j >= 0 && ((cur_pos - j) <= 12); --j) {
                        //bool smaller_d = g->d[arr[j]] > g->d[arr[sm_j]];
                        const int size_sm_j = c->ptn[arr[sm_j]];
                        if (c->ptn[size_sm_j] == 0)
                            break;
                        const bool smaller = (c->ptn[arr[j]] < size_sm_j);
                        //bool eq        = (c->ptn[arr[j]] == c->ptn[arr[sm_j]]);
                        sm_j = smaller ? j : sm_j;
                    }
                }

                // swap sm_j and j
                const int sm_col = arr[sm_j];
                arr[sm_j] = arr[cur_pos - 1];
                queue_pointer->set(arr[sm_j], sm_j);

                cur_pos--;
                queue_pointer->set(sm_col, -1);
                return sm_col;
            }

            void replace_cell(work_set_int *queue_pointer, int col_old, int col) {
                const int pos = queue_pointer->get(col_old);
                arr[pos] = col;
                assert(queue_pointer->get(col_old) != -1);
                queue_pointer->set(col_old, -1);
                queue_pointer->set(col, pos);
            }

            void reset(work_set_int *queue_pointer) {
                while (cur_pos > 0) {
                    cur_pos--;
                    queue_pointer->set(arr[cur_pos], -1);
                }
            }

            [[nodiscard]] bool empty() const {
                return (cur_pos == 0);
            }

            [[nodiscard]] int size() const {
                return cur_pos;
            }

        private:
            int *arr = nullptr;
            int arr_sz = -1;
            int cur_pos = -1;
            bool init = false;
        };

/**
 * \brief Color refinement and related algorithms
 *
 * Class that is used to preserve a workspace for color refinement, automorphism certification, and other miscellaneous
 * tasks. Once initialized to a certain size, the workspace can not be enlargened and methods can only be used for
 * graphs of the initial size, or smaller.
*/
        class refinement {
        public:
            /**
         * The color refinement algorithm. Refines a given coloring with respect to a given graph.
         * @param g The graph.
         * @param c The coloring to be refined.
         * @param init_color Initialize the worklist with a single color class (e.g., after individualization). The
         * default value -1 denotes that the worklist is initialized with all color classes of the coloring.
         * @param color_limit Integer which is used to stop refinement whenever the refined coloring reaches this number
         * of color classes. The default value -1 denotes that refinement is performed exhaustively.
         * @param split_hook Function pointer that is called whenever a color class is split. Return value can be used
             * to stop refinement early.
         * @param worklist_hook Function pointer that is called whenever a color class is considered for refinement. Return value
         * can be used to skip refinement of that color class.
         */
            void refine_coloring(sgraph *g, coloring *c, int init_color = -1, int color_limit = -1,
                                 const std::function<type_split_color_hook>    &split_hook = nullptr,
                                 const std::function<type_worklist_color_hook> &worklist_hook = nullptr) {
                singleton_hint.reset();
                assure_initialized(g);

                cell_todo.reset(&queue_pointer);

                if (init_color < 0) {
                    // initialize queue with all classes (except for largest one)
                    for (int i = 0; i < c->domain_size;) {
                        cell_todo.add_cell(&queue_pointer, i);
                        const int col_sz = c->ptn[i];
                        if (col_sz == 0) {
                            //singleton_hint.push_back(i);
                        }
                        i += col_sz + 1;

                    }
                } else {
                    const int col_sz = c->ptn[init_color];
                    assert(c->vertex_to_col[c->lab[init_color]] == init_color);
                    cell_todo.add_cell(&queue_pointer, init_color);
                    if (col_sz == 0) {
                        //singleton_hint.push_back(init_color);
                    }
                }

                while (!cell_todo.empty()) {
                    color_class_splits.reset();
                    const int next_color_class = cell_todo.next_cell(&queue_pointer, c);
                    const int next_color_class_sz = c->ptn[next_color_class] + 1;

                    if (worklist_hook && !worklist_hook(next_color_class, next_color_class_sz)) {
                        continue;
                    }

                    const bool very_dense = (g->d[c->lab[next_color_class]] > (g->v_size / (next_color_class_sz + 1)));

                    // this scheme is reverse-engineered from the color refinement in Traces by Adolfo Piperno
                    // we choose a separate algorithm depending on the size and density of the graph and/or color class
                    if (next_color_class_sz == 1 && !(g->dense && very_dense)) {
                        // singleton
                        refine_color_class_singleton(g, c, next_color_class,&color_class_splits);
                    } else if (g->dense) {
                        if (very_dense) { // dense-dense
                            refine_color_class_dense_dense(g, c, next_color_class,next_color_class_sz,
                                                           &color_class_splits);
                        } else { // dense-sparse
                            refine_color_class_dense(g, c, next_color_class, next_color_class_sz,
                                                     &color_class_splits);
                        }
                    } else { // sparse
                        refine_color_class_sparse(g, c, next_color_class, next_color_class_sz,
                                                  &color_class_splits);
                    }

                    // add all new classes except for the first, largest one
                    int latest_old_class = -1;
                    bool skipped_largest = false;
                    bool early_out = false;

                    // color class splits are sorted in reverse
                    // the old color class will always come last
                    while (!color_class_splits.empty()) {
                        int old_class = color_class_splits.last()->first.first;
                        int new_class = color_class_splits.last()->first.second;
                        bool is_largest = color_class_splits.last()->second;

                        c->cells += (old_class != new_class);
                        int class_size = c->ptn[new_class];

                        if (split_hook && !split_hook(old_class, new_class, c->ptn[new_class] + 1)) {
                            early_out = true;
                        }

                        if (latest_old_class != old_class) {
                            latest_old_class = old_class;
                            skipped_largest = false;
                        }

                        // management code for skipping largest class resulting from old_class
                        color_class_splits.pop_back();
                        //int new_class_sz = c->ptn[new_class] + 1;

                        if (skipped_largest || !is_largest) {
                            cell_todo.add_cell(&queue_pointer, new_class);
                            //if(new_class_sz == 1)
                            //singleton_hint.push_back(new_class);
                        } else {
                            skipped_largest = true;

                            // since old color class will always appear last, the queue pointer of old color class is still valid!
                            int i = queue_pointer.get(old_class);
                            if (i >= 0) {
                                cell_todo.replace_cell(&queue_pointer, old_class, new_class);
                                //if(new_class_sz == 1)
                                //    singleton_hint.push_back(new_class);
                            }
                        }
                    }

                    if (early_out) {
                        color_class_splits.reset();
                        cell_todo.reset(&queue_pointer);
                        break;
                    }

                    // detection if coloring is discrete
                    if (c->cells == g->v_size) {
                        color_class_splits.reset();
                        cell_todo.reset(&queue_pointer);
                        break;
                    }

                    // partition is at least as large as the one of target invariant, can skip to the end of the entire refinement
                    if (c->cells == color_limit) {
                        color_class_splits.reset();
                        cell_todo.reset(&queue_pointer);
                        break;
                    }
                }
            }

            /**
         * Individualizes a vertex in a coloring.
         * @param c Coloring in which the vertex is individualized.
         * @param v Vertex to be individualized.
         * @param split_hook Function to be called whenever a color class is split. Return value is not used.
         * @return
         */
            static int
            individualize_vertex(coloring *c, int v, const std::function<type_split_color_hook> &split_hook = nullptr) {
                const int color = c->vertex_to_col[v];
                const int pos = c->vertex_to_lab[v];

                int color_class_size = c->ptn[color];

                assert(color_class_size > 0);

                const int vertex_at_pos = c->lab[color + color_class_size];
                c->lab[pos] = vertex_at_pos;
                c->vertex_to_lab[vertex_at_pos] = pos;

                c->lab[color + color_class_size] = v;
                c->vertex_to_lab[v] = color + color_class_size;
                c->vertex_to_col[v] = color + color_class_size;

                c->ptn[color] -= 1;
                c->ptn[color + color_class_size] = 0;
                c->ptn[color + color_class_size - 1] = 0;
                c->cells += 1;

                if (split_hook) {
                    split_hook(color, color + color_class_size, 1);
                    split_hook(color, color, c->ptn[color] + 1);
                }

                return color + color_class_size;
            }

            // color refinement that does not produce an isomorphism-invariant partitioning, but uses more optimization
            // techniques -- meant to be used as the first refinement in automorphism computation
            bool refine_coloring_first(sgraph *g, coloring *c, int init_color_class = -1) {
                assure_initialized(g);
                singleton_hint.reset();

                cell_todo.reset(&queue_pointer);

                if (init_color_class < 0) {
                    for (int i = 0; i < c->domain_size;) {
                        cell_todo.add_cell(&queue_pointer, i);
                        const int col_sz = c->ptn[i];
                        if (col_sz == 0) {
                            singleton_hint.push_back(i);
                        }
                        i += col_sz + 1;
                    }
                } else {
                    cell_todo.add_cell(&queue_pointer, init_color_class);
                }

                while (!cell_todo.empty()) {
                    color_class_splits.reset();
                    const int next_color_class = cell_todo.next_cell(&queue_pointer, c, &singleton_hint);
                    const int next_color_class_sz = c->ptn[next_color_class] + 1;
                    //bool dense_dense = false;
                    //if(g->d[c->lab[next_color_class]] == 0) {
                    //   continue;
                    //}
                    //if (g->d[c->lab[next_color_class]] > 5) {
                    //    dense_dense = (g->d[c->lab[next_color_class]] > (g->v_size / (next_color_class_sz)));
                    //}
                    const bool very_dense = (g->d[c->lab[next_color_class]] > (g->v_size / (next_color_class_sz + 1)));

                    if (next_color_class_sz == 1 && !(g->dense && very_dense)) {
                        // singleton
                        refine_color_class_singleton_first(g, c, next_color_class,&color_class_splits);
                    } else if (g->dense) {
                        if (very_dense) { // dense-dense
                            refine_color_class_dense_dense_first(g, c, next_color_class, next_color_class_sz,
                                                                 &color_class_splits);
                        } else { // dense-sparse
                            refine_color_class_dense_first(g, c, next_color_class, next_color_class_sz,
                                                           &color_class_splits);
                        }
                    } else { // sparse
                        refine_color_class_sparse_first(g, c, next_color_class, next_color_class_sz,
                                                        &color_class_splits);
                    }

                    // add all new classes except for the first, largest one
                    int latest_old_class = -1;
                    bool skipped_largest = false;

                    // color class splits are sorted in reverse
                    // the old color class will always come last
                    while (!color_class_splits.empty()) {
                        const int old_class = color_class_splits.last()->first.first;
                        const int new_class = color_class_splits.last()->first.second;
                        const bool is_largest = color_class_splits.last()->second;

                        c->cells += (old_class != new_class);
                        const int class_size = c->ptn[new_class] + 1;

#ifndef NDEBUG // debug code
                        if(color_class_splits.empty()) {
int actual_cells = 0;
for (int i = 0; i < c->domain_size;) {
actual_cells += 1;
i += c->ptn[i] + 1;
}

assert(c->cells == actual_cells);
}
#endif

                        if (c->cells == g->v_size) {
                            color_class_splits.reset();
                            cell_todo.reset(&queue_pointer);
                            return true;
                        }

                        if (latest_old_class != old_class) {
                            latest_old_class = old_class;
                            skipped_largest = false;
                        }

                        color_class_splits.pop_back();
                        const int new_class_sz = c->ptn[new_class] + 1;


                        if (skipped_largest || !is_largest) {
                            cell_todo.add_cell(&queue_pointer, new_class);
                            if (new_class_sz == 1)
                                singleton_hint.push_back(new_class);
                        } else {
                            skipped_largest = true;

                            // since old color class will always appear last, the queue pointer of old color class is still valid!
                            int i = queue_pointer.get(old_class);
                            if (i >= 0) {
                                cell_todo.replace_cell(&queue_pointer, old_class, new_class);
                                if (new_class_sz == 1)
                                    singleton_hint.push_back(new_class);
                            }
                        }
                    }
                }

                // assert(assert_is_equitable(g, c));
                return true;
            }

            // certify an automorphism on a graph
            bool certify_automorphism(sgraph *g, const int *p) {
                int i, found;

                assure_initialized(g);

                for (i = 0; i < g->v_size; ++i) {
                    const int image_i = p[i];
                    if (image_i == i)
                        continue;

                    scratch_set.reset();
                    // automorphism must preserve neighbours
                    found = 0;
                    for (int j = g->v[i]; j < g->v[i] + g->d[i]; ++j) {
                        const int vertex_j = g->e[j];
                        const int image_j = p[vertex_j];
                        scratch_set.set(image_j);
                        found += 1;
                    }
                    for (int j = g->v[image_i]; j < g->v[image_i] + g->d[image_i]; ++j) {
                        const int vertex_j = g->e[j];
                        if (!scratch_set.get(vertex_j)) {
                            return false;
                        }
                        scratch_set.unset(vertex_j);
                        found -= 1;
                    }
                    if (found != 0) {
                        return false;
                    }
                }

                return true;
            }

            // certify an automorphism on a graph
            bool certify_automorphism(sgraph *g, const int *colmap, const int *p) {
                int i, found;

                assure_initialized(g);

                for (i = 0; i < g->v_size; ++i) {
                    const int image_i = p[i];
                    if (image_i == i)
                        continue;
                    if (colmap[i] != colmap[image_i]) // colors must be equal
                        return false;

                    scratch_set.reset();
                    // automorphism must preserve neighbours
                    found = 0;
                    for (int j = g->v[i]; j < g->v[i] + g->d[i]; ++j) {
                        const int vertex_j = g->e[j];
                        const int image_j = p[vertex_j];
                        scratch_set.set(image_j);
                        found += 1;
                    }
                    for (int j = g->v[image_i]; j < g->v[image_i] + g->d[image_i]; ++j) {
                        const int vertex_j = g->e[j];
                        if (!scratch_set.get(vertex_j)) {
                            return false;
                        }
                        scratch_set.unset(vertex_j);
                        found -= 1;
                    }
                    if (found != 0) {
                        return false;
                    }
                }

                return true;
            }

            // certify an automorphism on a graph, sparse
            bool certify_automorphism_sparse(const sgraph *g, const int *p, int supp, const int *supp_arr) {
                int i, found;

                assure_initialized(g);

                for (int f = 0; f < supp; ++f) {
                    i = supp_arr[f];
                    const int image_i = p[i];
                    scratch_set.reset();
                    // automorphism must preserve neighbours
                    found = 0;
                    for (int j = g->v[i]; j < g->v[i] + g->d[i]; ++j) {
                        const int vertex_j = g->e[j];
                        const int image_j = p[vertex_j];
                        scratch_set.set(image_j);
                        found += 1;
                    }
                    for (int j = g->v[image_i]; j < g->v[image_i] + g->d[image_i]; ++j) {
                        const int vertex_j = g->e[j];
                        if (!scratch_set.get(vertex_j)) {
                            scratch_set.reset();
                            return false;
                        }
                        found -= 1;
                    }
                    if (found != 0) {
                        scratch_set.reset();
                        return false;
                    }
                }
                scratch_set.reset();
                return true;
            }

            // certify an automorphism on a graph, sparse
            bool certify_automorphism_sparse(const sgraph *g, const int *colmap, const int *p, int supp,
                                             const int *supp_arr) {
                int i, found;

                assure_initialized(g);

                for (int f = 0; f < supp; ++f) {
                    i = supp_arr[f];
                    const int image_i = p[i];
                    if (colmap[i] != colmap[image_i]) // colors must be equal
                        return false;

                    scratch_set.reset();
                    // automorphism must preserve neighbours
                    found = 0;
                    for (int j = g->v[i]; j < g->v[i] + g->d[i]; ++j) {
                        const int vertex_j = g->e[j];
                        const int image_j = p[vertex_j];
                        if (colmap[vertex_j] != colmap[image_j])
                            return false;
                        scratch_set.set(image_j);
                        found += 1;
                    }
                    for (int j = g->v[image_i]; j < g->v[image_i] + g->d[image_i]; ++j) {
                        const int vertex_j = g->e[j];
                        if (!scratch_set.get(vertex_j)) {
                            return false;
                        }
                        scratch_set.unset(vertex_j);
                        found -= 1;
                    }
                    if (found != 0) {
                        return false;
                    }
                }

                return true;
            }

            // certify an automorphism on a graph, sparse, report on which vertex failed
            std::pair<bool, int>
            certify_automorphism_sparse_report_fail(const dejavu::sgraph *g, const int *colmap,
                                                    const int *p, int supp, const int *supp_arr) {
                int i, found;

                assure_initialized(g);

                //for(i = 0; i < g->v_size; ++i) {
                for (int f = 0; f < supp; ++f) {
                    i = supp_arr[f];
                    const int image_i = p[i];
                    if (colmap[i] != colmap[image_i]) // colors must be equal
                        return {false, -1};

                    scratch_set.reset();
                    // automorphism must preserve neighbours
                    found = 0;
                    for (int j = g->v[i]; j < g->v[i] + g->d[i]; ++j) {
                        const int vertex_j = g->e[j];
                        const int image_j = p[vertex_j];
                        if (colmap[vertex_j] != colmap[image_j])
                            return {false, i};
                        scratch_set.set(image_j);
                        found += 1;
                    }
                    for (int j = g->v[image_i]; j < g->v[image_i] + g->d[image_i]; ++j) {
                        const int vertex_j = g->e[j];
                        if (!scratch_set.get(vertex_j)) {
                            return {false, i};
                        }
                        scratch_set.unset(vertex_j);
                        found -= 1;
                    }
                    if (found != 0) {
                        return {false, i};
                    }
                }

                return {true, -1};
            }

            // certify an automorphism on a graph, sparse, report on which vertex failed
            std::tuple<bool, int, int>
            certify_automorphism_sparse_report_fail_resume(const sgraph *g, const int*, const int *p, int supp,
                                                           const int *supp_arr, int pos_start) {
                assure_initialized(g);

                //for(i = 0; i < g->v_size; ++i) {
                for (int f = pos_start; f < supp; ++f) {
                    const int i = supp_arr[f];
                    const int image_i = p[i];
                    assert(image_i != i);

                    scratch_set.reset();
                    // automorphism must preserve neighbours
                    const int i_deg = g->d[i];
                    for (int j = g->v[i]; j < g->v[i] + i_deg; ++j) {
                        const int vertex_j = g->e[j];
                        const int image_j = p[vertex_j];
                        scratch_set.set(image_j);
                    }
                    const int image_i_deg = g->d[image_i];
                    int fail = 0;
                    for (int j = g->v[image_i]; j < g->v[image_i] + image_i_deg; ++j) {
                        const int vertex_j = g->e[j];
                        fail += !scratch_set.get(vertex_j);
                        scratch_set.unset(vertex_j);
                    }

                    if (fail != 0 || i_deg != image_i_deg) {
                        return {false, i, f};
                    }
                }

                return {true, -1, supp};
            }

            // certify an automorphism, for a single vertex
            bool check_single_failure(const sgraph *g, const int *, const int *p, int failure) {
                assure_initialized(g);

                const int i = failure;
                const int image_i = p[i];

                scratch_set.reset();
                // automorphism must preserve neighbours
                const int i_deg = g->d[i];
                for (int j = g->v[i]; j < g->v[i] + i_deg; ++j) {
                    const int vertex_j = g->e[j];
                    const int image_j = p[vertex_j];
                    scratch_set.set(image_j);
                }
                const int image_i_deg = g->d[image_i];
                int fail = 0;
                for (int j = g->v[image_i]; j < g->v[image_i] + image_i_deg; ++j) {
                    const int vertex_j = g->e[j];
                    fail += !scratch_set.get(vertex_j);
                    scratch_set.unset(vertex_j);
                }

                if (fail != 0 || i_deg != image_i_deg) {
                    return false;
                }

                return true;
            }

            ~refinement() {
                if (domain_size >= 0)
                    delete[] workspace_int;
            }

        private:
            int domain_size = -1;
            work_set_int queue_pointer;
            cell_worklist cell_todo;
            mark_set scratch_set;
            work_list_t<int> vertex_worklist;
            work_set_t<int> color_vertices_considered; // todo this should use a different datastructure, with n space and not 2n
            work_set_t<int> neighbours;
            work_set_t<int> neighbour_sizes;
            work_list_t<int> singleton_hint;
            work_list_t<int> old_color_classes;
            work_list_pair_bool color_class_splits;

            int *scratch = nullptr;
            int *workspace_int = nullptr;

            void assure_initialized(const sgraph *g) {
                if (g->v_size > domain_size) {
                    const int n = g->v_size;

                    // reducing contention on heap allocator through bulk allocation...
                    if(workspace_int) delete[] workspace_int;
                    workspace_int = new int[n];

                    vertex_worklist.allocate(n * 2);
                    //singletons.initialize(n);
                    singleton_hint.allocate(n);
                    old_color_classes.allocate(n);
                    neighbours.initialize(n);
                    neighbour_sizes.initialize(n);
                    queue_pointer.initialize(n);
                    color_vertices_considered.initialize(n);

                    scratch = (int *) workspace_int;
                    //scratch_set.initialize_from_array(workspace_int + n, n);
                    scratch_set.initialize(n);

                    color_class_splits.allocate(n);
                    cell_todo.initialize(n * 2);

                    memset(scratch, 0, n * sizeof(int));
                    domain_size = n;
                }
            }

            bool refine_color_class_sparse(sgraph *g, coloring *c, int color_class, int class_size,
                                           work_list_pair_bool *report_splits) {
                // for all vertices of the color class...
                bool comp, mark_as_largest;
                int i, j, cc, end_cc, largest_color_class_size, acc;
                int *vertex_to_lab = c->vertex_to_lab;
                int *lab = c->lab;
                int *ptn = c->ptn;
                int *vertex_to_col = c->vertex_to_col;

                cc = color_class; // iterate over color class
                comp = true;

                old_color_classes.reset();
                neighbours.reset();
                color_vertices_considered.reset();

                end_cc = color_class + class_size;
                while (cc < end_cc) { // increment value of neighbours of vc by 1
                    const int vc = lab[cc];
                    const int pe = g->v[vc];
                    const int end_i = pe + g->d[vc];
                    for (i = pe; i < end_i; i++) {
                        const int v = g->e[i];
                        const int col = vertex_to_col[v];
                        if (ptn[col] == 0) {
                            continue;
                        }
                        neighbours.inc_nr(v);
                        if (neighbours.get(v) == 0) {
                            color_vertices_considered.inc_nr(col);
                            assert(col + color_vertices_considered.get(col) < g->v_size);
                            scratch[col + color_vertices_considered.get(col)] = v; // hit vertices
                                                                                  // TODO could consolidate or checkerboard color_vertices_considered and scratch?
                            if (color_vertices_considered.get(col) == 0) { // TODO use color_vertices_considered[col] == 0?
                                old_color_classes.push_back(col);
                            }
                        }
                    }
                    cc += 1;
                }

                // sort split color classes
                old_color_classes.sort();

                // split color classes according to neighbour count
                while (!old_color_classes.empty()) {
                    const int _col = old_color_classes.pop_back();
                    const int _col_sz = ptn[_col] + 1;
                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    for (i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
                        const int v = scratch[_col + i];
                        const int index = neighbours.get(v) + 1;
                        if (neighbour_sizes.inc(index) == 0)
                            vertex_worklist.push_back(index);
                    }

                    vertex_worklist.sort();

                    // enrich neighbour_sizes to accumulative counting array
                    acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int k = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(k) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(k, val + acc);
                            acc += val;
                            const int _ncol = _col + _col_sz - (neighbour_sizes.get(k));
                            if (_ncol != _col)
                                ptn[_ncol] = -1;
                        }
                    }

                    const int vcount = color_vertices_considered.get(_col);

                    vertex_worklist.reset();
                    j = 0;
                    color_vertices_considered.set(_col, -1);

                    // determine colors and rearrange vertices
                    while (j < vcount + 1) {
                        const int v = scratch[_col + j];
                        ++j;
                        if ((neighbours.get(v) == -1))
                            continue;
                        const int v_new_color = _col + _col_sz - neighbour_sizes.get(neighbours.get(v) + 1);
                        neighbours.set(v, -1);
                        if (v_new_color == _col)
                            continue;

                        const int lab_pt = v_new_color + ptn[v_new_color] + 1;
                        ptn[v_new_color] += 1;
                        ptn[_col] -= (_col != v_new_color);

                        const int vertex_old_pos = vertex_to_lab[v];
                        const int vertex_at_pos = lab[lab_pt];
                        lab[vertex_old_pos] = vertex_at_pos;
                        vertex_to_lab[vertex_at_pos] = vertex_old_pos;
                        lab[lab_pt] = v;
                        vertex_to_col[v] = v_new_color;
                        vertex_to_lab[v] = lab_pt;
                    }

                    // add new colors to worklist
                    largest_color_class_size = -1;
                    for (i = _col; i < _col + _col_sz;) {
                        //assert(i >= 0 && i < ptn_sz);
                        assert(ptn[i] + 1 > 0);
                        mark_as_largest = largest_color_class_size < (ptn[i] + 1);
                        largest_color_class_size = mark_as_largest ? (ptn[i] + 1) : largest_color_class_size;

                        report_splits->push_back(std::pair<std::pair<int, int>, bool>(
                                std::pair<int, int>(_col, i), mark_as_largest));

                        i += ptn[i] + 1;
                    }
                }

                neighbour_sizes.reset();
                vertex_worklist.reset();

                return comp;
            }

            bool refine_color_class_dense(sgraph *g, coloring *c, int color_class, int class_size,
                                          work_list_pair_bool *report_splits) {
                bool comp;
                int i, cc, acc, largest_color_class_size, pos;
                cc = color_class; // iterate over color class
                comp = true;

                neighbours.reset_hard();
                scratch_set.reset();
                old_color_classes.reset();
                vertex_worklist.reset();

                const int end_cc = color_class + class_size;

                // for all vertices of the color class...
                while (cc < end_cc) { // increment value of neighbours of vc by 1
                    const int vc = c->lab[cc];
                    const int pe = g->v[vc];
                    const int end_i = pe + g->d[vc];
                    for (i = pe; i < end_i; i++) {
                        const int v = g->e[i];
                        const int col = c->vertex_to_col[v];
                        if (c->ptn[col] > 0) {
                            neighbours.inc(v); // want to use reset variant?
                            if (!scratch_set.get(col)) {
                                scratch_set.set(col);
                                old_color_classes.push_back(col);
                            }
                        } /*else {
                            if (config.CONFIG_IR_FULL_INVARIANT)
                                vertex_worklist.push_back(col);
                            else {
                                singleton_inv += MASH4(col);
                            }
                        }*/
                    }
                    cc += 1;
                }

                // write singletons
                //comp = I->write_top_and_compare(g->v_size * 3 + old_color_classes.cur_pos) && comp;

                /*if(config.CONFIG_IR_FULL_INVARIANT) {
                vertex_worklist.sort();
                while (!vertex_worklist.empty()) {
                    const int col = vertex_worklist.pop_back();
                    //comp = I->write_top_and_compare(g->v_size * 9 + col) && comp;
                }
            } else {
                comp = I->write_top_and_compare(singleton_inv) && comp;
            }*/

                old_color_classes.sort();

                // for every cell to be split...
                while (!old_color_classes.empty()) {
                    const int col = old_color_classes.pop_back();
                    const int col_sz = c->ptn[col] + 1;

                    if (col_sz == 1)
                        continue;

                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    const int first_index = neighbours.get(c->lab[col]) + 1;
                    vertex_worklist.push_back(first_index);
                    int total = 0;
                    for (i = col; i < col + col_sz; ++i) {
                        const int v = c->lab[i];
                        int index = neighbours.get(v) + 1;
                        if (index == first_index)
                            continue;
                        total += 1;
                        if (neighbour_sizes.inc(index) == 0)
                            vertex_worklist.push_back(index);
                    }

                    neighbour_sizes.inc(first_index);
                    neighbour_sizes.set(first_index, col_sz - total - 1);

                    //comp = I->write_top_and_compare(vertex_worklist.cur_pos) && comp;

                    if (vertex_worklist.cur_pos == 1) {
                        // no split
                        //const int v_degree = neighbours.get(c->lab[col]);
                        //comp = I->write_top_and_compare(col + g->v_size * v_degree) && comp;
                        continue;
                    }

                    vertex_worklist.sort();
                    // enrich neighbour_sizes to accumulative counting array
                    acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int j = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(j) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(j, val + acc);
                            acc += val;
                            const int _col = col + col_sz - (neighbour_sizes.get(j));
                            c->ptn[_col] = -1; // this is val - 1, actually...
                        }
                    }

                    // copy cell for rearranging
                    memcpy(scratch, c->lab + col, col_sz * sizeof(int));
                    //vertex_worklist.cur_pos = col_sz;
                    pos = col_sz;

                    // determine colors and rearrange
                    // split color classes according to count in counting array
                    while (pos > 0) {
                        const int v = scratch[--pos];
                        const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
                        // i could immediately determine size of v_new_color here and write invariant before rearrange?
                        assert(v_new_color >= col);
                        assert(v_new_color < col + col_sz);

                        c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                        c->vertex_to_col[v] = v_new_color;
                        c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                        c->ptn[v_new_color] += 1;
                    }

                    largest_color_class_size = -1;

                    // determine largest class to throw away and finish (fourth iteration)
                    for (i = col; i < col + col_sz;) {
                        assert(i >= 0 && i < c->domain_size);
                        assert(c->ptn[i] + 1 > 0);
                        const int v_color = i;
                        const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
                        largest_color_class_size = mark_as_largest ? c->ptn[i] + 1 : largest_color_class_size;

                        report_splits->push_back(std::pair<std::pair<int, int>, bool>(
                                std::pair<int, int>(col, i), mark_as_largest));
                        if (v_color != 0)
                            c->ptn[v_color - 1] = 0;
                        i += c->ptn[i] + 1;
                    }
                }

                return comp;
            }

            bool refine_color_class_dense_dense(sgraph *g, coloring *c,
                                                int color_class, int class_size,
                                                work_list_pair_bool *report_splits) {
                bool comp;
                int i, j, acc, cc, largest_color_class_size;
                cc = color_class; // iterate over color class
                comp = true;

                neighbours.reset();
                old_color_classes.reset();

                const int end_cc = color_class + class_size;
                while (cc < end_cc) { // increment value of neighbours of vc by 1
                    const int vc = c->lab[cc];
                    const int pe = g->v[vc];
                    const int end_i = pe + g->d[vc];
                    for (i = pe; i < end_i; i++) {
                        const int v = g->e[i];
                        neighbours.inc_nr(v);
                    }
                    cc += 1;
                }

                // in dense-dense we dont sort, instead we just iterate over all cells (which are already sorted)
                old_color_classes.sort();
                // for every cell...
                for (j = 0; j < g->v_size;) {
                    const int col = j;
                    j += c->ptn[j] + 1;
                    const int col_sz = c->ptn[col] + 1;
                    if (col_sz == 1) {
                        const int v_degree = neighbours.get(c->lab[col]) + 1;
                        if (v_degree == -1) // not connected! skip...
                            continue;
                        //comp = I->write_top_and_compare(col + v_degree * g->v_size) && comp;
                        continue;
                    }

                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    const int first_index = neighbours.get(c->lab[col]) + 1;
                    vertex_worklist.push_back(first_index);
                    int total = 0;
                    for (i = col; i < col + col_sz; ++i) {
                        const int v = c->lab[i];
                        const int index = neighbours.get(v) + 1;
                        if (index == first_index)
                            continue;
                        total += 1;
                        if (neighbour_sizes.inc(index) == 0)
                            vertex_worklist.push_back(index);
                    }

                    neighbour_sizes.inc(first_index);
                    neighbour_sizes.set(first_index, col_sz - total - 1);

                    if (vertex_worklist.cur_pos == 1) {
                        // no split
                        continue;
                    }

                    vertex_worklist.sort();
                    // enrich neighbour_sizes to accumulative counting array
                    acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int k = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(k) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(k, val + acc);
                            acc += val;
                            const int _col = col + col_sz - (neighbour_sizes.get(k));
                            c->ptn[_col] = -1; // this is val - 1, actually...
                        }
                    }

                    vertex_worklist.reset();

                    // copy cell for rearranging
                    memcpy(vertex_worklist.get_array(), c->lab + col, col_sz * sizeof(int));
                    vertex_worklist.cur_pos = col_sz;

                    // determine colors and rearrange
                    // split color classes according to count in counting array
                    while (!vertex_worklist.empty()) {
                        const int v = vertex_worklist.pop_back();
                        const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
                        // i could immediately determine size of v_new_color here and write invariant before rearrange?
                        assert(v_new_color >= col);
                        assert(v_new_color < col + col_sz);

                        c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                        c->vertex_to_col[v] = v_new_color;
                        c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                        c->ptn[v_new_color] += 1;
                    }

                    largest_color_class_size = -1;

                    // determine largest class to throw away and finish (fourth iteration)
                    for (i = col; i < col + col_sz;) {
                        assert(i >= 0 && i < c->domain_size);
                        assert(c->ptn[i] + 1 > 0);
                        const int v_color = i;
                        const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
                        largest_color_class_size = mark_as_largest ? c->ptn[i] + 1 : largest_color_class_size;

                        report_splits->push_back(std::pair<std::pair<int, int>, bool>(
                                std::pair<int, int>(col, i), mark_as_largest));
                        if (v_color != 0)
                            c->ptn[v_color - 1] = 0;
                        i += c->ptn[i] + 1;
                    }
                }

                neighbours.reset_hard();
                return comp;
            }

            void refine_color_class_singleton(sgraph *g, coloring *c, int color_class,
                                              work_list_pair_bool *report_splits) {
                int i, cc, deg1_write_pos, deg1_read_pos;
                cc = color_class; // iterate over color class

                neighbours.reset();
                vertex_worklist.reset();
                old_color_classes.reset();

                const int vc = c->lab[cc];
                const int pe = g->v[vc];
                const int end_i = pe + g->d[vc];
                for (i = pe; i < end_i; i++) {
                    const int v = g->e[i];
                    const int col = c->vertex_to_col[v];

                    if (c->ptn[col] == 0) {
                        // if full invariant, sort -- else use a hash value
                        //singleton_inv += MASH5(col);
                        continue;
                    }

                    if(neighbours.get(col) == -1) {
                        old_color_classes.push_back(col);
                        // neighbours acts as the write position for degree 1 vertices of col in the scratchpad
                        neighbours.set(col, col);
                    }
                    // write vertex to scratchpad for later use
                    scratch[neighbours.get(col)] = v;
                    neighbours.inc_nr(col); // we reset neighbours later, since old_color_classes for reset
                }

                //singleton_inv += old_color_classes.cur_pos;
                //if(add_hook) (write_to_trace)(singleton_inv);

                old_color_classes.sort();

                // write invariant first...
                //for(i = 0; i < old_color_classes.cur_pos && comp; ++i) {
                //comp = comp && I->write_top_and_compare(old_color_classes.arr[i]); // color class
                //if(write_to_trace) (add_hook)(neighbours.get(old_color_classes[i]));
                // contains information about color degree (= 1)
                //}

                while (!old_color_classes.empty()) {
                    const int deg0_col = old_color_classes.pop_back();
                    const int deg1_col_sz = neighbours.get(deg0_col) - deg0_col;
                    const int deg0_col_sz = (c->ptn[deg0_col] + 1) - deg1_col_sz;
                    const int deg1_col = deg0_col + deg0_col_sz;

                    assert(c->vertex_to_col[c->lab[deg0_col]] == deg0_col);

                    // no split? done...
                    if (deg0_col == deg1_col) {
                        neighbours.set(deg1_col, -1);
                        assert(c->vertex_to_col[c->lab[deg0_col]] == deg0_col);
                        continue;
                    }

                    assert(deg0_col_sz + deg1_col_sz - 1 == c->ptn[deg0_col]);

                    // set ptn
                    c->ptn[deg0_col] = deg0_col_sz - 1;
                    c->ptn[deg1_col] = deg1_col_sz - 1;
                    c->ptn[deg1_col - 1] = 0;

                    deg1_write_pos = deg1_col;
                    deg1_read_pos = neighbours.get(deg0_col) - 1;

                    //c->vertex_to_col[c->lab[deg1_col]] = deg1_col;

                    // rearrange vertices of deg1 to the back of deg0 color
                    assert(deg1_read_pos >= deg0_col);

                    while (deg1_read_pos >= deg0_col) {
                        const int v = scratch[deg1_read_pos];
                        const int vertex_at_pos = c->lab[deg1_write_pos];
                        const int lab_pos = c->vertex_to_lab[v];

                        c->lab[deg1_write_pos] = v;
                        c->vertex_to_lab[v] = deg1_write_pos;
                        c->vertex_to_col[v] = deg1_col;

                        c->lab[lab_pos] = vertex_at_pos;
                        c->vertex_to_lab[vertex_at_pos] = lab_pos;

                        deg1_write_pos++;
                        deg1_read_pos--;
                    }

                    assert(c->vertex_to_col[c->lab[deg0_col]] == deg0_col);
                    assert(c->vertex_to_col[c->lab[deg1_col]] == deg1_col);

                    // add new classes to report_splits
                    const bool leq = deg1_col_sz > deg0_col_sz;
                    report_splits->push_back(std::pair<std::pair<int, int>, bool>(
                            std::pair<int, int>(deg0_col, deg0_col), !leq));
                    report_splits->push_back(std::pair<std::pair<int, int>, bool>(
                            std::pair<int, int>(deg0_col, deg1_col), leq));

                    // reset neighbours count to -1
                    neighbours.set(deg0_col, -1);
                }
            }

            bool refine_color_class_singleton_first(sgraph *g, coloring *c, int color_class,
                                                    work_list_pair_bool *report_splits) {
                bool comp;
                int i, cc, deg1_write_pos, deg1_read_pos;
                cc = color_class; // iterate over color class
                comp = true;

                neighbours.reset();
                scratch_set.reset();
                vertex_worklist.reset();
                old_color_classes.reset();

                const int vc = c->lab[cc];
                const int pe = g->v[vc];
                const int end_i = pe + g->d[vc];
                for (i = pe; i < end_i; i++) {
                    const int v = g->e[i];
                    const int col = c->vertex_to_col[v];

                    if (c->ptn[col] == 0)
                        continue;

                    if (!scratch_set.get(col)) {
                        scratch_set.set(col);
                        old_color_classes.push_back(col);
                        // neighbours acts as the write position for degree 1 vertices of col
                        // in the singleton scratchpad
                        neighbours.set(col, col);
                    }
                    // write vertex to singleton scratchpad for later use
                    scratch[neighbours.get(col)] = v;
                    neighbours.inc_nr(col); // we reset neighbours later...
                    // old_color_classes can be used as reset information
                }

                while (!old_color_classes.empty()) {
                    const int deg0_col = old_color_classes.pop_back();
                    const int deg1_col_sz = neighbours.get(deg0_col) - deg0_col;
                    const int deg0_col_sz = (c->ptn[deg0_col] + 1) - deg1_col_sz;
                    const int deg1_col = deg0_col + deg0_col_sz;

                    // no split? done...
                    if (deg0_col == deg1_col) {
                        neighbours.set(deg1_col, -1);
                        continue;
                    }

                    // set ptn
                    c->ptn[deg0_col] = deg0_col_sz - 1;
                    c->ptn[deg1_col] = deg1_col_sz - 1;
                    c->ptn[deg1_col - 1] = 0;

                    deg1_write_pos = deg1_col;
                    deg1_read_pos = neighbours.get(deg0_col) - 1;

                    // rearrange vertices of deg1 to the back of deg0 color
                    while (deg1_read_pos >= deg0_col) {
                        const int v = scratch[deg1_read_pos];
                        const int vertex_at_pos = c->lab[deg1_write_pos];
                        const int lab_pos = c->vertex_to_lab[v];

                        c->lab[deg1_write_pos] = v;
                        c->vertex_to_lab[v] = deg1_write_pos;
                        c->vertex_to_col[v] = deg1_col;
                        c->lab[lab_pos] = vertex_at_pos;
                        c->vertex_to_lab[vertex_at_pos] = lab_pos;

                        deg1_write_pos++;
                        deg1_read_pos--;
                    }

                    // add new classes to report_splits
                    const bool leq = deg1_col_sz > deg0_col_sz;
                    report_splits->push_back(std::pair<std::pair<int, int>, bool>(
                            std::pair<int, int>(deg0_col, deg0_col), !leq));
                    report_splits->push_back(std::pair<std::pair<int, int>, bool>(
                            std::pair<int, int>(deg0_col, deg1_col), leq));

                    // reset neighbours count to -1
                    neighbours.set(deg0_col, -1);
                }

                return comp;
            }

            bool refine_color_class_dense_first(sgraph *g, coloring *c,
                                                int color_class, int class_size,
                                                work_list_pair_bool *report_splits) {
                int i, cc, acc, largest_color_class_size, pos;
                cc = color_class; // iterate over color class

                neighbours.reset();
                scratch_set.reset();
                old_color_classes.reset();

                const int end_cc = color_class + class_size;
                while (cc < end_cc) { // increment value of neighbours of vc by 1
                    const int vc = c->lab[cc];
                    const int pe = g->v[vc];
                    const int end_i = pe + g->d[vc];

                    for (i = pe; i < end_i; i++) {
                        const int v = g->e[i];
                        const int col = c->vertex_to_col[v];
                        if (c->ptn[col] == 0)
                            continue;

                        neighbours.inc(v);
                        if (!scratch_set.get(col)) {
                            scratch_set.set(col);
                            old_color_classes.push_back(col);
                        }
                    }
                    cc += 1;
                }

                while (!old_color_classes.empty()) {
                    const int col = old_color_classes.pop_back();
                    const int col_sz = c->ptn[col] + 1;

                    if (col_sz == 1) {
                        const int v_degree = neighbours.get(c->lab[col]) + 1;
                        if (v_degree == -1)
                            continue;
                    }

                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    const int first_index = neighbours.get(c->lab[col]) + 1;
                    vertex_worklist.push_back(first_index);
                    int total = 0;
                    for (i = col; i < col + col_sz; ++i) {
                        const int v = c->lab[i];
                        const int index = neighbours.get(v) + 1;
                        if (index == first_index)
                            continue;
                        total += 1;
                        if (neighbour_sizes.inc(index) == 0)
                            vertex_worklist.push_back(index);
                    }

                    neighbour_sizes.inc(first_index);
                    neighbour_sizes.set(first_index, col_sz - total - 1);

                    if (vertex_worklist.cur_pos == 1) {
                        continue;
                    }

                    // enrich neighbour_sizes to accumulative counting array
                    acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int k = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(k) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(k, val + acc);
                            acc += val;
                            const int _col = col + col_sz - (neighbour_sizes.get(k));
                            c->ptn[_col] = -1; // this is val - 1, actually...
                        }
                    }

                    // copy cell for rearranging
                    memcpy(scratch, c->lab + col, col_sz * sizeof(int));
                    pos = col_sz;

                    // determine colors and rearrange
                    // split color classes according to count in counting array
                    while (pos > 0) {
                        const int v = scratch[--pos];
                        const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
                        // i could immediately determine size of v_new_color here and write invariant before rearrange?
                        assert(v_new_color >= col);
                        assert(v_new_color < col + col_sz);

                        c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                        c->vertex_to_col[v] = v_new_color;
                        c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                        c->ptn[v_new_color] += 1;
                    }

                    largest_color_class_size = -1;

                    // determine largest class to throw away and finish (fourth iteration)
                    for (i = col; i < col + col_sz;) {
                        assert(i >= 0 && i < c->domain_size);
                        assert(c->ptn[i] + 1 > 0);
                        const int v_color = i;
                        const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
                        largest_color_class_size = mark_as_largest ? c->ptn[i] + 1 : largest_color_class_size;

                        report_splits->push_back(std::pair<std::pair<int, int>, bool>(
                                std::pair<int, int>(col, i), mark_as_largest));
                        if (v_color != 0)
                            c->ptn[v_color - 1] = 0;
                        i += c->ptn[i] + 1;
                    }
                }

                return true;
            }

            bool refine_color_class_dense_dense_first(sgraph *g, coloring *c,
                                                      int color_class, int class_size,
                                                      work_list_pair_bool *report_splits) {
                // for all vertices of the color class...
                int i, j, cc, acc, largest_color_class_size, pos;
                cc = color_class; // iterate over color class

                neighbours.reset();
                old_color_classes.reset();

                const int end_cc = color_class + class_size;
                while (cc < end_cc) { // increment value of neighbours of vc by 1
                    const int vc = c->lab[cc];
                    const int pe = g->v[vc];
                    const int end_i = pe + g->d[vc];
                    for (i = pe; i < end_i; i++) {
                        const int v = g->e[i];
                        neighbours.inc_nr(v);
                    }
                    cc += 1;
                }

                // dont sort, just iterate over all cells
                for (j = 0; j < g->v_size;) {
                    const int col = j;
                    const int c_ptn = c->ptn[j];
                    j += c_ptn + 1;
                    const int col_sz = c_ptn + 1;

                    if (col_sz == 1)
                        continue;

                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    const int first_index = neighbours.get(c->lab[col]) + 1;
                    vertex_worklist.push_back(first_index);
                    int total = 0;
                    for (i = col; i < col + col_sz; ++i) {
                        const int v = c->lab[i];
                        const int index = neighbours.get(v) + 1;
                        if (index == first_index)
                            continue;
                        total += 1;
                        if (neighbour_sizes.inc(index) == 0)
                            vertex_worklist.push_back(index);
                    }

                    neighbour_sizes.inc(first_index);
                    neighbour_sizes.set(first_index, col_sz - total - 1);

                    if (vertex_worklist.cur_pos == 1) {
                        continue;
                    }

                    // enrich neighbour_sizes to accumulative counting array
                    acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int k = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(k) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(k, val + acc);
                            acc += val;
                            const int _col = col + col_sz - (neighbour_sizes.get(k));
                            c->ptn[_col] = -1;
                        }
                    }

                    // copy cell for rearranging
                    memcpy(scratch, c->lab + col, col_sz * sizeof(int));
                    pos = col_sz;

                    // determine colors and rearrange
                    // split color classes according to count in counting array
                    while (pos > 0) {
                        const int v = scratch[--pos];
                        const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
                        assert(v_new_color >= col);
                        assert(v_new_color < col + col_sz);

                        c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
                        c->vertex_to_col[v] = v_new_color;
                        c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
                        c->ptn[v_new_color] += 1;
                    }

                    largest_color_class_size = -1;

                    // determine largest class to throw away and finish (fourth iteration)
                    for (i = col; i < col + col_sz;) {
                        assert(i >= 0 && i < c->domain_size);
                        assert(c->ptn[i] + 1 > 0);
                        const int v_color = i;
                        const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
                        largest_color_class_size = mark_as_largest ? c->ptn[i] + 1 : largest_color_class_size;

                        report_splits->push_back(std::pair<std::pair<int, int>, bool>(
                                std::pair<int, int>(col, i), mark_as_largest));
                        if (v_color != 0)
                            c->ptn[v_color - 1] = 0;
                        i += c->ptn[i] + 1;
                    }
                }

                neighbours.reset_hard();
                return true;
            }

            bool refine_color_class_sparse_first(sgraph *g, coloring *c,
                                                 int color_class, int class_size,
                                                 work_list_pair_bool *report_splits) {
                bool comp;
                int v_new_color, cc, largest_color_class_size, acc;

                cc = color_class; // iterate over color class
                comp = true;

                scratch_set.reset();
                old_color_classes.reset();
                neighbours.reset();
                color_vertices_considered.reset();

                const int end_cc = color_class + class_size;
                while (cc < end_cc) { // increment value of neighbours of vc by 1
                    const int vc = c->lab[cc];
                    const int pe = g->v[vc];
                    const int end_i = pe + g->d[vc];
                    for (int i = pe; i < end_i; i++) {
                        const int v = g->e[i];
                        const int col = c->vertex_to_col[v];

                        if (c->ptn[col] == 0)
                            continue;

                        neighbours.inc_nr(v);
                        if (neighbours.get(v) == 0) {
                            color_vertices_considered.inc_nr(col);
                            assert(col + color_vertices_considered.get(col) < g->v_size);
                            scratch[col + color_vertices_considered.get(col)] = v; // hit vertices
                            if (!scratch_set.get(col)) {
                                old_color_classes.push_back(col);
                                scratch_set.set(col);
                            }
                        }
                    }
                    cc += 1;
                }

                // split color classes according to neighbour count
                while (!old_color_classes.empty()) {
                    const int _col = old_color_classes.pop_back();
                    const int _col_sz = c->ptn[_col] + 1;

                    neighbour_sizes.reset();
                    vertex_worklist.reset();

                    for (int i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
                        const int v = scratch[_col + i];
                        int index = neighbours.get(v) + 1;
                        if (neighbour_sizes.inc(index) == 0)
                            vertex_worklist.push_back(index);
                    }

                    // enrich neighbour_sizes to accumulative counting array
                    acc = 0;
                    while (!vertex_worklist.empty()) {
                        const int i = vertex_worklist.pop_back();
                        const int val = neighbour_sizes.get(i) + 1;
                        if (val >= 1) {
                            neighbour_sizes.set(i, val + acc);
                            acc += val;
                        }
                    }

                    // determine colors
                    vertex_worklist.reset();

                    for (int i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
                        const int v = scratch[_col + i];
                        if (neighbours.get(v) == -1) {
                            v_new_color = _col;
                            assert(false);
                        } else {
                            assert((neighbours.get(v) + 1) > 0);
                            v_new_color = _col + _col_sz - neighbour_sizes.get(neighbours.get(v) + 1);
                        }
                        if (v_new_color != _col) {
                            vertex_worklist.push_back(v_new_color);
                            vertex_worklist.push_back(v);
                            c->ptn[v_new_color] = -1;
                        }
                        neighbours.set(v, -1);
                    }

                    color_vertices_considered.set(_col, -1);

                    // rearrange vertices
                    while (!vertex_worklist.empty()) {
                        const int v = vertex_worklist.pop_back();
                        const int v_new_color2 = vertex_worklist.pop_back();

                        const int vertex_old_pos = c->vertex_to_lab[v];
                        const int vertex_at_pos = c->lab[v_new_color2 + c->ptn[v_new_color2] + 1];
                        c->lab[vertex_old_pos] = vertex_at_pos;
                        c->vertex_to_lab[vertex_at_pos] = vertex_old_pos;

                        c->lab[v_new_color2 + c->ptn[v_new_color2] + 1] = v;
                        c->vertex_to_col[v] = v_new_color2;
                        c->vertex_to_lab[v] = v_new_color2 + c->ptn[v_new_color2] + 1;
                        c->ptn[v_new_color2] += 1;

                        if (_col != v_new_color2) {
                            assert(v_new_color2 > _col);
                            c->ptn[_col] -= 1;
                        } else
                            assert(false);
                    }

                    // add new colors to worklist
                    largest_color_class_size = -1;
                    int i;
                    for (i = _col; i < _col + _col_sz;) {
                        assert(i >= 0 && i < c->domain_size);
                        assert(c->ptn[i] + 1 > 0);
                        const bool mark_as_largest = largest_color_class_size < (c->ptn[i] + 1);
                        largest_color_class_size = mark_as_largest ? (c->ptn[i] + 1) : largest_color_class_size;

                        report_splits->push_back(std::pair<std::pair<int, int>, bool>(
                                std::pair<int, int>(_col, i), mark_as_largest));
                        if (i != 0)
                            c->ptn[i - 1] = 0;
                        i += c->ptn[i] + 1;
                    }

                    if (i != 0)
                        c->ptn[i - 1] = 0;
                }

                neighbour_sizes.reset();
                vertex_worklist.reset();

                return comp;
            }
        };
    }
}

#endif //DEJAVU_REFINEMENT_H
