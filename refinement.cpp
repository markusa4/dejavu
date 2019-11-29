#include "refinement.h"
#include "configuration.h"
#include "utility.h"
#include <list>
#include <tuple>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <cstring>

struct pairhash {
public:
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &x) const
    {
        return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
    }
};

bool refinement::certify_automorphism(sgraph *g, bijection* p) {
    assert(p->map_sz == g->v_size);

    for(int i = 0; i < g->v_size; ++i) {
        int image_i = p->map_vertex(i);
        if(image_i == i)
            continue;
        if(g->d[i] != g->d[image_i]) // degrees must be equal
            return false;

        scratch_set.reset();
        // automorphism must preserve neighbours
        int found = 0;
        for(int j = g->v[i]; j < g->v[i] + g->d[i]; ++j) {
            int vertex_j = g->e[j];
            int image_j  = p->map_vertex(vertex_j);
            scratch_set.set(image_j);
            found += 1;
        }
        for(int j = g->v[image_i]; j < g->v[image_i] + g->d[image_i]; ++j) {
            int vertex_j = g->e[j];
            if(!scratch_set.get(vertex_j)) {
                return false;
            }
            scratch_set.unset(vertex_j);
            found -= 1;
        }
        if(found != 0) {
            return false;
        }
    }

    return true;
}

bool refinement::refine_coloring(sgraph *g, coloring *c, change_tracker *changes, invariant* I, int init_color_class, bool track_changes, strategy_metrics* m) {
    bool comp = true;
    assure_initialized(g);

    cell_todo.reset(&queue_pointer);

    if(init_color_class < 0) {
        // initialize queue with all classes (except for largest one)
        for (int i = 0; i < c->ptn_sz;) {
            cell_todo.add_cell(&queue_pointer, i, c->ptn[i] + 1);
            i += c->ptn[i] + 1;
        }
    } else {
        cell_todo.add_cell(&queue_pointer, init_color_class, c->ptn[init_color_class] + 1);
    }
    int its = 0;

    while(!cell_todo.empty()) {
        its += 1;
        color_class_splits.reset();
        std::pair<int, int> next_color_class = cell_todo.next_cell(&queue_pointer);

        if(m)
            m->color_refinement_cost += next_color_class.second;

        // write color class and size to invariant
        comp = comp && I->write_top_and_compare(next_color_class.first);
        comp = comp && I->write_top_and_compare(next_color_class.second);

        bool dense_dense = (g->d[c->lab[next_color_class.first]] > (g->v_size / (next_color_class.second + 1)));


        if(next_color_class.second == 1 && !(config.CONFIG_IR_DENSE && dense_dense)) {
            // SINGLETON
            comp = comp && refine_color_class_singleton(g, c, next_color_class.first, next_color_class.second,
                                                        &color_class_splits, I);
        } else if(config.CONFIG_IR_DENSE) {
            if(dense_dense) { // DENSE-DENSE
                if(g->max_degree > 127) {
                    comp = comp && refine_color_class_dense_dense(g, c, next_color_class.first, next_color_class.second,
                                                                  &color_class_splits, I);
                } else {
                    comp = comp && refine_color_class_dense_dense127(g, c, next_color_class.first, next_color_class.second,
                                                                  &color_class_splits, I);
                }
            } else { // DENSE-SPARSE
                comp = comp && refine_color_class_dense(g, c, next_color_class.first, next_color_class.second,
                                               &color_class_splits, I);
            }
        } else { // SPARSE
            if(g->max_degree > 127) {
                comp = comp &&
                       refine_color_class_sparse(g, c, next_color_class.first, next_color_class.second,
                                                 &color_class_splits, I);
            } else {
                comp = comp &&
                       refine_color_class_sparse127(g, c, next_color_class.first, next_color_class.second,
                                                 &color_class_splits, I);
            }
        }
        if(!comp)
            break;

        // add all new classes except for the first, largest one
        int skip = 0;
        //for(auto cc = color_class_splits.begin(); cc != color_class_splits.end(); ++cc) {

        int  latest_old_class = -1;
        bool skipped_largest = false;

        // color class splits are sorted in reverse
        // the old color class will always come last
        while(!color_class_splits.empty()) {
            int  old_class  = color_class_splits.last()->first.first;
            int  new_class  = color_class_splits.last()->first.second;
            bool is_largest = color_class_splits.last()->second;

            if(track_changes)
                changes->track(old_class);

            if(latest_old_class != old_class) {
                latest_old_class = old_class;
                skipped_largest = false;
            }

            color_class_splits.pop_back();
            int new_class_sz = c->ptn[new_class] + 1;


            if(skipped_largest || !is_largest) {
                cell_todo.add_cell(&queue_pointer, new_class, new_class_sz);
            } else {
                skipped_largest = true;
                skip += 1;

                // since old color class will always appear last, the queue pointer of old color class is still valid!
                int i = queue_pointer.get(old_class);
                if(i >= 0) {
                    cell_todo.replace_cell(&queue_pointer, old_class, new_class, new_class_sz);
                }
            }
        }
        if(!comp) break;
    }

    assert(comp?assert_is_equitable(g, c):true);

    return comp;
}

void refinement::assure_initialized(sgraph *g) {
    if(!initialized) {
        color_workset.initialize(g->v_size);
        vertex_worklist.initialize(g->v_size);
        scratch_set.initialize(g->v_size);
        singletons.initialize(g->v_size);
        color_worklist_color.initialize(g->v_size);
        color_class_splits.initialize(g->v_size);
        queue_pointer.initialize(g->v_size);
        degrees_worklist.initialize(g->v_size);
        neighbours.initialize(g->v_size);
        neighbours127.initialize(g->v_size);
        color_vertices_considered.initialize(g->v_size);
        neighbour_sizes.initialize(g->v_size);
        old_color_classes.initialize(g->v_size);
        cell_todo.initialize(g->v_size * 2);
        scratch = new int[g->v_size];
        memset(scratch, 0, g->v_size * sizeof(int));
        initialized = true;
    }
}


bool refinement::refine_coloring_first(sgraph *g, coloring *c, int init_color_class) {
    bool comp = true;
    assure_initialized(g);

    cell_todo.reset(&queue_pointer);

    if(init_color_class < 0) {
        for (int i = 0; i < c->ptn_sz;) {
            cell_todo.add_cell(&queue_pointer, i, c->ptn[i] + 1);
            i += c->ptn[i] + 1;
        }
    } else {
        cell_todo.add_cell(&queue_pointer, init_color_class, c->ptn[init_color_class] + 1);
    }
    int its = 0;

    while(!cell_todo.empty()) {
        its += 1;
        color_class_splits.reset();
        std::pair<int, int> next_color_class = cell_todo.next_cell(&queue_pointer);

        bool dense_dense = (g->d[c->lab[next_color_class.first]] > (g->v_size / (next_color_class.second + 1)));


        if(next_color_class.second == 1 && !(config.CONFIG_IR_DENSE && dense_dense)) {
            // SINGLETON
            comp = comp && refine_color_class_singleton_first(g, c, next_color_class.first, next_color_class.second,
                                                        &color_class_splits);
        } else if(config.CONFIG_IR_DENSE) {
            if(dense_dense) { // DENSE-DENSE
                comp = comp && refine_color_class_dense_dense_first(g, c, next_color_class.first, next_color_class.second,
                                                              &color_class_splits);
            } else { // DENSE-SPARSE
                comp = comp && refine_color_class_dense_first(g, c, next_color_class.first, next_color_class.second,
                                                        &color_class_splits);
            }
        } else { // SPARSE
            comp = comp &&
                    refine_color_class_sparse_first(g, c, next_color_class.first, next_color_class.second, &color_class_splits);
        }
        if(!comp)
            break;

        // add all new classes except for the first, largest one
        int skip = 0;

        int  latest_old_class = -1;
        bool skipped_largest = false;

        // color class splits are sorted in reverse
        // the old color class will always come last
        while(!color_class_splits.empty()) {
            int  old_class  = color_class_splits.last()->first.first;
            int  new_class  = color_class_splits.last()->first.second;
            bool is_largest = color_class_splits.last()->second;

            if(latest_old_class != old_class) {
                latest_old_class = old_class;
                skipped_largest = false;
            }

            color_class_splits.pop_back();
            int new_class_sz = c->ptn[new_class] + 1;


            if(skipped_largest || !is_largest) {
                cell_todo.add_cell(&queue_pointer, new_class, new_class_sz);
            } else {
                skipped_largest = true;
                skip += 1;

                // since old color class will always appear last, the queue pointer of old color class is still valid!
                int i = queue_pointer.get(old_class);
                if(i >= 0) {
                    cell_todo.replace_cell(&queue_pointer, old_class, new_class, new_class_sz);
                }
            }
        }
        if(!comp) break;
    }

    assert(c->check());

    //std::cout << "its: " << its << std::endl;

    return comp;
}

bool refinement::refine_color_class_sparse(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist, invariant* I) {
    // for all vertices of the color class...
    bool comp, mark_as_largest;
    int i, cc, vc, pe, end_i, end_cc, v, col, v_new_color, vertex_old_pos, vertex_at_pos, v_degree,
                     largest_color_class_size, acc_in, _col, _col_sz, acc, val, __col, singleton_inv, total;

    cc = color_class; // iterate over color class
    comp = true;

    singletons.reset();
    color_workset.reset();
    old_color_classes.reset();
    neighbours.reset();
    color_vertices_considered.reset();

    end_cc = color_class + class_size;
    acc_in = 0;
    singleton_inv = 0;
    while (cc < end_cc) { // increment value of neighbours of vc by 1
        vc = c->lab[cc];
        pe = g->v[vc];
        end_i = pe + g->d[vc];
        for (i = pe; i < end_i; i++) {
            v   = g->e[i];
            col = c->vertex_to_col[v];
            neighbours.inc_nr(v);
            if(neighbours.get(v) == 0) {
                color_vertices_considered.inc_nr(col);
                assert(col + color_vertices_considered.get(col) < g->v_size);
                scratch[col + color_vertices_considered.get(col)] = v; // hit vertices
                if (!color_workset.get(col)) {
                    if (c->ptn[col] == 0)
                        singletons.push_back(col);
                    else
                        old_color_classes.push_back(col);
                    acc_in += (col + 1);
                    color_workset.set_nr(col);
                }
            }
        }
        cc += 1;
    }
    // write invariant for singleton color classes
    // could also hash this...
    singletons.sort();
    while (!singletons.empty()) {
        col = singletons.pop_back();
        color_workset.unset(col);
        v_degree = neighbours.get(c->lab[col]);
        neighbours.set(c->lab[col], -1);
        color_vertices_considered.set(col, -1);
        comp = comp && I->write_top_and_compare(-col);
        comp = comp && I->write_top_and_compare(v_degree);
    }
    comp = comp && I->write_top_and_compare(INT32_MAX - 9);
    comp = comp && I->write_top_and_compare(acc_in);

    // early out before sorting color classes
    if(!comp) {
        while (!old_color_classes.empty()) {
            _col = old_color_classes.pop_back();
            color_workset.unset(_col);
            for(i = 0; i < color_vertices_considered.get(_col) + 1; ++i)
                neighbours.set(scratch[_col + i], -1);
            color_vertices_considered.set(_col, -1);
        }
        return comp;
    }

    // sort split color classes
    old_color_classes.sort();

    // split color classes according to neighbour count
    while(!old_color_classes.empty()) {
        _col    = old_color_classes.pop_back();
        _col_sz = c->ptn[_col] + 1;
        //std::cout << "it " << _col << ", " << _col_sz << std::endl;

        color_workset.unset(_col);
        neighbour_sizes.reset();
        degrees_worklist.reset();

        total = 0;
        for(i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
            v = scratch[_col + i];
            int index = neighbours.get(v) + 1;
            total += 1;
            if(neighbour_sizes.inc(index) == 0)
                degrees_worklist.push_back(index);
        }

        degrees_worklist.sort();

        // enrich neighbour_sizes to accumulative counting array
        acc = 0;
        comp = comp && I->write_top_and_compare(INT32_MAX - 5);
        while(!degrees_worklist.empty()) {
            i   = degrees_worklist.pop_back();
            val = neighbour_sizes.get(i) + 1;
            if(val >= 1) {
                neighbour_sizes.set(i, val + acc);
                acc += val;
                __col = _col + _col_sz - (neighbour_sizes.get(i));
                v_degree = i;
                comp = comp && I->write_top_and_compare(-__col);
                comp = comp && I->write_top_and_compare(v_degree);
                comp = comp && I->write_top_and_compare(val + 1);
            }
        }

        // determine colors
        color_worklist_color.reset();

        for(i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
            v = scratch[_col + i];
            if (neighbours.get(v) == -1) {
                v_new_color = _col; assert(false);
            } else {
                assert((neighbours.get(v) + 1) > 0);
                v_new_color = _col + _col_sz - neighbour_sizes.get(neighbours.get(v) + 1);
            }
            if (v_new_color != _col) {
                color_worklist_color.push_back(std::pair<int, int>(v, v_new_color));
                c->ptn[v_new_color] = -1;
            }
            neighbours.set(v, -1);
        }

        color_vertices_considered.set(_col, -1);

        // rearrange vertices
        while (!color_worklist_color.empty()) {
            v           = color_worklist_color.last()->first;
            v_new_color = color_worklist_color.last()->second;
            color_worklist_color.pop_back();

            vertex_old_pos = c->vertex_to_lab[v];
            vertex_at_pos  = c->lab[v_new_color + c->ptn[v_new_color] + 1];
            c->lab[vertex_old_pos]          = vertex_at_pos;
            c->vertex_to_lab[vertex_at_pos] = vertex_old_pos;

            c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
            c->vertex_to_col[v] = v_new_color;
            c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
            c->ptn[v_new_color] += 1;

            if (_col != v_new_color) {
                assert(v_new_color > _col);
                c->ptn[_col] -= 1;
            } else assert(false);
        }

        // add new colors to worklist
        largest_color_class_size = -1;
        for(i = _col; i < _col + _col_sz;) {
            assert(i >= 0 && i < c->ptn_sz);
            assert(c->ptn[i] + 1 > 0);
            mark_as_largest = largest_color_class_size < (c->ptn[i] + 1);
            largest_color_class_size = mark_as_largest?(c->ptn[i] + 1):largest_color_class_size;

            color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(std::pair<int, int>(_col, i), mark_as_largest));
            if(i != 0)
                c->ptn[i - 1] = 0;
            i += c->ptn[i] + 1;
        }

        if(i != 0)
            c->ptn[i - 1] = 0;

        // early out, but later
        if(!comp) {
            while (!old_color_classes.empty()) {
                _col = old_color_classes.pop_back();
                color_workset.unset(_col);
                for(i = 0; i < color_vertices_considered.get(_col) + 1; ++i)
                    neighbours.set(scratch[_col + i], -1);
                color_vertices_considered.set(_col, -1);
            }
            return comp;
        }
    }

    neighbour_sizes.reset();
    degrees_worklist.reset();
    color_worklist_color.reset();
    color_workset.reset();

    return comp;
}

bool refinement::refine_color_class_sparse127(sgraph *g, coloring *c, const int color_class, const int class_size, work_list_pair_bool* color_class_split_worklist, invariant* I) {
    // for all vertices of the color class...
    bool comp, mark_as_largest;
    int i, cc, end_cc,
            largest_color_class_size, acc_in, acc, singleton_inv, total;

    cc = color_class; // iterate over color class
    comp = true;

    singletons.reset();
    color_workset.reset();
    old_color_classes.reset();
    neighbours127.reset();
    color_vertices_considered.reset();

    end_cc = color_class + class_size;
    acc_in = 0;
    singleton_inv = 0;
    while (cc < end_cc) { // increment value of neighbours of vc by 1
        const int vc = c->lab[cc];
        const int pe = g->v[vc];
        const int end_i = pe + g->d[vc];
        for (i = pe; i < end_i; i++) {
            const int v   = g->e[i];
            const int col = c->vertex_to_col[v];
            neighbours127.inc_nr(v);
            if(neighbours127.get(v) == 0) {
                color_vertices_considered.inc_nr(col);
                assert(col + color_vertices_considered.get(col) < g->v_size);
                scratch[col + color_vertices_considered.get(col)] = v; // hit vertices
                if (!color_workset.get(col)) {
                    if (c->ptn[col] == 0)
                        singletons.push_back(col);
                    else
                        old_color_classes.push_back(col);
                    acc_in += (col + 1);
                    color_workset.set_nr(col);
                }
            }
        }
        cc += 1;
    }

    // write invariant for singleton color classes
    // could also hash this...
    singletons.sort();
    while (!singletons.empty()) {
        const int col = singletons.pop_back();
        color_workset.unset(col);
        const int vt = c->lab[col];
        const int v_degree = neighbours127.get(vt);
        neighbours127.set(vt, -1);
        color_vertices_considered.set(col, -1);
        comp = comp && I->write_top_and_compare(-col);
        comp = comp && I->write_top_and_compare(v_degree);
    }
    comp = comp && I->write_top_and_compare(INT32_MAX - 9);
    comp = comp && I->write_top_and_compare(acc_in);

    // early out before sorting color classes
    if(!comp) {
        while (!old_color_classes.empty()) {
            const int _col = old_color_classes.pop_back();
            color_workset.unset(_col);
            for(i = 0; i < color_vertices_considered.get(_col) + 1; ++i)
                neighbours127.set(scratch[_col + i], -1);
            color_vertices_considered.set(_col, -1);
        }
        return comp;
    }

    // sort split color classes
    old_color_classes.sort();

    // split color classes according to neighbour count
    while(!old_color_classes.empty()) {
        const int _col    = old_color_classes.pop_back();
        const int _col_sz = c->ptn[_col] + 1;
        //std::cout << "it " << _col << ", " << _col_sz << std::endl;

        color_workset.unset(_col);
        neighbour_sizes.reset();
        degrees_worklist.reset();

        total = 0;
        for(i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
            const int v = scratch[_col + i];
            int index = neighbours127.get(v) + 1;
            total += 1;
            if(neighbour_sizes.inc(index) == 0)
                degrees_worklist.push_back(index);
        }

        degrees_worklist.sort();

        // enrich neighbour_sizes to accumulative counting array
        acc = 0;
        comp = comp && I->write_top_and_compare(INT32_MAX - 5);
        while(!degrees_worklist.empty()) {
            const int i   = degrees_worklist.pop_back();
            const int val = neighbour_sizes.get(i) + 1;
            if(val >= 1) {
                neighbour_sizes.set(i, val + acc);
                acc += val;
                const int __col = _col + _col_sz - (neighbour_sizes.get(i));
                const int v_degree = i;
                comp = comp && I->write_top_and_compare(-__col);
                comp = comp && I->write_top_and_compare(v_degree);
                comp = comp && I->write_top_and_compare(val + 1);
            }
        }

        // determine colors
        color_worklist_color.reset();

        for(i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
            const int v = scratch[_col + i];
            int v_new_color;
            if (neighbours127.get(v) == -1) {
                v_new_color = _col; assert(false);
            } else {
                assert((neighbours127.get(v) + 1) > 0);
                v_new_color = _col + _col_sz - neighbour_sizes.get(neighbours127.get(v) + 1);
            }
            if (v_new_color != _col) {
                color_worklist_color.push_back(std::pair<int, int>(v, v_new_color));
                c->ptn[v_new_color] = -1;
            }
            neighbours127.set(v, -1);
        }

        color_vertices_considered.set(_col, -1);

        // rearrange vertices
        while (!color_worklist_color.empty()) {
            const int v           = color_worklist_color.last()->first;
            const int v_new_color = color_worklist_color.last()->second;
            color_worklist_color.pop_back();

            const int vertex_old_pos = c->vertex_to_lab[v];
            const int vertex_at_pos  = c->lab[v_new_color + c->ptn[v_new_color] + 1];
            c->lab[vertex_old_pos]          = vertex_at_pos;
            c->vertex_to_lab[vertex_at_pos] = vertex_old_pos;

            c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
            c->vertex_to_col[v] = v_new_color;
            c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
            c->ptn[v_new_color] += 1;

            if (_col != v_new_color) {
                assert(v_new_color > _col);
                c->ptn[_col] -= 1;
            } else assert(false);
        }

        // add new colors to worklist
        largest_color_class_size = -1;
        for(i = _col; i < _col + _col_sz;) {
            assert(i >= 0 && i < c->ptn_sz);
            assert(c->ptn[i] + 1 > 0);
            mark_as_largest = largest_color_class_size < (c->ptn[i] + 1);
            largest_color_class_size = mark_as_largest?(c->ptn[i] + 1):largest_color_class_size;

            color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(std::pair<int, int>(_col, i), mark_as_largest));
            if(i != 0)
                c->ptn[i - 1] = 0;
            i += c->ptn[i] + 1;
        }

        if(i != 0)
            c->ptn[i - 1] = 0;

        // early out, but later
        if(!comp) {
            while (!old_color_classes.empty()) {
                const int __col = old_color_classes.pop_back();
                color_workset.unset(__col);
                for(i = 0; i < color_vertices_considered.get(__col) + 1; ++i)
                    neighbours127.set(scratch[__col + i], -1);
                color_vertices_considered.set(__col, -1);
            }
            return comp;
        }
    }

    neighbour_sizes.reset();
    degrees_worklist.reset();
    color_worklist_color.reset();
    color_workset.reset();

    return comp;
}

bool refinement::refine_color_class_sparse_first(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist) {
    bool comp, mark_as_largest;
    int i, cc, vc, pe, end_i, end_cc, v, col, v_new_color, vertex_old_pos, vertex_at_pos, v_degree,
            largest_color_class_size, acc_in, _col, _col_sz, acc, val, __col;

    cc = color_class; // iterate over color class
    comp = true;

    singletons.reset();
    color_workset.reset();
    old_color_classes.reset();
    neighbours.reset();
    color_vertices_considered.reset();

    end_cc = color_class + class_size;
    acc_in = 0;
    while (cc < end_cc) { // increment value of neighbours of vc by 1
        vc = c->lab[cc];
        pe = g->v[vc];
        end_i = pe + g->d[vc];
        for (i = pe; i < end_i; i++) {
            v   = g->e[i];
            col = c->vertex_to_col[v];
            neighbours.inc_nr(v);
            if(neighbours.get(v) == 0) {
                color_vertices_considered.inc_nr(col);
                assert(col + color_vertices_considered.get(col) < g->v_size);
                scratch[col + color_vertices_considered.get(col)] = v; // hit vertices
                if (!color_workset.get(col)) {
                    if (c->ptn[col] == 0)
                        singletons.push_back(col);
                    else
                        old_color_classes.push_back(col);
                    acc_in += (col + 1);
                    color_workset.set(col);
                }
            }
        }
        cc += 1;
    }

    while(!singletons.empty()) {
        singletons.pop_back();
        neighbours.set(c->lab[col], -1);
        color_vertices_considered.set(col, -1);
    }

    // split color classes according to neighbour count
    while(!old_color_classes.empty()) {
        _col    = old_color_classes.pop_back();
        _col_sz = c->ptn[_col] + 1;
        //std::cout << "it " << _col << ", " << _col_sz << std::endl;

        neighbour_sizes.reset();
        degrees_worklist.reset();

        int total = 0;
        for(i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
            v = scratch[_col + i];
            int index = neighbours.get(v) + 1;
            total += 1;
            if(neighbour_sizes.inc(index) == 0)
                degrees_worklist.push_back(index);
        }

        degrees_worklist.sort();

        // enrich neighbour_sizes to accumulative counting array
        acc = 0;
        while(!degrees_worklist.empty()) {
            i   = degrees_worklist.pop_back();
            val = neighbour_sizes.get(i) + 1;
            if(val >= 1) {
                neighbour_sizes.set(i, val + acc);
                acc += val;
                __col = _col + _col_sz - (neighbour_sizes.get(i));
                v_degree = i;
            }
        }

        // determine colors
        color_worklist_color.reset();

        for(i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
            v = scratch[_col + i];
            if (neighbours.get(v) == -1) {
                v_new_color = _col; assert(false);
            } else {
                assert((neighbours.get(v) + 1) > 0);
                v_new_color = _col + _col_sz - neighbour_sizes.get(neighbours.get(v) + 1);
            }
            if (v_new_color != _col) {
                color_worklist_color.push_back(std::pair<int, int>(v, v_new_color));
                c->ptn[v_new_color] = -1;
            }
            neighbours.set(v, -1);
        }

        color_vertices_considered.set(_col, -1);

        // rearrange vertices
        while (!color_worklist_color.empty()) {
            v           = color_worklist_color.last()->first;
            v_new_color = color_worklist_color.last()->second;
            color_worklist_color.pop_back();

            vertex_old_pos = c->vertex_to_lab[v];
            vertex_at_pos  = c->lab[v_new_color + c->ptn[v_new_color] + 1];
            c->lab[vertex_old_pos]          = vertex_at_pos;
            c->vertex_to_lab[vertex_at_pos] = vertex_old_pos;

            c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
            c->vertex_to_col[v] = v_new_color;
            c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
            c->ptn[v_new_color] += 1;

            if (_col != v_new_color) {
                assert(v_new_color > _col);
                c->ptn[_col] -= 1;
            } else assert(false);
        }

        // add new colors to worklist
        largest_color_class_size = -1;
        for(i = _col; i < _col + _col_sz;) {
            assert(i >= 0 && i < c->ptn_sz);
            assert(c->ptn[i] + 1 > 0);
            mark_as_largest = largest_color_class_size < (c->ptn[i] + 1);
            largest_color_class_size = mark_as_largest?(c->ptn[i] + 1):largest_color_class_size;

            color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(std::pair<int, int>(_col, i), mark_as_largest));
            if(i != 0)
                c->ptn[i - 1] = 0;
            i += c->ptn[i] + 1;
        }

        if(i != 0)
            c->ptn[i - 1] = 0;
    }

    neighbour_sizes.reset();
    degrees_worklist.reset();
    color_worklist_color.reset();
    color_workset.reset();

    return comp;
}

bool refinement::refine_color_class_singleton(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist, invariant* I) {
    bool comp;
    int i, cc, deg1_write_pos, deg1_read_pos, singleton_inv;
    cc = color_class; // iterate over color class
    comp = true;

    neighbours.reset();
    color_workset.reset();
    vertex_worklist.reset();
    old_color_classes.reset();

    singleton_inv = 0;
    const int vc = c->lab[cc];
    const int pe = g->v[vc];
    const int end_i = pe + g->d[vc];
    for (i = pe; i < end_i; i++) {
        const int v   = g->e[i];
        const int col = c->vertex_to_col[v];

        if(c->ptn[col] == 0) {
            // if full invariant, sort -- else use a hash value
            if(config.CONFIG_IR_FULL_INVARIANT)
                vertex_worklist.push_back(col); // treat singletons in separate list (more efficient sorting)
            else
                singleton_inv += (col + 1) * (23524361 - col * 3);
            continue;
        }

        if(!color_workset.get(col)) {
            color_workset.set(col);
            old_color_classes.push_back(col);
            // neighbours acts as the write position for degree 1 vertices of col in the scratchpad
            neighbours.set(col, col);
        }
        // write vertex to scratchpad for later use
        scratch[neighbours.get(col)] = v;
        neighbours.inc_nr(col); // we reset neighbours later, since old_color_classes for reset
    }

    //comp = comp && I->write_top_and_compare(INT32_MAX - 3);
    if(!config.CONFIG_IR_FULL_INVARIANT)
        comp = comp && I->write_top_and_compare(g->v_size * 3 + singleton_inv);
    comp = comp && I->write_top_and_compare(g->v_size * 8 + old_color_classes.cur_pos);
    if(!comp) {
        while(!old_color_classes.empty())
            neighbours.set(old_color_classes.pop_back(), -1);
        return comp;
    }

    old_color_classes.sort();

    // write invariant first...
    for(i = 0; i < old_color_classes.cur_pos && comp; ++i) {
        //comp = comp && I->write_top_and_compare(old_color_classes.arr[i]); // color class
        comp = comp && I->write_top_and_compare(g->v_size * 14 + neighbours.get(old_color_classes.arr[i])); // size
        // contains information about color degree (= 1)
    }

    // sort and write down singletons in invariant
    if(comp && config.CONFIG_IR_FULL_INVARIANT)
        vertex_worklist.sort();

    for(i = 0; i < vertex_worklist.cur_pos && comp; ++i) {
        comp = comp && I->write_top_and_compare(g->v_size * 11 + vertex_worklist.arr[i]); // size
        // should contain information about color degree
    }

    if(!comp) {
        while(!old_color_classes.empty())
            neighbours.set(old_color_classes.pop_back(), -1);
        return comp;
    }

    while(!old_color_classes.empty()) {
        const int deg0_col    = old_color_classes.pop_back();
        const int deg1_col_sz = neighbours.get(deg0_col) - deg0_col;
        const int deg0_col_sz = (c->ptn[deg0_col] + 1) - deg1_col_sz;
        const int deg1_col    = deg0_col + deg0_col_sz;

        // no split? done...
        if(deg0_col == deg1_col) {
            neighbours.set(deg1_col, -1);
            continue;
        }

        // set ptn
        c->ptn[deg0_col]     = deg0_col_sz - 1;
        c->ptn[deg1_col]     = deg1_col_sz - 1;
        c->ptn[deg1_col - 1] = 0;

        deg1_write_pos = deg1_col;
        deg1_read_pos  = neighbours.get(deg0_col) - 1;

        // rearrange vertices of deg1 to the back of deg0 color
        while(deg1_read_pos >= deg0_col) {
            const int v             = scratch[deg1_read_pos];
            const int vertex_at_pos = c->lab[deg1_write_pos];
            const int lab_pos       = c->vertex_to_lab[v];

            c->lab[deg1_write_pos]          = v;
            c->vertex_to_lab[v]             = deg1_write_pos;
            c->vertex_to_col[v]             = deg1_col;
            c->lab[lab_pos]                 = vertex_at_pos;
            c->vertex_to_lab[vertex_at_pos] = lab_pos;

            deg1_write_pos++;
            deg1_read_pos--;
        }

        // add new classes to color_class_split_worklist
        const bool leq = deg1_col_sz > deg0_col_sz;
        color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(std::pair<int, int>(deg0_col, deg1_col), leq));
        color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(std::pair<int, int>(deg0_col, deg0_col), !leq));

        // reset neighbours count to -1
        neighbours.set(deg0_col, -1);
    }

    return comp;
}

bool refinement::refine_color_class_dense(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist, invariant* I) {
    bool comp;
    int i, cc, acc, largest_color_class_size, singleton_inv;
    cc = color_class; // iterate over color class
    comp = true;

    neighbours.reset_hard();
    color_workset.reset();
    old_color_classes.reset();
    vertex_worklist.reset();

    singleton_inv = 0;
    const int end_cc = color_class + class_size;

    // for all vertices of the color class...
    while (cc < end_cc) { // increment value of neighbours of vc by 1
        const int vc = c->lab[cc];
        const int pe = g->v[vc];
        const int end_i = pe + g->d[vc];
        for (i = pe; i < end_i; i++) {
            const int v   = g->e[i];
            const int col = c->vertex_to_col[v];
            if(c->ptn[col] > 0) {
                neighbours.inc_nr(v);
                if (!color_workset.get(col)) {
                    color_workset.set_nr(col);
                    old_color_classes.push_back(col);
                }
            } else {
                if(config.CONFIG_IR_FULL_INVARIANT)
                    vertex_worklist.push_back(col);
                else
                    singleton_inv += (col + 1) * (23524361 - col * 3);
            }
        }
        cc += 1;
    }

    color_workset.reset_hard();

    // write singletons

    //comp = comp && I->write_top_and_compare(INT32_MAX - 3);
    comp = comp && I->write_top_and_compare(g->v_size * 3 + old_color_classes.cur_pos);
    if(!comp) {vertex_worklist.reset(); return comp;}

    if(config.CONFIG_IR_FULL_INVARIANT) {
        vertex_worklist.sort();
        while (!vertex_worklist.empty()) {
            const int col = vertex_worklist.pop_back();
            comp = comp && I->write_top_and_compare(g->v_size * 9 + col);
        }
        //comp = comp && I->write_top_and_compare(INT32_MAX - 9);
    } else {
        comp = comp && I->write_top_and_compare(singleton_inv);
    }
    if(!comp) {return comp;}

    old_color_classes.sort();

    // for every cell to be split...
    while(!old_color_classes.empty()) {
        const int col = old_color_classes.pop_back();
        //j     += c->ptn[j] + 1;
        const int col_sz = c->ptn[col] + 1;

        if(col_sz == 1) {
            const int v_degree = neighbours.get(c->lab[col]) + 1;
            if(v_degree == -1) {
                continue;
            }
            //comp = comp && I->write_top_and_compare(INT32_MAX - 2);
            comp = comp && I->write_top_and_compare(- g->v_size * 2 - col);
            comp = comp && I->write_top_and_compare(col + v_degree);
            comp = comp && I->write_top_and_compare(col + c->ptn[col] + 1);
            if(!comp) return comp;
            continue;
        }

        neighbour_sizes.reset();
        degrees_worklist.reset();

        const int first_index = neighbours.get(c->lab[col]) + 1;
        degrees_worklist.push_back(first_index);
        int total = 0;
        for(i = col; i < col + col_sz; ++i) {
            const int v = c->lab[i];
            int index = neighbours.get(v) + 1;
            if(index == first_index)
                continue;
            total += 1;
            if(neighbour_sizes.inc(index) == 0)
                degrees_worklist.push_back(index);
        }

        neighbour_sizes.inc(first_index);
        neighbour_sizes.set(first_index, col_sz - total - 1);

        comp = comp && I->write_top_and_compare(degrees_worklist.cur_pos);
        if(!comp) return comp;

        if(degrees_worklist.cur_pos == 1) {
            // no split
            const int v_degree = neighbours.get(c->lab[col]);
            //comp = comp && I->write_top_and_compare(INT32_MAX - 4);
            comp = comp && I->write_top_and_compare(- g->v_size * 4 - col);
            comp = comp && I->write_top_and_compare(col + v_degree);
            comp = comp && I->write_top_and_compare(col + c->ptn[col] + 1);
            if(!comp) return comp;
            continue;
        }

        degrees_worklist.sort();
        // enrich neighbour_sizes to accumulative counting array
        acc = 0;
        //comp = comp && I->write_top_and_compare(INT32_MAX - 5);
        while(!degrees_worklist.empty()) {
            const int i   = degrees_worklist.pop_back();
            const int val = neighbour_sizes.get(i) + 1;
            if(val >= 1) {
                neighbour_sizes.set(i, val + acc);
                acc += val;
                const int _col = col + col_sz - (neighbour_sizes.get(i)); // use new color here to remove color_workset?
                c->ptn[_col] = -1; // this is val - 1, actually...

                const int v_degree = i;
                comp = comp && I->write_top_and_compare(-g->v_size * 5 - _col);
                comp = comp && I->write_top_and_compare(_col + v_degree);
                comp = comp && I->write_top_and_compare(_col + val + 1);
                if(!comp) return comp;
            }
        }

        // copy cell for rearranging
        memcpy(vertex_worklist.arr, c->lab + col, col_sz * sizeof(int));
        vertex_worklist.cur_pos = col_sz;

        // determine colors and rearrange
        // split color classes according to count in counting array
        while(!vertex_worklist.empty()) {
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
        for(i = col; i < col + col_sz;) {
            assert(i >= 0 && i < c->ptn_sz);
            assert(c->ptn[i] + 1 > 0);
            const int v_color  = i;
            const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
            largest_color_class_size = mark_as_largest?c->ptn[i] + 1:largest_color_class_size;

            color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(std::pair<int, int>(col, i), mark_as_largest));
            if(v_color != 0)
                c->ptn[v_color - 1] = 0;
            i += c->ptn[i] + 1;
        }
    }

    return comp;
}

bool refinement::refine_color_class_dense_dense(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist, invariant* I) {
    bool comp;
    int i, j, acc, cc, largest_color_class_size;
    cc = color_class; // iterate over color class
    comp = true;

    neighbours.reset_hard();
    color_workset.reset();
    old_color_classes.reset();

    const int end_cc = color_class + class_size;
    while (cc < end_cc) { // increment value of neighbours of vc by 1
        const int vc = c->lab[cc];
        const int pe = g->v[vc];
        const int end_i = pe + g->d[vc];
        for (i = pe; i < end_i; i++) {
            const int v   = g->e[i];
            neighbours.inc_nr(v);
        }
        cc += 1;
    }

    // in dense-dense we dont sort, instead we just iterate over all cells (which are already sorted)
    old_color_classes.sort();
    // for every cell...
    for(j = 0; j < g->v_size;) {
        const int col    = j;
        j     += c->ptn[j] + 1;
        const int col_sz = c->ptn[col] + 1;

        if(col_sz == 1) {
            const int v_degree = neighbours.get(c->lab[col]) + 1;
            if(v_degree == -1) // not connected! skip...
                continue;
            //comp = comp && I->write_top_and_compare(INT32_MAX - 2);
            //comp = comp && I->write_top_and_compare(-g->v_size * 11 - col);
            comp = comp && I->write_top_and_compare(col + v_degree * g->v_size);
            //comp = comp && I->write_top_and_compare(g->v_size + col + c->ptn[col] + 1);
            if(!comp) return comp;
            continue;
        }

        neighbour_sizes.reset();
        degrees_worklist.reset();

        const int first_index = neighbours.get(c->lab[col]) + 1;
        degrees_worklist.push_back(first_index);
        int total = 0;
        for(i = col; i < col + col_sz; ++i) {
            const int v = c->lab[i];
            const int index = neighbours.get(v) + 1;
            if(index == first_index)
                continue;
            total += 1;
            if(neighbour_sizes.inc(index) == 0)
                degrees_worklist.push_back(index);
        }

        neighbour_sizes.inc(first_index);
        neighbour_sizes.set(first_index, col_sz - total - 1);

        comp = comp && I->write_top_and_compare(g->v_size * 12 + degrees_worklist.cur_pos);
        if(!comp) return comp;

        if(degrees_worklist.cur_pos == 1) {
            // no split
            const int v_degree = neighbours.get(c->lab[col]);
            //comp = comp && I->write_top_and_compare(INT32_MAX - 4);
            //comp = comp && I->write_top_and_compare(-g->v_size * 10 - col);
            comp = comp && I->write_top_and_compare(col + v_degree * g->v_size);
            comp = comp && I->write_top_and_compare(col + c->ptn[col] + 1);
            if(!comp) return comp;
            continue;
        }

        degrees_worklist.sort();
        // enrich neighbour_sizes to accumulative counting array
        acc = 0;
        //comp = comp && I->write_top_and_compare(INT32_MAX - 5);
        while(!degrees_worklist.empty()) {
            const int i   = degrees_worklist.pop_back();
            const int val = neighbour_sizes.get(i) + 1;
            if(val >= 1) {
                neighbour_sizes.set(i, val + acc);
                acc += val;
                const int _col = col + col_sz - (neighbour_sizes.get(i)); // use new color here to remove color_workset?
                c->ptn[_col] = -1; // this is val - 1, actually...

                const int v_degree = i;
                comp = comp && I->write_top_and_compare(-g->v_size * 5 - _col);
                comp = comp && I->write_top_and_compare(_col + v_degree);
                comp = comp && I->write_top_and_compare(_col + val + 1);
                if(!comp) return comp;
            }
        }

        // copy cell for rearranging
        memcpy(vertex_worklist.arr, c->lab + col, col_sz * sizeof(int));
        vertex_worklist.cur_pos = col_sz;

        // determine colors and rearrange
        // split color classes according to count in counting array
        while(!vertex_worklist.empty()) {
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
        for(i = col; i < col + col_sz;) {
            assert(i >= 0 && i < c->ptn_sz);
            assert(c->ptn[i] + 1 > 0);
            const int v_color  = i;
            const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
            largest_color_class_size = mark_as_largest?c->ptn[i] + 1:largest_color_class_size;

            color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(std::pair<int, int>(col, i), mark_as_largest));
            if(v_color != 0)
                c->ptn[v_color - 1] = 0;
            i += c->ptn[i] + 1;
        }
    }

    return comp;
}

bool refinement::refine_color_class_dense_dense127(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist, invariant* I) {
    bool comp;
    int i, j, acc, cc, largest_color_class_size;
    cc = color_class; // iterate over color class
    comp = true;

    neighbours127.reset_hard();
    color_workset.reset();
    old_color_classes.reset();

    const int end_cc = color_class + class_size;
    while (cc < end_cc) { // increment value of neighbours of vc by 1
        const int vc = c->lab[cc];
        const int pe = g->v[vc];
        const int end_i = pe + g->d[vc];
        for (i = pe; i < end_i; i++) {
            const int v   = g->e[i];
            neighbours127.inc_nr(v);
        }
        cc += 1;
    }

    // in dense-dense we dont sort, instead we just iterate over all cells (which are already sorted)
    old_color_classes.sort();
    // for every cell...
    for(j = 0; j < g->v_size;) {
        const int col    = j;
        j     += c->ptn[j] + 1;
        const int col_sz = c->ptn[col] + 1;

        if(col_sz == 1) {
            const int v_degree = neighbours127.get(c->lab[col]) + 1;
            if(v_degree == -1) // not connected! skip...
                continue;
            //comp = comp && I->write_top_and_compare(INT32_MAX - 2);
            //comp = comp && I->write_top_and_compare(-g->v_size * 11 - col);
            comp = comp && I->write_top_and_compare(col + v_degree * g->v_size);
            //comp = comp && I->write_top_and_compare(g->v_size + col + c->ptn[col] + 1);
            if(!comp) return comp;
            continue;
        }

        neighbour_sizes.reset();
        degrees_worklist.reset();

        const int first_index = neighbours127.get(c->lab[col]) + 1;
        degrees_worklist.push_back(first_index);
        int total = 0;
        for(i = col; i < col + col_sz; ++i) {
            const int v = c->lab[i];
            const int index = neighbours127.get(v) + 1;
            if(index == first_index)
                continue;
            total += 1;
            if(neighbour_sizes.inc(index) == 0)
                degrees_worklist.push_back(index);
        }

        neighbour_sizes.inc(first_index);
        neighbour_sizes.set(first_index, col_sz - total - 1);

        comp = comp && I->write_top_and_compare(g->v_size * 12 + degrees_worklist.cur_pos);
        if(!comp) return comp;

        if(degrees_worklist.cur_pos == 1) {
            // no split
            const int v_degree = neighbours127.get(c->lab[col]);
            //comp = comp && I->write_top_and_compare(INT32_MAX - 4);
            //comp = comp && I->write_top_and_compare(-g->v_size * 10 - col);
            comp = comp && I->write_top_and_compare(col + v_degree * g->v_size);
            comp = comp && I->write_top_and_compare(col + c->ptn[col] + 1);
            if(!comp) return comp;
            continue;
        }

        degrees_worklist.sort();
        // enrich neighbour_sizes to accumulative counting array
        acc = 0;
        //comp = comp && I->write_top_and_compare(INT32_MAX - 5);
        while(!degrees_worklist.empty()) {
            const int i   = degrees_worklist.pop_back();
            const int val = neighbour_sizes.get(i) + 1;
            if(val >= 1) {
                neighbour_sizes.set(i, val + acc);
                acc += val;
                const int _col = col + col_sz - (neighbour_sizes.get(i)); // use new color here to remove color_workset?
                c->ptn[_col] = -1; // this is val - 1, actually...

                const int v_degree = i;
                comp = comp && I->write_top_and_compare(-g->v_size * 5 - _col);
                comp = comp && I->write_top_and_compare(_col + v_degree);
                comp = comp && I->write_top_and_compare(_col + val + 1);
                if(!comp) return comp;
            }
        }

        // copy cell for rearranging
        memcpy(vertex_worklist.arr, c->lab + col, col_sz * sizeof(int));
        vertex_worklist.cur_pos = col_sz;

        // determine colors and rearrange
        // split color classes according to count in counting array
        while(!vertex_worklist.empty()) {
            const int v = vertex_worklist.pop_back();
            const int v_new_color = col + col_sz - (neighbour_sizes.get(neighbours127.get(v) + 1));
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
        for(i = col; i < col + col_sz;) {
            assert(i >= 0 && i < c->ptn_sz);
            assert(c->ptn[i] + 1 > 0);
            const int v_color  = i;
            const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
            largest_color_class_size = mark_as_largest?c->ptn[i] + 1:largest_color_class_size;

            color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(std::pair<int, int>(col, i), mark_as_largest));
            if(v_color != 0)
                c->ptn[v_color - 1] = 0;
            i += c->ptn[i] + 1;
        }
    }

    return comp;
}

bool refinement::refine_color_class_singleton_first(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist) {
    // for all vertices of the color class...
    bool comp, mark_as_largest;
    int i, cc, vc, pe, end_i, v, col, deg1_col, deg1_col_sz, deg0_col, deg0_col_sz, deg1_write_pos, deg1_read_pos,
            vertex_at_pos, lab_pos;
    cc = color_class; // iterate over color class
    comp = true;

    neighbours.reset();
    color_workset.reset();
    vertex_worklist.reset();
    old_color_classes.reset();

    vc = c->lab[cc];
    pe = g->v[vc];
    end_i = pe + g->d[vc];
    for (i = pe; i < end_i; i++) {
        v   = g->e[i];
        col = c->vertex_to_col[v];

        if(c->ptn[col] == 0) {
            //vertex_worklist.push_back(col); // treat singletons in separate list
            continue;
        }

        if(!color_workset.get(col)) {
            color_workset.set(col);
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

    //old_color_classes.sort();

    // sort and write down singletons in invariant

    while(!old_color_classes.empty()) {
        deg0_col    = old_color_classes.pop_back();
        deg1_col_sz = neighbours.get(deg0_col) - deg0_col;
        deg0_col_sz = (c->ptn[deg0_col] + 1) - deg1_col_sz;
        deg1_col    = deg0_col + deg0_col_sz;

        // no split? done...
        if(deg0_col == deg1_col) {
            neighbours.set(deg1_col, -1);
            continue;
        }

        // set ptn
        c->ptn[deg0_col]     = deg0_col_sz - 1;
        c->ptn[deg1_col]     = deg1_col_sz - 1;
        c->ptn[deg1_col - 1] = 0;

        deg1_write_pos = deg1_col;
        deg1_read_pos  = neighbours.get(deg0_col) - 1;

        // rearrange vertices of deg1 to the back of deg0 color
        while(deg1_read_pos >= deg0_col) {
            v             = scratch[deg1_read_pos];
            vertex_at_pos = c->lab[deg1_write_pos];
            lab_pos       = c->vertex_to_lab[v];

            c->lab[deg1_write_pos]          = v;
            c->vertex_to_lab[v]             = deg1_write_pos;
            c->vertex_to_col[v]             = deg1_col;
            c->lab[lab_pos]                 = vertex_at_pos;
            c->vertex_to_lab[vertex_at_pos] = lab_pos;

            deg1_write_pos++;
            deg1_read_pos--;
        }

        // add new classes to color_class_split_worklist
        bool leq = deg1_col_sz > deg0_col_sz;
        color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(std::pair<int, int>(deg0_col, deg1_col), leq));
        color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(std::pair<int, int>(deg0_col, deg0_col), !leq));

        // reset neighbours count to -1
        neighbours.set(deg0_col, -1);
    }

    return comp;
}

bool refinement::refine_color_class_dense_first(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist) {
    // for all vertices of the color class...
    bool mark_as_largest;
    int i, cc, vc, pe, end_i, end_cc, v, col, _col, acc, val, col_sz, v_new_color, v_color, v_degree, largest_color_class_size;
    cc = color_class; // iterate over color class

    neighbours.reset_hard();
    color_workset.reset();
    old_color_classes.reset();

    end_cc = color_class + class_size;
    while (cc < end_cc) { // increment value of neighbours of vc by 1
        vc = c->lab[cc];
        pe = g->v[vc];
        end_i = pe + g->d[vc];
        for (i = pe; i < end_i; i++) {
            v   = g->e[i];
            neighbours.inc_nr(v);
            col = c->vertex_to_col[v];
            if(!color_workset.get(col)) {
                color_workset.set_nr(col);
                old_color_classes.push_back(col);
            }
        }
        cc += 1;
    }

    color_workset.reset_hard();

    // DENSE-DENSE: dont sort, just iterate over all cells
    while(!old_color_classes.empty()) {
        col = old_color_classes.pop_back();
        //j     += c->ptn[j] + 1;
        col_sz = c->ptn[col] + 1;

        if(col_sz == 1) {
            v_degree = neighbours.get(c->lab[col]) + 1;
            if(v_degree == -1)
                continue;
        }

        neighbour_sizes.reset();
        degrees_worklist.reset();

        int first_index = neighbours.get(c->lab[col]) + 1;
        degrees_worklist.push_back(first_index);
        int total = 0;
        for(i = col; i < col + col_sz; ++i) {
            v = c->lab[i];
            int index = neighbours.get(v) + 1;
            if(index == first_index)
                continue;
            total += 1;
            if(neighbour_sizes.inc(index) == 0)
                degrees_worklist.push_back(index);
        }

        neighbour_sizes.inc(first_index);
        neighbour_sizes.set(first_index, col_sz - total - 1);

        if(degrees_worklist.cur_pos == 1) {
            continue;
        }

        degrees_worklist.sort();
        // enrich neighbour_sizes to accumulative counting array
        acc = 0;
        while(!degrees_worklist.empty()) {
            i   = degrees_worklist.pop_back();
            val = neighbour_sizes.get(i) + 1;
            if(val >= 1) {
                neighbour_sizes.set(i, val + acc);
                acc += val;
                _col = col + col_sz - (neighbour_sizes.get(i)); // use new color here to remove color_workset?
                c->ptn[_col] = -1; // this is val - 1, actually...

                v_degree = i;
            }
        }

        // copy cell for rearranging
        memcpy(vertex_worklist.arr, c->lab + col, col_sz * sizeof(int));
        vertex_worklist.cur_pos = col_sz;

        // determine colors and rearrange
        // split color classes according to count in counting array
        while(!vertex_worklist.empty()) {
            v = vertex_worklist.pop_back();
            v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
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
        for(i = col; i < col + col_sz;) {
            assert(i >= 0 && i < c->ptn_sz);
            assert(c->ptn[i] + 1 > 0);
            v_color  = i;
            mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
            largest_color_class_size = mark_as_largest?c->ptn[i] + 1:largest_color_class_size;

            color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(std::pair<int, int>(col, i), mark_as_largest));
            if(v_color != 0)
                c->ptn[v_color - 1] = 0;
            i += c->ptn[i] + 1;
        }
    }

    return true;
}

bool refinement::refine_color_class_dense_dense_first(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist) {
    // for all vertices of the color class...
    bool mark_as_largest;
    int i, j, cc, vc, pe, end_i, end_cc, v, col, _col, acc, val, col_sz, v_new_color, v_color, v_degree, largest_color_class_size;
    cc = color_class; // iterate over color class

    neighbours.reset_hard();
    color_workset.reset();
    old_color_classes.reset();

    end_cc = color_class + class_size;
    while (cc < end_cc) { // increment value of neighbours of vc by 1
        vc = c->lab[cc];
        pe = g->v[vc];
        end_i = pe + g->d[vc];
        for (i = pe; i < end_i; i++) {
            v   = g->e[i];
            neighbours.inc_nr(v);
        }
        cc += 1;
    }

    // dont sort, just iterate over all cells
    for(j = 0; j < g->v_size;) {
        col    = j;
        j     += c->ptn[j] + 1;
        col_sz = c->ptn[col] + 1;

        if(col_sz == 1)
            continue;

        neighbour_sizes.reset();
        degrees_worklist.reset();

        int first_index = neighbours.get(c->lab[col]) + 1;
        degrees_worklist.push_back(first_index);
        int total = 0;
        for(i = col; i < col + col_sz; ++i) {
            v = c->lab[i];
            int index = neighbours.get(v) + 1;
            if(index == first_index)
                continue;
            total += 1;
            if(neighbour_sizes.inc(index) == 0)
                degrees_worklist.push_back(index);
        }

        neighbour_sizes.inc(first_index);
        neighbour_sizes.set(first_index, col_sz - total - 1);

        if(degrees_worklist.cur_pos == 1) {
            continue;
        }

        degrees_worklist.sort();
        // enrich neighbour_sizes to accumulative counting array
        acc = 0;
        while(!degrees_worklist.empty()) {
            i   = degrees_worklist.pop_back();
            val = neighbour_sizes.get(i) + 1;
            if(val >= 1) {
                neighbour_sizes.set(i, val + acc);
                acc += val;
                _col = col + col_sz - (neighbour_sizes.get(i)); // use new color here to remove color_workset?
                c->ptn[_col] = -1; // this is val - 1, actually...
            }
        }

        // copy cell for rearranging
        memcpy(vertex_worklist.arr, c->lab + col, col_sz * sizeof(int));
        vertex_worklist.cur_pos = col_sz;

        // determine colors and rearrange
        // split color classes according to count in counting array
        while(!vertex_worklist.empty()) {
            v = vertex_worklist.pop_back();
            v_new_color = col + col_sz - (neighbour_sizes.get(neighbours.get(v) + 1));
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
        for(i = col; i < col + col_sz;) {
            assert(i >= 0 && i < c->ptn_sz);
            assert(c->ptn[i] + 1 > 0);
            v_color  = i;
            mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
            largest_color_class_size = mark_as_largest?c->ptn[i] + 1:largest_color_class_size;

            color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(std::pair<int, int>(col, i), mark_as_largest));
            if(v_color != 0)
                c->ptn[v_color - 1] = 0;
            i += c->ptn[i] + 1;
        }
    }

    return true;
}

int refinement::individualize_vertex(coloring *c, int v) {
    //assert(initialized);

    int color = c->vertex_to_col[v];
    int pos   = c->vertex_to_lab[v];
    //int pos   = color;
    //while(c->lab[pos] != v) pos += 1;

    int color_class_size = c->ptn[color];

    c->color_choices.emplace_back(std::pair<int, int>(color, color_class_size + 1));

    assert(color_class_size > 0);
    int new_color = color + color_class_size;

    int vertex_at_pos = c->lab[color + color_class_size];
    c->lab[pos] = vertex_at_pos;
    c->vertex_to_lab[vertex_at_pos] = pos;

    c->lab[color + color_class_size] = v;
    c->vertex_to_lab[v] = color + color_class_size;
    c->vertex_to_col[v] = color + color_class_size;

    c->ptn[color] -= 1;
    c->ptn[color + color_class_size - 1] = 0;
    return color + color_class_size;
}


bool refinement::assert_is_equitable(sgraph *g, coloring *c) {
    std::vector<int> neighbour_col_canon;
    std::vector<int> neighbour_col;
    bool first_of_color = true;
    bool test = true;

    for(int i = 0; i < g->v_size; ++i) {
        neighbour_col.clear();
        int v = c->lab[i];
        for(int j = g->v[v]; j < g->v[v] + g->d[v]; ++j) {
            neighbour_col.push_back(c->vertex_to_col[g->e[j]]);
        }
        std::sort(neighbour_col.begin(), neighbour_col.end());
        if(first_of_color) {
            first_of_color = false;
            neighbour_col_canon.swap(neighbour_col);
        } else {
            for (auto j = 0; j < neighbour_col.size(); ++j) {
                test = test && (neighbour_col[j] == neighbour_col_canon[j]);
            }
        }

        if(c->ptn[i] == 0) {
            first_of_color = true;
        }
    }
    int expect0 = 0;
    bool expectsize = true;
    //std::cout << "part";
    for(int i = 0; i < c->lab_sz; ++i) {
        if(expectsize) {
            expectsize = false;
            expect0 = c->ptn[i];
        }
        if(expect0 == 0) {
            //std::cout << i << " ";
            test = test && (c->ptn[i] == 0);
            expectsize = true;
        } else {
            test = test && (c->ptn[i] > 0);
        }
        expect0 -= 1;
    }
    return test;
}

refinement::~refinement() {
    //if(initialized)
       // delete[] p;
}

void refinement::undo_changes_with_reference(change_tracker *changes, coloring *changed_c, coloring *old_c) {
    color_workset.reset();
    while(!changes->empty()) {
        int v, i;
        int c      = changes->pop();
        int c_size = old_c->ptn[c] + 1;
        if(c_size > 1 && !color_workset.get(c)) {
            color_workset.set(c);
            int changed_size = changed_c->ptn[c];
            changed_c->ptn[c] = old_c->ptn[c];
            memcpy(changed_c->ptn + c + changed_size, old_c->ptn + c + changed_size, c_size*sizeof(int));
            for(i = changed_size; i < c_size; ++i) {
                   v = changed_c->lab[c + i];
                   changed_c->vertex_to_col[v] = c;
            }
        }
    }
    color_workset.reset();
}

void work_set::initialize(int size) {
    s = new bool[size];
    //s.reserve(size);

    //for(int i = 0; i < size; ++i)
    //    s.push_back(false);
    memset(s, false, size * sizeof(bool));

    reset_queue.initialize(size);
    init = true;
    sz = size;
}

void work_set::set(int index) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    if(!s[index]) {
        s[index] = true;
        reset_queue.push(index);
    }
}

void work_set::unset(int index) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    s[index] = false;
}

void work_set::set_nr(int index) {
    s[index] = true;
}

bool work_set::get(int index) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    return s[index];
}

void work_set::reset() {
    while(!reset_queue.empty()) {
        int index = reset_queue.pop();
        assert(init);
        assert(index >= 0);
        assert(index < sz);
        s[index] = false;
    }
}

void work_set::reset_hard() {
    reset_queue.pos = 0;
    //s.
    // std::fill(s.begin(), s.end(), false);
    memset(s, false, sz*sizeof(bool));
}

work_set::~work_set() {
    if(init)
        delete[] s;
}

void work_set::reset_soft() {
    reset_queue.reset();
}

void work_set_int::initialize(int size) {
    s = new int[size];
    reset_queue.initialize(size);

    memset(s, -1, size * sizeof(int));

    init = true;
    sz = size;
}

void work_set_int::set(int index, int value) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    s[index] = value;
}

int work_set_int::inc(int index) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    if(s[index]++ == -1)
        reset_queue.push(index);
    return s[index];
}

int work_set_int::get(int index) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    return s[index];
}

void work_set_int::reset() {
    while(!reset_queue.empty())
        s[reset_queue.pop()] = -1;
}

work_set_int::~work_set_int() {
    if(init)
        delete[] s;
}

void work_set_int::inc_nr(int index) {
    assert(index >= 0 && index < sz);
    ++s[index];
}

void work_set_int::reset_hard() {
    memset(s, -1, sz*sizeof(int));
    reset_queue.pos = 0;
}

void work_set_char::initialize(int size) {
    s = new char[size];
    reset_queue.initialize(size);

    memset(s, -1, size * sizeof(char));

    init = true;
    sz = size;
}

void work_set_char::set(int index, char value) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    s[index] = value;
}

int work_set_char::inc(int index) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    if(s[index]++ == -1)
        reset_queue.push(index);
    return s[index];
}

char work_set_char::get(int index) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    return s[index];
}

void work_set_char::reset() {
    while(!reset_queue.empty())
        s[reset_queue.pop()] = -1;
}

work_set_char::~work_set_char() {
    if(init)
        delete[] s;
}

void work_set_char::inc_nr(int index) {
    assert(index >= 0 && index < sz);
    ++s[index];
}

void work_set_char::reset_hard() {
    memset(s, -1, sz*sizeof(char));
    reset_queue.pos = 0;
}

void work_list::initialize(int size) {
    arr     = new int[size];
    arr_sz  = size;
    cur_pos = 0;
}

void work_list::push_back(int value) {
    assert(cur_pos >= 0 && cur_pos < arr_sz);
    arr[cur_pos] = value;
    cur_pos += 1;
}

int work_list::pop_back() {
    //assert(cur_pos >= 0 && cur_pos < arr_sz);
    return arr[--cur_pos];
}

void work_list::reset() {
    cur_pos = 0;
}

bool work_list::empty() {
    return cur_pos == 0;
}

inline void sort_int(int* arr, int sz) {
#define min(x, y) (x<y?x:y)
#define max(x, y) (x<y?y:x)
#define SWAP(x,y) { const int a = min(arr[x], arr[y]); \
                    const int b = max(arr[x], arr[y]); \
                    arr[x] = a; arr[y] = b; }
    switch(sz) {
        case 0:
        case 1:
            break;
        case 2:
            SWAP(0, 1);
            break;
        case 3:
            SWAP(0, 1);
            SWAP(0, 2);
            SWAP(1, 2);
            break;
        case 4:
            SWAP(0, 1);
            SWAP(2, 3);
            SWAP(0, 2);
            SWAP(1, 3);
            SWAP(1,2);
            break;
        case 5:
            SWAP(1, 2);
            SWAP(0, 2);
            SWAP(0, 1);
            SWAP(3, 4);
            SWAP(1, 4);
            SWAP(0, 3);
            SWAP(1, 3);
            SWAP(2, 4);
            SWAP(2, 3);
            break;
        case 6:
            SWAP(1, 2);
            SWAP(4, 5);
            SWAP(0, 2);
            SWAP(3, 5);
            SWAP(0, 1);
            SWAP(3, 4);
            SWAP(1, 4);
            SWAP(0, 3);
            SWAP(2, 5);
            SWAP(1, 3);
            SWAP(2, 4);
            SWAP(2, 3);
            break;
        default:
            std::sort(arr, arr + sz);
    }
#undef SWAP
#undef min
#undef max
}

void work_list::sort() {
        sort_int(arr, cur_pos);
        //std::sort(arr, arr + cur_pos);
}

work_list::~work_list() {
    if(arr_sz >= 0) {
        delete[] arr;
    }
}

void work_queue::initialize(int size) {
    assert(!init);
    sz = size;
    pos = 0;
    queue = new int[size];
    init = true;
}

void work_queue::push(int val) {
    assert(init);
    assert(pos != sz);
    queue[pos] = val;
    pos++;
}

int work_queue::pop() {
    assert(init);
    assert(pos > 0);
    pos--;
    return queue[pos];
}

bool work_queue::empty() {
    assert(init);
    return (pos == 0);
}

work_queue::~work_queue() {
    if(init) {
        delete[] queue;
    }
}

void work_queue::reset() {
    pos = 0;
}

void work_list_pair::push_back(std::pair<int, int> value) {
    arr[arr_sz] = value;
    arr_sz += 1;
}

std::pair<int, int> *work_list_pair::last() {
    return &arr[arr_sz - 1];
}

void work_list_pair::initialize(int size) {
    init = true;
    arr = new std::pair<int, int>[size];
    arr_sz = 0;
}

void work_list_pair::initialize_from_array(std::pair<int, int>* p, int size) {
    arr = p;
    arr_sz = 0;
}

work_list_pair::~work_list_pair() {
    if(init)
        delete[] arr;
}

void work_list_pair::pop_back() {
    arr_sz -= 1;
}

void work_list_pair::sort() {
    if(arr_sz > 1) {
        std::sort(arr, arr + arr_sz);
    }
}

bool work_list_pair::empty() {
    return arr_sz == 0;
}

void work_list_pair::reset() {
    arr_sz = 0;
}

void work_list_pair_bool::push_back(std::pair<std::pair<int, int>, bool> value) {
    assert(arr_sz >= 0);
    arr[arr_sz] = value;
    arr_sz += 1;
}

std::pair<std::pair<int, int>, bool>* work_list_pair_bool::last() {
    return &arr[arr_sz - 1];
}

void work_list_pair_bool::initialize(int size) {
    init = true;
    arr = new std::pair<std::pair<int, int>, bool>[size];
    arr_sz = 0;
}

work_list_pair_bool::~work_list_pair_bool() {
    if(init)
        delete[] arr;
}

void work_list_pair_bool::pop_back() {
    arr_sz -= 1;
}

void work_list_pair_bool::sort() {
    if(arr_sz > 1) {
        std::sort(arr, arr + arr_sz);
    }
}

bool work_list_pair_bool::empty() {
    return arr_sz == 0;
}

void work_list_pair_bool::reset() {
    arr_sz = 0;
}

void ring_pair::initialize(int size) {
    arr = new std::pair<int, int>[size];
    arr_sz = size;
    back_pos  = 0;
    front_pos = 0;
    init = true;
}

void ring_pair::initialize_from_array(std::pair<int, int>* p, int size) {
    arr = p;
    arr_sz    = size;
    back_pos  = 0;
    front_pos = 0;
}

void ring_pair::push_back(std::pair<int, int> value) {
    arr[back_pos] = value;
    back_pos = (back_pos + 1) % arr_sz;
}

std::pair<int, int> *ring_pair::front() {
    return &arr[front_pos];
}

void ring_pair::pop() {
    front_pos = (front_pos + 1) % arr_sz;
}

bool ring_pair::empty() {
    return (front_pos == back_pos);
}

ring_pair::~ring_pair() {
    if(init)
        delete[] arr;
}

void ring_pair::reset() {
    front_pos = back_pos;
}

void change_tracker::initialize(int _limit, int domain_size) {
    l.initialize(_limit);
    s.initialize(domain_size);
    limit = _limit;
    sz = 0;
}

void change_tracker::reset() {
    l.reset();
    s.reset();
    overflow = false;
    sz       = 0;
}

void change_tracker::track(int oldcolor) {
    if(overflow) return;
    bool check = s.get(oldcolor);
    if(!check && sz < limit - 1) {
        s.set(oldcolor);
        l.push_back(oldcolor);
        sz += 1;
    } else if (!check) {
        overflow = true;
    }
}

int change_tracker::pop() {
    return l.pop_back();
}

bool change_tracker::empty() {
    return l.empty();
}

bool change_tracker::did_overflow() {
    return overflow;
}

void cell_worklist::initialize(int domain_size) {
    arr = new std::pair<int, int>[domain_size];
    arr_sz = domain_size;
    init = true;
    cur_pos = 0;
}

int cell_worklist::add_cell(work_set_int *queue_pointer, int col, int col_sz) {
    assert(init);
    assert(cur_pos >= 0 && cur_pos < arr_sz - 1);
    queue_pointer->set(col, cur_pos);
    arr[cur_pos] = std::pair<int, int>(col, col_sz);
    cur_pos++;
    return 0;
}

std::pair<int, int> cell_worklist::next_cell(work_set_int *queue_pointer) {
    // look at first 12 positions and scramble if possible
    int sm_j = cur_pos - 1;
    for(int j = cur_pos - 1; j >= 0 && ((cur_pos - j) <= 14); --j) {
        if(arr[j].second < arr[sm_j].second)
            sm_j = j;
        if(arr[sm_j].second == 1)
            break;
    }

    // swap sm_j and j
    std::pair<int, int> sm_col = arr[sm_j];
    arr[sm_j] = arr[cur_pos - 1];
    queue_pointer->set(arr[sm_j].first, sm_j);

    cur_pos--;
    queue_pointer->set(sm_col.first, -1);
    return sm_col;
}

void cell_worklist::replace_cell(work_set_int *queue_pointer, int col_old, int col, int col_sz) {
    int pos = queue_pointer->get(col_old);
    arr[pos].first  = col;
    arr[pos].second = col_sz;
    assert(queue_pointer->get(col_old) != -1);
    queue_pointer->set(col_old, -1);
    queue_pointer->set(col, pos);
}

bool cell_worklist::empty() {
    return (cur_pos == 0);
}

void cell_worklist::reset(work_set_int *queue_pointer) {
    while(cur_pos > 0) {
        cur_pos--;
        queue_pointer->set(arr[cur_pos].first, -1);
    }
}

