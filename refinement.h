#ifndef DEJAVU_REFINEMENT_H
#define DEJAVU_REFINEMENT_H

#include "coloring.h"
#include "sgraph.h"
#include "invariant.h"
#include "configuration.h"
#include "utility.h"
#include <list>
#include <iostream>
#include <cstring>

class mark_set {
    int mark = 0;
    int *s;
    int sz;
    bool init = false;
public:
    void initialize(int size) {
        s = new int[size];
        sz = size;
        init = true;
        memset(s, mark, sz * sizeof(int));
        reset();
    }
    void initialize_from_array(int* arr, int size) {
        s  = arr;
        sz = size;
        init = false;
        memset(s, mark, sz * sizeof(int));
        reset();
    }
    bool get(int pos) {
        return s[pos] == mark;
    }
    void set(int pos) {
        s[pos] = mark;
    }
    void unset(int pos) {
        s[pos] = mark - 1;
    }
    void reset() {
        if(mark == -1) {
            memset(s, mark, sz * sizeof(int));
        }
        ++mark;
    }
    ~mark_set() {
        if(init)
            delete[] s;
    }
};

template<class T>
void sort_temp(T* arr, int sz) {
#define min(x, y) (x<y?x:y)
#define max(x, y) (x<y?y:x)
#define SWAP(x,y) { const T a = min(arr[x], arr[y]); \
                    const T b = max(arr[x], arr[y]); \
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
            SWAP(0, 1);
            SWAP(2, 3);
            SWAP(0, 2);
            SWAP(1, 4);
            SWAP(0, 1);
            SWAP(2, 3);
            SWAP(1, 2);
            SWAP(3, 4);
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


template<class T>
class work_list_temp {
public:
    void initialize(int size) {
        arr     = new T[size];
        arr_sz  = size;
        init     = true;
        cur_pos = 0;
    }

    void initialize_from_array(T* arr, int size) {
        this->arr = arr;
        arr_sz    = size;
        init      = false;
        cur_pos   = 0;
    }

    void push_back(T value) {
        assert(cur_pos >= 0 && cur_pos < arr_sz);
        arr[cur_pos] = value;
        cur_pos += 1;
    }

    T pop_back() {
        return arr[--cur_pos];
    }

    T* last() {
        return &arr[cur_pos - 1];
    }

    bool empty() {
        return cur_pos == 0;
    }

    void reset() {
        cur_pos = 0;
    }

    ~work_list_temp() {
        if(init) {
            delete[] arr;
        }
    }

    void sort() {
        sort_temp<T>(arr, cur_pos);
    }

    int  cur_pos = 0;
    T*   arr  = nullptr;
    bool init = false;
private:
    int arr_sz = -1;
};

typedef work_list_temp<std::pair<std::pair<int, int>, bool>> work_list_pair_bool;
typedef work_list_temp<int> work_list;

class work_queue {
public:
    void initialize(int size);
    void initialize_from_array(int* arr, int size);
    void push(int val);
    int  pop();
    void reset();
    bool empty();
    ~work_queue();
    int* queue;
    int  pos;
private:
    bool init = false;
    int  sz;
};

class work_set {
public:
    void initialize(int size);
    void set(int index);
    bool get(int index);
    void reset();
    ~work_set();
    void reset_hard();
    void set_nr(int index);
    void unset(int index);
    void reset_soft();

private:
    work_queue reset_queue;
    bool init = false;
    bool* s;
    int sz;
};

template<class T>
class alignas(16) work_set_temp {
public:
    void initialize(int size) {
        s = new T[size];
        reset_queue.initialize(size);

        memset(s, -1, size * sizeof(T));

        init = true;
        sz = size;
    }

    void initialize_from_array(T* arr, int size) {
        s = arr;
        reset_queue.initialize_from_array(arr + size, size);

        memset(s, -1, size * sizeof(T));

        init = false;
        sz = size;
    }

    void set(int index, T value) {
        assert(index >= 0);
        assert(index < sz);
        s[index] = value;
    }

    T   get(int index) {
        assert(index >= 0);
        assert(index < sz);
        return s[index];
    }

    void reset() {
        while(!reset_queue.empty())
            s[reset_queue.pop()] = -1;
    }

    void reset_hard() {
        memset(s, -1, sz*sizeof(T));
        reset_queue.pos = 0;
    }

    T   inc(int index) {
        assert(index >= 0);
        assert(index < sz);
        if(s[index]++ == -1)
            reset_queue.push(index);
        return s[index];
    }

    void inc_nr(int index) {
        assert(index >= 0 && index < sz);
        ++s[index];
    }

    ~work_set_temp() {
        if(init)
            delete[] s;
    }

    work_queue reset_queue;
private:
    bool  init = false;
    T*    s;
    int   sz;
};

typedef work_set_temp<int>  work_set_int;
typedef work_set_temp<char> work_set_char;

class ring_pair {
public:
    void initialize(int size);
    void push_back(std::pair<int, int> value);
    std::pair<int, int>* front();
    void pop();
    bool empty();
    void reset();
    ~ring_pair();
    std::pair<int, int>* arr;
    bool init = false;
    int arr_sz = -1;
    int front_pos = -1;
    int back_pos = -1;

    void initialize_from_array(std::pair<int, int> *p, int size);
};

class change_tracker {
    work_list l;
    work_set   s;
    int    limit   = -1;
    int       sz   = -1;
    bool  overflow = false;
public:
    void initialize(int limit, int domain_size);
    void reset();
    void track(int oldcolor);
    int  pop();
    bool empty();
    bool did_overflow();
};

template<class vertex_type>
class cell_worklist_temp {
public:
    void initialize(int domain_size) {
        arr = new int[domain_size];
        arr_sz = domain_size;
        init = true;
        cur_pos = 0;
    }

    int add_cell(work_set_int* queue_pointer, int col) {
        assert(init);
        assert(cur_pos >= 0 && cur_pos < arr_sz - 1);
        queue_pointer->set(col, cur_pos);
        arr[cur_pos] = col;
        cur_pos++;
        return 0;
    }

    int next_cell(work_set_int* queue_pointer, coloring_temp<vertex_type>* c) {
        // look at first 12 positions and scramble if possible
        int sm_j = cur_pos - 1;
        for(int j = cur_pos - 1; j >= 0 && ((cur_pos - j) <= 12); --j) {
            if(c->ptn[arr[j]] < c->ptn[arr[sm_j]]) {
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

    void replace_cell(work_set_int* queue_pointer, int col_old, int col) {
        const int pos = queue_pointer->get(col_old);
        arr[pos] = col;
        assert(queue_pointer->get(col_old) != -1);
        queue_pointer->set(col_old, -1);
        queue_pointer->set(col, pos);
    }

    void reset(work_set_int* queue_pointer) {
        while(cur_pos > 0) {
            cur_pos--;
            queue_pointer->set(arr[cur_pos], -1);
        }
    }

    bool empty() {
        return (cur_pos == 0);
    }

private:
    int* arr = nullptr;
    int  arr_sz  = -1;
    int  cur_pos = -1;
    bool init    = false;
};

template<class vertex_type, class degree_type, class edge_type>
class refinement_temp {
public:
    bool refine_coloring(sgraph_temp<vertex_type, degree_type, edge_type> *g, coloring_temp<vertex_type> *c,
                         invariant *I, int init_color_class, strategy_metrics *m, int cell_early) {
        bool comp = true;
        assure_initialized(g);

        cell_todo.reset(&queue_pointer);

        if(init_color_class < 0) {
            // initialize queue with all classes (except for largest one)
            for (int i = 0; i < c->ptn_sz;) {
                cell_todo.add_cell(&queue_pointer, i);
                i += c->ptn[i] + 1;
            }
        } else {
            cell_todo.add_cell(&queue_pointer, init_color_class);
        }
        int its = 0;

        while(!cell_todo.empty()) {
            its += 1;
            color_class_splits.reset();
            const int next_color_class    = cell_todo.next_cell(&queue_pointer, c);
            const int next_color_class_sz = c->ptn[next_color_class] + 1;

            if(m)
                m->color_refinement_cost += next_color_class_sz;

            bool dense_dense = (g->d[c->lab[next_color_class]] > (g->v_size / (next_color_class_sz + 1)));

            if(next_color_class_sz == 1 && !(config.CONFIG_IR_DENSE && dense_dense)) {
                // singleton
                comp = comp && refine_color_class_singleton(g, c, next_color_class, next_color_class_sz,
                                                            &color_class_splits, I);
            } else if(config.CONFIG_IR_DENSE) {
                if(dense_dense) { // dense-dense
                    comp = comp && refine_color_class_dense_dense(g, c, next_color_class, next_color_class_sz,
                                                                  &color_class_splits, I);
                } else { // dense-sparse
                    comp = comp && refine_color_class_dense(g, c, next_color_class, next_color_class_sz,
                                                            &color_class_splits, I);
                }
            } else { // sparse
                comp = comp && refine_color_class_sparse(g, c, next_color_class, next_color_class_sz, &color_class_splits, I);
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

                c->cells += (old_class != new_class);
                const int class_size = c->ptn[new_class];
                c->smallest_cell_lower_bound = (class_size < c->smallest_cell_lower_bound)?
                                               class_size:c->smallest_cell_lower_bound;

#ifndef NDEBUG // debug code
                if(color_class_splits.empty()) {
                int actual_cells = 0;
                for (int i = 0; i < c->ptn_sz;) {
                    actual_cells += 1;
                    i += c->ptn[i] + 1;
                }

                assert(c->cells == actual_cells);
            }
#endif

                if(c->cells == g->v_size) {
                    color_class_splits.reset();
                    cell_todo.reset(&queue_pointer);
                    I->write_cells(c->cells);
                    comp = comp && I->write_top_and_compare(ENDREF_MARK);
                    return comp;
                }

                // fast forward coloring
                if(c->cells == cell_early) {
                    I->fast_forward(ENDREF_MARK);
                    color_class_splits.reset();
                    cell_todo.reset(&queue_pointer);
                    return comp;
                }

                if(latest_old_class != old_class) {
                    latest_old_class = old_class;
                    skipped_largest = false;
                }

                color_class_splits.pop_back();
                int new_class_sz = c->ptn[new_class] + 1;


                if(skipped_largest || !is_largest) {
                    cell_todo.add_cell(&queue_pointer, new_class);
                } else {
                    skipped_largest = true;
                    skip += 1;

                    // since old color class will always appear last, the queue pointer of old color class is still valid!
                    int i = queue_pointer.get(old_class);
                    if(i >= 0) {
                        cell_todo.replace_cell(&queue_pointer, old_class, new_class);
                    }
                }
            }
            if(!comp) break;
        }

        I->write_cells(c->cells);
        comp = comp && I->write_top_and_compare(ENDREF_MARK);
        assert((comp && !I->never_fail)?assert_is_equitable(g, c):true);

        return comp;
    }

    int  individualize_vertex(coloring_temp<vertex_type>* c, int v) {
        const int color = c->vertex_to_col[v];
        const int pos   = c->vertex_to_lab[v];

        int color_class_size = c->ptn[color];

        assert(color_class_size > 0);

        const int vertex_at_pos = c->lab[color + color_class_size];
        c->lab[pos] = vertex_at_pos;
        c->vertex_to_lab[vertex_at_pos] = pos;

        c->lab[color + color_class_size] = v;
        c->vertex_to_lab[v] = color + color_class_size;
        c->vertex_to_col[v] = color + color_class_size;

        c->ptn[color] -= 1;
        c->ptn[color + color_class_size - 1] = 0;
        c->cells += 1;
        return color + color_class_size;
    }

    bool refine_coloring_first(sgraph_temp<vertex_type, degree_type, edge_type>  *g, coloring_temp<vertex_type> *c,
                               int init_color_class) {
        assure_initialized(g);

        cell_todo.reset(&queue_pointer);

        if(init_color_class < 0) {
            for (int i = 0; i < c->ptn_sz;) {
                cell_todo.add_cell(&queue_pointer, i);
                i += c->ptn[i] + 1;
            }
        } else {
            cell_todo.add_cell(&queue_pointer, init_color_class);
        }
        int its = 0;

        while(!cell_todo.empty()) {
            its += 1;
            color_class_splits.reset();
            const int next_color_class = cell_todo.next_cell(&queue_pointer, c);
            const int next_color_class_sz = c->ptn[next_color_class] + 1;
            const bool dense_dense = (g->d[c->lab[next_color_class]] > (g->v_size / (next_color_class_sz + 1)));

            if(next_color_class_sz == 1 && !(config.CONFIG_IR_DENSE && dense_dense)) {
                // singleton
                refine_color_class_singleton_first(g, c, next_color_class, next_color_class_sz,
                                                                  &color_class_splits);
            } else if(config.CONFIG_IR_DENSE) {
                if(dense_dense) { // dense-dense
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
            int skip = 0;

            int  latest_old_class = -1;
            bool skipped_largest = false;

            // color class splits are sorted in reverse
            // the old color class will always come last
            while(!color_class_splits.empty()) {
                const int  old_class  = color_class_splits.last()->first.first;
                const int  new_class  = color_class_splits.last()->first.second;
                const bool is_largest = color_class_splits.last()->second;

                c->cells += (old_class != new_class);
                const int class_size             = c->ptn[new_class] + 1;
                const int class_size_non_trivial = (class_size == 1)?INT32_MAX:class_size;
                c->smallest_cell_lower_bound = (class_size < c->smallest_cell_lower_bound)?
                        class_size:c->smallest_cell_lower_bound;

#ifndef NDEBUG // debug code
                if(color_class_splits.empty()) {
                int actual_cells = 0;
                for (int i = 0; i < c->ptn_sz;) {
                    actual_cells += 1;
                    i += c->ptn[i] + 1;
                }

                assert(c->cells == actual_cells);
            }
#endif

                if(c->cells == g->v_size) {
                    color_class_splits.reset();
                    cell_todo.reset(&queue_pointer);
                    assert(c->check());
                    assert(assert_is_equitable(g, c));
                    return true;
                }

                if(latest_old_class != old_class) {
                    latest_old_class = old_class;
                    skipped_largest = false;
                }

                color_class_splits.pop_back();
                const int new_class_sz = c->ptn[new_class] + 1;


                if(skipped_largest || !is_largest) {
                    cell_todo.add_cell(&queue_pointer, new_class);
                } else {
                    skipped_largest = true;
                    skip += 1;

                    // since old color class will always appear last, the queue pointer of old color class is still valid!
                    int i = queue_pointer.get(old_class);
                    if(i >= 0) {
                        cell_todo.replace_cell(&queue_pointer, old_class, new_class);
                    }
                }
            }
        }

        assert(c->check());
        assert(assert_is_equitable(g, c));
        return true;
    }

    bool certify_automorphism(sgraph_temp<vertex_type, degree_type, edge_type>  *g, bijection_temp<vertex_type> *p) {
        assert(p->map_sz == g->v_size);
        int i, found;

        for(i = 0; i < g->v_size; ++i) {
            const int image_i = p->map_vertex(i);
            if(image_i == i)
                continue;
            if(g->d[i] != g->d[image_i]) // degrees must be equal
                return false;

            scratch_set.reset();
            // automorphism must preserve neighbours
            found = 0;
            for(int j = g->v[i]; j < g->v[i] + g->d[i]; ++j) {
                const int vertex_j = g->e[j];
                const int image_j  = p->map_vertex(vertex_j);
                scratch_set.set(image_j);
                found += 1;
            }
            for(int j = g->v[image_i]; j < g->v[image_i] + g->d[image_i]; ++j) {
                const int vertex_j = g->e[j];
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

    bool assert_is_equitable(sgraph_temp<vertex_type, degree_type, edge_type>  *g, coloring_temp<vertex_type> *c) {
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
                for (size_t j = 0; j < neighbour_col.size(); ++j) {
                    test = test && (neighbour_col[j] == neighbour_col_canon[j]);
                }
            }

            if(c->ptn[i] == 0) {
                first_of_color = true;
            }
        }
        int expect0 = 0;
        bool expectsize = true;
        for(int i = 0; i < c->lab_sz; ++i) {
            if(expectsize) {
                expectsize = false;
                expect0 = c->ptn[i];
            }
            if(expect0 == 0) {
                test = test && (c->ptn[i] == 0);
                expectsize = true;
            } else {
                test = test && (c->ptn[i] > 0);
            }
            expect0 -= 1;
        }
        return test;
    }

    ~refinement_temp() {
        if(initialized)
            delete[] workspace_int;
    }

private:
    bool initialized = false;
    work_set_int  queue_pointer;
    cell_worklist_temp<vertex_type> cell_todo;
    mark_set      scratch_set;
    work_list_temp<vertex_type> vertex_worklist;
    work_set_temp<vertex_type>  color_vertices_considered;
    work_set_temp<vertex_type>  neighbours; // degree type instead?
    work_set_temp<vertex_type>  neighbour_sizes;
    work_list_temp<vertex_type> singletons;
    work_list_temp<vertex_type> old_color_classes;
    work_list_pair_bool color_class_splits;

    vertex_type* scratch;
    int* workspace_int;

    void assure_initialized(sgraph_temp<vertex_type, degree_type, edge_type>  *g) {
        if(!initialized) {
            const int n = g->v_size;

            // reducing contention on heap allocator through bulk allocation...
            workspace_int = new int[n * 2];

            vertex_worklist.initialize(n * 2);
            singletons.initialize(n);
            old_color_classes.initialize(n);
            neighbours.initialize(n);
            neighbour_sizes.initialize(n);
            queue_pointer.initialize(n);
            color_vertices_considered.initialize(n);

            scratch = (vertex_type*) workspace_int;
            scratch_set.initialize_from_array(workspace_int + n, n);

            color_class_splits.initialize(n);
            cell_todo.initialize(n * 2);

            memset(scratch, 0, n * sizeof(vertex_type));
            initialized = true;
        }
    }

    bool refine_color_class_sparse(sgraph_temp<vertex_type, degree_type, edge_type>  *g, coloring_temp<vertex_type> *c,
                                   int color_class, int class_size,
                                   work_list_pair_bool* color_class_split_worklist, invariant* I) {
        // for all vertices of the color class...
        bool comp, mark_as_largest;
        int i, j, cc, end_cc, largest_color_class_size, acc_in, singleton_inv1, singleton_inv2, acc, total;

        cc = color_class; // iterate over color class
        comp = true;

        singletons.reset();
        scratch_set.reset();
        old_color_classes.reset();
        neighbours.reset();
        color_vertices_considered.reset();

        end_cc = color_class + class_size;
        acc_in = 0;
        singleton_inv1 = 0;
        singleton_inv2 = 0;
        while (cc < end_cc) { // increment value of neighbours of vc by 1
            const int vc = c->lab[cc];
            const int pe = g->v[vc];
            const int end_i = pe + g->d[vc];
            for (i = pe; i < end_i; i++) {
                const int v   = g->e[i];
                const int col = c->vertex_to_col[v];
                if (c->ptn[col] == 0) {
                    singleton_inv1 += (col + 1) * (423733 - (col + 1));
                    singleton_inv2 += (col + 3) * (723732 - (col + 2));
                    continue;
                }
                neighbours.inc_nr(v);
                if(neighbours.get(v) == 0) {
                    color_vertices_considered.inc_nr(col);
                    assert(col + color_vertices_considered.get(col) < g->v_size);
                    scratch[col + color_vertices_considered.get(col)] = v; // hit vertices
                    if (!scratch_set.get(col)) {
                        old_color_classes.push_back(col);
                        acc_in += (col + 1) * (423233 - (col + 1));
                        scratch_set.set(col);
                    }
                }
            }
            cc += 1;
        }

        // write invariant for singleton color classes
        comp = comp && I->write_top_and_compare(singleton_inv1);
        comp = comp && I->write_top_and_compare(singleton_inv2);
        comp = comp && I->write_top_and_compare(-acc_in);

        // early out before sorting color classes
        if(!comp) {
            while (!old_color_classes.empty()) {
                const int _col = old_color_classes.pop_back();
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
            const int _col    = old_color_classes.pop_back();
            const int _col_sz = c->ptn[_col] + 1;
            neighbour_sizes.reset();
            vertex_worklist.reset();

            total = 0;
            for(i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
                const int v = scratch[_col + i];
                int index = neighbours.get(v) + 1;
                total += 1;
                if(neighbour_sizes.inc(index) == 0)
                    vertex_worklist.push_back(index);
            }

            vertex_worklist.sort();

            // enrich neighbour_sizes to accumulative counting array
            acc = 0;
            while(!vertex_worklist.empty()) {
                const int i   = vertex_worklist.pop_back();
                const int val = neighbour_sizes.get(i) + 1;
                if(val >= 1) {
                    neighbour_sizes.set(i, val + acc);
                    acc += val;
                    const int __col = _col + _col_sz - (neighbour_sizes.get(i));
                    const int v_degree = i;
                    if(__col != _col)
                        c->ptn[__col] = -1;
                    comp = comp && I->write_top_and_compare(__col + v_degree * g->v_size);
                    comp = comp && I->write_top_and_compare(g->v_size * 7 + val + 1);
                }
            }

            const int vcount = color_vertices_considered.get(_col);

            // early out
            if(!comp) {
                bool hard_reset = false;
                if(2 * vcount > g->v_size) {
                    neighbours.reset_hard();
                    hard_reset = true;
                } else {
                    j = 0;
                    while (j < vcount + 1) {
                        const int v = scratch[_col + j];
                        neighbours.set(v, -1);
                        ++j;
                    }
                }
                color_vertices_considered.set(_col, -1);
                while (!old_color_classes.empty()) {
                    const int __col = old_color_classes.pop_back();
                    if(!hard_reset) {
                        for (i = 0; i < color_vertices_considered.get(__col) + 1; ++i)
                            neighbours.set(scratch[__col + i], -1);
                    }
                    color_vertices_considered.set(__col, -1);
                }
                neighbour_sizes.reset();
                vertex_worklist.reset();
                color_vertices_considered.reset();
                return comp;
            }
            //

            vertex_worklist.reset();
            j = 0;
            color_vertices_considered.set(_col, -1);

            // determine colors and rearrange vertices
            while (j < vcount + 1) {
                const int v = scratch[_col + j];
                ++j;
                if((neighbours.get(v) == -1))
                    continue;
                const int v_new_color = _col + _col_sz - neighbour_sizes.get(neighbours.get(v) + 1);
                neighbours.set(v, -1);
                if(v_new_color == _col)
                    continue;
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

                color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                        std::pair<int, int>(_col, i), mark_as_largest));
                if(i != 0)
                    c->ptn[i - 1] = 0;
                i += c->ptn[i] + 1;
            }

            if(i != 0)
                c->ptn[i - 1] = 0;
        }

        neighbour_sizes.reset();
        vertex_worklist.reset();

        return comp;
    }

    bool refine_color_class_dense(sgraph_temp<vertex_type, degree_type, edge_type>  *g, coloring_temp<vertex_type> *c,
                                  int color_class, int class_size,
                                  work_list_pair_bool* color_class_split_worklist, invariant* I) {
        bool comp;
        int i, cc, acc, largest_color_class_size, singleton_inv, pos;
        cc = color_class; // iterate over color class
        comp = true;

        neighbours.reset_hard();
        scratch_set.reset();
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
                    neighbours.inc(v); // want to use reset variant?
                    if (!scratch_set.get(col)) {
                        scratch_set.set(col);
                        old_color_classes.push_back(col);
                    }
                } else {
                    if(config.CONFIG_IR_FULL_INVARIANT)
                        vertex_worklist.push_back(col);
                    else {
                        singleton_inv += (col + 1) * (23524361 - col * 3);
                    }
                }
            }
            cc += 1;
        }

        // write singletons

        comp = comp && I->write_top_and_compare(g->v_size * 3 + old_color_classes.cur_pos);

        if(!comp) return comp;

        if(config.CONFIG_IR_FULL_INVARIANT) {
            vertex_worklist.sort();
            while (!vertex_worklist.empty()) {
                const int col = vertex_worklist.pop_back();
                comp = comp && I->write_top_and_compare(g->v_size * 9 + col);
            }
        } else {
            comp = comp && I->write_top_and_compare(singleton_inv);
        }

        if(!comp) return comp;

        old_color_classes.sort();

        // for every cell to be split...
        while(!old_color_classes.empty()) {
            const int col = old_color_classes.pop_back();
            const int col_sz = c->ptn[col] + 1;

            if(col_sz == 1)
                continue;

            neighbour_sizes.reset();
            vertex_worklist.reset();

            const int first_index = neighbours.get(c->lab[col]) + 1;
            vertex_worklist.push_back(first_index);
            int total = 0;
            for(i = col; i < col + col_sz; ++i) {
                const int v = c->lab[i];
                int index = neighbours.get(v) + 1;
                if(index == first_index)
                    continue;
                total += 1;
                if(neighbour_sizes.inc(index) == 0)
                    vertex_worklist.push_back(index);
            }

            neighbour_sizes.inc(first_index);
            neighbour_sizes.set(first_index, col_sz - total - 1);

            comp = comp && I->write_top_and_compare(vertex_worklist.cur_pos);
            if(!comp) return comp;

            if(vertex_worklist.cur_pos == 1) {
                // no split
                const int v_degree = neighbours.get(c->lab[col]);
                comp = comp && I->write_top_and_compare(col + g->v_size * v_degree);
                continue;
            }

            vertex_worklist.sort();
            // enrich neighbour_sizes to accumulative counting array
            acc = 0;
            while(!vertex_worklist.empty()) {
                const int i   = vertex_worklist.pop_back();
                const int val = neighbour_sizes.get(i) + 1;
                if(val >= 1) {
                    neighbour_sizes.set(i, val + acc);
                    acc += val;
                    const int _col = col + col_sz - (neighbour_sizes.get(i));
                    c->ptn[_col] = -1; // this is val - 1, actually...

                    const int v_degree = i;
                    comp = comp && I->write_top_and_compare(_col + g->v_size * v_degree);
                    comp = comp && I->write_top_and_compare(_col + val + 1);
                    if(!comp) return comp;
                }
            }

            // copy cell for rearranging
            memcpy(scratch, c->lab + col, col_sz * sizeof(vertex_type));
            //vertex_worklist.cur_pos = col_sz;
            pos = col_sz;

            // determine colors and rearrange
            // split color classes according to count in counting array
            while(pos > 0) {
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
            for(i = col; i < col + col_sz;) {
                assert(i >= 0 && i < c->ptn_sz);
                assert(c->ptn[i] + 1 > 0);
                const int v_color  = i;
                const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
                largest_color_class_size = mark_as_largest?c->ptn[i] + 1:largest_color_class_size;

                color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                        std::pair<int, int>(col, i), mark_as_largest));
                if(v_color != 0)
                    c->ptn[v_color - 1] = 0;
                i += c->ptn[i] + 1;
            }
        }

        return comp;
    }

    bool refine_color_class_dense_dense(sgraph_temp<vertex_type, degree_type, edge_type>  *g, coloring_temp<vertex_type> *c,
                                        int color_class, int class_size,
                                        work_list_pair_bool* color_class_split_worklist, invariant* I) {
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
                comp = comp && I->write_top_and_compare(col + v_degree * g->v_size);
                if(!comp) {neighbours.reset_hard(); return comp;}
                continue;
            }

            neighbour_sizes.reset();
            vertex_worklist.reset();

            const int first_index = neighbours.get(c->lab[col]) + 1;
            vertex_worklist.push_back(first_index);
            int total = 0;
            for(i = col; i < col + col_sz; ++i) {
                const int v = c->lab[i];
                const int index = neighbours.get(v) + 1;
                if(index == first_index)
                    continue;
                total += 1;
                if(neighbour_sizes.inc(index) == 0)
                    vertex_worklist.push_back(index);
            }

            neighbour_sizes.inc(first_index);
            neighbour_sizes.set(first_index, col_sz - total - 1);

            comp = comp && I->write_top_and_compare(g->v_size * 12 + vertex_worklist.cur_pos);
            if(!comp) {neighbours.reset_hard(); return comp;}

            if(vertex_worklist.cur_pos == 1) {
                // no split
                const int v_degree = neighbours.get(c->lab[col]);
                //comp = comp && I->write_top_and_compare(-g->v_size * 10 - col);
                comp = comp && I->write_top_and_compare(col + v_degree * g->v_size);
                comp = comp && I->write_top_and_compare(col + c->ptn[col] + 1);
                if(!comp) {neighbours.reset_hard(); return comp;}
                continue;
            }

            vertex_worklist.sort();
            // enrich neighbour_sizes to accumulative counting array
            acc = 0;
            while(!vertex_worklist.empty()) {
                const int i   = vertex_worklist.pop_back();
                const int val = neighbour_sizes.get(i) + 1;
                if(val >= 1) {
                    neighbour_sizes.set(i, val + acc);
                    acc += val;
                    const int _col = col + col_sz - (neighbour_sizes.get(i));
                    c->ptn[_col] = -1; // this is val - 1, actually...

                    const int v_degree = i;
                    comp = comp && I->write_top_and_compare(-g->v_size * 5 - _col);
                    comp = comp && I->write_top_and_compare(_col + v_degree);
                    comp = comp && I->write_top_and_compare(_col + val + 1);
                    if(!comp) {neighbours.reset_hard(); return comp;}
                }
            }

            vertex_worklist.reset();

            // copy cell for rearranging
            memcpy(vertex_worklist.arr, c->lab + col, col_sz * sizeof(vertex_type));
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

                color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                        std::pair<int, int>(col, i), mark_as_largest));
                if(v_color != 0)
                    c->ptn[v_color - 1] = 0;
                i += c->ptn[i] + 1;
            }
        }

        neighbours.reset_hard();
        return comp;
    }

    bool refine_color_class_singleton(sgraph_temp<vertex_type, degree_type, edge_type>  *g, coloring_temp<vertex_type> *c,
                                      int color_class, int class_size,
                                      work_list_pair_bool *color_class_split_worklist, invariant *I) {
        bool comp;
        int i, cc, deg1_write_pos, deg1_read_pos, singleton_inv;
        cc = color_class; // iterate over color class
        comp = true;

        neighbours.reset();
        scratch_set.reset();
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

            if(!scratch_set.get(col)) {
                scratch_set.set(col);
                old_color_classes.push_back(col);
                // neighbours acts as the write position for degree 1 vertices of col in the scratchpad
                neighbours.set(col, col);
            }
            // write vertex to scratchpad for later use
            scratch[neighbours.get(col)] = v;
            neighbours.inc_nr(col); // we reset neighbours later, since old_color_classes for reset
        }

        singleton_inv += old_color_classes.cur_pos;
        if(!config.CONFIG_IR_FULL_INVARIANT)
            comp = comp && I->write_top_and_compare(g->v_size * 3 + singleton_inv);

        if(!comp) {
            while(!old_color_classes.empty())
                neighbours.set(old_color_classes.pop_back(), -1);
            return comp;
        }

        old_color_classes.sort();

        // write invariant first...
        for(i = 0; i < old_color_classes.cur_pos && comp; ++i) {
            //comp = comp && I->write_top_and_compare(old_color_classes.arr[i]); // color class
            comp = comp && I->write_top_and_compare(g->v_size * 14 + neighbours.get(old_color_classes.arr[i]));
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
            color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                    std::pair<int, int>(deg0_col, deg0_col), !leq));
            color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                    std::pair<int, int>(deg0_col, deg1_col), leq));

            // reset neighbours count to -1
            neighbours.set(deg0_col, -1);
        }

        return comp;
    }

    bool refine_color_class_singleton_first(sgraph_temp<vertex_type, degree_type, edge_type>  *g, coloring_temp<vertex_type> *c,
                                            int color_class, int class_size,
                                            work_list_pair_bool *color_class_split_worklist) {
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
            const int v   = g->e[i];
            const int col = c->vertex_to_col[v];

            if(c->ptn[col] == 0)
                continue;

            if(!scratch_set.get(col)) {
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
            color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                    std::pair<int, int>(deg0_col, deg0_col), !leq));
            color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                    std::pair<int, int>(deg0_col, deg1_col), leq));

            // reset neighbours count to -1
            neighbours.set(deg0_col, -1);
        }

        return comp;
    }

    bool refine_color_class_dense_first(sgraph_temp<vertex_type, degree_type, edge_type>  *g, coloring_temp<vertex_type> *c,
                                        int color_class, int class_size,
                                        work_list_pair_bool *color_class_split_worklist) {
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
                const int v   = g->e[i];
                const int col = c->vertex_to_col[v];
                if(c->ptn[col] == 0)
                    continue;

                neighbours.inc(v);
                if(!scratch_set.get(col)) {
                    scratch_set.set(col);
                    old_color_classes.push_back(col);
                }
            }
            cc += 1;
        }

        while(!old_color_classes.empty()) {
            const int col = old_color_classes.pop_back();
            const int col_sz = c->ptn[col] + 1;

            if(col_sz == 1) {
                const int v_degree = neighbours.get(c->lab[col]) + 1;
                if(v_degree == -1)
                    continue;
            }

            neighbour_sizes.reset();
            vertex_worklist.reset();

            const int first_index = neighbours.get(c->lab[col]) + 1;
            vertex_worklist.push_back(first_index);
            int total = 0;
            for(i = col; i < col + col_sz; ++i) {
                const int v = c->lab[i];
                const int index = neighbours.get(v) + 1;
                if(index == first_index)
                    continue;
                total += 1;
                if(neighbour_sizes.inc(index) == 0)
                    vertex_worklist.push_back(index);
            }

            neighbour_sizes.inc(first_index);
            neighbour_sizes.set(first_index, col_sz - total - 1);

            if(vertex_worklist.cur_pos == 1) {
                continue;
            }

            // enrich neighbour_sizes to accumulative counting array
            acc = 0;
            while(!vertex_worklist.empty()) {
                const int i   = vertex_worklist.pop_back();
                const int val = neighbour_sizes.get(i) + 1;
                if(val >= 1) {
                    neighbour_sizes.set(i, val + acc);
                    acc += val;
                    const int _col = col + col_sz - (neighbour_sizes.get(i));
                    c->ptn[_col] = -1; // this is val - 1, actually...
                }
            }

            // copy cell for rearranging
            memcpy(scratch, c->lab + col, col_sz * sizeof(vertex_type));
            pos = col_sz;

            // determine colors and rearrange
            // split color classes according to count in counting array
            while(pos > 0) {
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
            for(i = col; i < col + col_sz;) {
                assert(i >= 0 && i < c->ptn_sz);
                assert(c->ptn[i] + 1 > 0);
                const int  v_color  = i;
                const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
                largest_color_class_size = mark_as_largest?c->ptn[i] + 1:largest_color_class_size;

                color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                        std::pair<int, int>(col, i), mark_as_largest));
                if(v_color != 0)
                    c->ptn[v_color - 1] = 0;
                i += c->ptn[i] + 1;
            }
        }

        return true;
    }

    bool refine_color_class_dense_dense_first(sgraph_temp<vertex_type, degree_type, edge_type>  *g, coloring_temp<vertex_type> *c,
                                              int color_class, int class_size,
                                              work_list_pair_bool *color_class_split_worklist) {
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
                const int v   = g->e[i];
                neighbours.inc_nr(v);
            }
            cc += 1;
        }

        // dont sort, just iterate over all cells
        for(j = 0; j < g->v_size;) {
            const int col    = j;
            const int c_ptn = c->ptn[j];
            j += c_ptn + 1;
            const int col_sz = c_ptn + 1;

            if(col_sz == 1)
                continue;

            neighbour_sizes.reset();
            vertex_worklist.reset();

            const int first_index = neighbours.get(c->lab[col]) + 1;
            vertex_worklist.push_back(first_index);
            int total = 0;
            for(i = col; i < col + col_sz; ++i) {
                const int v = c->lab[i];
                const int  index = neighbours.get(v) + 1;
                if(index == first_index)
                    continue;
                total += 1;
                if(neighbour_sizes.inc(index) == 0)
                    vertex_worklist.push_back(index);
            }

            neighbour_sizes.inc(first_index);
            neighbour_sizes.set(first_index, col_sz - total - 1);

            if(vertex_worklist.cur_pos == 1) {
                continue;
            }

            // enrich neighbour_sizes to accumulative counting array
            acc = 0;
            while(!vertex_worklist.empty()) {
                const int i   = vertex_worklist.pop_back();
                const int val = neighbour_sizes.get(i) + 1;
                if(val >= 1) {
                    neighbour_sizes.set(i, val + acc);
                    acc += val;
                    const int _col = col + col_sz - (neighbour_sizes.get(i));
                    c->ptn[_col] = -1;
                }
            }

            // copy cell for rearranging
            memcpy(scratch, c->lab + col, col_sz * sizeof(vertex_type));
            pos = col_sz;

            // determine colors and rearrange
            // split color classes according to count in counting array
            while(pos > 0) {
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
            for(i = col; i < col + col_sz;) {
                assert(i >= 0 && i < c->ptn_sz);
                assert(c->ptn[i] + 1 > 0);
                const int v_color  = i;
                const bool mark_as_largest = largest_color_class_size < c->ptn[i] + 1;
                largest_color_class_size = mark_as_largest?c->ptn[i] + 1:largest_color_class_size;

                color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                        std::pair<int, int>(col, i), mark_as_largest));
                if(v_color != 0)
                    c->ptn[v_color - 1] = 0;
                i += c->ptn[i] + 1;
            }
        }

        neighbours.reset_hard();
        return true;
    }

    bool refine_color_class_sparse_first(sgraph_temp<vertex_type, degree_type, edge_type>  *g, coloring_temp<vertex_type> *c,
                                         int color_class, int class_size,
                                         work_list_pair_bool* color_class_split_worklist) {
        bool comp;
        int i, cc, v_new_color, largest_color_class_size, acc;

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
            for (i = pe; i < end_i; i++) {
                const int v   = g->e[i];
                const int col = c->vertex_to_col[v];

                if (c->ptn[col] == 0)
                    continue;

                neighbours.inc_nr(v);
                if(neighbours.get(v) == 0) {
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
        while(!old_color_classes.empty()) {
            const int _col    = old_color_classes.pop_back();
            const int _col_sz = c->ptn[_col] + 1;

            neighbour_sizes.reset();
            vertex_worklist.reset();

            int total = 0;
            for(i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
                const int v = scratch[_col + i];
                int index = neighbours.get(v) + 1;
                total += 1;
                if(neighbour_sizes.inc(index) == 0)
                    vertex_worklist.push_back(index);
            }

            // enrich neighbour_sizes to accumulative counting array
            acc = 0;
            while(!vertex_worklist.empty()) {
                const int i   = vertex_worklist.pop_back();
                const int val = neighbour_sizes.get(i) + 1;
                if(val >= 1) {
                    neighbour_sizes.set(i, val + acc);
                    acc += val;
                }
            }

            // determine colors
            vertex_worklist.reset();

            for(i = 0; i < color_vertices_considered.get(_col) + 1; ++i) {
                const int v = scratch[_col + i];
                if (neighbours.get(v) == -1) {
                    v_new_color = _col; assert(false);
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
                const int v           = vertex_worklist.pop_back();
                const int v_new_color = vertex_worklist.pop_back();

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
                const bool mark_as_largest = largest_color_class_size < (c->ptn[i] + 1);
                largest_color_class_size = mark_as_largest?(c->ptn[i] + 1):largest_color_class_size;

                color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(
                        std::pair<int, int>(_col, i), mark_as_largest));
                if(i != 0)
                    c->ptn[i - 1] = 0;
                i += c->ptn[i] + 1;
            }

            if(i != 0)
                c->ptn[i - 1] = 0;
        }

        neighbour_sizes.reset();
        vertex_worklist.reset();

        return comp;
    }
};

// typedef refinement_temp<int, int, int> refinement;

#endif //DEJAVU_REFINEMENT_H
