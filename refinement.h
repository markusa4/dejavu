#ifndef DEJAVU_REFINEMENT_H
#define DEJAVU_REFINEMENT_H

#include "coloring.h"
#include "sgraph.h"
#include "invariant.h"
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

    int  cur_pos;
    T*   arr;
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

class cell_worklist {
public:
    void initialize(int domain_size);
    int add_cell(work_set_int* queue_pointer, int col);
    int next_cell(work_set_int* queue_pointer, coloring* c);
    void replace_cell(work_set_int* queue_pointer, int col_old, int col);
    void reset(work_set_int* queue_pointer);
    bool empty();

private:
    int* arr = nullptr;
    int  arr_sz  = -1;
    int  cur_pos = -1;
    bool init    = false;
};

class refinement {
public:
    bool refine_coloring(sgraph *g, coloring *c, change_tracker *, invariant *I, int init_color_class,
                         bool track_changes, strategy_metrics *m);

    int  individualize_vertex(coloring* c, int v);

    bool refine_coloring_first(sgraph *g, coloring *c, int init_color_class);

    bool certify_automorphism(sgraph *g, bijection *p);

    bool assert_is_equitable(sgraph *g, coloring *c);

    ~refinement();

private:
    bool initialized = false;
    work_set_int  queue_pointer;
    cell_worklist cell_todo;
    mark_set      scratch_set;
    work_list     vertex_worklist;
    work_set_int  color_vertices_considered;
    work_set_int  neighbours;
    work_set_int  neighbour_sizes;
    work_list     singletons;
    work_list     old_color_classes;
    work_list_pair_bool color_class_splits;

    int* scratch;
    int* workspace_int;

    void assure_initialized(sgraph *g);

    bool refine_color_class_sparse(sgraph *g, coloring *c, int color_class, int class_size,
                                   work_list_pair_bool* color_class_split_worklist, invariant* I);

    bool refine_color_class_dense(sgraph *g, coloring *c, int color_class, int class_size,
                                  work_list_pair_bool* color_class_split_worklist, invariant* I);

    bool refine_color_class_dense_dense(sgraph *g, coloring *c, int color_class, int class_size,
                                        work_list_pair_bool* color_class_split_worklist, invariant* I);

    bool refine_color_class_singleton(sgraph *g, coloring *c, int color_class, int class_size,
                                      work_list_pair_bool *color_class_split_worklist, invariant *I);

    bool refine_color_class_singleton_first(sgraph *g, coloring *c, int color_class, int class_size,
                                            work_list_pair_bool *color_class_split_worklist);

    bool refine_color_class_dense_first(sgraph *g, coloring *c, int color_class, int class_size,
                                        work_list_pair_bool *color_class_split_worklist);

    bool refine_color_class_dense_dense_first(sgraph *g, coloring *c, int color_class, int class_size,
                                              work_list_pair_bool *color_class_split_worklist);

    bool refine_color_class_sparse_first(sgraph *g, coloring *c, int color_class, int class_size,
                                                     work_list_pair_bool* color_class_split_worklist);
};




#endif //DEJAVU_REFINEMENT_H
