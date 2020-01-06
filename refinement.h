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

class work_queue {
public:
    void initialize(int size);
    void initialize_from_array(int* arr, int size);
    void push(int val);
    int pop();
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
   // std::vector<bool> s;
    bool init = false;
    bool* s;

    int sz;
};

class work_list {
public:
    void initialize(int size);
    void initialize_from_array(int* arr, int size);
    void push_back(int index);
    int pop_back();
    bool empty();
    void reset();
    ~work_list();
    void sort();
    int  cur_pos;
    int* arr;
    bool init = false;
private:

    int arr_sz = -1;
};


class alignas(16) work_set_int {
public:
    void initialize(int size);
    void initialize_from_array(int* arr, int size);
    void set(int index, int value);
    int  get(int index);
    void reset();
    void reset_hard();
    int inc(int index);
    void inc_nr(int index);
    ~work_set_int();
    work_queue reset_queue;
private:
    bool  init = false;
    int*  s;
    int   sz;
};

class alignas(16) work_set_char {
public:
    void initialize(int size);
    void set(int index, char value);
    char  get(int index);
    void reset();
    void reset_hard();
    int inc(int index);
    void inc_nr(int index);
    ~work_set_char();
    work_queue reset_queue;
private:
    bool  init = false;
    char*  s;
    int   sz;
};

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


class work_list_pair_bool {
public:
    void initialize(int size);
    void push_back(std::pair<std::pair<int, int>, bool> value);
    std::pair<std::pair<int, int>, bool>* last();
    void pop_back();
    void sort();
    bool empty();
    void reset();
    ~work_list_pair_bool();
private:
    std::pair<std::pair<int, int>, bool>* arr;
    bool init = false;
    int arr_sz = -1;
    int cur_pos;
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
    bool refine_coloring(sgraph* g, coloring* c, change_tracker *, invariant* I, int init_color_class, bool track_changes, strategy_metrics* m);
    int  individualize_vertex(coloring* c, int v);
    bool refine_coloring_first(sgraph *g, coloring *c, int init_color_class);
    bool certify_automorphism(sgraph *g, bijection *p);
    bool assert_is_equitable(sgraph *g, coloring *c);
    void undo_changes_with_reference(change_tracker* changes, coloring* c, coloring* old_c);
    ~refinement();

private:
    bool initialized = false;
    work_set_int queue_pointer;
    cell_worklist  cell_todo;
    mark_set scratch_set;
    work_list vertex_worklist;
    work_set_int color_vertices_considered;
    work_set_int neighbours;
    work_set_char neighbours127;
    work_set_int neighbour_sizes;
    work_list singletons;
    work_list_pair_bool color_class_splits;
    work_list old_color_classes;
    int* scratch;

    int* workspace_int;

    bool refine_color_class_sparse(sgraph *g, coloring *c, int color_class, int class_size,
                                   work_list_pair_bool* color_class_split_worklist, invariant* I);

    bool refine_color_class_sparse127(sgraph *g, coloring *c, int color_class, int class_size,
                                   work_list_pair_bool* color_class_split_worklist, invariant* I);

    bool refine_color_class_dense(sgraph *g, coloring *c, int color_class, int class_size,
                                  work_list_pair_bool* color_class_split_worklist, invariant* I);

    bool refine_color_class_dense_dense(sgraph *g, coloring *c, int color_class, int class_size,
                                        work_list_pair_bool* color_class_split_worklist, invariant* I);

    bool refine_color_class_dense_dense127(sgraph *g, coloring *c, int color_class, int class_size,
                                           work_list_pair_bool *color_class_split_worklist, invariant *I);

    bool refine_color_class_singleton(sgraph *g, coloring *c, int color_class, int class_size,
                                      work_list_pair_bool *color_class_split_worklist, invariant *I);

    void assure_initialized(sgraph *g);

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
