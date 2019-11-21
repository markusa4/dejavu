//
// Created by markus on 19.09.19.
//

#ifndef BRUTUS_REFINEMENT_H
#define BRUTUS_REFINEMENT_H


#include <queue>
#include <set>
#include "coloring.h"
#include "sgraph.h"
#include "invariant.h"
#include "utility.h"
#include <list>
#include <iostream>
#include <bitset>
#include <mutex>
#include <cstring>

static std::mutex malloc_lock;

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

class cumulative_counting {
public:
    void initialize(int size, int maxdegree);
    void reset();
    void increment(int index, int color);
    void increment_r(int index, int color);
    int get_size(int index,  int color);
    int get_count(int index);
    ~cumulative_counting();
private:
    int* count;
    bool init = false;
    //std::vector<int>* sizes;
    int sz;
    int* sizes;
    int m;
    work_queue reset_queue;
    work_queue reset_queue_sizes;
};

class work_set {
public:
    void initialize(int size);
    void set(int index);
    bool get(int index);
    void reset();
    ~work_set();
    // std::bitset< s;
    void reset_hard();
    void set_nr(int index);
    void initialize_from_array(bool *p, int size);
    void unset(int index);
    void reset_soft();
private:
    work_queue reset_queue;
    //std::vector<bool> s;
    bool init = false;
    bool* s;

    int sz;
};

class work_set_int {
public:
    void initialize(int size);
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

class work_list {
public:
    void initialize(int size);
    void push_back(int index);
    int pop_back();
    bool empty();
    void reset();
    ~work_list();
    void sort();
    int cur_pos;
    int* arr;
private:

    int arr_sz = -1;
};

class work_list_pair {
public:
    void initialize(int size);
    void push_back(std::pair<int, int> value);
    std::pair<int, int>* last();
    void pop_back();
    void sort();
    bool empty();
    void reset();
    ~work_list_pair();
    std::pair<int, int>* arr;
    bool init = false;
    int arr_sz = -1;

    void initialize_from_array(int *p, int size);
    void initialize_from_array(std::pair<int, int> *p, int size);
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
    int add_cell(work_set_int* queue_pointer, int col, int col_sz);
    std::pair<int, int> next_cell(work_set_int* queue_pointer);
    void replace_cell(work_set_int* queue_pointer, int col_old, int col, int col_sz);
    void reset(work_set_int* queue_pointer);
    bool empty();

private:
    std::pair<int, int>* arr = nullptr;
    int  arr_sz  = -1;
    int  cur_pos = -1;
    bool init    = false;
};


class refinement {
public:
    bool refine_coloring(sgraph* g, coloring* c, change_tracker *, invariant* I, int init_color_class, bool track_changes, strategy_metrics* m);
    int  individualize_vertex(coloring* c, int v);
    void undo_individualize_vertex(sgraph *g, coloring *c, int v);
    bool refine_color_class_cumulative(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist, invariant* I);
    void undo_refine_color_class(sgraph *g, coloring *c, std::list<std::pair<int, int>> *changes);
    bool refine_coloring_first(sgraph *g, coloring *c, int init_color_class);
    bool refine_color_class_dense(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist, invariant* I);
    bool refine_color_class_dense_dense(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist, invariant* I);
    bool assert_is_equitable(sgraph *g, coloring *c);
    void undo_changes_with_reference(change_tracker* changes, coloring* c, coloring* old_c);
    bool old_refine_color_class_first(sgraph *g, coloring *c, int color_class, int class_size,
                                      work_list_pair_bool *color_class_split_worklist);
    bool old_refine_coloring_first(sgraph *g, coloring *c, int init_color_class);
    ~refinement();

    bool initialized = false;
    bool counting_initialized = false;
    work_set_int queue_pointer;
    cell_worklist cell_todo;
    work_list color_worklist_vertex;
    work_list_pair color_worklist_color;
    work_set    color_workset;
    mark_set    scratch_set;
    work_list vertex_worklist;
    work_list degrees_worklist;
    ring_pair worklist_color_classes;
    work_set_int color_vertices_considered;
    work_set_int neighbours;
    work_set_int neighbour_sizes;
    std::pair<int, int>* p;
    work_list      singletons;
    work_list_pair old_color_classes;
    cumulative_counting counting_array;
    work_list_pair_bool color_class_splits;
    work_list dense_old_color_classes;
    int* scratch;

    bool refine_color_class_sparse(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist, invariant* I);

    bool refine_color_class_first(sgraph *g, coloring *c, int color_class, int class_size,
                                  work_list_pair_bool *color_class_split_worklist);

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

    void assure_initialized_accumulate_counting(sgraph *g);

    bool certify_automorphism(sgraph *g, bijection *p);
};




#endif //BRUTUS_REFINEMENT_H
