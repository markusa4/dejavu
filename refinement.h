//
// Created by markus on 19.09.19.
//

#ifndef BRUTUS_REFINEMENT_H
#define BRUTUS_REFINEMENT_H


#include <queue>
#include <set>
#include "coloring.h"
#include "sgraph.h"
#include "invariant_acc.h"
#include <list>
#include <iostream>
#include <bitset>
#include <mutex>

static std::mutex malloc_lock;

class work_queue {
public:
    void initialize(int size);
    void push(int val);
    int pop();
    bool empty();
    ~work_queue();
private:
    int* queue;
    bool init = false;
    int  sz;
    int  pos;
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
    void initialize_from_array(bool *p, int size);
private:
    work_queue reset_queue;
    //std::vector<bool> s;
    bool init = false;
    bool* s;
    int sz;
};

class work_list {
public:
    void initialize(int size);
    void push_back(int index);
    int pop_back();
    bool empty();
    void reset();
    ~work_list();
private:
    int* arr;
    int arr_sz = -1;
    int cur_pos;
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

class refinement {
public:
    bool refine_coloring(sgraph* g, coloring* c, std::list<std::pair<int, int>> *changes, invariant* I, int init_color_class, bool track_changes);
    int  individualize_vertex(coloring* c, int v);
    void undo_individualize_vertex(sgraph *g, coloring *c, int v);
    bool refine_color_class(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist, invariant* I);
    void undo_refine_color_class(sgraph *g, coloring *c, std::list<std::pair<int, int>> *changes);
    void complete_colorclass_invariant(sgraph *g, coloring *c, invariant_acc *I);
    bool refine_coloring_first(sgraph *g, coloring *c, int init_color_class);
    bool assert_is_equitable(sgraph *g, coloring *c);
    ~refinement();
private:
    bool initialized = false;
    cumulative_counting counting_array;
    work_set  color_workset;
    work_list color_worklist_vertex;
    work_list_pair color_worklist_color;
    work_list_pair vertex_worklist;
    work_list_pair old_color_classes;
    work_list_pair_bool color_class_splits;
    ring_pair worklist_color_classes;

    std::pair<int, int>* p;

    bool refine_color_class_first(sgraph *g, coloring *c, int color_class, int class_size,
                                  work_list_pair_bool *color_class_split_worklist);
};


#endif //BRUTUS_REFINEMENT_H
