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
#include <list>

class cumulative_counting {
public:
    void initialize(int size, coloring *c);
    void reset();
    void increment(int index);
    int get_size(int index);
    const int operator[](const size_t index);
private:
    coloring* c;
    std::vector<int> count;
    std::vector<std::vector<int>> sizes;
    std::queue<int>  reset_queue;
    std::vector<int> col_list;
    std::vector<int> col_list_short;

    void write_color_degrees(invariant *I);
};

class refinement {
public:
    void refine_coloring(sgraph* g, coloring* c, std::set<std::pair<int, int>> *changes, invariant* I);
    void individualize_vertex(sgraph* g, coloring* c, int v);
    void undo_individualize_vertex(sgraph *g, coloring *c, int v);
    void refine_color_class(sgraph *g, coloring *c, int color_class, int class_size, std::list<std::pair<int, int>> *color_class_split_worklist, invariant* I);
    void undo_refine_color_class(sgraph *g, coloring *c, std::set<std::pair<int, int>> *changes);
private:
    cumulative_counting counting_array;
};


#endif //BRUTUS_REFINEMENT_H
