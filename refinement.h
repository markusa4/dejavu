//
// Created by markus on 19.09.19.
//

#ifndef BRUTUS_REFINEMENT_H
#define BRUTUS_REFINEMENT_H


#include <queue>
#include <set>
#include "coloring.h"
#include "graph.h"
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
};

class refinement {
public:
    void refine_coloring(graph* g, coloring* c, std::set<std::pair<int, int>> *changes);
    void individualize_vertex(graph* g, coloring* c, int v);
    void undo_individualize_vertex(graph *g, coloring *c, int v);
    void refine_color_class(graph *g, coloring *c, int color_class, int class_size, std::list<std::pair<int, int>> *color_class_split_worklist);
    void undo_refine_color_class(graph *g, coloring *c, std::set<std::pair<int, int>> *changes);
private:
    cumulative_counting counting_array;
};


#endif //BRUTUS_REFINEMENT_H
