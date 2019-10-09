//
// Created by markus on 09.10.19.
//

#ifndef DEJAVU_REFINEMENT_BUCKET_H
#define DEJAVU_REFINEMENT_BUCKET_H


#include <list>
#include "sgraph.h"
#include "invariant.h"
#include "refinement.h"
#include "nauty/nausparse.h"
#include "nauty_refine.h"

class refinement_bucket {
public:
    bool refine_coloring(sgraph* g, coloring_bucket* c, invariant* I, int level);
    void individualize_vertex(sgraph* g, coloring_bucket* c, int cell, int v, int level);
    void undo_individualize_vertex(sgraph *g, coloring *c, int v);
    bool refine_color_class(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair* color_class_split_worklist, invariant* I, int* largest_color_class_index);
    void undo_refine_color_class(sgraph *g, coloring *c, std::list<std::pair<int, int>> *changes);
    void complete_colorclass_invariant(sgraph *g, coloring *c, invariant_acc *I);
    bool assert_is_equitable(sgraph *g, coloring *c);
    ~refinement_bucket();
    void initialize_active(sgraph *g);
private:
    bool has_nauty_graph = false;
    sparsegraph* sg;
    int m;
    int n;
    size_t active_sz=0;

    workspace W;
    set*   active;

    void assure_graph(sgraph *g);
};


#endif //DEJAVU_REFINEMENT_BUCKET_H
