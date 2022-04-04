#ifndef DEJAVU_PREP_H
#define DEJAVU_PREP_H

#include "sgraph.h"
#include "refinement.h"
#include "schreier_shared.h"

class sparse_automorphism {
private:
    int cnt_supp = 0;
    int support  = 0;
    std::vector<int> perm_from;
    std::vector<int> perm_to;
public:
    sparse_automorphism(int supp) {
        support = supp;
        perm_from.reserve(supp);
        perm_to.reserve(supp);
    }
    void add_map(int from, int to) {
        assert(cnt_supp < support);
        perm_from.push_back(from);
        perm_to.push_back(to);
        ++cnt_supp;
    }
};

class preprocessor {
private:
    std::vector<int> orig_graph_map;

    void red_singleton_only_refinement() {

    }

    void res_singleton_only_refinement() {

    }

    void red_low_deg(sgraph* g, coloring<int>* c) {
        mark_set del;
        del.initialize(g->v_size);

        work_list_t<int> worklist_deg0;
        mark_set deg0;
        work_list_t<int> worklist_deg1;
        mark_set deg1;
        work_list_t<int> worklist_deg2;
        mark_set deg2;

        worklist_deg0.initialize(g->v_size);
        deg0.initialize(g->v_size);
        worklist_deg1.initialize(g->v_size);
        deg1.initialize(g->v_size);
        worklist_deg2.initialize(g->v_size);
        deg2.initialize(g->v_size);

        for(int i = 0; i < g->v_size; ++i) {
            if(g->d[i] == 0) {
                worklist_deg0.push_back(i);
                deg0.set(i);
            }
            if(g->d[i] == 1) {
                worklist_deg1.push_back(i);
                deg1.set(i);
            }
            if(g->d[i] == 2) {
                worklist_deg2.push_back(i);
                deg2.set(i);
            }
        }

        //while(!worklist_deg1.empty()) {
            // reverse-split c according to deg1 vertices, such that they form own color classes, collect these color classes
            // for each of these colors X, split parents according to number of X vertices
            // remove the X vertices from graph, save to which parent vertices of X were connected (for each parent save list of vertices) -- just keep the worklist?
            // restoration: each generator that permutes parents is extended to also permute X_1 -> X_2
            // restoration: add generators for symmetrical group for vertices of each parent
        //}
    }

    void res_low_deg() {

    }
public:
    void reduce(sgraph* g, int* colmap) {
        coloring<int> c;
        // singleton-only refinement, then cut graph
        g->initialize_coloring(&c, colmap);
        red_low_deg(g, &c);
    }

    void restore(sgraph* g, automorphism_info* a) {

    }
};

#endif //DEJAVU_PREP_H
