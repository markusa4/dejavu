//
// Created by markus on 09.10.19.
//

#include <assert.h>
#include "refinement_bucket.h"
#include "nauty_refine.h"
extern "C" {
#include "nauty/nausparse.h";
}

#define DYNALLSTAT_NOSTATIC(type,name,name_sz) \
	 TLS_ATTR type *name; TLS_ATTR size_t name_sz=0


 void refinement_bucket::assure_graph(sgraph *g) {
     if(!has_nauty_graph) {
         sg = new sparsegraph;
         SG_INIT(*sg);
         SG_ALLOC(*sg, g->v.size(), g->e.size(), "malloc");
         sg->nv = g->v.size();
         sg->nde = g->e.size();
         for (int i = 0; i < g->v.size(); ++i) {
             sg->v[i] = g->v[i];
             sg->d[i] = g->d[i];
         }
         for (int i = 0; i < g->e.size(); ++i) {
             sg->e[i] = g->e[i];
         }
         n = sg->nv;
         m = SETWORDSNEEDED(n);
         has_nauty_graph = true;
         DYNALLOC1(set,active,active_sz,m,"nauty");
     }
 }

 void refinement_bucket::initialize_active(sgraph *g) {
     assure_graph(g);
     EMPTYSET(active,m);
     ADDELEMENT(active, 0);
 }

bool refinement_bucket::refine_coloring(sgraph *g, coloring_bucket *c, invariant *I, int level) {
    assure_graph(g);

    int code = 0;
    int count;
    //extern void refine_sg(graph*,int*,int*,int,int*,int*,set*,int*,int,int);
    bool comp = dynamic_refine_sg((graph*)sg, c->lab, c->ptn, level, &c->numcells, &count, active, &code, m, n, I, &W);
    return comp && I->write_top_and_compare(code);
}

void refinement_bucket::individualize_vertex(sgraph *g, coloring_bucket *c, int cell, int v, int level) {
    assure_graph(g);
    int i,prev,next;

    EMPTYSET(active,m);
    ADDELEMENT(active, cell);

    i = cell;
    prev = v;

    do {
        assert(i < c->lab_sz);
        next = c->lab[i];
        c->lab[i++] = prev;
        prev = next;
    } while(prev != v);

    c->ptn[cell] = level;
}

refinement_bucket::~refinement_bucket() {
     if(has_nauty_graph) {
         SG_FREE(*sg);
         delete sg;
     }
}
