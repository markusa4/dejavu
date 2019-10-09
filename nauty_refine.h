//
// Created by markus on 09.10.19.
//

#ifndef DEJAVU_NAUTY_REFINE_H
#define DEJAVU_NAUTY_REFINE_H

#include "invariant.h"

extern "C" {
#include "nauty/nausparse.h"
}

#define DYNALLSTAT_NOSTATIC(type,name,name_sz) \
	 TLS_ATTR type *name=NULL; TLS_ATTR size_t name_sz=0

struct workspace {
    DYNALLSTAT_NOSTATIC(int,work1,work1_sz);
    DYNALLSTAT_NOSTATIC(int,work2,work2_sz);
    DYNALLSTAT_NOSTATIC(int,work3,work3_sz);
    DYNALLSTAT_NOSTATIC(int,work4,work4_sz);
    DYNALLSTAT_NOSTATIC(set,snwork,snwork_sz);
    DYNALLSTAT_NOSTATIC(short,vmark1,vmark1_sz);
    DYNALLSTAT_NOSTATIC(short,vmark2,vmark2_sz);
    TLS_ATTR short vmark1_val = 32000;
    TLS_ATTR short vmark2_val = 32000;
};

bool dynamic_refine_sg(graph *g, int *lab, int *ptn, int level, int *numcells,
                   int *count, set *active, int *code, int m, int n, invariant* I, workspace* w);

#endif //DEJAVU_NAUTY_REFINE_H
