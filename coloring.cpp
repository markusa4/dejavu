//
// Created by markus on 19.09.19.
//

#include "coloring.h"

void coloring::rewrite_ptn(coloring *c) {
    for(int i = 0; i < c->ptn_sz; i += 1) {
        ptn[i] = c->ptn[i];
    }
}

coloring::~coloring() {
    if(init) {
        delete[] ptn;
        delete[] lab;
    }
}
