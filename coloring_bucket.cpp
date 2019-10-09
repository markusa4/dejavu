//
// Created by markus on 09.10.19.
//

#include <iostream>
#include "coloring_bucket.h"

coloring_bucket::~coloring_bucket() {
    if(init) {
        delete[] lab;
        delete[] ptn;
    }
}

coloring_bucket::coloring_bucket() {
    init = false;
}

void coloring_bucket::copy(coloring_bucket *c) {
    if(init) {
        delete[] lab;
        delete[] ptn;
    }

    lab = new int[c->lab_sz];
    ptn = new int[c->ptn_sz];

    for(int i = 0; i < c->lab_sz; ++i) {
        lab[i] = c->lab[i];
        ptn[i] = c->ptn[i];
    }
    lab_sz = c->lab_sz;
    ptn_sz = c->ptn_sz;
    numcells = c->numcells;
    init = true;
}

void coloring_bucket::print() {
    if(init) {
        for (int i = 0; i < ptn_sz; ++i) {
            std::cout << ptn[i] << " ";
        }
        std::cout << std::endl;
    }
}
