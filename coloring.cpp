#include "coloring.h"
#include <cstring>

coloring::~coloring() {
    if(init) {
        delete[] ptn;
        delete[] lab;
        delete[] vertex_to_lab;
        delete[] vertex_to_col;
    }
}

void coloring::copy(coloring *c) {
    if(init) {
        if(lab_sz != c->lab_sz || ptn_sz != c->ptn_sz) {
            delete[] lab;
            delete[] ptn;
            delete[] vertex_to_lab;
            delete[] vertex_to_col;
            init = false;
        } else {
             memcpy(ptn, c->ptn, c->ptn_sz*sizeof(int));
             memcpy(vertex_to_col, c->vertex_to_col, c->lab_sz*sizeof(int));
             return;
        }
    }

    if(!init) {
        lab = new int[c->lab_sz];
        ptn = new int[c->ptn_sz];
        vertex_to_col = new int[c->lab_sz];
        vertex_to_lab = new int[c->lab_sz];
    }

    color_choices = c->color_choices;

    memcpy(lab, c->lab, c->lab_sz*sizeof(int));
    memcpy(ptn, c->ptn, c->ptn_sz*sizeof(int));
    memcpy(vertex_to_col, c->vertex_to_col, c->lab_sz*sizeof(int));
    memcpy(vertex_to_lab, c->vertex_to_lab, c->lab_sz*sizeof(int));

    lab_sz = c->lab_sz;
    ptn_sz = c->ptn_sz;

    init = true;
}

void coloring::copy_force(coloring *c) {
    if(init) {
        if(lab_sz != c->lab_sz || ptn_sz != c->ptn_sz) {
            delete[] lab;
            delete[] ptn;
            delete[] vertex_to_lab;
            delete[] vertex_to_col;
            init = false;
        }
    }

    if(!init) {
        lab = new int[c->lab_sz];
        ptn = new int[c->ptn_sz];
        vertex_to_col = new int[c->lab_sz];
        vertex_to_lab = new int[c->lab_sz];
    }

    color_choices = c->color_choices;

    memcpy(lab, c->lab, c->lab_sz*sizeof(int));
    memcpy(ptn, c->ptn, c->ptn_sz*sizeof(int));
    memcpy(vertex_to_col, c->vertex_to_col, c->lab_sz*sizeof(int));
    memcpy(vertex_to_lab, c->vertex_to_lab, c->lab_sz*sizeof(int));

    lab_sz = c->lab_sz;
    ptn_sz = c->ptn_sz;

    init = true;
}

void coloring::initialize(int domain_size) {
    lab = new int[domain_size];
    ptn = new int[domain_size];
    vertex_to_col = new int[domain_size];
    vertex_to_lab = new int[domain_size];

    init = true;

    lab_sz = domain_size;
    ptn_sz = domain_size;

    color_choices.clear();
    color_choices.reserve(16);
}

bool coloring::check() {
    bool comp = true;

    for(int i = 0; i < lab_sz;++i) {
        comp = comp && (lab[i] >= 0 && lab[i] < lab_sz);
        comp = comp && (lab[vertex_to_lab[i]] == i);
    }

    // assert proper ptn
    for(int i = 0; i < lab_sz;) {
        if(i != 0)
            comp = comp && (ptn[i - 1] == 0);
        i += ptn[i] + 1;
    }

    return comp;
}