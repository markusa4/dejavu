#include "coloring.h"
#include <cstring>
#include <iostream>

std::pair<int*, int*> coloring_bulk_allocator(int domain_size) {
    thread_local int* bulk_domain = nullptr;
    thread_local int bulk_domain_sz = -1, bulk_domain_cnt = -1;

    if(bulk_domain_sz < 0) {
        bulk_domain     = new int[20 * domain_size + 1];
        bulk_domain[0]  = 1;
        bulk_domain_sz  = 20 * domain_size + 1;
        bulk_domain_cnt = 1;
    }

    bulk_domain_cnt += domain_size;
    if(bulk_domain_cnt == bulk_domain_sz)
        bulk_domain_sz = -1;
    else
        bulk_domain[0]  += 1;


    return std::pair<int*, int*>(bulk_domain, bulk_domain + bulk_domain_cnt - domain_size);
}

void coloring_bulk_deallocator(int* bulk_domain) {
    if(--bulk_domain[0] == 0)
        delete[] bulk_domain;
}

coloring::~coloring() {
    if(init) {
        dealloc();
    }
}

void coloring::dealloc() {
    if(!efficient_alloc) {
        delete[] ptn;
        delete[] lab;
        delete[] vertex_to_lab;
        delete[] vertex_to_col;
    } else {
        coloring_bulk_deallocator(bulk_alloc);
    }
};

void coloring::copy(coloring *c) {
    if(init) {
        if(lab_sz != c->lab_sz || ptn_sz != c->ptn_sz) {
            dealloc();
            init = false;
        } else {
             memcpy(ptn, c->ptn, c->ptn_sz*sizeof(int));
             memcpy(vertex_to_col, c->vertex_to_col, c->ptn_sz*sizeof(int));
             return;
        }
    }

    if(!init) {
        std::pair<int*, int*> alloc = coloring_bulk_allocator(c->lab_sz * 4);
        bulk_alloc = alloc.first;
        bulk_pt    = alloc.second;

        lab           = bulk_pt;
        ptn           = lab + c->lab_sz;
        vertex_to_col = lab + c->lab_sz * 2;
        vertex_to_lab = lab + c->lab_sz * 3;
        efficient_alloc = true;
    }

    // color_choices = c->color_choices;

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
            dealloc();
            init = false;
        }
    }

    if(!init) {
        std::pair<int*, int*> alloc = coloring_bulk_allocator(c->lab_sz * 4);
        bulk_alloc = alloc.first;
        bulk_pt    = alloc.second;

        lab           = bulk_pt;
        ptn           = lab + c->lab_sz;
        vertex_to_col = lab + c->lab_sz * 2;
        vertex_to_lab = lab + c->lab_sz * 3;
        efficient_alloc = true;
    }

    // color_choices = c->color_choices;

    memcpy(lab, c->lab, c->lab_sz*sizeof(int));
    memcpy(ptn, c->ptn, c->ptn_sz*sizeof(int));
    memcpy(vertex_to_col, c->vertex_to_col, c->lab_sz*sizeof(int));
    memcpy(vertex_to_lab, c->vertex_to_lab, c->lab_sz*sizeof(int));

    lab_sz = c->lab_sz;
    ptn_sz = c->ptn_sz;

    init = true;
}

void coloring::initialize(int domain_size) {
    std::pair<int*, int*> alloc = coloring_bulk_allocator(domain_size * 4);
    bulk_alloc = alloc.first;
    bulk_pt    = alloc.second;

    lab           = bulk_pt;
    ptn           = lab + domain_size;
    vertex_to_col = lab + domain_size * 2;
    vertex_to_lab = lab + domain_size * 3;
    efficient_alloc = true;
    init = true;

    lab_sz = domain_size;
    ptn_sz = domain_size;

    // color_choices.clear();
    // color_choices.reserve(16);
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