// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_COLORING_H
#define DEJAVU_COLORING_H

#include <vector>
#include <cassert>
#include "utility.h"

/**
 * \brief Vertex coloring for a graph
 *
 * Stores a vertex coloring for a graph with \a domain_size many vertices. Has mappings to and from colorings. The
 * datastructure mostly follows the design used by Traces, with some slight adjustments. The class mainly exposes
 * raw arrays, to be used by the solver internally.
 */
class coloring {
public:
    int* lab = nullptr;
    int* ptn = nullptr;
    int* vertex_to_col = nullptr;
    int* vertex_to_lab = nullptr;

    int cells = 1;
    int domain_size = 0;

private:
    int* alloc_pt = nullptr;

public:

    void alloc(int sz) {
        assert(sz >= 0);

        if (alloc_pt) dealloc();
        alloc_pt = (int*) malloc(sz * 4 * sizeof(int));
        assert(alloc_pt != nullptr);

        lab = alloc_pt;
        ptn = lab + sz;
        vertex_to_col = lab + sz * 2;
        vertex_to_lab = lab + sz * 3;

        domain_size = sz;
    }

    void dealloc() {
        if (alloc_pt) free(alloc_pt);
        lab           = nullptr;
        ptn           = nullptr;
        vertex_to_col = nullptr;
        vertex_to_lab = nullptr;
    };
    ~coloring() {
        if(alloc_pt) {
            dealloc();
        }
    }

    void copy_ptn(coloring *c) const {
        assert(alloc_pt);
        assert(c->alloc_pt);
        memcpy(ptn, c->ptn, c->domain_size*sizeof(int));
    }

    void copy_from_ir_ancestor(coloring *c) {
        if(alloc_pt) {
            if(domain_size != c->domain_size) {
                dealloc();
            } else {
                cells = c->cells;
                memcpy(vertex_to_col, c->vertex_to_col, c->domain_size * sizeof(int));
                return;
            }
        }

        if(!alloc_pt) alloc(c->domain_size);

        if(c->cells > c->domain_size / 4) {
            memcpy(ptn, c->ptn, c->domain_size * sizeof(int));
        } else {
            for (int i = 0; i < c->domain_size;) {
                const int rd = c->ptn[i];
                ptn[i] = rd;
                i += rd + 1;
            }
        }
        memcpy(lab, c->lab, c->domain_size*sizeof(int));
        memcpy(vertex_to_col, c->vertex_to_col, c->domain_size*sizeof(int));
        memcpy(vertex_to_lab, c->vertex_to_lab, c->domain_size*sizeof(int));

        domain_size = c->domain_size;
        cells = c->cells;
    }

    /**
     * Copies the given coloring into this coloring.
     * 
     * @param c The coloring to copy from.
     */
    void copy_any(coloring *c) {
        if(alloc_pt) {
            if(domain_size != c->domain_size) dealloc();
        }

        if(!alloc_pt) {
            alloc(c->domain_size);
        }

        if(c->cells > c->domain_size / 4) {
            memcpy(ptn, c->ptn, c->domain_size * sizeof(int));
        } else {
            for (int i = 0; i < c->domain_size;) {
                const int rd = c->ptn[i];
                ptn[i] = rd;
                i += rd + 1;
            }
        }
        memcpy(lab, c->lab, c->domain_size*sizeof(int));
        memcpy(vertex_to_col, c->vertex_to_col, c->domain_size*sizeof(int));
        memcpy(vertex_to_lab, c->vertex_to_lab, c->domain_size*sizeof(int));

        domain_size = c->domain_size;
        cells = c->cells;
    }

    void initialize(int new_domain_size) {
        alloc(new_domain_size);
    }

    [[maybe_unused]] void check() const {
        bool comp = true;

        for(int i = 0; i < domain_size;++i) {
            comp = comp && (lab[i] >= 0 && lab[i] < domain_size);
            comp = comp && (lab[vertex_to_lab[i]] == i);
        }

        [[maybe_unused]] int last_col = -1;
        int counter  = 1;
        for (int i = 0; i < domain_size; ++i) {
            --counter;
            if(counter == 0) {
                counter = ptn[i] + 1;
                last_col = i;
                assert(ptn[i] >= 0 && ptn[i] < domain_size);
            } else {
                assert(vertex_to_col[lab[i]] == last_col);
            }
        }

        for(int i = 0; i < domain_size;) {
            assert(vertex_to_col[lab[i]] == i);
            i += ptn[i] + 1;
        }
        assert(comp);
    }
};

#endif //DEJAVU_COLORING_H
