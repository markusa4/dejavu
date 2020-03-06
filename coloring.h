#ifndef DEJAVU_COLORING_H
#define DEJAVU_COLORING_H


#include <vector>
#include <cstring>
#include <iostream>

template <class vertex_type>
std::pair<vertex_type*, vertex_type*> coloring_bulk_allocator(int domain_size) {
    thread_local vertex_type* bulk_domain = nullptr;
    thread_local int  buffer_const = 20;
    thread_local int  bulk_domain_sz = -1, bulk_domain_cnt = -1;

    if(bulk_domain_sz < 0) {
        bulk_domain     = new vertex_type[buffer_const * domain_size + 1];
        bulk_domain[0]  = 1;
        bulk_domain_sz  = buffer_const * domain_size + 1;
        bulk_domain_cnt = 1;
        buffer_const *= 2;
    }

    bulk_domain_cnt += domain_size;
    if(bulk_domain_cnt == bulk_domain_sz)
        bulk_domain_sz = -1;
    else
        bulk_domain[0]  += 1;


    return std::pair<vertex_type*, vertex_type*>(bulk_domain, bulk_domain + bulk_domain_cnt - domain_size);
}

template <class vertex_type>
void coloring_bulk_deallocator(vertex_type* bulk_domain) {
    if(--bulk_domain[0] == 0)
        delete[] bulk_domain;
}

template <class vertex_type>
class alignas(16) coloring_temp {
public:
    vertex_type* bulk_alloc;
    vertex_type* bulk_pt;

    vertex_type* lab;
    vertex_type* ptn;

    int lab_sz;
    int ptn_sz;
    bool init = false;
    bool efficient_alloc = false;
    vertex_type* vertex_to_col;
    vertex_type* vertex_to_lab;

    int cells = 1;
    int smallest_cell_lower_bound = INT32_MAX;

    ~coloring_temp() {
        if(init) {
            dealloc();
        }
    }

    void dealloc() {
        if(!efficient_alloc) {
            delete[] ptn;
            delete[] lab;
            delete[] vertex_to_lab;
            delete[] vertex_to_col;
        } else {
            coloring_bulk_deallocator<vertex_type>(bulk_alloc);
        }
    };

    void copy(coloring_temp<vertex_type> *c) {
        if(init) {
        if(lab_sz != c->lab_sz || ptn_sz != c->ptn_sz) {
            dealloc();
            init = false;
        } else {
            cells = c->cells;
            if(!efficient_alloc || !c->efficient_alloc) {
                for(int i = 0; i < c->ptn_sz;) {
                    const vertex_type rd = c->ptn[i];
                    ptn[i] = rd;
                    i += rd +1;
                }
                // memcpy(ptn, c->ptn, c->ptn_sz * sizeof(vertex_type));
                memcpy(vertex_to_col, c->vertex_to_col, c->ptn_sz * sizeof(vertex_type));
            } else {
                for(int i = 0; i < c->ptn_sz;) {
                    const vertex_type rd = c->ptn[i];
                    ptn[i] = rd;
                    i += rd +1;
                }
                // memcpy(ptn, c->ptn, c->ptn_sz * sizeof(vertex_type));
                memcpy(vertex_to_col, c->vertex_to_col, c->ptn_sz * sizeof(vertex_type));
            }
            return;
        }
        }

        if(!init) {
            std::pair<vertex_type*, vertex_type*> alloc = coloring_bulk_allocator<vertex_type>(c->lab_sz * 4);
            bulk_alloc = alloc.first;
            bulk_pt    = alloc.second;

            lab           = bulk_pt;
            ptn           = lab + c->lab_sz;
            vertex_to_col = lab + c->lab_sz * 2;
            vertex_to_lab = lab + c->lab_sz * 3;
            efficient_alloc = true;
        }

        if(c->cells > c->ptn_sz / 4) {
            memcpy(ptn, c->ptn, c->ptn_sz * sizeof(vertex_type));
        } else {
            for (int i = 0; i < c->ptn_sz;) {
                const vertex_type rd = c->ptn[i];
                ptn[i] = rd;
                i += rd + 1;
            }
        }
        // memcpy(ptn, c->ptn, c->ptn_sz * sizeof(vertex_type));
        memcpy(lab, c->lab, c->lab_sz*sizeof(vertex_type));
        memcpy(vertex_to_col, c->vertex_to_col, c->lab_sz*sizeof(vertex_type));
        memcpy(vertex_to_lab, c->vertex_to_lab, c->lab_sz*sizeof(vertex_type));

        lab_sz = c->lab_sz;
        ptn_sz = c->ptn_sz;

        cells = c->cells;
        init = true;
    }

    void copy_force(coloring_temp<vertex_type> *c) {
        if(init) {
            if(lab_sz != c->lab_sz || ptn_sz != c->ptn_sz) {
                dealloc();
                init = false;
            }
        }

        if(!init) {
        std::pair<vertex_type*, vertex_type*> alloc = coloring_bulk_allocator<vertex_type>(c->lab_sz * 4);
        bulk_alloc = alloc.first;
        bulk_pt    = alloc.second;

        lab           = bulk_pt;
        ptn           = lab + c->lab_sz;
        vertex_to_col = lab + c->lab_sz * 2;
        vertex_to_lab = lab + c->lab_sz * 3;
        efficient_alloc = true;
        }

        if(c->cells > c->ptn_sz / 4) {
            memcpy(ptn, c->ptn, c->ptn_sz * sizeof(vertex_type));
        } else {
            for (int i = 0; i < c->ptn_sz;) {
                const vertex_type rd = c->ptn[i];
                ptn[i] = rd;
                i += rd + 1;
            }
        }
        memcpy(lab, c->lab, c->lab_sz*sizeof(vertex_type));
        memcpy(vertex_to_col, c->vertex_to_col, c->lab_sz*sizeof(vertex_type));
        memcpy(vertex_to_lab, c->vertex_to_lab, c->lab_sz*sizeof(vertex_type));

        lab_sz = c->lab_sz;
        ptn_sz = c->ptn_sz;

        cells = c->cells;

        init = true;
    }

    void initialize(int domain_size) {
        std::pair<vertex_type*, vertex_type*> alloc = coloring_bulk_allocator<vertex_type>(domain_size * 4);
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

    bool check() {
        bool comp = true;

        for(int i = 0; i < lab_sz;++i) {
            comp = comp && (lab[i] >= 0 && lab[i] < lab_sz);
            comp = comp && (lab[vertex_to_lab[i]] == i);
        }

        // assert proper ptn
        /*for(int i = 0; i < lab_sz;) {
            if(i != 0)
                comp = comp && (ptn[i - 1] == 0);
            i += ptn[i] + 1;
        }*/

        return comp;
    }
};

// typedef coloring_temp<int> coloring;

#endif //DEJAVU_COLORING_H
