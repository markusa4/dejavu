#ifndef DEJAVU_COLORING_H
#define DEJAVU_COLORING_H


#include <vector>
#include <cstring>
#include <iostream>
#include "utility.h"

extern thread_local bool bulk_domain_reset; // ToDo: parallelization messes this up, can create dangling pointer

template<class vertex_t>
class garbage_collector {
    static std::vector<vertex_t*> junk;
    static std::unique_ptr<std::mutex> lock;
public:
    static void free_trash() {
        lock->lock();
        while(!junk.empty()) {
            delete[] junk.back();
            junk.pop_back();
        }
        lock->unlock();
    }
    static void add_trash(vertex_t* arr) {
        lock->lock();
        junk.push_back(arr);
        lock->unlock();
    }
};

template<class vertex_t> std::unique_ptr<std::mutex> garbage_collector<vertex_t>::lock = std::make_unique<std::mutex>();
template<class vertex_t> std::vector<vertex_t*> garbage_collector<vertex_t>::junk = std::vector<vertex_t*>();

template<class vertex_t>
class coloring_allocator {
public:
    static std::pair<vertex_t *, vertex_t *> coloring_bulk_allocator(int domain_size) {
        thread_local vertex_t *bulk_domain = nullptr;
        thread_local int buffer_const = 20;
        thread_local int bulk_domain_sz = -1, bulk_domain_cnt = -1;

        if (bulk_domain_sz < 0 || bulk_domain_reset) {
            if (bulk_domain_reset)
                buffer_const = 20;
            bulk_domain_reset = false;
            bulk_domain = new vertex_t[buffer_const * domain_size + 1];
            garbage_collector<vertex_t>::add_trash(bulk_domain);
            bulk_domain[0] = 0;
            bulk_domain_sz = buffer_const * domain_size + 1;
            bulk_domain_cnt = 1;
            buffer_const *= 2;
        }

        bulk_domain_cnt += domain_size;
        bulk_domain[0] += 1;
        if (bulk_domain_cnt >= bulk_domain_sz)
            bulk_domain_sz = -1;

        return std::pair<vertex_t *, vertex_t *>(bulk_domain, bulk_domain + bulk_domain_cnt - domain_size);
    }

    static void coloring_bulk_deallocator(vertex_t *bulk_domain) {
        bulk_domain_reset = true;
        /*if (--bulk_domain[0] == 0) {
            std::cout << "removeB " << bulk_domain << std::endl;
            bulk_domain_reset = true;
            delete[] bulk_domain;
        }*/
    }
};

template <class vertex_t>
class coloring {
public:
    vertex_t* bulk_alloc;
    vertex_t* bulk_pt;

    vertex_t* lab;
    vertex_t* ptn;

    int lab_sz;
    int ptn_sz;
    bool init = false;
    bool efficient_alloc = false;
    vertex_t* vertex_to_col;
    vertex_t* vertex_to_lab;

    int cells = 1;
    int smallest_cell_lower_bound = INT32_MAX;

    ~coloring() {
        if(init) {
            dealloc();
        }
    }

    void alloc(int sz) {
        if(!init) {
            std::pair<vertex_t*, vertex_t*> alloc = coloring_allocator<vertex_t>::coloring_bulk_allocator(sz * 4);
            bulk_alloc = alloc.first;
            bulk_pt    = alloc.second;

            lab           = bulk_pt;
            ptn           = lab + sz;
            vertex_to_col = lab + sz * 2;
            vertex_to_lab = lab + sz * 3;
            efficient_alloc = true;
            init = true;

            lab_sz = sz;
            ptn_sz = sz;
        }
        /*if(!init) {
            lab           = new vertex_t[sz];
            ptn           = new vertex_t[sz];
            vertex_to_col = new vertex_t[sz];
            vertex_to_lab = new vertex_t[sz];
            efficient_alloc = false;
            init = true;

            lab_sz = sz;
            ptn_sz = sz;
        }*/
    }

    void dealloc() {
        if(!efficient_alloc) {
            delete[] ptn;
            delete[] lab;
            delete[] vertex_to_lab;
            delete[] vertex_to_col;
        } else {
            coloring_allocator<vertex_t>::coloring_bulk_deallocator(bulk_alloc);
        }
    };

    void copy(coloring<vertex_t> *c) {
        if(init) {
            if(lab_sz != c->lab_sz || ptn_sz != c->ptn_sz) {
                dealloc();
                init = false;
            } else {
                cells = c->cells;
                if(!efficient_alloc || !c->efficient_alloc) {
                    for(int i = 0; i < c->ptn_sz;) {
                        const vertex_t rd = c->ptn[i];
                        ptn[i] = rd;
                        i += rd +1;
                    }
                    memcpy(vertex_to_col, c->vertex_to_col, c->ptn_sz * sizeof(vertex_t));
                } else {
                    for(int i = 0; i < c->ptn_sz;) {
                        const vertex_t rd = c->ptn[i];
                        ptn[i] = rd;
                        i += rd +1;
                    }
                    memcpy(vertex_to_col, c->vertex_to_col, c->ptn_sz * sizeof(vertex_t));
                }
                return;
            }
        }

        if(!init) {
            alloc(c->lab_sz);
        }

        if(c->cells > c->ptn_sz / 4) {
            memcpy(ptn, c->ptn, c->ptn_sz * sizeof(vertex_t));
        } else {
            for (int i = 0; i < c->ptn_sz;) {
                const vertex_t rd = c->ptn[i];
                ptn[i] = rd;
                i += rd + 1;
            }
        }
        memcpy(lab, c->lab, c->lab_sz*sizeof(vertex_t));
        memcpy(vertex_to_col, c->vertex_to_col, c->lab_sz*sizeof(vertex_t));
        memcpy(vertex_to_lab, c->vertex_to_lab, c->lab_sz*sizeof(vertex_t));

        lab_sz = c->lab_sz;
        ptn_sz = c->ptn_sz;

        cells = c->cells;
        smallest_cell_lower_bound = c->smallest_cell_lower_bound;
        init = true;
    }

    void copy_force(coloring<vertex_t> *c) {
        if(init) {
            if(lab_sz != c->lab_sz || ptn_sz != c->ptn_sz) {
                dealloc();
                init = false;
            }
        }

        if(!init) {
            alloc(c->lab_sz);
        }

        if(c->cells > c->ptn_sz / 4) {
            memcpy(ptn, c->ptn, c->ptn_sz * sizeof(vertex_t));
        } else {
            for (int i = 0; i < c->ptn_sz;) {
                const vertex_t rd = c->ptn[i];
                ptn[i] = rd;
                i += rd + 1;
            }
        }
        memcpy(lab, c->lab, c->lab_sz*sizeof(vertex_t));
        memcpy(vertex_to_col, c->vertex_to_col, c->lab_sz*sizeof(vertex_t));
        memcpy(vertex_to_lab, c->vertex_to_lab, c->lab_sz*sizeof(vertex_t));

        lab_sz = c->lab_sz;
        ptn_sz = c->ptn_sz;

        cells = c->cells;
        smallest_cell_lower_bound = c->smallest_cell_lower_bound;
        init = true;
    }

    void initialize(int domain_size) {
        alloc(domain_size);
    }

    bool check() {
        bool comp = true;

        for(int i = 0; i < lab_sz;++i) {
            comp = comp && (lab[i] >= 0 && lab[i] < lab_sz);
            comp = comp && (lab[vertex_to_lab[i]] == i);
        }
        return comp;
    }
};

#endif //DEJAVU_COLORING_H
