#ifndef DEJAVU_BIJECTION_H
#define DEJAVU_BIJECTION_H

#include "coloring.h"
#include <algorithm>
#include <chrono>
#include <random>
#include <iostream>


template<class vertex_t>
class bijection {
public:
    bool init = false;
    bool mark = false;
    vertex_t* map;
    int map_sz;
    bool non_uniform  = false;
    bool foreign_base = false;
    bool certified    = false;

    int map_vertex(int v) {
        return map[v];
    }

    void copy(bijection* p) {
        if(init)
            delete[] map;
        init = p->init;
        mark = p->mark;
        map_sz = p->map_sz;
        non_uniform = p->non_uniform;
        foreign_base = p->foreign_base;
        certified = p->certified;
        memcpy(map, p->map, p->map_sz * sizeof(vertex_t));
    }

    void print() {
        for(int i = 0; i < map_sz; ++i)
            std::cout << map[i] << " ";
        std::cout << std::endl;
    }

    void read_from_coloring(coloring<vertex_t> *c) {
        if(!init) {
            //    delete[] map;
            //}
            map = new vertex_t[c->lab_sz];
        }
        init = true;
        map_sz = c->lab_sz;
        for(int i = 0; i < c->lab_sz; ++i) {
            map[i] = c->lab[i];
        }
    }

    void inverse() {
        thread_local vertex_t* switch_map;
        thread_local bool         switch_map_init;

        if(!switch_map_init) {
            switch_map_init = true;
            switch_map      = new vertex_t[map_sz];
        }

        vertex_t* swap = map;
        map = switch_map;
        switch_map = swap;

        for(int i = 0; i < map_sz; ++i) {
            map[switch_map[i]] = i;
        }
    }

    void compose(bijection<vertex_t>* p) {
        for(int i = 0; i < map_sz; ++i) {
            map[i] = p->map[map[i]];
        }
    }

    void not_deletable() {
        init = false;
    }


    void deletable() {
        init = true;
    }

    bijection() {
        init = false;
    }

    ~bijection() {
        if(init)
            delete[] map;
    }

    static void random_bijection(bijection<vertex_t>* p, int n, unsigned seed) {
        p->map = new vertex_t[n];
        p->init = true;
        p->map_sz = n;
        for(int i = 0; i < n; ++i) {
            p->map[i] = i;
        }
        std::default_random_engine re = std::default_random_engine(seed);
        std::shuffle(p->map, p->map + p->map_sz, re);
    }
};

#endif //DEJAVU_BIJECTION_H
