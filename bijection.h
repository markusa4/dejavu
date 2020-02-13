#ifndef DEJAVU_BIJECTION_H
#define DEJAVU_BIJECTION_H

#include "coloring.h"
#include <algorithm>
#include <chrono>
#include <random>
#include <iostream>

extern thread_local int* switch_map;
extern thread_local bool switch_map_init;

template<class vertex_type>
class bijection_temp {
public:
    bool init = false;
    bool mark = false;
    vertex_type* map;
    int map_sz;
    bool non_uniform  = false;
    bool foreign_base = false;
    bool certified    = false;

    int map_vertex(int v) {
        return map[v];
    }

    void print() {
        for(int i = 0; i < map_sz; ++i)
            std::cout << map[i] << " ";
        std::cout << std::endl;
    }

    void read_from_coloring(coloring *c) {
        if(!init) {
            //    delete[] map;
            //}
            map = new vertex_type[c->lab_sz];
        }
        init = true;
        map_sz = c->lab_sz;
        for(int i = 0; i < c->lab_sz; ++i) {
            map[i] = c->lab[i];
        }
    }

    void inverse() {
        if(!switch_map_init) {
            switch_map_init = true;
            switch_map      = new vertex_type[map_sz];
        }

        int* swap = map;
        map = switch_map;
        switch_map = swap;

        for(int i = 0; i < map_sz; ++i) {
            map[switch_map[i]] = i;
        }
    }

    void compose(bijection_temp<vertex_type>* p) {
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

    bijection_temp() {
        init = false;
    }

    ~bijection_temp() {
        if(init)
            delete[] map;
    }

    static void random_bijection(bijection_temp<vertex_type>* p, int n, unsigned seed) {
        p->map = new vertex_type[n];
        p->init = true;
        p->map_sz = n;
        for(int i = 0; i < n; ++i) {
            p->map[i] = i;
        }
        std::default_random_engine re = std::default_random_engine(seed);
        std::shuffle(p->map, p->map + p->map_sz, re);
    }
};

typedef  bijection_temp<int> bijection;


#endif //DEJAVU_BIJECTION_H
