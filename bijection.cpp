#include <algorithm>
#include <chrono>
#include <random>
#include <iostream>
#include "bijection.h"

int bijection::map_vertex(int v) {
    return map[v];
}

void bijection::print() {
    for(int i = 0; i < map_sz; ++i)
        std::cout << map[i] << " ";
    std::cout << std::endl;
}

void bijection::read_from_coloring(coloring *c) {
    if(init) {
        delete[] map;
    }
    map = new int[c->lab_sz];
    init = true;
    map_sz = c->lab_sz;
    for(int i = 0; i < c->lab_sz; ++i) {
        map[i] = c->lab[i];
    }
}

void bijection::inverse() {
    int* old_map = map;
    map = new int[map_sz];
    for(int i = 0; i < map_sz; ++i) {
        map[old_map[i]] = i;
    }
    delete[] old_map;
}

void bijection::compose(bijection* p) {
    for(int i = 0; i < map_sz; ++i) {
        map[i] = p->map[map[i]];
    }
}

void bijection::not_deletable() {
    init = false;
}


void bijection::deletable() {
    init = true;
}

bijection::bijection() {
    init = false;
}

bijection::~bijection() {
    if(init)
        delete[] map;
}

void bijection::random_bijection(bijection* p, int n) {
    p->map = new int[n];
    p->init = true;
    p->map_sz = n;
    for(int i = 0; i < n; ++i) {
        p->map[i] = i;
    }
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine re = std::default_random_engine(seed);
    std::shuffle(p->map, p->map + p->map_sz, re);
}