//
// Created by markus on 19.09.19.
//

#include <algorithm>
#include <chrono>
#include <random>
#include <iostream>
#include "bijection.h"
#include "coloring_bucket.h"

int bijection::map_vertex(int v) {
    return map[v];
}

void bijection::print() {
    for(int i = 0; i < map.size(); ++i)
        std::cout << map[i] << " ";
    std::cout << std::endl;
}

void bijection::read_from_coloring(coloring *c) {
    map.clear();
    map.reserve(c->lab_sz);
    for(int i = 0; i < c->lab_sz; ++i) {
        map.push_back(c->lab[i]);
    }
}

void bijection::inverse() {
    std::vector<int> old_map = map;
    for(int i = 0; i < map.size(); ++i) {
        map[old_map[i]] = i;
    }
}

void bijection::compose(bijection p) {
    for(int i = 0; i < map.size(); ++i) {
        map[i] = p.map[map[i]];
    }
}

bijection::bijection() {
    map = std::vector<int>();
}

bijection::~bijection() {

}

void bijection::random_bijection(bijection* p, int n) {
    p->map = std::vector<int>();
    p->map.reserve(n);
    for(int i = 0; i < n; ++i) {
        p->map.push_back(i);
    }
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine re = std::default_random_engine(seed);
    std::shuffle(p->map.begin(), p->map.end(), re);
}

void bijection::read_from_coloring_bucket(coloring_bucket *c) {
    map.clear();
    map.reserve(c->lab_sz);
    for(int i = 0; i < c->lab_sz; ++i) {
        map.push_back(c->lab[i]);
    }
}
