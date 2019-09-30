//
// Created by markus on 19.09.19.
//

#include <algorithm>
#include <chrono>
#include <random>
#include "bijection.h"

int bijection::map_vertex(int v) {
    return map[v];
}

void bijection::read_from_coloring(coloring *c) {
    map.clear();
    map.reserve(c->lab.size());
    for(int i = 0; i < c->lab.size(); ++i) {
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

bijection bijection::random_bijection(int n) {
    bijection p;
    p.map = std::vector<int>();
    p.map.reserve(n);
    for(int i = 0; i < n; ++i) {
        p.map.push_back(i);
    }
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine re = std::default_random_engine(seed);
    std::shuffle(p.map.begin(), p.map.end(), re);
    return p;
}
