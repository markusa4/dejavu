//
// Created by markus on 19.09.19.
//

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
