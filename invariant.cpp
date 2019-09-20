#include <assert.h>
#include <iostream>
#include "invariant.h"

std::vector<int> invariant::top_level() {
    assert(vec_invariant.size() > 0);
    return vec_invariant[vec_invariant.size() - 1];
}

int invariant::top_is_geq(std::vector<int> *other) {
    for(int i = 0; i < other->size(); ++i) {
        if(i >= vec_invariant[vec_invariant.size() - 1].size())
            return -1;
        if((*other)[i] < vec_invariant[vec_invariant.size() - 1][i]) {
            return 1;
        }
    }
    if(vec_invariant[vec_invariant.size() - 1].size() > other->size())
        return 1;
    assert(vec_invariant[vec_invariant.size() - 1].size() == other->size());
    return 0;
}

void invariant::pop_level() {
    vec_invariant.pop_back();
}

void invariant::push_level() {
    vec_invariant.push_back(std::vector<int>());
}

void invariant::write_top(int i) {
    vec_invariant[vec_invariant.size() - 1].push_back(i);
}

void invariant::print() {
    for(int i = 0; i < vec_invariant.size(); ++i) {
        for(int j = 0; j < vec_invariant[i].size(); ++j) {
            std::cout << vec_invariant[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int invariant::current_level() {
    return vec_invariant.size() - 1;
}

std::vector<int>* invariant::get_level(int i) {
    return &vec_invariant[i];
}
