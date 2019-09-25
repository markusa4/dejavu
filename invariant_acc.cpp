//
// Created by markus on 24.09.19.
//

#include "invariant_acc.h"
#include <assert.h>
#include <iostream>

std::vector<int> invariant_acc::top_level() {
    assert(vec_invariant.size() > 0);
    return vec_invariant[vec_invariant.size() - 1];
}

int invariant_acc::top_is_geq(std::vector<int> *other) {
    for(int i = 0; i < other->size(); ++i) {
        if(i >= vec_invariant[vec_invariant.size() - 1].size())
            return -1;
        if((*other)[i] < vec_invariant[vec_invariant.size() - 1][i]) {
            return 1;
        }
        if((*other)[i] > vec_invariant[vec_invariant.size() - 1][i]) {
            return -1;
        }
    }
    if(vec_invariant[vec_invariant.size() - 1].size() > other->size())
        return 1;
    assert(vec_invariant[vec_invariant.size() - 1].size() == other->size());
    return 0;
}

bool invariant_acc::top_is_eq(std::vector<int> *other) {
    std::vector<int>* top_vector = &vec_invariant[vec_invariant.size() - 1];
    if(other->size() != top_vector->size())
        return false;
    for(int i = 0; i < other->size(); ++i) {
        if((*other)[i] != (*top_vector)[i])
            return false;
    }
    return true;
}

void invariant_acc::pop_level() {
    vec_invariant.pop_back();
}

void invariant_acc::push_level() {
    vec_invariant.emplace_back(std::vector<int>());
}

void invariant_acc::write_top(int i) {
    vec_invariant[vec_invariant.size() - 1][0] = vec_invariant[vec_invariant.size() - 1][0] * i * 2654435761 % INT32_MAX;
}

void invariant_acc::print() {
    for(int i = 0; i < vec_invariant.size(); ++i) {
        for(int j = 0; j < vec_invariant[i].size(); ++j) {
            std::cout << vec_invariant[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int invariant_acc::current_level() {
    return vec_invariant.size() - 1;
}

std::vector<int>* invariant_acc::get_level(int i) {
    return &vec_invariant[i];
}
