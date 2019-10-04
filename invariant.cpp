#include <assert.h>
#include <iostream>
#include "invariant.h"

std::vector<int> invariant::top_level() {
    assert(vec_invariant.size() > 0);
    return vec_invariant[vec_invariant.size() - 1];
}

bool invariant::level_is_eq(invariant* other, int level) {
    std::vector<int>* my_vec = &vec_invariant[level];
    std::vector<int>* other_vec = other->get_level(level);

    if(other_vec->size() != my_vec->size())
        return false;
    for(int i = 0; i < other_vec->size(); ++i) {
        if((*other_vec)[i] != (*my_vec)[i])
            return false;
    }
    return true;
}

int invariant::top_is_geq(std::vector<int> *other) {
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

bool invariant::top_is_eq(std::vector<int> *other) {
    std::vector<int>* top_vector = &vec_invariant[vec_invariant.size() - 1];
    if(other->size() != top_vector->size())
        return false;
    for(int i = 0; i < other->size(); ++i) {
        if((*other)[i] != (*top_vector)[i])
            return false;
    }
    return true;
}

void invariant::pop_level() {
    vec_invariant.pop_back();
}

void invariant::push_level() {
    vec_invariant.emplace_back(std::vector<int>());
}

void invariant::write_top(int i) {
    vec_invariant[vec_invariant.size() - 1].push_back(i);
}

bool invariant::write_top_and_compare(int i) {
    int pos1 = vec_invariant.size() - 1;
    vec_invariant[pos1].push_back(i);
    int pos2 = vec_invariant[pos1].size() - 1;
    if(has_compare) {
        return vec_invariant[pos1][pos2] == (*compareI->get_level(pos1))[pos2];
    } else {
        return true;
    }
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

void invariant::set_compare_invariant(invariant* I) {
    has_compare = true;
    compareI = I;
}