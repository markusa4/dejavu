#include <assert.h>
#include <iostream>
#include "invariant.h"

std::vector<int> invariant::top_level() {
    assert(vec_invariant.size() > 0);
    return vec_invariant[vec_invariant.size() - 1];
}

void invariant::top_is_geq(std::vector<int> *other) {

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
