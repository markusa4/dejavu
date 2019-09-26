//
// Created by markus on 24.09.19.
//

#ifndef DEJAVU_INVARIANT_ACC_H
#define DEJAVU_INVARIANT_ACC_H


#include <vector>
#include "invariant.h"

class invariant_acc {
    std::vector<std::vector<int>> vec_invariant;
public:
    void push_level();
    void pop_level();
    int current_level();
    std::vector<int>* get_level(int i);
    std::vector<int> top_level();
    int top_is_geq(std::vector<int> *other);
    void write_top(int i);
    void print();

    bool top_is_eq(std::vector<int> *other);

    bool level_is_eq(invariant_acc *other, int level);
};


#endif //DEJAVU_INVARIANT_ACC_H
