//
// Created by markus on 23.09.19.
//

#include <assert.h>
#include <iostream>
#include <algorithm>
#include <random>
#include "trail.h"

void trail::push_op_r(std::list<std::pair<int, int>>* color_class_changes) {
    trail_operation.push(OP_R);
    trail_color_class_changes.push(std::list<std::pair<int, int>>());
    trail_color_class_changes.top().swap(*color_class_changes);
}

void trail::push_op_i(std::list<int>* individualizaiton_todo, int v) {
    trail_operation.push(OP_I);
    trail_op_i_class.push(std::list<int>());
    trail_op_i_class.top().swap(*individualizaiton_todo);
    push_op_i_v(v);
}

trail::trail(int domain_size) {
    trail_operation.push(OP_END);
    this->domain_size = domain_size;
    ipath = new int[domain_size];
    ipos = 0;
}

ir_operation trail::last_op() {
    return trail_operation.top();
}

std::list<std::pair<int, int>>& trail::top_op_r() {
    return trail_color_class_changes.top();
}

void trail::pop_op_r() {
    trail_operation.pop();
    trail_color_class_changes.pop();
}

std::list<int>& trail::top_op_i_class() {
    return trail_op_i_class.top();
}
void trail::shuffle_top_i_class(std::default_random_engine* re) {
    //std::shuffle(trail_op_i_class.top().begin(), trail_op_i_class.top().end(), *re);
}

void trail::pop_op_i_class() {
    trail_operation.pop();
    trail_op_i_class.pop();
}

int trail::top_op_i_v() {
    return trail_op_i_v.top();
}

void trail::pop_op_i_v() {
    ipos -= 1;
    trail_op_i_v.pop();
    path.pop_back();
}

void trail::push_op_i_v(int v) {
    assert(ipos >= 0 && ipos < domain_size);
    ipath[ipos] = v;
    ipos += 1;
    trail_op_i_v.push(v);
    path.push_back(v);
}

void trail::print_path() {
    for(int i = 0; i < ipos; ++i) {
        std::cout << ipath[i] << " ";
    }
    std::cout << std::endl;
}

void trail::free_path() {
    delete[] ipath;
}

int trail::size() {
    return trail_operation.size();
}