//
// Created by markus on 23.09.19.
//

#include <assert.h>
#include "trail.h"

void trail::push_op_r(std::set<std::pair<int, int>>* color_class_changes) {
    trail_operation.push(OP_R);
    trail_color_class_changes.push(std::set<std::pair<int, int>>());
    trail_color_class_changes.top().swap(*color_class_changes);
}

void trail::push_op_i(std::deque<int>* individualizaiton_todo, int v) {
    trail_operation.push(OP_I);
    trail_op_i_class.push(std::deque<int>());
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

std::set<std::pair<int, int>>& trail::top_op_r() {
    return trail_color_class_changes.top();
}

void trail::pop_op_r() {
    trail_operation.pop();
    trail_color_class_changes.pop();
}

std::deque<int>& trail::top_op_i_class() {
    return trail_op_i_class.top();
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
}

void trail::push_op_i_v(int v) {
    assert(ipos >= 0 && ipos < domain_size);
    ipath[ipos] = v;
    ipos += 1;
    trail_op_i_v.push(v);
}

void trail::free_path() {
    delete[] ipath;
}
