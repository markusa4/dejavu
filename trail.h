//
// Created by markus on 23.09.19.
//

#ifndef DEJAVU_TRAIL_H
#define DEJAVU_TRAIL_H


#include <set>
#include <stack>
#include "bijection.h"
#include "sgraph.h"
#include "ir_tools.h"

enum ir_operation {
    OP_I, OP_R, OP_END
};

class trail {
    std::stack<ir_operation> trail_operation;
    std::stack<std::set<std::pair<int, int>>> trail_color_class_changes;
    std::stack<std::deque<int>> trail_op_i_class;
    std::stack<int> trail_op_i_v;
    int domain_size;
    int* ipath;
    int ipos;
public:
    trail(int domain_size);
    void push_op_r(std::set<std::pair<int, int>> *color_class_changes);
    void push_op_i(std::deque<int> *individualizaiton_todo, int v);
    ir_operation last_op();
    std::set<std::pair<int, int>>& top_op_r();
    std::deque<int>& top_op_i_class();
    int top_op_i_v();
    void pop_op_r();
    void free_path();
    void pop_op_i_class();
    void pop_op_i_v();
    void push_op_i_v(int v);
};

#endif //DEJAVU_TRAIL_H
