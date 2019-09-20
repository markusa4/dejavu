#ifndef BRUTUS_IR_TOOLS_H
#define BRUTUS_IR_TOOLS_H


#include <set>
#include <stack>
#include "bijection.h"
#include "sgraph.h"

class ir_tools {
public:
    void label_graph(sgraph* g, bijection* canon_p);
};

enum ir_operation {
    OP_I, OP_R, OP_END
};

class trail {
    std::stack<ir_operation> trail_operation;
    std::stack<std::set<std::pair<int, int>>> trail_color_class_changes;
    std::stack<std::stack<int>> trail_op_i_class;
    std::stack<int> trail_op_i_v;
    int domain_size;
    int* ipath;
    int ipos;
public:
    trail(int domain_size);
    void push_op_r(std::set<std::pair<int, int>> *color_class_changes);
    void push_op_i(std::stack<int> *individualizaiton_todo, int v);
    ir_operation last_op();
    std::set<std::pair<int, int>>& top_op_r();
    std::stack<int>& top_op_i_class();
    int top_op_i_v();
    void pop_op_r();
    void free_path();
    void pop_op_i_class();
    void pop_op_i_v();
    void push_op_i_v(int v);
};


#endif //BRUTUS_IR_TOOLS_H
