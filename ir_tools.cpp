//
// Created by markus on 20.09.19.
//

#include <stack>
#include <iostream>
#include <assert.h>
#include "ir_tools.h"
#include "refinement.h"
#include "selector.h"

void ir_tools::label_graph(graph *g, bijection *canon_p) {
    coloring c;
    g->initialize_coloring(&c);

    refinement R;
    selector S;
    trail T;

    std::set<std::pair<int, int>> changes;
    R.refine_coloring(g, &c, &changes);
    T.push_op_r(&changes);
    bool backtrack = false;

    while(T.last_op() != OP_END) {
        int s;
        if(!backtrack) {
            s = S.select_color(g, &c);
            if (s == -1) {
                std::cout << "Discrete coloring found." << std::endl;
                backtrack = true;
            }
        }

        if(T.last_op() == OP_I && !backtrack) {
            changes.clear();
            R.refine_coloring(g, &c, &changes);
            T.push_op_r(&changes);
        } else if(T.last_op() == OP_R  && !backtrack) {
            // collect all elements of color s
            std::stack<int> color_s;
            int i = s;
            while((i == 0) || c.ptn[i - 1] != 0) {
                color_s.push(c.lab[i]);
                i += 1;
            }
            int v = color_s.top();
            color_s.pop();
            assert(color_s.size() > 0);
            // individualize first element
            R.individualize_vertex(g, &c, v);
            T.push_op_i(&color_s, v);
        } else if(T.last_op() == OP_I && backtrack){
            assert(backtrack);
            // backtrack until we can do new operation
            int v = T.top_op_i_v();
            R.undo_individualize_vertex(g, &c, v);
            T.pop_op_i_v();

            // undo individualization
            if(T.top_op_i_class().empty()) {
                T.pop_op_i_class();
                continue;
            } else {
                int v = T.top_op_i_class().top();
                T.top_op_i_class().pop();
                T.push_op_i_v(v);
                R.individualize_vertex(g, &c, v);
                backtrack = false;
            }
        }  else if(T.last_op() == OP_R && backtrack){
            R.undo_refine_color_class(g, &c, &T.top_op_r());
            T.pop_op_r();
        }
    }
}

void trail::push_op_r(std::set<std::pair<int, int>>* color_class_changes) {
    trail_operation.push(OP_R);
    trail_color_class_changes.push(std::set<std::pair<int, int>>());
    trail_color_class_changes.top().swap(*color_class_changes);
}

void trail::push_op_i(std::stack<int>* individualizaiton_todo, int v) {
    trail_operation.push(OP_I);
    trail_op_i_class.push(std::stack<int>());
    trail_op_i_class.top().swap(*individualizaiton_todo);
    trail_op_i_v.push(v);
}

trail::trail() {
    trail_operation.push(OP_END);
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

std::stack<int>& trail::top_op_i_class() {
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
    trail_op_i_v.pop();
}

void trail::push_op_i_v(int v) {
    trail_op_i_v.push(v);
}