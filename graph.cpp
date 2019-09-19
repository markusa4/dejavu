#include <algorithm>
#include <assert.h>
#include <iostream>
#include "graph.h"


// initialize a coloring of this graph, partitioning degrees of vertices
void graph::initialize_coloring(coloring *c) {
    c->lab.reserve(this->v.size());
    c->ptn.reserve(this->v.size());
    c->vertex_to_col.reserve(this->v.size());
    c->vertex_to_lab.reserve(this->v.size());
    for(int i = 0; i < v.size(); i++) {
        c->vertex_to_col.push_back(-1);
        c->vertex_to_lab.push_back(-1);
        c->lab.push_back(i);
        c->ptn.push_back(1);
    }

    std::sort(c->lab.begin(), c->lab.end(), vertexComparator(*this));

    int cells = 0;
    int last_new_cell   = 0;
    for(int i = 0; i < c->lab.size(); i++) {
        c->vertex_to_col[c->lab[i]] = last_new_cell;
        c->vertex_to_lab[c->lab[i]] = i;
        if(i + 1 == c->lab.size()) {
            cells += 1;
            c->ptn[last_new_cell] = i - last_new_cell;
            c->ptn[i] = 0;
            break;
        }
        assert(this->d[c->lab[i]] <= this->d[c->lab[i + 1]]);
        if(this->d[c->lab[i]] < this->d[c->lab[i + 1]]) {
            c->ptn[i] = 0;
            cells += 1;
            c->ptn[last_new_cell] = i - last_new_cell;
            last_new_cell = i + 1;
            continue;
        }
    }

    std::cout << "Cells: " << cells << std::endl;
}