//
// Created by markus on 19.09.19.
//

#include "refinement.h"
#include "invariant_acc.h"
#include <list>
#include <set>
#include <iostream>
#include <assert.h>
#include <unordered_map>
#include <algorithm>

void refinement::refine_coloring(sgraph *g, coloring *c, std::list<std::pair<int, int>> *changes, invariant* I, std::list<int>* init_color_class, bool track_changes) {
    //std::cout << "Refining..." << std::endl;
    if(!initialized) {
        counting_array.initialize(g->v.size(), c);
        vertex_workset.initialize(g->v.size());
        color_worklset.initialize(g->v.size());
        vertex_worklist.initialize(g->v.size());
        color_worklist_vertex.initialize(g->v.size());
        color_worklist_color.initialize(g->v.size());
        initialized = true;
        largest_color_class_index = new int[c->lab.size()];
    }
    counting_array.set_coloring(c);
    std::queue<std::pair<int, int>> worklist_color_classes;

    if(init_color_class->empty()) {
        // initialize queue with all classes
        for (int i = 0; i < c->ptn.size();) {
            worklist_color_classes.push(std::pair<int, int>(i, c->ptn[i] + 1));
            i += c->ptn[i] + 1;
        }
    } else {
        for(auto it = init_color_class->begin(); it != init_color_class->end(); ++it) {
            worklist_color_classes.push(std::pair<int, int>(*it, c->ptn[*it] + 1));
        }
    }

    while(!worklist_color_classes.empty()) {
        std::list<std::pair<int, int>> color_class_splits;
        std::pair<int, int> next_color_class = worklist_color_classes.front();
        worklist_color_classes.pop();

       //std::cout << "Refining color class " << next_color_class.first << ", size: " << next_color_class.second << std::endl;
        // write color class and size to invariant
        I->write_top(-1);
        I->write_top(next_color_class.first);
        I->write_top(next_color_class.second);

        refine_color_class(g, c, next_color_class.first, next_color_class.second, &color_class_splits, I, largest_color_class_index);

        // add all new classes except for the first, largest one
        int skip = 0;
        for(auto cc = color_class_splits.begin(); cc != color_class_splits.end(); ++cc) {
            if(track_changes)
                changes->push_back(*cc);
            int old_class    = cc->first;
            int new_class    = cc->second;
            int new_class_sz = c->ptn[new_class] + 1;
            assert(largest_color_class_index[old_class] != -1);
            if(largest_color_class_index[old_class] != new_class ) {
                worklist_color_classes.push(std::pair<int, int>(new_class, new_class_sz));
            } else {
                skip += 1;
            }
        }
        assert(color_class_splits.size() > 0? skip > 0: true);
    }
    //assert(assert_is_equitable(g, c));
}

void refinement::refine_color_class(sgraph *g, coloring *c, int color_class, int class_size, std::list<std::pair<int, int>> *color_class_split_worklist, invariant* I, int* largest_color_class_index) {
    // for all vertices of the color class...
    // ToDo: can replace worklists with fixed size arrays
    //std::list<std::pair<int, int>> color_set_worklist;
    std::set<int> new_colors;
    std::list<std::pair<int, int>> old_color_classes;

    int cc = color_class; // iterate over color class
    while (cc < color_class + class_size) { // increment value of neighbours of vc by 1
        int vc = c->lab[cc];
        int pe = g->v[vc];
        for (int i = pe; i < pe + g->d[vc]; i++) {
            // v is a neighbour of vc
            int v = g->e[i];
            if (!vertex_workset.get(v)) { // <- ToDo: think about this: && c->ptn[c->vertex_to_col[v]] > 0
                vertex_workset.set(v);
                vertex_worklist.push_back(v);
            }
            counting_array.increment(v);
        }
        cc += 1;
    }

    // split color classes according to count in counting array
    // care: does not respect lab order
    while(!vertex_worklist.empty()) {
        int v = vertex_worklist.pop_back();
        int v_class_size = c->ptn[c->vertex_to_col[v]] + 1;
        int v_old_color  = c->vertex_to_col[v];
        int v_new_color;
        if (counting_array.get_count(v) == 0) {
            v_new_color = v_old_color;
        } else {
            v_new_color = v_old_color + v_class_size - counting_array.get_size(v);
        }
        if (!color_worklset.get(v_old_color)) {
            assert(c->ptn[v_old_color] >= 0);
            old_color_classes.emplace_back(std::pair<int, int>(v_old_color, v_class_size));
            color_worklset.set(v_old_color);
        }
        if (v_new_color != v_old_color) {
            color_worklist_vertex.push_back(v);
            color_worklist_color.push_back(v_new_color);
            //color_set_worklist.emplace_back(std::pair<int, int>(v, v_new_color));
            new_colors.insert(v_new_color);
            c->ptn[v_new_color] = -1;
        }
    }

    while(!color_worklist_vertex.empty()) {
        assert(!color_worklist_color.empty());
        int vertex = color_worklist_vertex.pop_back();
        int color  = color_worklist_color.pop_back();
        int old_color = c->vertex_to_col[vertex];

        int vertex_old_pos = c->vertex_to_lab[vertex];
        int vertex_at_pos = c->lab[color + c->ptn[color] + 1];
        c->lab[vertex_old_pos] = vertex_at_pos;
        c->vertex_to_lab[vertex_at_pos] = vertex_old_pos;

        c->lab[color + c->ptn[color] + 1] = vertex;
        c->vertex_to_col[vertex] = color;
        c->vertex_to_lab[vertex] = color + c->ptn[color] + 1;
        c->ptn[color] += 1;

        if (old_color != color) {
            assert(color > old_color);
            c->ptn[old_color] -= 1;
        }
    }

    for(auto new_color = new_colors.begin(); new_color != new_colors.end(); new_color++) {
        if(*new_color != 0) {
            c->ptn[*new_color - 1] = 0;
        }
    }

    old_color_classes.sort();

    for(auto it = old_color_classes.begin(); it != old_color_classes.end(); ++it) {
        int largest_color_class      = -1;
        int largest_color_class_size = -1;

        for(int i = it->first; i < it->first + it->second;){
            assert(i >= 0 && i < c->ptn.size());
            assert(c->ptn[i] + 1 > 0);
            int v_color  = i;
            int v_degree = counting_array.get_count(c->lab[i]);
            I->write_top(-40);
            I->write_top(v_color);
            I->write_top(v_degree);
            I->write_top(c->ptn[i] + 1);
            if(!(i == it->first && c->ptn[i] + 1== it->second)) {
                if(largest_color_class_size < c->ptn[i] + 1) {
                    largest_color_class_size = c->ptn[i] + 1;
                    largest_color_class = i;
                }
                color_class_split_worklist->push_front(std::pair<int, int>(it->first, i));
            } else {
                break;
            }
            i += c->ptn[i] + 1;
        }

        largest_color_class_index[it->first] = largest_color_class;
    }
    vertex_worklist.reset();
    vertex_workset.reset();
    color_worklset.reset();
    color_worklist_color.reset();
    color_worklist_vertex.reset();
    counting_array.reset();
}

void refinement::individualize_vertex(sgraph *g, coloring *c, int v) {
    int color = c->vertex_to_col[v];
    int pos   = c->vertex_to_lab[v];
    int color_class_size = c->ptn[color];
    assert(color_class_size > 0);
    int new_color = color + color_class_size;

    int vertex_at_pos = c->lab[color + color_class_size];
    c->lab[pos] = vertex_at_pos;
    c->vertex_to_lab[vertex_at_pos] = pos;

    c->lab[color + color_class_size] = v;
    c->vertex_to_lab[v] = color + color_class_size;
    c->vertex_to_col[v] = color + color_class_size;

    c->ptn[color] -= 1;
    c->ptn[color + color_class_size - 1] = 0;
}

void refinement::undo_individualize_vertex(sgraph *g, coloring *c, int v) {
    int color = c->vertex_to_col[v];
    int pos   = c->vertex_to_lab[v];
    assert(color == pos);
    assert(pos > 0);
    assert(c->ptn[pos] == 0);
    int new_color = c->vertex_to_col[c->lab[pos - 1]];

    c->ptn[new_color] += 1;
    c->vertex_to_col[v] = new_color;
    assert(c->ptn[new_color] > 0);
    if(c->ptn[pos - 1] == 0) {
        c->ptn[pos - 1] = 1;
    }
}

void refinement::undo_refine_color_class(sgraph *g, coloring *c, std::list<std::pair<int, int>> *changes) {
    std::queue<int> reset_ends;
    for(auto p = changes->begin(); p != changes->end(); ++p) {
        int old_color = c->vertex_to_col[c->lab[p->first]];
        int new_color = p->second;
        //std::cout << new_color << " back to " << old_color << std::endl;
        if(old_color < new_color && c->vertex_to_col[c->lab[new_color]] == new_color) {
            reset_ends.push(new_color - 1);
            for (int j = new_color; j <= new_color + c->ptn[new_color]; ++j) {
                c->vertex_to_col[c->lab[j]] = old_color;
                c->ptn[old_color] += 1;
            }
        }
    }

    while(!reset_ends.empty()) {
        int i = reset_ends.front();
        reset_ends.pop();
        if(i >= 0 && c->ptn[i] == 0) {
            c->ptn[i] = 1;
        }
    }

    // sanity check
    int expect0 = 0;
    bool expectsize = true;
    for(int i = 0; i < c->lab.size(); ++i) {
        if(expectsize) {
            expectsize = false;
            expect0 = c->ptn[i];
        }
        if(expect0 == 0) {
            assert(c->ptn[i] == 0);
            expectsize = true;
        } else {
            assert(c->ptn[i] > 0);
        }
        expect0 -= 1;
    }
}

void refinement::complete_colorclass_invariant(sgraph *g, coloring *c, invariant_acc *I) {
    for(int i = 0; i < g->v.size(); ++i) {
        std::list<int> neighbour_col;
        int v = c->lab[i];
        I->write_top(-5);
        I->write_top(c->vertex_to_col[v]);
        for(int j = g->v[v]; j < g->v[v] + g->d[v]; ++j) {
            neighbour_col.push_back(c->vertex_to_col[g->e[j]]);
        }
        neighbour_col.sort();
        for(auto n = neighbour_col.begin(); n != neighbour_col.end(); ++n) {
            I->write_top(*n);
        }
    }
}

bool refinement::assert_is_equitable(sgraph *g, coloring *c) {
    std::vector<int> neighbour_col_canon;
    std::vector<int> neighbour_col;
    bool first_of_color = true;
    bool test = true;

    for(int i = 0; i < g->v.size(); ++i) {
        neighbour_col.clear();
        int v = c->lab[i];
        for(int j = g->v[v]; j < g->v[v] + g->d[v]; ++j) {
            neighbour_col.push_back(c->vertex_to_col[g->e[j]]);
        }
        std::sort(neighbour_col.begin(), neighbour_col.end());
        if(first_of_color) {
            first_of_color = false;
            neighbour_col_canon.swap(neighbour_col);
        } else {
            for (auto j = 0; j < neighbour_col.size(); ++j) {
                test = test && (neighbour_col[j] == neighbour_col_canon[j]);
            }
        }

        if(c->ptn[i] == 0) {
            first_of_color = true;
        }
    }
    int expect0 = 0;
    bool expectsize = true;
    //std::cout << "part";
    for(int i = 0; i < c->lab.size(); ++i) {
        if(expectsize) {
            expectsize = false;
            expect0 = c->ptn[i];
        }
        if(expect0 == 0) {
            //std::cout << i << " ";
            test = test && (c->ptn[i] == 0);
            expectsize = true;
        } else {
            test = test && (c->ptn[i] > 0);
        }
        expect0 -= 1;
    }
    return test;
}

refinement::~refinement() {
    if(initialized) {
        delete[] largest_color_class_index;
    }

}

void cumulative_counting::initialize(int size, coloring *c) {
    this->c = c;
    for(int i = 0; i < size; ++i) {
        this->count.push_back(0);
        this->sizes.emplace_back(std::vector<int>());
        this->sizes[this->sizes.size() - 1].push_back(-1);
        //this->sizes[this->sizes.size() - 1].reserve(16);
    }
    reset();
}

void cumulative_counting::reset() {
    while(!reset_queue.empty()) {
        int index = reset_queue.front();
        reset_queue.pop();
        count[index] = 0;
    }
    while(!reset_queue_sizes.empty()) {
        int index = reset_queue_sizes.front();
        reset_queue_sizes.pop();
        sizes[index].clear();
        sizes[index].push_back(-1);
        //sizes[index].reserve(16);
    }
}

void cumulative_counting::increment(int index) {
    if(count[index] == 0) {
        reset_queue.push(index);
    }
    count[index] += 1;
    assert(c->vertex_to_col[index] == c->vertex_to_col[c->lab[c->vertex_to_col[index]]]);
    int b_col_index = c->vertex_to_col[index];
    if(sizes[b_col_index].size() <= 1) {
        assert(sizes[b_col_index].size() == 1);
        reset_queue_sizes.push(b_col_index);
    }
    if(sizes[b_col_index].size() <= count[index]) {
        sizes[b_col_index].push_back(1);
    } else {
        sizes[b_col_index][count[index]] += 1;
    }
}

int cumulative_counting::get_size(int index) {
    assert(count[index] > 0);
    assert(sizes[c->vertex_to_col[index]].size() > 1);
    return sizes[c->vertex_to_col[index]][count[index]];
}

int cumulative_counting::get_count(int index) {
    assert(index < count.size());
    return count[index];
}

void cumulative_counting::set_coloring(coloring *c) {
    this->c = c;
}

void work_set::initialize(int size) {
    for(int i = 0; i < size; i++) {
        s.push_back(false);
    }
}

void work_set::set(int index) {
    if(s[index] == false) {
        s[index] = true;
        reset_queue.push(index);
    }

}

bool work_set::get(int index) {
    return s[index];
}

void work_set::reset() {
    while(!reset_queue.empty()) {
        int index = reset_queue.front();
        reset_queue.pop();
        s[index] = false;
    }
}

void work_list::initialize(int size) {
    arr = new int[size];
    arr_sz = size;
    cur_pos = 0;
}

void work_list::push_back(int value) {
    assert(cur_pos >= 0 && cur_pos < arr_sz);
    arr[cur_pos] = value;
    cur_pos += 1;
}

int work_list::pop_back() {
    cur_pos -= 1;
    assert(cur_pos >= 0 && cur_pos < arr_sz);
    return arr[cur_pos];
}

void work_list::reset() {
    cur_pos = 0;
}

bool work_list::empty() {
    return cur_pos == 0;
}

work_list::~work_list() {
    if(arr_sz >= 0) {
        delete[] arr;
    }
}
