//
// Created by markus on 19.09.19.
//

#include "refinement.h"
#include "invariant_acc.h"
#include <list>
#include <set>
#include <tuple>
#include <iostream>
#include <assert.h>
#include <unordered_map>
#include <algorithm>

struct pairhash {
public:
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &x) const
    {
        return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
    }
};

bool refinement::refine_coloring(sgraph *g, coloring *c, std::list<std::pair<int, int>> *changes, invariant* I, int init_color_class, bool track_changes) {
    //std::cout << "Refining..." << std::endl;
    bool comp = true;
    //std::unordered_map<std::pair<int, int>, std::pair<int, int>, pairhash> reduce_class;
    if(!initialized) {
        counting_array.initialize(g->v_size, g->max_degree);
        color_workset.initialize(g->v_size);
        vertex_worklist.initialize(g->v_size);
        old_color_classes.initialize(g->v_size);
        color_worklist_vertex.initialize(g->v_size);
        color_worklist_color.initialize(g->v_size);
        color_class_splits.initialize(g->v_size);
        initialized = true;
        worklist_color_classes.initialize(g->v_size * 2);
    }

    worklist_color_classes.reset();

    if(init_color_class < 0) {
        // initialize queue with all classes (except for largest one)
        for (int i = 0; i < c->ptn_sz;) {
                worklist_color_classes.push_back(std::pair<int, int>(i, c->ptn[i] + 1));
            i += c->ptn[i] + 1;
        }
    } else {
        worklist_color_classes.push_back(std::pair<int, int>(init_color_class, c->ptn[init_color_class] + 1));
    }

    while(!worklist_color_classes.empty()) {
        //color_class_splits.clear(); // <- 2%
        color_class_splits.reset();
        std::pair<int, int> next_color_class = *worklist_color_classes.front();
        worklist_color_classes.pop();

       //std::cout << "Refining color class " << next_color_class.first << ", size: " << next_color_class.second << std::endl;
        // write color class and size to invariant
        comp = comp && I->write_top_and_compare(-1);
        comp = comp && I->write_top_and_compare(next_color_class.first);
        comp = comp && I->write_top_and_compare(next_color_class.second);
        comp = comp && refine_color_class(g, c, next_color_class.first, next_color_class.second, &color_class_splits, I);

        // add all new classes except for the first, largest one
        int skip = 0;
        //for(auto cc = color_class_splits.begin(); cc != color_class_splits.end(); ++cc) {

        int  latest_old_class = -1;
        bool skipped_largest = false;

        while(!color_class_splits.empty()) {
            if(track_changes)
                changes->push_back(color_class_splits.last()->first);
            int old_class    = color_class_splits.last()->first.first;
            int new_class    = color_class_splits.last()->first.second;
            bool is_largest  = color_class_splits.last()->second;

            if(latest_old_class != old_class) {
                latest_old_class = old_class;
                skipped_largest = false;
            }

            color_class_splits.pop_back();
            int new_class_sz = c->ptn[new_class] + 1;
            if(skipped_largest || !is_largest) {
                worklist_color_classes.push_back(std::pair<int, int>(new_class, new_class_sz));
            } else {
                skipped_largest = true;
                skip += 1;
            }
        }
        if(!comp) break;
        //assert(color_class_splits.size() > 0? skip > 0: true);
    }
    //assert(assert_is_equitable(g, c));
    return comp;
}

//__attribute__((optimize("unroll-loops")))
bool refinement::refine_color_class(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist, invariant* I) {
    // for all vertices of the color class...
    bool comp = true;
    int cc = color_class; // iterate over color class

    while (cc < color_class + class_size) { // increment value of neighbours of vc by 1
        int vc = c->lab[cc];
        int pe = g->v[vc];
        int end_i = pe + g->d[vc];
        for (int i = pe; i < end_i; i++) {
            // v is a neighbour of vc
            int v   = g->e[i];
            int col = c->vertex_to_col[v];
            if(counting_array.get_count(v) == 0) {
                vertex_worklist.push_back(std::pair<int, int>(v, col));
                counting_array.increment_r(v, col);
                if (!color_workset.get(col)) {
                    old_color_classes.push_back(std::pair<int, int>(col, c->ptn[col] + 1));
                    color_workset.set(col);
                }
            } else {
                counting_array.increment(v, col);
            }
        }
        cc += 1;
    }

    // split color classes according to count in counting array
    // care: does not respect lab order
    // ToDo: consider catching the case where vertices old_color does not appear (color degree = 0)
    // ToDo: then we can immediately assign instead of swapping
    // Problem: how do we get all "other vertices" of the color, though?

    while(!vertex_worklist.empty()) {
        auto w = vertex_worklist.last();
        int v            = w->first;
        int v_old_color  = w->second;
        vertex_worklist.pop_back();

        int v_class_size = c->ptn[v_old_color] + 1;
        int v_new_color;
        if (counting_array.get_count(v) == 0) {
            v_new_color = v_old_color;
        } else {
            v_new_color = v_old_color + v_class_size - counting_array.get_size(v, v_old_color);
        }

        if (v_new_color != v_old_color) {
            color_worklist_vertex.push_back(v);
            color_worklist_color.push_back(std::pair<int, int>(v_old_color, v_new_color));
            c->ptn[v_new_color] = -1;
        }
    }

    color_workset.reset();

    while(!color_worklist_vertex.empty()) {
        assert(!color_worklist_color.empty());
        int vertex    = color_worklist_vertex.pop_back();
        int color     = color_worklist_color.last()->second;
        int old_color = color_worklist_color.last()->first;
        color_worklist_color.pop_back();

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

    //std::sort(old_color_classes_.begin(), old_color_classes_.end());
    old_color_classes.sort();

    //for(auto it = old_color_classes_.begin(); it != old_color_classes_.end(); ++it) {
    while(!old_color_classes.empty()) {
        int fst = old_color_classes.last()->first;
        int snd = old_color_classes.last()->second;
        old_color_classes.pop_back();
        int largest_color_class_size = -1;

        for(int i = fst; i < fst + snd;){
            assert(i >= 0 && i < c->ptn_sz);
            assert(c->ptn[i] + 1 > 0);
            int v_color  = i;
            int v_degree = counting_array.get_count(c->lab[i]);
            comp = comp && I->write_top_and_compare(-v_color);
            comp = comp && I->write_top_and_compare(v_degree);
            comp = comp && I->write_top_and_compare(c->ptn[i] + 1);
            bool mark_as_largest = false;
            if(!(i == fst && c->ptn[i] + 1== snd)) {
                if(largest_color_class_size < c->ptn[i] + 1) {
                    mark_as_largest = true;
                    largest_color_class_size = c->ptn[i] + 1;
                }
                color_class_split_worklist->push_back(std::pair<std::pair<int, int>, bool>(std::pair<int, int>(fst, i), mark_as_largest));
            } else {
                break;
            }

            /* */
            if(v_color != 0) {
                c->ptn[v_color - 1] = 0;
            }
            /* */

            i += c->ptn[i] + 1;
        }
    }

    counting_array.reset();
    old_color_classes.reset();

    return comp;
}

int refinement::individualize_vertex(sgraph *g, coloring *c, int v) {
    //assert(initialized);

    int color = c->vertex_to_col[v];
    int pos   = c->vertex_to_lab[v];
    //int pos   = color;
    //while(c->lab[pos] != v) pos += 1;

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
    return color + color_class_size;

}

void refinement::undo_individualize_vertex(sgraph *g, coloring *c, int v) {
    assert(false);
    /*int color = c->vertex_to_col[v];
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
    }*/
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
    for(int i = 0; i < c->lab_sz; ++i) {
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
    for(int i = 0; i < g->v_size; ++i) {
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

    for(int i = 0; i < g->v_size; ++i) {
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
    for(int i = 0; i < c->lab_sz; ++i) {
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
}

void cumulative_counting::initialize(int size, int maxdegree) {
    m     = maxdegree + 1;
    sizes = new int[size * m];
    count = new int[size];
    init = true;

    for(int i = 0; i < size; ++i) {
        count[i]       = 0;
        this->sizes[i * m] = 0;
    }
    reset_queue.initialize(size);
    reset_queue_sizes.initialize(size);

    reset();
}

void cumulative_counting::reset() {
    while(!reset_queue.empty()) {
        int index = reset_queue.pop();
        count[index] = 0;
    }
    while(!reset_queue_sizes.empty()) {
        int index = reset_queue_sizes.pop();
        sizes[index] = 0;
        //sizes[index].push_back(-1);
    }
}

void inline cumulative_counting::increment(int index, int color) {
    //if(count[index] == 0) {
    count[index] += 1;
    int in = color * m;
    int sz = sizes[in];
    if(sz < count[index]) {
        if(sz <= 1) {
            reset_queue_sizes.push(in);
        }
        sizes[in + count[index]] = 1;
        sizes[in] += 1;
    } else {
        sizes[in + count[index]] += 1;
    }
}

void inline cumulative_counting::increment_r(int index, int color) {
    //if(count[index] == 0) {
    count[index] += 1;
    reset_queue.push(index);
    int in = color * m;
    int sz = sizes[in];
    if(sz < count[index]) {
        if(sz <= 1) {
            reset_queue_sizes.push(in);
        }
        sizes[in + count[index]] = 1;
        sizes[in] += 1;
    } else {
        sizes[in + count[index]] += 1;
    }
}


int cumulative_counting::get_size(int index, int color) {
    assert(count[index] > 0);
    //assert(sizes[c->vertex_to_col[index]].size() > 1);
    return sizes[color * m + count[index]];
}

int cumulative_counting::get_count(int index) {
    //assert(index < count);
    return count[index];
}

cumulative_counting::~cumulative_counting() {
    if(init) {
        delete[] count;
        delete[] sizes;
    }
}

void work_set::initialize(int size) {
    s = new bool[size];
    //s.reserve(size);
    for(int i = 0; i < size; i++) {
        s[i] = false;
        //s.push_back(false);
    }
    reset_queue.initialize(size);
    init = true;
}

void work_set::set(int index) {
    if(!s[index]) {
        s[index] = true;
        reset_queue.push(index);
    }

}

bool work_set::get(int index) {
    return s[index];
}

void work_set::reset() {
    while(!reset_queue.empty()) {
        int index = reset_queue.pop();
        s[index] = false;
    }
}

work_set::~work_set() {
    if(init)
        delete[] s;
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

void work_queue::initialize(int size) {
    assert(!init);
    sz = size;
    pos = 0;
    queue = new int[size];
    init = true;
}

void work_queue::push(int val) {
    assert(init);
    assert(pos != sz);
    queue[pos] = val;
    pos++;
}

int work_queue::pop() {
    assert(init);
    assert(pos > 0);
    pos--;
    return queue[pos];
}

bool work_queue::empty() {
    assert(init);
    return (pos == 0);
}

work_queue::~work_queue() {
    if(init) {
        delete[] queue;
    }
}

void work_list_pair::push_back(std::pair<int, int> value) {
    arr[arr_sz] = value;
    arr_sz += 1;
}

std::pair<int, int> *work_list_pair::last() {
    return &arr[arr_sz - 1];
}

void work_list_pair::initialize(int size) {
    init = true;
    arr = new std::pair<int, int>[size];
    arr_sz = 0;
}

work_list_pair::~work_list_pair() {
    if(init)
        delete[] arr;
}

void work_list_pair::pop_back() {
    arr_sz -= 1;
}

void work_list_pair::sort() {
    if(arr_sz > 1) {
        std::sort(arr, arr + arr_sz);
    }
}

bool work_list_pair::empty() {
    return arr_sz == 0;
}

void work_list_pair::reset() {
    arr_sz = 0;
}

void work_list_pair_bool::push_back(std::pair<std::pair<int, int>, bool> value) {
    arr[arr_sz] = value;
    arr_sz += 1;
}

std::pair<std::pair<int, int>, bool>* work_list_pair_bool::last() {
    return &arr[arr_sz - 1];
}

void work_list_pair_bool::initialize(int size) {
    init = true;
    arr = new std::pair<std::pair<int, int>, bool>[size];
    arr_sz = 0;
}

work_list_pair_bool::~work_list_pair_bool() {
    if(init)
        delete[] arr;
}

void work_list_pair_bool::pop_back() {
    arr_sz -= 1;
}

void work_list_pair_bool::sort() {
    if(arr_sz > 1) {
        std::sort(arr, arr + arr_sz);
    }
}

bool work_list_pair_bool::empty() {
    return arr_sz == 0;
}

void work_list_pair_bool::reset() {
    arr_sz = 0;
}

void ring_pair::initialize(int size) {
    arr = new std::pair<int, int>[size];
    arr_sz = size;
    back_pos  = 0;
    front_pos = 0;
    init = true;
}

void ring_pair::push_back(std::pair<int, int> value) {
    arr[back_pos] = value;
    back_pos = (back_pos + 1) % arr_sz;
}

std::pair<int, int> *ring_pair::front() {
    return &arr[front_pos];
}

void ring_pair::pop() {
    front_pos = (front_pos + 1) % arr_sz;
}

bool ring_pair::empty() {
    return (front_pos == back_pos);
}

ring_pair::~ring_pair() {
    if(init)
        delete[] arr;
}

void ring_pair::reset() {
    front_pos = back_pos;
}
