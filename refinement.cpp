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
#include <cstring>

struct pairhash {
public:
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &x) const
    {
        return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
    }
};

bool refinement::refine_coloring(sgraph *g, coloring *c, change_tracker *changes, invariant* I, int init_color_class, bool track_changes) {
    //std::cout << "Refining..." << std::endl;
    bool comp = true;
    //std::unordered_map<std::pair<int, int>, std::pair<int, int>, pairhash> reduce_class;
    if(!initialized) {
        p = new std::pair<int, int>[g->v_size * 5];
        counting_array.initialize(g->v_size, g->max_degree);
        color_workset.initialize(g->v_size);
        vertex_worklist.initialize(g->v_size);
        old_color_classes.initialize_from_array(p, g->v_size);
        color_worklist_vertex.initialize(g->v_size);
        color_worklist_color.initialize_from_array(p + (g->v_size), g->v_size);
        color_class_splits.initialize(g->v_size);
        queue_pointer.initialize(g->v_size);
        worklist_color_classes.initialize_from_array(p + (g->v_size) * 3, g->v_size * 2);
        initialized = true;
    }

    worklist_color_classes.reset();

    if(init_color_class < 0) {
        // initialize queue with all classes (except for largest one)
        for (int i = 0; i < c->ptn_sz;) {
            queue_pointer.set(i, worklist_color_classes.back_pos);
            worklist_color_classes.push_back(std::pair<int, int>(i, c->ptn[i] + 1));
            assert(worklist_color_classes.arr[queue_pointer.get(i)].first == i);
            i += c->ptn[i] + 1;
        }
    } else {
        queue_pointer.set(init_color_class, worklist_color_classes.back_pos);
        worklist_color_classes.push_back(std::pair<int, int>(init_color_class, c->ptn[init_color_class] + 1));
        assert(worklist_color_classes.arr[queue_pointer.get(init_color_class)].first == init_color_class);
    }
    int its = 0;

    while(!worklist_color_classes.empty()) {
        its += 1;
        //color_class_splits.clear(); // <- 2%
        color_class_splits.reset();
        std::pair<int, int> next_color_class = *worklist_color_classes.front();
        assert(queue_pointer.get(next_color_class.first) == worklist_color_classes.front_pos);
        worklist_color_classes.pop();
        queue_pointer.set(next_color_class.first, -1);

       //std::cout << "Refining color class " << next_color_class.first << ", size: " << next_color_class.second << std::endl;
        // write color class and size to invariant
        //comp = comp && I->write_top_and_compare(-1);
        comp = comp && I->write_top_and_compare(next_color_class.first);
        comp = comp && I->write_top_and_compare(next_color_class.second);
        comp = comp && refine_color_class(g, c, next_color_class.first, next_color_class.second, &color_class_splits, I);

        if(!comp)
            break;

        // add all new classes except for the first, largest one
        int skip = 0;
        //for(auto cc = color_class_splits.begin(); cc != color_class_splits.end(); ++cc) {

        int  latest_old_class = -1;
        bool skipped_largest = false;

        // color class splits are sorted in reverse
        // the old color class will always come last
        while(!color_class_splits.empty()) {
            int  old_class  = color_class_splits.last()->first.first;
            int  new_class  = color_class_splits.last()->first.second;
            bool is_largest = color_class_splits.last()->second;

            if(track_changes)
                changes->track(old_class);

            if(latest_old_class != old_class) {
                latest_old_class = old_class;
                skipped_largest = false;
            }

            color_class_splits.pop_back();
            int new_class_sz = c->ptn[new_class] + 1;


            if(skipped_largest || !is_largest) {
                queue_pointer.set(new_class, worklist_color_classes.back_pos);
                worklist_color_classes.push_back(std::pair<int, int>(new_class, new_class_sz));
            } else {
                skipped_largest = true;
                skip += 1;

                // since old color class will always appear last, the queue pointer of old color class is still valid!
                int i = queue_pointer.get(old_class);
                if(i >= 0) {
                    worklist_color_classes.arr[i].first  = new_class;
                    worklist_color_classes.arr[i].second = new_class_sz;
                    queue_pointer.set(old_class, -1);
                    queue_pointer.set(new_class, i);
                }
            }
        }
        if(!comp) break;
        //assert(color_class_splits.size() > 0? skip > 0: true);
    }
    //assert(assert_is_equitable(g, c));
    //std::cout << its << std::endl;
    // reset queue pointers
    while(!worklist_color_classes.empty()) {
        std::pair<int, int> next_color_class = *worklist_color_classes.front();
        worklist_color_classes.pop();
        queue_pointer.set(next_color_class.first, -1);
    }

    return comp;
}

bool refinement::refine_coloring_first(sgraph *g, coloring *c, int init_color_class) {
    //std::cout << "Refining..." << std::endl;
    bool comp = true;
    //std::unordered_map<std::pair<int, int>, std::pair<int, int>, pairhash> reduce_class;
    if(!initialized) {
        p = new std::pair<int, int>[g->v_size * 5];
        counting_array.initialize(g->v_size, g->max_degree);
        color_workset.initialize(g->v_size);
        vertex_worklist.initialize(g->v_size);
        old_color_classes.initialize_from_array(p, g->v_size);
        color_worklist_vertex.initialize(g->v_size);
        color_worklist_color.initialize_from_array(p + (g->v_size), g->v_size);
        color_class_splits.initialize(g->v_size);
        queue_pointer.initialize(g->v_size);
        worklist_color_classes.initialize_from_array(p + (g->v_size) * 3, g->v_size * 2);
        initialized = true;
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
        // ToDo: rearrange worklist here

        //color_class_splits.clear(); // <- 2%
        color_class_splits.reset();
        std::pair<int, int> next_color_class = *worklist_color_classes.front();
        worklist_color_classes.pop();

        refine_color_class_first(g, c, next_color_class.first, next_color_class.second, &color_class_splits);

        // add all new classes except for the first, largest one
        int skip = 0;
        //for(auto cc = color_class_splits.begin(); cc != color_class_splits.end(); ++cc) {

        int  latest_old_class = -1;
        bool skipped_largest = false;

        while(!color_class_splits.empty()) {
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
    }
    //assert(assert_is_equitable(g, c));
    return comp;
}

//__attribute__((optimize("unroll-loops")))
bool refinement::refine_color_class(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist, invariant* I) {
    // for all vertices of the color class...
    thread_local bool comp, mark_as_largest;
    thread_local int i, cc, vc, pe, end_i, end_cc, v, col, v_old_color, v_class_size, v_new_color, vertex_old_pos, vertex_at_pos, color, v_color, v_degree, fst, snd, largest_color_class_size;
    cc = color_class; // iterate over color class
    comp = true;

    old_color_classes.reset();
    counting_array.reset();

    end_cc = color_class + class_size;
    while (cc < end_cc) { // increment value of neighbours of vc by 1
        vc = c->lab[cc];
        pe = g->v[vc];
        end_i = pe + g->d[vc];
        for (i = pe; i < end_i; i++) {
            // v is a neighbour of vc
            v   = g->e[i];
            col = c->vertex_to_col[v];
            if(counting_array.get_count(v) == 0) {
                vertex_worklist.push_back(v);
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
        v            = vertex_worklist.pop_back();
        v_old_color  = c->vertex_to_col[v];

        v_class_size = c->ptn[v_old_color] + 1;
        v_new_color;
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
        v    = color_worklist_vertex.pop_back();
        color     = color_worklist_color.last()->second;
        v_old_color = color_worklist_color.last()->first;
        color_worklist_color.pop_back();

        vertex_old_pos = c->vertex_to_lab[v];
        vertex_at_pos = c->lab[color + c->ptn[color] + 1];
        c->lab[vertex_old_pos] = vertex_at_pos;
        c->vertex_to_lab[vertex_at_pos] = vertex_old_pos;

        c->lab[color + c->ptn[color] + 1] = v;
        c->vertex_to_col[v] = color;
        c->vertex_to_lab[v] = color + c->ptn[color] + 1;
        c->ptn[color] += 1;

        if (v_old_color != color) {
           // assert(color > old_color);
            c->ptn[v_old_color] -= 1;
        }
    }

    //std::sort(old_color_classes_.begin(), old_color_classes_.end());
    old_color_classes.sort();

    //for(auto it = old_color_classes_.begin(); it != old_color_classes_.end(); ++it) {
    while(!old_color_classes.empty()) {
        fst = old_color_classes.last()->first;
        snd = old_color_classes.last()->second;
        old_color_classes.pop_back();
        largest_color_class_size = -1;

        for(i = fst; i < fst + snd;){
            assert(i >= 0 && i < c->ptn_sz);
            assert(c->ptn[i] + 1 > 0);
            v_color  = i;
            v_degree = counting_array.get_count(c->lab[i]);
            comp = comp && I->write_top_and_compare(-v_color);
            comp = comp && I->write_top_and_compare(v_degree);
            comp = comp && I->write_top_and_compare(c->ptn[i] + 1);

            mark_as_largest = false;
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

        //if(!comp)
        //    break;
    }

    return comp;
}

bool refinement::refine_color_class_first(sgraph *g, coloring *c, int color_class, int class_size, work_list_pair_bool* color_class_split_worklist) {
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
                vertex_worklist.push_back(v);
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

    while(!vertex_worklist.empty()) {
        int v            = vertex_worklist.pop_back();
        int v_old_color  = c->vertex_to_col[v];

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

int refinement::individualize_vertex(coloring *c, int v) {
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
    if(initialized)
        delete[] p;
}

void refinement::undo_changes_with_reference(change_tracker *changes, coloring *changed_c, coloring *old_c) {
    color_workset.reset();
    while(!changes->empty()) {
        int v, i;
        int c      = changes->pop();
        int c_size = old_c->ptn[c] + 1;
        if(c_size > 1 && !color_workset.get(c)) {
            color_workset.set(c);
            int changed_size = changed_c->ptn[c];
            changed_c->ptn[c] = old_c->ptn[c];
            memcpy(changed_c->ptn + c + changed_size, old_c->ptn + c + changed_size, c_size*sizeof(int));
            for(i = changed_size; i < c_size; ++i) {
                   v = changed_c->lab[c + i];
                   changed_c->vertex_to_col[v] = c;
            }
        }
    }
    color_workset.reset();
}

void cumulative_counting::initialize(int size, int maxdegree) {
    m     = maxdegree + 1;
    sz    = size;
    sizes = new int[size * m + 64];
    count = new int[size + 64];
    init = true;

    for(int i = 0; i < size; ++i) {
        count[i]           = 0;
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
    thread_local int in;
    thread_local int _sz;
    assert(count[index] > 0);
    assert(index >= 0 && index < sz);
    //if(count[index] == 0) {
    count[index]++;
    assert(count[index] <= m);

    in  = color * m;
    _sz = sizes[in];
    if(_sz < count[index]) {
        if(_sz == 0)
            reset_queue_sizes.push(in);
        sizes[in + count[index]] = 1;
        sizes[in]++;
    } else {
        sizes[in + count[index]]++;
    }
}

void inline cumulative_counting::increment_r(int index, int color) {
    thread_local int in;
    thread_local int _sz;
    //if(count[index] == 0) {
    assert(index >= 0 && index < sz);
    assert(count[index] == 0);
    count[index]++;
    assert(count[index] <= m);
    reset_queue.push(index);
    in = color * m;
    _sz = sizes[in];
    if(_sz < count[index]) {
        if(_sz == 0) {
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
    assert(index >= 0 && index < sz);
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
    sz = size;
}

void work_set::initialize_from_array(bool* p, int size) {
    s = p;
    //s.reserve(size);
    for(int i = 0; i < size; i++) {
        s[i] = false;
        //s.push_back(false);
    }
    reset_queue.initialize(size);
    //init = true;
}

void work_set::set(int index) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    if(!s[index]) {
        s[index] = true;
        reset_queue.push(index);
    }

}

bool work_set::get(int index) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    return s[index];
}

void work_set::reset() {
    while(!reset_queue.empty()) {
        int index = reset_queue.pop();
        assert(init);
        assert(index >= 0);
        assert(index < sz);
        s[index] = false;
    }
}

work_set::~work_set() {
    if(init)
        delete[] s;
}

void work_set_int::initialize(int size) {
    s = new int[size];
    for(int i = 0; i < size; i++) {
        s[i] = -1;
    }
    init = true;
    sz = size;
}

void work_set_int::set(int index, int value) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    s[index] = value;
}

int work_set_int::get(int index) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    return s[index];
}

void work_set_int::reset() {
    assert(false);
}

work_set_int::~work_set_int() {
    if(init)
        delete[] s;
}

void work_list::initialize(int size) {
    arr = new int[size + 64];
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
    queue = new int[size + 64];
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

void work_queue::reset() {
    pos = 0;
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
    arr = new std::pair<int, int>[size + 64];
    arr_sz = 0;
}

void work_list_pair::initialize_from_array(std::pair<int, int>* p, int size) {
    arr = p;
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
    assert(arr_sz >= 0);
    arr[arr_sz] = value;
    arr_sz += 1;
}

std::pair<std::pair<int, int>, bool>* work_list_pair_bool::last() {
    return &arr[arr_sz - 1];
}

void work_list_pair_bool::initialize(int size) {
    init = true;
    arr = new std::pair<std::pair<int, int>, bool>[size + 64];
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
    arr = new std::pair<int, int>[size + 64];
    arr_sz = size;
    back_pos  = 0;
    front_pos = 0;
    init = true;
}

void ring_pair::initialize_from_array(std::pair<int, int>* p, int size) {
    arr = p;
    arr_sz    = size;
    back_pos  = 0;
    front_pos = 0;
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

void change_tracker::initialize(int _limit, int domain_size) {
    l.initialize(_limit);
    s.initialize(domain_size);
    limit = _limit;
    sz = 0;
}

void change_tracker::reset() {
    l.reset();
    s.reset();
    overflow = false;
    sz       = 0;
}

void change_tracker::track(int oldcolor) {
    if(overflow) return;
    bool check = s.get(oldcolor);
    if(!check && sz < limit - 1) {
        s.set(oldcolor);
        l.push_back(oldcolor);
        sz += 1;
    } else if (!check) {
        overflow = true;
    }
}

int change_tracker::pop() {
    return l.pop_back();
}

bool change_tracker::empty() {
    return l.empty();
}

bool change_tracker::did_overflow() {
    return overflow;
}
