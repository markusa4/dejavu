#include "refinement.h"
#include "utility.h"
#include <list>
#include <assert.h>
#include <algorithm>
#include <cstring>

void work_set::initialize(int size) {
    s = new bool[size];
    //s.reserve(size);

    //for(int i = 0; i < size; ++i)
    //    s.push_back(false);
    memset(s, false, size * sizeof(bool));

    reset_queue.initialize(size);
    init = true;
    sz = size;
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

void work_set::unset(int index) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    s[index] = false;
}

void work_set::set_nr(int index) {
    s[index] = true;
}

bool work_set::get(int index) {
    assert(init);
    assert(index >= 0);
    assert(index < sz);
    return s[index];
}

void work_set::reset() {
    while(!reset_queue.empty()) {
        const int index = reset_queue.pop();
        assert(init);
        assert(index >= 0);
        assert(index < sz);
        s[index] = false;
    }
}

void work_set::reset_hard() {
    reset_queue.pos = 0;
    memset(s, false, sz*sizeof(bool));
}

work_set::~work_set() {
    if(init)
        delete[] s;
}

void work_set::reset_soft() {
    reset_queue.reset();
}

void work_queue::initialize(int size) {
    assert(!init);
    sz = size;
    pos = 0;
    queue = new int[size];
    init = true;
}

void work_queue::push(int val) {
    //assert(init);
    assert(pos != sz);
    queue[pos] = val;
    pos++;
}

int work_queue::pop() {
    //assert(init);
    assert(pos > 0);
    pos--;
    return queue[pos];
}

bool work_queue::empty() {
    //assert(init);
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

void work_queue::initialize_from_array(int* arr, int size) {
    assert(!init);
    sz = size;
    pos = 0;
    queue = arr;
    init = false;
}

void ring_pair::initialize(int size) {
    arr = new std::pair<int, int>[size];
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

void cell_worklist::initialize(int domain_size) {
    arr = new int[domain_size];
    arr_sz = domain_size;
    init = true;
    cur_pos = 0;
}

int cell_worklist::add_cell(work_set_int *queue_pointer, int col) {
    assert(init);
    assert(cur_pos >= 0 && cur_pos < arr_sz - 1);
    queue_pointer->set(col, cur_pos);
    arr[cur_pos] = col;
    cur_pos++;
    return 0;
}

int cell_worklist::next_cell(work_set_int *queue_pointer, coloring* c) {
    // look at first 12 positions and scramble if possible
    int sm_j = cur_pos - 1;
    for(int j = cur_pos - 1; j >= 0 && ((cur_pos - j) <= 12); --j) {
        if(c->ptn[arr[j]] < c->ptn[arr[sm_j]]) {
            sm_j = j;
            if (c->ptn[arr[sm_j]] == 0)
                break;
        }
    }

    // swap sm_j and j
    const int sm_col = arr[sm_j];
    arr[sm_j] = arr[cur_pos - 1];
    queue_pointer->set(arr[sm_j], sm_j);

    cur_pos--;
    queue_pointer->set(sm_col, -1);
    return sm_col;
}

void cell_worklist::replace_cell(work_set_int *queue_pointer, int col_old, int col) {
    const int pos = queue_pointer->get(col_old);
    arr[pos] = col;
    assert(queue_pointer->get(col_old) != -1);
    queue_pointer->set(col_old, -1);
    queue_pointer->set(col, pos);
}

bool cell_worklist::empty() {
    return (cur_pos == 0);
}

void cell_worklist::reset(work_set_int *queue_pointer) {
    while(cur_pos > 0) {
        cur_pos--;
        queue_pointer->set(arr[cur_pos], -1);
    }
}

