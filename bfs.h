//
// Created by markus on 18/10/2019.
//

#ifndef DEJAVU_BFS_H
#define DEJAVU_BFS_H


#include "concurrentqueue.h"
#include "coloring.h"
#include "invariant.h"

class bfs_element {
public:
    // coloring and invariant for the specified path / base
    coloring*  c = nullptr;
    invariant* I = nullptr;
    int* base = nullptr;
    int  base_sz = -1;

    // position of element in level_states of workspace
    int level = -1;
    int id = -1;

    // probability weight for probing in level
    int weight = -1;

    // memory management
    bool init_c = false;
    bool init_I = false;
    bool init_base = false;
    ~bfs_element();
};

class bfs_workspace {
public:
    // level array that keeps a queue with tasks
    // bfs_element* is a pointer to the state where int has to be individualized
    moodycamel::ConcurrentQueue<std::tuple<bfs_element*, int>>* bfs_level_todo;

    // commit finished elements (or nullptr if element was deleted)
    // integer determines how many todos were added for this element
    moodycamel::ConcurrentQueue<std::pair<bfs_element*, int>>* bfs_level_finished_elements;

    // elements of all levels
    // delete elements of level once work of level + 1 is fully commited, then it is safe
    bfs_element*** level_states;
    int*          level_sizes;
    int*          level_expecting_finished;
    bool done = false;

    // current level where work has to be done, and where experimential paths should probe
    std::atomic_int current_level;

    // the level which was determined to be reached via BFS
    std::atomic_int target_level;

    int chunk_size = 4;
};

class bfs {
public:
    bfs_workspace BW;
    bfs();
    void initialize(bfs_element* root_node, int init_c, int domain_size, int base_size);
    void work_queues();
};


#endif //DEJAVU_BFS_H