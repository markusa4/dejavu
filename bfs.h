//
// Created by markus on 18/10/2019.
//

#ifndef DEJAVU_BFS_H
#define DEJAVU_BFS_H


#include "concurrentqueue.h"
#include "coloring.h"
#include "invariant.h"

struct bfs_element {
    coloring*  c;
    invariant* I;
    int* base = nullptr;
    int level = -1;
    int id = -1;
    int weight = -1;
};

struct bfs_workspace {
    moodycamel::ConcurrentQueue<std::tuple<int, int>>* bfs_level_todo;
    moodycamel::ConcurrentQueue<bfs_element*>*         bfs_level_finished_elements;
    moodycamel::ConcurrentQueue<int>*                  bfs_commit_finished_work;

    bfs_element** level_states;
    int* level_sizes;
    std::atomic_int current_level;
};


#endif //DEJAVU_BFS_H
