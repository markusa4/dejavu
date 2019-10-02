//
// Created by markus on 01.10.19.
//

#ifndef DEJAVU_PIPELINE_GROUP_H
#define DEJAVU_PIPELINE_GROUP_H


#include "my_schreier.h"
#include "bijection.h"
#include "concurrentqueue.h"

class pipeline_group {
public:
    // input and pipeline management
    moodycamel::ConcurrentQueue<bijection> automorphisms;
    std::vector<moodycamel::ConcurrentQueue<filterstate>> pipeline_queues;
    moodycamel::ConcurrentQueue<std::pair<bool, bool>> sift_results;
    std::vector<int> intervals;
    std::vector<std::thread> work_threads;
    int stages;

    // group structures
    int domain_size;
    int base_size;
    int* b;
    int added;
    mschreier *gp;
    mpermnode *gens;

    pipeline_group(int domain_size, bijection* base_points, int stages);
    ~pipeline_group();
    bool add_permutation(bijection* p);
    void print_group_size();
    void pipeline_stage(int n, bool* done);
    void launch_pipeline_threads(bool *done);
    void join_threads();
};


#endif //DEJAVU_PIPELINE_GROUP_H
