//
// Created by markus on 01.10.19.
//

#ifndef DEJAVU_PIPELINE_GROUP_H
#define DEJAVU_PIPELINE_GROUP_H


#include "pipeline_schreier.h"
#include "bijection.h"
#include "concurrentqueue.h"
#include "blockingconcurrentqueue.h"

class pipeline_group {
public:
    // input and pipeline management
    moodycamel::ConcurrentQueue<bijection> automorphisms = moodycamel::ConcurrentQueue<bijection>(50, 4, 4);
    std::vector<moodycamel::ConcurrentQueue<filterstate>> pipeline_queues;
    moodycamel::ConcurrentQueue<std::pair<bool, bool>> sift_results = moodycamel::ConcurrentQueue<std::pair<bool, bool>>(20, 1, 1);
    std::vector<int> intervals;
    std::vector<std::thread> work_threads;
    moodycamel::ConsumerToken ctoken = moodycamel::ConsumerToken(automorphisms);
    int stages;

    // group structures
    int domain_size;
    int  base_size;
    int* b;
    int added;
    mschreier *gp;
    mpermnode *gens;

    pipeline_group();
    void initialize(int domain_size, bijection *base_points, int stages);
    ~pipeline_group();
    bool add_permutation(bijection* p, int* idle_ms, bool* done);
    void print_group_size();
    void pipeline_stage(int n, bool* done, bool* done_fast);
    void launch_pipeline_threads(bool *done, bool* done_fast);
    void join_threads();

    void determine_stages();
};


#endif //DEJAVU_PIPELINE_GROUP_H
