//
// Created by markus on 24/10/2019.
//

#ifndef DEJAVU_DIY_GROUP_H
#define DEJAVU_DIY_GROUP_H


#include "pipeline_schreier.h"
#include "bijection.h"
#include "concurrentqueue.h"
#include "blockingconcurrentqueue.h"
#include "utility.h"

struct auto_workspace;

class diy_group {
public:
    // input and pipeline management
    moodycamel::ConcurrentQueue<std::pair<sift_type, bool>> sift_results        = moodycamel::ConcurrentQueue<std::pair<sift_type, bool>>(20, 16, 1);
    moodycamel::ConsumerToken                          sift_results_ctoken = moodycamel::ConsumerToken(sift_results);

    // group structures
    int  domain_size;
    int  base_size;
    int* b;
    int  added;
    std::atomic<int> _ack_done_shared_group;

    std::pair<int,int>* dequeue_space;
    int dequeue_space_sz;

    // abort criterion
    int abort_counter = 0;
    int random_abort_counter = 0;
    int non_uniform_abort_counter = 0;
    int gens_added = 0;

    std::atomic_int  shared_group_todo;
    std::atomic_bool shared_group_trigger;

    mschreier *gp;
    mpermnode *gens;

    diy_group(int domain_size);
    void initialize(int domain_size, bijection *base_points);
    ~diy_group();
    bool add_permutation(bijection* p, int* idle_ms, bool* done);
    void print_group_size();
    void manage_results(shared_switches* switches);

    void wait_for_ack_done_shared(int n);

    void ack_done_shared();
};



#endif //DEJAVU_DIY_GROUP_H
