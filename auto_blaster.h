//
// Created by markus on 23.09.19.
//

#ifndef DEJAVU_AUTO_BLASTER_H
#define DEJAVU_AUTO_BLASTER_H


#include <random>
#include "sgraph.h"
#include "invariant.h"
#include "concurrentqueue.h"
#include "invariant_acc.h"
#include "refinement_bucket.h"
#include "selector.h"
#include "bfs.h"


typedef std::vector<moodycamel::ConcurrentQueue<std::tuple<int, int>>> com_pad;
class pipeline_group;

struct alignas(64) auto_workspace {
    refinement R;
    selector S;
    coloring c;
    invariant I;

    coloring* work_c;
    invariant* work_I;

    work_set first_level_fail;
    work_set first_level_succ;
    int first_level_sz = 0;
    int first_level = 1;
    int base_size = 0;
    int first_level_succ_point = -1;
    int skiplevels = 1;

    int first_skiplevel = 1;
    coloring  skip_c;
    invariant skip_I;

    std::tuple<int, int>* dequeue_space;
    int dequeue_space_sz = 0;

    std::tuple<int, int>* enqueue_space;
    int enqueue_space_sz = 0;

    int measure1 = 0;
    int measure2 = 0;

    // shared state
    moodycamel::ConsumerToken* ctok;
    moodycamel::ProducerToken* ptok;
    std::vector<moodycamel::ProducerToken*> ptoks;

    pipeline_group* G;

    coloring* start_c;
    invariant start_I;
    com_pad* communicator_pad;
    int communicator_id;
    int** shared_orbit;

    // bfs workspace
    bfs* BW;
    std::tuple<bfs_element*, int>* todo_dequeue;
    std::pair<bfs_element *, int>* finished_elements;
    std::tuple<bfs_element *, int>* todo_elements;
    bfs_element* prev_bfs_element = nullptr;
    bool init_bfs = false;
};

#include "pipeline_group.h"

class auto_blaster {
public:
    void sample_pipelined(sgraph *g, bool master, shared_switches* switches, pipeline_group* G, coloring* start_c, bijection* canon_leaf, invariant* canon_I,
                          com_pad* communicator_pad, int communicator_id, int** shared_orbit, bfs* bwork);

private:
    void find_automorphism_prob(auto_workspace *w, sgraph *g, bool compare, invariant *canon_I,
                                bijection *canon_leaf, bijection *automorphism, int *restarts,
                                shared_switches *switches, int selector_seed);

    void fast_automorphism_non_uniform(sgraph *g, bool compare, invariant *canon_I, bijection *canon_leaf,
                                       bijection *automorphism, int *restarts,
                                       bool *done, int selector_seed, auto_workspace *w);

    void find_automorphism_from_bfs(auto_workspace *w, sgraph *g, bool compare, invariant *canon_I,
                                                  bijection *canon_leaf, bijection *automorphism, int *restarts,
                                                  shared_switches *switches, int selector_seed);

    bool proceed_state(auto_workspace* w, sgraph* g, coloring* c, invariant* I, int v);

    bool bfs_chunk(sgraph *g, invariant *canon_I, bijection *canon_leaf, bool *done,
              int selector_seed,
              auto_workspace *w);
};


#endif //DEJAVU_AUTO_BLASTER_H
