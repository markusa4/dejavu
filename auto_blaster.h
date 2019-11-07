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
#include "pipeline_schreier.h"
#include "diy_group.h"


typedef std::vector<moodycamel::ConcurrentQueue<std::pair<int, int>>> com_pad;
class pipeline_group;

struct alignas(64) auto_workspace {
    refinement R;
    selector   S;
    coloring   c;
    invariant  I;

    coloring*  work_c;
    invariant* work_I;

    work_set first_level_fail;
    work_set first_level_succ;
    int first_level_sz = 0;
    int first_level = 1;
    int base_size = 0;
    int first_level_succ_point = -1;
    int skiplevels = 1;

    int        first_skiplevel = 1;
    coloring   skip_c;
    invariant  skip_I;
    mschreier* skip_schreier_level;

    std::pair<int, int>* dequeue_space;
    int dequeue_space_sz = 0;

    std::pair<int, int>* enqueue_space;
    int enqueue_space_sz = 0;

    int measure1 = 0;
    int measure2 = 0;

    int* my_base_points;
    int  my_base_points_sz;
    bool is_foreign_base;

    // shared state
    moodycamel::ConsumerToken* ctok;
    moodycamel::ProducerToken* ptok;
    std::vector<moodycamel::ProducerToken*> ptoks;

    pipeline_group* G;
    diy_group*      G_;

    coloring* start_c;
    invariant start_I;
    com_pad*  communicator_pad;
    int       communicator_id;

    // shared orbit and generators
    int**       shared_orbit;
    int**       shared_orbit_weights;
    mpermnode** shared_generators;
    int*        shared_generators_size;

    //
    work_set  orbit_considered;
    work_list orbit_vertex_worklist;
    work_list orbit;
    int canonical_v;
    mpermnode** generator_fix_base;
    int         generator_fix_base_size;

    // bfs workspace
    bfs* BW;
    std::pair<bfs_element*, int>* todo_dequeue;
    int todo_deque_sz        = -1;
    std::pair<bfs_element *, int>* finished_elements;
    int finished_elements_sz = -1;
    std::pair<bfs_element *, int>* todo_elements;
    int todo_elements_sz     = -1;
    bfs_element* prev_bfs_element = nullptr;
    change_tracker changes;
    bool init_bfs = false;
};

#include "pipeline_group.h"

class auto_blaster {
public:
    void sample_pipelined(sgraph *g, bool master, shared_switches* switches, pipeline_group* G, coloring* start_c, bijection* canon_leaf, invariant* canon_I,
                          com_pad* communicator_pad, int communicator_id, int** shared_orbit, bfs* bwork, mpermnode** gens, int* shared_group_identity);

    void sample_shared(sgraph *g_, bool master, shared_switches *switches, diy_group *G, coloring *start_c,
                       bijection **canon_leaf, invariant **canon_I, com_pad *communicator_pad, int communicator_id,
                       int **shared_orbit, int** shared_orbit_weights, bfs *bwork, mpermnode **gens, int *shared_group_size);

private:
    void find_automorphism_prob(auto_workspace *w, sgraph *g, bool compare, invariant *canon_I,
                                bijection *canon_leaf, bijection *automorphism, int *restarts,
                                shared_switches *switches, int selector_seed);

    void fast_automorphism_non_uniform(sgraph *g, bool compare, invariant *canon_I, bijection *canon_leaf,
                                       bijection *automorphism, int *restarts,
                                       bool *done, int selector_seed, auto_workspace *w, int tolerance);

    void find_automorphism_from_bfs(auto_workspace *w, sgraph *g, bool compare, invariant *canon_I,
                                                  bijection *canon_leaf, bijection *automorphism, int *restarts,
                                                  shared_switches *switches, int selector_seed);

    bool proceed_state(auto_workspace* w, sgraph* g, coloring* c, invariant* I, int v, change_tracker* changes);

    bool bfs_chunk(sgraph *g, invariant *canon_I, bijection *canon_leaf, bool *done,
              int selector_seed,
              auto_workspace *w);

    bool get_orbit(auto_workspace *w, int *base, int base_sz, int v, int v_base, work_list *orbit, bool reuse_generators);

    void fast_automorphism_non_uniform_from_bfs(auto_workspace *w, sgraph *g, bool compare, invariant *canon_I,
                                                bijection *canon_leaf, bijection *automorphism, int *restarts,
                                                bool *done,
                                                int selector_seed);

    void reset_skiplevels(auto_workspace *w);

    void bfs_reduce_tree(auto_workspace *w);
};


#endif //DEJAVU_AUTO_BLASTER_H
