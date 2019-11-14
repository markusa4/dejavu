//
// Created by markus on 23.09.19.
//

#ifndef DEJAVU_DEJAVU_H
#define DEJAVU_DEJAVU_H


#include <random>
#include "sgraph.h"
#include "invariant.h"
#include "concurrentqueue.h"
#include "selector.h"
#include "bfs.h"
#include "schreier_shared.h"
#include "group_diy.h"

struct alignas(64) dejavu_workspace {
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

    group_diy*      G_;

    coloring* start_c;
    invariant start_I;
    int       communicator_id;

    // shared orbit and generators
    int**       shared_orbit;
    int**       shared_orbit_weights;
    mpermnode** shared_generators;
    int*        shared_generators_size;
    int         generator_fix_base_alloc = -1;

    //
    work_set  orbit_considered;
    work_list orbit_vertex_worklist;
    work_list orbit;
    int canonical_v;
    mpermnode** generator_fix_base;
    int         generator_fix_base_size;

    // bfs workspace
    bfs* BW;
    std::tuple<bfs_element*, int, int>* todo_dequeue;
    int todo_deque_sz        = -1;
    std::pair<bfs_element *, int>* finished_elements;
    int finished_elements_sz = -1;
    std::pair<bfs_element *, int>* todo_elements;
    int todo_elements_sz     = -1;
    bfs_element* prev_bfs_element = nullptr;
    bool init_bfs = false;
};

class dejavu {
public:

    void sample_shared(sgraph *g_, bool master, shared_switches *switches, group_diy *G, coloring *start_c,
                       strategy* canon_strategy, int communicator_id,
                       int **shared_orbit, int** shared_orbit_weights, bfs *bwork, mpermnode **gens, int *shared_group_size);

private:
    void find_automorphism_prob(dejavu_workspace *w, sgraph *g, bool compare, invariant *canon_I,
                                bijection *canon_leaf, strategy* canon_strategy, bijection *automorphism, int *restarts,
                                shared_switches *switches, int selector_seed);

    void fast_automorphism_non_uniform(sgraph *g, bool compare, strategy* canon_strategy,
                                       bijection *automorphism, strategy_metrics *m,
                                       bool *done, int selector_seed, dejavu_workspace *w, int tolerance);

    void find_automorphism_from_bfs(dejavu_workspace *w, sgraph *g, bool compare, strategy* canon_strategy, bijection *automorphism, int *restarts,
                                    shared_switches *switches, int selector_seed);

    bool proceed_state(dejavu_workspace* w, sgraph* g, coloring* c, invariant* I, int v, change_tracker* changes);

    bool bfs_chunk(sgraph *g, strategy* canon_strategy, bool *done,
                   int selector_seed,
                   dejavu_workspace *w);

    bool get_orbit(dejavu_workspace *w, int *base, int base_sz, int v, int v_base, work_list *orbit, bool reuse_generators);

    void reset_skiplevels(dejavu_workspace *w);

    void bfs_reduce_tree(dejavu_workspace *w);

    void bfs_fill_queue(dejavu_workspace *w);

    void bfs_assure_init(dejavu_workspace *w);
};


#endif //DEJAVU_DEJAVU_H
