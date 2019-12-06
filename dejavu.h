#ifndef DEJAVU_DEJAVU_H
#define DEJAVU_DEJAVU_H


#include <random>
#include "sgraph.h"
#include "invariant.h"
#include "concurrentqueue.h"
#include "selector.h"
#include "bfs.h"
#include "schreier_shared.h"
#include "group_shared.h"
#include "schreier_sequential.h"

struct abort_code {
    abort_code()=default;
    abort_code(int reason):reason(reason){};
    int reason = 0;
};

struct base_cache {
    coloring*       base_cache_colorings;
    invariant*      base_cache_invariants;
    bool*           base_cache_done;
    std::atomic_int level_done;
};

struct alignas(64) dejavu_workspace {
    // workspace for normal search
    refinement R;
    selector   S;
    coloring   c;
    invariant  I;

    // workspace for bfs_workspace
    coloring*  work_c;
    invariant* work_I;

    // workspace for base aligned search
    int first_level = 1;
    int base_size = 0;
    int skiplevels = 1;
    int        first_skiplevel = 1;
    coloring   skip_c;
    invariant  skip_I;
    mschreier* skip_schreier_level;
    bool       skiplevel_is_uniform = false;

    int* my_base_points;
    int  my_base_points_sz;
    bool is_foreign_base;

    group_shared*      G;

    coloring* start_c;
    invariant start_I;

    // indicates which thread this is
    int id;

    // shared orbit and generators
    int**       shared_orbit;
    int**       shared_orbit_weights;
    mpermnode** shared_generators;
    int*        shared_generators_size;
    int         generator_fix_base_alloc = -1;

    // sequential, local group
    permnode*       sequential_gens;
    _schreierlevel* sequential_gp;
    bool            sequential_init = false;

    // deprecated workspace for simple orbit method
    work_set  orbit_considered;
    work_list orbit_vertex_worklist;
    work_list orbit;
    int canonical_v;
    mpermnode** generator_fix_base;
    int         generator_fix_base_size;

    // bfs_workspace workspace
    bfs_workspace* BW;
    std::tuple<bfs_element*, int, int>* todo_dequeue;
    int todo_deque_sz        = -1;
    std::pair<bfs_element *, int>* finished_elements;
    int finished_elements_sz = -1;
    std::pair<bfs_element *, int>* todo_elements;
    int todo_elements_sz     = -1;
    bfs_element* prev_bfs_element = nullptr;
    bool init_bfs = false;

    ~dejavu_workspace() {
        if(init_bfs) {
            delete[] todo_dequeue;
            delete[] todo_elements,
                    delete[] finished_elements;
            delete[] generator_fix_base;
        }
        if(sequential_init) {
            _freeschreier(&sequential_gp, &sequential_gens);
        }

        delete work_c;
        delete work_I;
        delete start_c;
    };
};

class dejavu {
public:
    void automorphisms(sgraph *g, mpermnode **gens);

private:
    void worker_thread(sgraph *g_, bool master, shared_workspace *switches, group_shared *G, coloring *start_c,
                       strategy* canon_strategy, int communicator_id,
                       int **shared_orbit, int** shared_orbit_weights, bfs_workspace *bwork, mpermnode **gens, int *shared_group_size);

    void find_first_leaf(dejavu_workspace *w, sgraph *g, bool compare, invariant *canon_I,
                         bijection *canon_leaf, strategy* canon_strategy, bijection *automorphism, int *restarts,
                         shared_workspace *switches, int selector_seed);

    abort_code base_aligned_search(dejavu_workspace *w, sgraph *g, strategy *canon_strategy, bijection *automorphism,
                             strategy_metrics *m, bool *done, shared_workspace *switches, int selector_seed);

    void reset_skiplevels(dejavu_workspace *w);

    abort_code uniform_from_bfs_search(dejavu_workspace *w, sgraph *g, bool compare, strategy* canon_strategy, bijection *automorphism, int *restarts,
                                       shared_workspace *switches, int selector_seed);

    bool proceed_state(dejavu_workspace* w, sgraph* g, coloring* c, invariant* I, int v, change_tracker* changes, strategy_metrics* m);

    bool get_orbit(dejavu_workspace *w, int *base, int base_sz, int v, int v_base, work_list *orbit, bool reuse_generators);

    bool bfs_chunk(dejavu_workspace *w, sgraph *g, strategy *canon_strategy, bool *done, int selector_seed);

    void bfs_reduce_tree(dejavu_workspace *w);

    void bfs_fill_queue(dejavu_workspace *w);

    void bfs_assure_init(dejavu_workspace *w);

    bool sequential_init_copy(dejavu_workspace *w);

    bool uniform_from_bfs_search_with_storage(dejavu_workspace *w, sgraph *g, shared_workspace* switches, bfs_element *elem, int selector_seed, strategy *strat,
                                              bijection *automorphism, bool look_close);
};


#endif //DEJAVU_DEJAVU_H
