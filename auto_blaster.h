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
#include "pipeline_group.h"
#include "refinement_bucket.h"

class auto_blaster {
    moodycamel::ConcurrentQueue<bijection> Q;
    invariant start_I;
    coloring start_c;
    coloring_bucket start_cb;
public:
    void sample(sgraph* g, bool master, bool* done);

    void
    find_automorphism_prob(sgraph *g, bool compare, invariant *canon_I, bijection *canon_leaf, bijection *automorphism,
                           std::default_random_engine *re, int *restarts, bool* done, int selector_seed, work_set* first_level_fail);

    void
    find_automorphism_prob_bucket(sgraph* g, bool compare, invariant* canon_I, bijection* canon_leaf,
                                                     bijection* automorphism, std::default_random_engine* re, int *restarts, bool *done, int selector_seed,  refinement_bucket* R);

    void sample_pipelined(sgraph *g, bool master, bool *done, pipeline_group* G);

    void
    find_automorphism_bt(sgraph *g, bool compare, invariant *canon_I, bijection *canon_leaf, bijection *automorphism,
                         std::default_random_engine *re, int *restarts, bool *done, int selector_seed);

    void sample_pipelined_bucket(sgraph *g, bool master, bool *done, pipeline_group *G);
};


#endif //DEJAVU_AUTO_BLASTER_H
