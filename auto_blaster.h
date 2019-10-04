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

class auto_blaster {
    moodycamel::ConcurrentQueue<bijection> Q;
    invariant start_I;
    coloring start_c;
public:
    void sample(sgraph* g, bool master, bool* done);

    void
    find_automorphism_prob(sgraph *g, bool compare, invariant *canon_I, bijection *canon_leaf, bijection *automorphism,
                           std::default_random_engine *re, int *restarts, bool* done);

    void sample_pipelined(sgraph *g, bool master, bool *done, pipeline_group* G);

    void
    find_automorphism_bt(sgraph *g, bool compare, invariant *canon_I, bijection *canon_leaf, bijection *automorphism,
                         std::default_random_engine *re, int *restarts, bool *done);
};


#endif //DEJAVU_AUTO_BLASTER_H
