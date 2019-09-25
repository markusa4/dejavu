//
// Created by markus on 23.09.19.
//

#ifndef DEJAVU_AUTO_BLASTER_H
#define DEJAVU_AUTO_BLASTER_H


#include <random>
#include "sgraph.h"
#include "invariant.h"
#include "concurrentqueue.h"

class auto_blaster {
    moodycamel::ConcurrentQueue<bijection> Q;
    invariant start_I;
    coloring start_c;
public:
    void find_automorphism(sgraph *g, bool compare, invariant* canon_I, bijection* canon_leaf, bijection* automorphism, std::default_random_engine* re);
    void sample(sgraph *g, bool master, bool* done);

    void
    find_automorphism_prob(sgraph *g, bool compare, invariant *canon_I, bijection *canon_leaf, bijection *automorphism,
                           std::default_random_engine *re, int *restarts);
};


#endif //DEJAVU_AUTO_BLASTER_H
