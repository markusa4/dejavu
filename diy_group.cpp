//
// Created by markus on 24/10/2019.
//

#include <iostream>
#include <tgmath.h>
#include "diy_group.h"
#include "configuration.h"

diy_group::diy_group(int domain_size) {
    this->domain_size = domain_size;
}

void diy_group::initialize(int domain_size, bijection *base_points) {
    mschreier_fails(-1);
    added = 0;
    shared_group_todo = 5;
    shared_group_trigger = false;

    _ack_done_shared_group = 0;

    dequeue_space    = new std::pair<int,int>[128];
    dequeue_space_sz = 128;

    b         = base_points->map;
    base_size = base_points->map_sz;

    mnewgroup(&gp, &gens, domain_size);
    mgetorbits(b, base_size, gp, &gens, domain_size);
}

void diy_group::manage_results(shared_switches *switches) {
    int num = sift_results.try_dequeue_bulk(dequeue_space, dequeue_space_sz);
    for(int j = 0; j < num; ++j) {
        switch(dequeue_space[j].first) {
            case SIFT_UNIFORM:
            if(!dequeue_space[j].second) {
                abort_counter += 1;
            } else {
                gens_added += 1;
                abort_counter = 0;
            }
            break;
            case SIFT_RANDOM:
            if(!dequeue_space[j].second) {
                random_abort_counter += 1;
            } else {
                gens_added += 1;
                random_abort_counter = 0;
            }
            break;
            case SIFT_NON_UNIFORM:
            if(!dequeue_space[j].second) {
                non_uniform_abort_counter += 1;
            } else {
                gens_added += 1;
                non_uniform_abort_counter = 0;
            }
            break;
        }
        if(non_uniform_abort_counter > (std::max(switches->tolerance - 1, 10))  && !switches->done_fast) {
            switches->done_fast = true;
            break;
        }
        if(abort_counter >= config.CONFIG_RAND_ABORT) {
            switches->done = true;
            break;
        }
    }
}

diy_group::~diy_group() {
    delete[] b;
    mfreeschreier(&gp, &gens);
}

bool diy_group::add_permutation(bijection *p, int *idle_ms, bool *done) {
    filterstate state;
    state.ingroup = false;
    state.counts_towards_abort = !p->non_uniform;
    state.level = -1;
    state.stype = p->non_uniform?SIFT_NON_UNIFORM:SIFT_UNIFORM;
    bool result;
    if(!p->foreign_base) {
        result = mfilterschreier_shared(gp, p->map, &gens, (state.ingroup ? TRUE : FALSE), domain_size + 1,
                                             domain_size, state.level + 1, domain_size + 1, &state, domain_size + 1);
        sift_results.enqueue(std::pair<sift_type, bool>(state.stype, result));
    } else {
        // finish sift, but return change according to sqrt(base) first levels such that we can switch to mor efficient base
        result = mfilterschreier_shared(gp, p->map, &gens, (state.ingroup ? TRUE : FALSE), domain_size + 1,
                                             domain_size, state.level + 1, sqrt(base_size + 1) + 1, &state, domain_size + 1); // sqrt(base_size + 1) + 1
        //sift_results.enqueue(std::pair<sift_type, bool>(state.stype, result));
    }
    return result;
}

void diy_group::ack_done_shared() {
    thread_local bool seen = false;
    if(!seen) {
        ++_ack_done_shared_group;
        seen = true;
    }
}

void diy_group::reset_ack_done_shared() {
    _ack_done_shared_group = 0;
}

void diy_group::wait_for_ack_done_shared(int n) {
    while(_ack_done_shared_group != n)
        continue;
    return;
}

void diy_group::print_group_size() {
    double grpsize1;
    int grpsize2;
    mgrouporder(b, base_size, gp,&gens, &grpsize1, &grpsize2, domain_size);
    std::cout << grpsize1 << " * 10^" << (grpsize2) << std::endl;
    //deleteunmarked(&gens);
    std::cout << "Generators: " << mschreier_gens(gens) << std::endl;
}

int diy_group::number_of_generators() {
    int k = 0;
    mpermnode *it = gens;
    if(it == NULL)
        return k;
    do {
        k += 1;
       // std::cout << it->next << std::endl;
        it = it->next;
    } while (it != gens);

    return k;
}
