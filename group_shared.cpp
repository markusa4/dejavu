#include <iostream>
#include <tgmath.h>
#include "group_shared.h"
#include "configuration.h"

group_shared::group_shared(int domain_size) {
    this->domain_size = domain_size;
}

void group_shared::initialize(int domain_size, bijection *base_points) {
    shared_schreier_fails(-1);
    added = 0;

    _ack_done_shared_group = 0;

    dequeue_space    = new std::pair<int,int>[128];
    dequeue_space_sz = 128;

    b         = base_points->map;
    base_size = base_points->map_sz;

    shared_newgroup(&gp, &gens, domain_size);
    shared_getorbits(b, base_size, gp, &gens, domain_size);
}

// master thread manages sifting results to determine abort criterion
void group_shared::manage_results(shared_workspace *switches) {
    int num = sift_results.try_dequeue_bulk(dequeue_space, dequeue_space_sz);
    for(int j = 0; j < num; ++j) {
        switch(dequeue_space[j].first) {
            case SIFT_UNIFORM:
            if(!dequeue_space[j].second) {
                abort_counter += 1;
            } else {
                gens_added += 1;
                if(abort_counter > 0) {
                    config.CONFIG_RAND_ABORT += 1;
                    std::cout << "d: " << config.CONFIG_RAND_ABORT << std::endl;
                }
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

        if(non_uniform_abort_counter > (std::min(switches->tolerance + 1, 5))  && !switches->done_fast) {
            // non-uniform done
            non_uniform_abort_counter = 0;
            switches->done_fast = true;
            break;
        }

        if(abort_counter >= config.CONFIG_RAND_ABORT) {
            // we are done
            switches->done = true;
            break;
        }
    }
}

group_shared::~group_shared() {
    delete[] b;
    delete[] dequeue_space;
    shared_freeschreier(&gp, &gens);
}

bool group_shared::add_permutation(bijection *p, int *idle_ms, bool *done) {
    filterstate state;
    state.ingroup = false;
    state.counts_towards_abort = !p->non_uniform;
    state.level = -1;
    state.stype = p->non_uniform?SIFT_NON_UNIFORM:SIFT_UNIFORM;
    bool result;
    if(!p->foreign_base || base_size < 10) {
        result = mfilterschreier_shared(gp, p->map, &gens, (state.ingroup ? TRUE : FALSE), domain_size + 1,
                                             domain_size, state.level + 1, domain_size + 1,
                                             &state, domain_size + 1);
        sift_results.enqueue(std::pair<sift_type, bool>(state.stype, result));
    } else {
        // finish sift, but return change according to sqrt(base) first levels
        // such that we can switch to more efficient base
        result = mfilterschreier_shared(gp, p->map, &gens, (state.ingroup ? TRUE : FALSE), domain_size + 1,
                                             domain_size, state.level + 1, sqrt(base_size + 1) + 1,
                                             &state, domain_size + 1);
    }
    return result;
}

void group_shared::ack_done_shared() {
    thread_local bool seen = false;
    if(!seen) {
        ++_ack_done_shared_group;
        seen = true;
    }
}

// sift some random elements to make schreier vectors more complete
void group_shared::sift_random() {
    bool result = true;
    int result_buffer = 1;
    while(result || (result_buffer > 0)) {
        if((result_buffer > 0) && !result) {
            --result_buffer;
        }
        random_element re;
        filterstate state;
        state.ingroup = true;
        state.counts_towards_abort = false;
        state.level = -1;
        state.stype = SIFT_RANDOM;
        bool generated = generate_random_element(gp, &gens, domain_size, &re);
        if(!generated) break;
        result = mfilterschreier_shared(gp, re.perm, &gens, TRUE, domain_size + 1,
                                            domain_size, state.level + 1, domain_size + 1,
                                            &state, domain_size + 1);

        if(result && result_buffer < 2) {
            result_buffer++;
        }
        delete[] re.perm;
    }
}

void group_shared::reset_ack_done_shared() {
    _ack_done_shared_group = 0;
}

void group_shared::wait_for_ack_done_shared(int n, bool* escape) {
    while(_ack_done_shared_group != n && !*escape)
        continue;
    return;
}

void group_shared::print_group_size() {
    double grpsize1;
    int grpsize2;
    shared_grouporder(b, base_size, gp, &gens, &grpsize1, &grpsize2, domain_size);
    std::cout << grpsize1 << " * 10^" << (grpsize2) << std::endl;
    //deleteunmarked(&gens);
    std::cout << "Generators: " << shared_schreier_gens(gens) << std::endl;
}

int group_shared::number_of_generators() {
    int k = 0;
    shared_permnode *it = gens;
    if(it == NULL)
        return k;
    do {
        k += 1;
       // std::cout << it->next << std::endl;
        it = it->next;
    } while (it != gens);

    return k;
}
