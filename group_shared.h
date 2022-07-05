#ifndef DEJAVU_GROUP_SHARED_H
#define DEJAVU_GROUP_SHARED_H


#include "schreier_shared.h"
#include "bijection.h"
#include "concurrentqueue.h"
#include "utility.h"
#include "configuration.h"

template <class vertex_t, class degree_t, class edge_t>
struct dejavu_workspace;

template <class vertex_t>
class group_shared {
public:
    // input and pipeline management
    moodycamel::ConcurrentQueue<std::pair<sift_type, bool>> sift_results
        = moodycamel::ConcurrentQueue<std::pair<sift_type, bool>>(20, 16, 1);
    moodycamel::ConsumerToken sift_results_ctoken = moodycamel::ConsumerToken(sift_results);

    // group structures
    int  domain_size;
    int  base_size;
    int* b;
    int  added;
    std::atomic<int> _ack_done_shared_group;

    std::pair<int,int>* dequeue_space;
    int dequeue_space_sz;

    // abort criterion
    int abort_counter             = 0;
    int random_abort_counter      = 0;
    int non_uniform_abort_counter = 0;
    int gens_added = 0;

    bool generators_persistent = false;

    shared_schreier *gp;
    shared_permnode *gens;

    group_shared(int domain_size) {
        this->domain_size = domain_size;
    }
    void initialize(int domain_size, bijection<vertex_t> *base_points) {
        shared_schreier_fails(-1);
        added = 0;

        _ack_done_shared_group = 0;

        dequeue_space    = new std::pair<int,int>[128];
        dequeue_space_sz = 128;

        b = new int[base_points->map_sz];
        base_points->copy_map(b);
        base_size = base_points->map_sz;

        shared_newgroup(&gp, &gens, domain_size);
        shared_getorbits(b, base_size, gp, &gens, domain_size);
    }

    ~group_shared() {
        delete[] b;
        delete[] dequeue_space;
        if(generators_persistent) {
            shared_freeschreier(&gp, nullptr);
            shared_schreier_freedyn();
        } else {
            shared_freeschreier(&gp, &gens);
            shared_schreier_freedyn();
        }
    }

    bool add_permutation(bijection<vertex_t>* p, int* idle_ms, bool* done) {
        int* map = p->extract_map();
        filterstate state;
        state.ingroup = false;
        state.counts_towards_abort = !p->non_uniform;
        state.level = -1;
        state.stype = p->non_uniform?SIFT_NON_UNIFORM:SIFT_UNIFORM;
        bool result;
        if(!p->foreign_base || base_size < 10) {
            result = mfilterschreier_shared(gp, map, &gens, (state.ingroup ? true : false), domain_size + 1,
                                            domain_size, state.level + 1, domain_size + 1,
                                            &state, domain_size + 1);
            sift_results.enqueue(std::pair<sift_type, bool>(state.stype, result)); // ToDo: replace sifting queue with counter!
        } else {
            // finish sift, but return change according to sqrt(base) first levels
            // such that we can switch to more efficient base
            result = mfilterschreier_shared(gp, map, &gens, (state.ingroup ? true : false), domain_size + 1,
                                            domain_size, state.level + 1, sqrt(base_size + 1) + 1,
                                            &state, domain_size + 1);
        }
        return result;
    }

    void print_group_size() {
        double grpsize1;
        int grpsize2;
        shared_grouporder(b, base_size, gp, &gens, &grpsize1, &grpsize2, domain_size);
        PRINT(grpsize1 << " * 10^" << (grpsize2));
        //deleteunmarked(&gens);
        PRINT("Generators: " << shared_schreier_gens(gens));
    }

    void print_group_size_stdout() {
        double grpsize1;
        int grpsize2;
        shared_grouporder(b, base_size, gp, &gens, &grpsize1, &grpsize2, domain_size);
        std::cout << grpsize1 << " * 10^" << (grpsize2) << std::endl;
        //deleteunmarked(&gens);
        PRINT("Generators: " << shared_schreier_gens(gens));
    }

    void compute_group_size(double* grpsize1, int* grpsize2) {
        shared_grouporder(b, base_size, gp, &gens, grpsize1, grpsize2, domain_size);
    }

    void manage_results(shared_workspace_auto<vertex_t>* switches) {
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
                            PRINT("[Dej] d: " << config.CONFIG_RAND_ABORT);
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
    int  number_of_generators() {
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

    void wait_for_ack_done_shared(int n, bool* escape) {
        while(_ack_done_shared_group != n && !*escape)
            continue;
        return;
    }

    void ack_done_shared(bool* seen) {
        if(!*seen) {
            ++_ack_done_shared_group;
            *seen = true;
        }
    }

    void reset_ack_done_shared() {
        _ack_done_shared_group = 0;
    }

    void sift_random() {
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
            result = mfilterschreier_shared(gp, re.perm, &gens, true, domain_size + 1,
                                            domain_size, state.level + 1, domain_size + 1,
                                            &state, domain_size + 1);
            if(result && result_buffer < 2) {
                result_buffer++;
            }
        }
    }
};



#endif //DEJAVU_GROUP_SHARED_H
