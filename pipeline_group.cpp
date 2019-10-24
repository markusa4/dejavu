//
// Created by markus on 01.10.19.
//

#include <iostream>
#include "pipeline_group.h"
#include "configuration.h"

extern std::mutex circ_mutex;

void pipeline_group::launch_pipeline_threads(shared_switches* switches, auto_workspace* w) {
    for(int i = 1; i < stages; ++i)
        work_threads.emplace_back(std::thread(&pipeline_group::pipeline_stage, this, i, switches, w));
}

void pipeline_group::join_threads() {
    while(!work_threads.empty()) {
        work_threads[work_threads.size()-1].join();
        work_threads.pop_back();
    }
}

void pipeline_group::pipeline_stage(int n, shared_switches* switches, auto_workspace* w) {
   bool* done = &switches->done;
   bool* done_fast = &switches->done_fast;
    int abort_counter = 0;
    int leafs_considered = 0;
    int random_abort_counter = 0;
    int bulk_deque_num = 8;
    int gens_added = 0;
    std::pair<bool, bool>* res = new std::pair<bool, bool>[bulk_deque_num];

    std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();

    if(base_size == 0)
        *done = true;

    bool is_first_stage = (n == 0);
    bool is_last_stage = (n == stages - 1);
    bool is_only_stage = is_first_stage && is_last_stage;
    int my_interval = intervals[n];

    int first_cell_size = -1;

    if(is_first_stage) {
        delete[] w->dequeue_space;
        delete[] w->enqueue_space;
        w->dequeue_space = new std::pair<int, int>[1024];
        w->dequeue_space_sz = 1024;
        w->enqueue_space = new std::pair<int, int>[2049];
        w->enqueue_space_sz = 2049;
        w->first_level = 1;
        int c = w->S.select_color_largest(w->start_c);
        first_cell_size = w->start_c->ptn[c] + 1;
        //std::cout << first_cell_size << std::endl;
    }
    int max_it = 0;

    while(!(*done)) {
        if(n == 0 && *done_fast && !(switches->done_shared_group)) {
            // copy gens and first orbit for shared use!
            circ_mutex.lock();
            *w->shared_orbit = new int[domain_size];
            memcpy(*w->shared_orbit, gp->orbits, domain_size * sizeof(int));



            if(gens_added > 0) {
                *w->shared_generators_size = 0;
                mpermnode *it = gens;
                do {
                    maddpermutation(w->shared_generators, it->p, domain_size);
                    *w->shared_generators_size += 1;
                    it = it->next;
                } while (it != gens);
            } else {
                *w->shared_generators_size = 0;
            }
            //std::cout << std::endl;
            circ_mutex.unlock();
            switches->done_shared_group.store(true);
            // need to determine target level!
        }

        if(n == 0)
            w->BW->work_queues();

        // work on pipeline_results and track
        filterstate state;
        bijection p;
        random_element re;
        int *_p = nullptr;
        bool is_random = false;

        // context switch
        if (is_first_stage) {
            // collect results and determine whether done
            int num = sift_results.try_dequeue_bulk(res, bulk_deque_num);
            for(int j = 0; j < num; ++j) {
                if(!res[j].first) {
                   // std::cout << leafs_considered << std::endl;
                    leafs_considered += 1;
                    if(!res[j].second) {
                        abort_counter += 1;
                    } else {
                        gens_added += 1;
                        abort_counter = 0;
                    }
                } else {
                    if(!res[j].second) {
                        random_abort_counter += 1;
                    } else {
                        gens_added += 1;
                        random_abort_counter = 0;
                    }
                }
                if(abort_counter >= config.CONFIG_RAND_ABORT) {
                    *done = true;
                    break;
                }
                //std::cout << abort_counter << "; " << random_abort_counter << std::endl;
                //d = sift_results.try_dequeue(res);
            }
            if(*done) {break;}

            // not done, so choose next element for the pipeline
            bool d = automorphisms.try_dequeue(ctoken, p);
            // automorphisms non-empty? take that element
            if(!d && random_abort_counter >= config.CONFIG_RAND_ABORT_RAND)  continue;
                //d = automorphisms.try_dequeue(ctoken, p);
            //}

            if(d) {
                double cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
                //std::cout << "[T] Automorphism arrived: " << cref / 1000000.0 << "ms" << std::endl;
            }

            if (d) {
                //std::cout << "automorphism in pipeline" << std::endl;
                _p = p.map;
                        /*new int[domain_size];
                for (int k = 0; k < p.map.size(); ++k) {
                    _p[k] = p.map[k];
                }*/
                is_random = false;
                state.ingroup = false;
                state.counts_towards_abort = !p.non_uniform;
                state.level = -1;
                random_abort_counter = 0;
            } else {
                //std::cout << "gen ran" << std::endl;
                // else generate random element and sift that
                bool gen = generate_random_element(gp, &gens, domain_size, &re);
                if (!gen) {continue;}
                _p = re.perm;
                is_random = true;
                state.ingroup = true;
                state.counts_towards_abort = true;
                state.level = -1;
            }
        } else {
            // take element from intermediate queue
            bool d = pipeline_queues[n].try_dequeue(state);
            while (!d && (!(*done))) {
                /*if(front_idle_ms % 10000 == 0) {
                    std::cout << "Pipeline(" << n << ") front idle " << front_idle_ms << ", "
                              << pipeline_queues[n].size_approx() << std::endl;
                }*/
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                //std::this_thread::yield();
                //front_idle_ms += 1;
                d = pipeline_queues[n].try_dequeue(state);
            }
            if(*done) {break;}
        }
        // perform filterschreier (partial)
        bool result = mfilterschreier_interval(gp, _p, &gens, (state.ingroup?TRUE:FALSE), domain_size + 1, domain_size, state.level + 1,
                                 my_interval, &state);

        if(is_first_stage) {
            if(is_random) {
                free_random_element(&re);
            } else {
                delete[] _p;
            }
        }

        if (is_last_stage) {
            shared_group_todo -= 1;
            if((!state.counts_towards_abort && !result && !state.ingroup && !(*done_fast)) || (shared_group_trigger && shared_group_todo < 0)) {
                //std::cout << "useless element" << result << std::endl;
                *done_fast = true;
               // std::cout << "random element: " << result << std::endl;
            }
            sift_results.enqueue(std::pair<bool, bool>(state.ingroup || !state.counts_towards_abort, result));
        } else {
            while(pipeline_queues[n + 1].size_approx() > 50 && (!(*done))) {
                //std::cout << "throttle" << n << std::endl;
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            pipeline_queues[n + 1].enqueue(state);
        }
    }
    if(n == 0 ){
        std::cout << "Leafs considered: " << leafs_considered << std::endl;
        std::cout << "Leftover: " << automorphisms.size_approx() << std::endl;
    }
    //std::cout << "Pipeline stage(" << n << ") idle: " << front_idle_ms << "ms / " << back_idle_ms << "ms" << std::endl;
}

void pipeline_group::skip_shared_group() {
    shared_group_todo.store(-1);
    shared_group_trigger.store(true);
}

bool pipeline_group::add_permutation(bijection *p, int* idle_ms, bool* done) {
    //std::cout << "enqueued" << std::endl;
    static thread_local moodycamel::ProducerToken ptoken = moodycamel::ProducerToken(automorphisms);

    /*if(!shared_group_trigger && p->non_uniform && automorphisms.size_approx() > 20) {
        shared_group_todo.store(20);
        shared_group_trigger.store(true);
        return false;
    }*/

    while(automorphisms.size_approx() > 100 && (!(*done))) {
        //std::cout << "throttle" << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        *idle_ms += 1;
    }
    automorphisms.enqueue(ptoken, *p);
    return true;
}

pipeline_group::pipeline_group(int domain_size) {
    this->domain_size = domain_size;
    first_level_points = moodycamel::ConcurrentQueue<std::tuple<int, int>>(domain_size, config.CONFIG_THREADS_REFINEMENT_WORKERS, config.CONFIG_THREADS_REFINEMENT_WORKERS);
}

void pipeline_group::initialize(int domain_size, bijection *base_points, int stages) {
    mschreier_fails(-1);
    added = 0;
    shared_group_todo = 5;
    shared_group_trigger = false;
    //this->base_size = 0;
    this->stages = stages;

    //std::cout << "Creating new pipeline_group... " << std::endl;
    mnewgroup(&gp, &gens, domain_size);
    mgetorbits(b, base_size, gp, &gens, domain_size);

    // create pipeline
    determine_stages();
}

// ToDo: do this better, more drastic towards start
// ToDo: remove unncessary stages (size = 0)
void pipeline_group::determine_stages() {
    std::cout << "Pipeline: ";
    int ded = 0;
    for(int i = 0; i < stages; ++i) {
        int stage_pos = (base_size / stages) * (i + 1) - (i + 1) * ((base_size / (stages* (i + 2))));
       // if(i == 0)
        //    stage_pos = stage_pos / 2;

        // skip stage if too small
        if((!intervals.empty()) &&
            (stage_pos - intervals[intervals.size() - 1] < config.CONFIG_THREADS_PIPELINE_STAGE_MIN)) {
            ded += 1;
            continue;
        }

        pipeline_queues.emplace_back(moodycamel::ConcurrentQueue<filterstate>(20));
        intervals.push_back(stage_pos);
        std::cout << "/" << intervals[intervals.size() - 1];
    }
    stages -= ded;
     std::cout << "(" <<  domain_size + 1 << ")" << std::endl;
    intervals[stages - 1] = domain_size + 1;
}

pipeline_group::~pipeline_group() {
    delete[] b;
    mfreeschreier(&gp, &gens);
}

void pipeline_group::print_group_size() {
    double grpsize1;
    int grpsize2;
    mgrouporder(b, base_size, gp,&gens, &grpsize1, &grpsize2, domain_size);
    std::cout << grpsize1 << " * 10^" << (grpsize2) << std::endl; // 10 ^ grpsize2 actually
    //deleteunmarked(&gens);
    std::cout << "Generators: " << mschreier_gens(gens) << std::endl;
}
