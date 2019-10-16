//
// Created by markus on 01.10.19.
//

#include <iostream>
#include "pipeline_group.h"
#include "configuration.h"

void pipeline_group::launch_pipeline_threads(bool* done, bool* done_fast) {
    for(int i = 1; i < stages; ++i) {
        work_threads.emplace_back(std::thread(&pipeline_group::pipeline_stage, this, i, done, done_fast));
    }
    //std::cout << "Pipeline workers (" << stages << ")" << std::endl;
}

void pipeline_group::join_threads() {
    while(!work_threads.empty()) {
        work_threads[work_threads.size()-1].join();
        work_threads.pop_back();
    }
}

void pipeline_group::pipeline_stage(int n, bool* done, bool* done_fast) {
   // int front_idle_ms = 0;
   // int back_idle_ms  = 0;
    int abort_counter = 0;
    int leafs_considered = 0;
    int random_abort_counter = 0;
    int bulk_deque_num = 8;
    std::pair<bool, bool>* res = new std::pair<bool, bool>[bulk_deque_num];

    std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();

    if(base_size == 0)
        *done = true;

    while(!(*done)) {
        // ToDo: load balancing
        // work on pipeline_results and track
        bool is_first_stage = (n == 0);
        bool is_last_stage = (n == stages - 1);
        bool is_only_stage = is_first_stage && is_last_stage;
        int my_interval = intervals[n];

        filterstate state;
        bijection p;
        random_element re;
        int *_p = nullptr;
        bool is_random = false;

        // context switch
        //std::this_thread::yield();
        if (is_first_stage) {
            // collect results and determine whether done
            int num = sift_results.try_dequeue_bulk(res, bulk_deque_num);
            for(int j = 0; j < num; ++j) {
                if(!res[j].first) {
                    leafs_considered += 1;
                    if(!res[j].second) {
                        abort_counter += 1;
                    } else {
                        abort_counter = 0;
                    }
                } else {
                    if(!res[j].second) {
                        random_abort_counter += 1;
                    } else {
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
            while(!d && random_abort_counter >= config.CONFIG_RAND_ABORT_RAND) {
                //std::this_thread::yield();
                //std::this_thread::sleep_for(std::chrono::milliseconds(1));
                //front_idle_ms += 1;
                d = automorphisms.try_dequeue(ctoken, p);
                /*if(front_idle_ms % 10000 == 0) {
                    std::cout << "Pipeline(" << n << ") front idle " << front_idle_ms << ", " << automorphisms.size_approx() << std::endl;
                }*/
            }

            if(leafs_considered == 0 && d) {
              //  double cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
               // std::cout << "First automorphism arrived: " << cref / 1000000.0 << "ms" << std::endl;
            }

            //std::cout << automorphisms.size_approx() << std::endl;

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
            if(!state.counts_towards_abort && !result && !state.ingroup && !(*done_fast)) {
                //std::cout << "useless element" << result << std::endl;
                *done_fast = true;
               // std::cout << "random element: " << result << std::endl;
            }
            sift_results.enqueue(std::pair<bool, bool>(state.ingroup || !state.counts_towards_abort, result));
        } else {
            while(pipeline_queues[n + 1].size_approx() > 50 && (!(*done))) {
                /*if(back_idle_ms % 10000 == 0) {
                    std::cout << "Pipeline(" << n << ") back idle " << back_idle_ms << ", "
                              << pipeline_queues[n + 1].size_approx() << std::endl;
                }*/
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                //std::this_thread::yield();
                //back_idle_ms += 1;
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

bool pipeline_group::add_permutation(bijection *p, int* idle_ms, bool* done) {
    //std::cout << "enqueued" << std::endl;
    static thread_local moodycamel::ProducerToken ptoken = moodycamel::ProducerToken(automorphisms);

    while(automorphisms.size_approx() > 100 && (!(*done))) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        *idle_ms += 1;
    }
    automorphisms.enqueue(ptoken, *p);
    return true;
}

pipeline_group::pipeline_group() {
}

void pipeline_group::initialize(int domain_size, bijection *base_points, int stages) {
    mschreier_fails(-1);
    added = 0;
    this->domain_size = domain_size;
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
   // std::cout << "Pipeline: ";
    int ded = 0;
    for(int i = 0; i < stages; ++i) {
        int stage_pos = (base_size / stages) * (i + 1) - (i + 1) * ((base_size / (stages* (i + 2))));

        // skip stage if too small
        if((!intervals.empty()) &&
            (stage_pos - intervals[intervals.size() - 1] < config.CONFIG_THREADS_PIPELINE_STAGE_MIN)) {
            ded += 1;
            continue;
        }

        pipeline_queues.emplace_back(moodycamel::ConcurrentQueue<filterstate>(20));
        intervals.push_back(stage_pos);
       // std::cout << "/" << intervals[intervals.size() - 1];
    }
    stages -= ded;
   // std::cout << "(" <<  domain_size + 1 << ")" << std::endl;
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
