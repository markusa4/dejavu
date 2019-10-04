//
// Created by markus on 01.10.19.
//

#include <iostream>
#include "pipeline_group.h"

void pipeline_group::launch_pipeline_threads(bool* done) {
    for(int i = 1; i < stages; ++i) {
        std::cout << "Launching pipeline worker (" << i << ")..." << std::endl;
        work_threads.emplace_back(std::thread(&pipeline_group::pipeline_stage, this, i, done));
    }
}

void pipeline_group::join_threads() {
    while(!work_threads.empty()) {
        work_threads[work_threads.size()-1].join();
        work_threads.pop_back();
    }
}

void pipeline_group::pipeline_stage(int n, bool* done) {
    int front_idle_ms = 0;
    int back_idle_ms  = 0;
    int abort_counter = 0;
    int random_abort_counter = 0;

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

        if (is_first_stage) {
            // collect results and determine whether done
            std::pair<bool, bool> res;
            bool d = sift_results.try_dequeue(res);
            while(d) {
                if(!res.first) {
                    if(!res.second) {
                        abort_counter += 1;
                    } else {
                        abort_counter = 0;
                    }
                } else {
                    if(!res.second) {
                        random_abort_counter += 1;
                    } else {
                        random_abort_counter = 0;
                    }
                }
                if(abort_counter >= 5) {
                    *done = true;
                    break;
                }
                //std::cout << abort_counter << "; " << random_abort_counter << std::endl;
                d = sift_results.try_dequeue(res);
            }
            if(*done) {break;}

            // not done, so choose next element for the pipeline
            d = automorphisms.try_dequeue(p);
            // automorphisms non-empty? take that element
            while(!d && random_abort_counter >= 5) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                front_idle_ms += 1;
                d = automorphisms.try_dequeue(p);
            }

            if (d) {
                //std::cout << "automorphism in pipeline" << std::endl;
                _p = new int[domain_size];
                for (int k = 0; k < p.map.size(); ++k) {
                    _p[k] = p.map[k];
                }
                is_random = false;
                state.ingroup = false;
                state.level = -1;
                random_abort_counter = 0;
            } else {
                // else generate random element and sift that
                bool gen = generate_random_element(gp, &gens, domain_size, &re);
                if (!gen) {continue;}
                _p = re.perm;
                is_random = true;
                state.ingroup = true;
                state.level = -1;
            }
        } else {
            // take element from intermediate queue
            bool d = pipeline_queues[n].try_dequeue(state);
            while (!d && (!(*done))) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                front_idle_ms += 1;
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
            sift_results.enqueue(std::pair<bool, bool>(state.ingroup, result));
        } else {
            while(pipeline_queues[n + 1].size_approx() > 50) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                back_idle_ms += 1;
            }
            pipeline_queues[n + 1].enqueue(state);
        }
    }
    std::cout << "Pipeline stage(" << n << ") idle: " << front_idle_ms << "ms / " << back_idle_ms << "ms" << std::endl;
}

bool pipeline_group::add_permutation(bijection *p, int* idle_ms) {
    //std::cout << "enqueued" << std::endl;
    while(automorphisms.size_approx() > 50) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        *idle_ms += 1;
    }
    automorphisms.enqueue(*p);
    return true;
}

pipeline_group::pipeline_group() {
}

void pipeline_group::initialize(int domain_size, bijection *base_points, int stages) {
    mschreier_fails(1);
    added = 0;
    this->domain_size = domain_size;
    this->base_size = 0;
    this->stages = stages;

    std::cout << "Creating new pipeline_group... " << std::endl;
    mnewgroup(&gp, &gens, domain_size);
    b = new int[domain_size];
    for(int i = 0; i < base_points->map.size(); ++i) {
        b[i] = base_points->map[i];
        base_size += 1;
    }

    mgetorbits(b, base_size, gp, &gens, domain_size);

    // create pipeline
    determine_stages();
}

void pipeline_group::determine_stages() {
    std::cout << "Pipeline: ";
    for(int i = 0; i < stages; ++i) {
        pipeline_queues.emplace_back(moodycamel::ConcurrentQueue<filterstate>());
        intervals.push_back((base_size / stages) * (i + 1) - ((base_size / (stages* 2))));
        std::cout << "/" << intervals[intervals.size() - 1];
    }
    std::cout << "(" <<  domain_size + 1 << ")" << std::endl;
    intervals[stages - 1] = domain_size + 1;
}

pipeline_group::~pipeline_group() {
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
