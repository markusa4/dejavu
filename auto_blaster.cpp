//
// Created by markus on 23.09.19.
//
//#define NDEBUG

#include "auto_blaster.h"
#include <stack>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <random>
#include <chrono>
#include "ir_tools.h"
#include "trail.h"
#include "refinement.h"
#include "selector.h"
#include "invariant.h"
#include "configuration.h"
#include "group.h"
#include "concurrentqueue.h"
#include "invariant_acc.h"
#include <pthread.h>

void auto_blaster::find_automorphism_prob(sgraph* g, bool compare, invariant* canon_I, bijection* canon_leaf, bijection* automorphism, std::default_random_engine* re, int *restarts, bool *done) {
    bool backtrack = true;
    std::set<std::pair<int, int>> changes;
    refinement R;
    selector S;
    coloring c;
    invariant I;
    *restarts -= 1;
    ir_operation last_op;

    while (true) {
        if (backtrack) {
            // initialize a search state
            if(*done) {
                return;
            }
            *restarts += 1;
            c = start_c;
            I = start_I; // invariant, hopefully becomes complete in leafs such that automorphisms can be found
            last_op = OP_R;
            backtrack = false;
        }

        int s;
        if (!backtrack) {
            s = S.select_color(g, &c);
            if (s == -1) {
                if (compare) {
                    assert(I.level_is_eq(canon_I, I.current_level()));
                    I.push_level();
                    //R.complete_colorclass_invariant(g, &c, &I);
                    if (!I.level_is_eq(canon_I, I.current_level())) {
                        assert(false);
                        I.pop_level();
                        backtrack = true;
                        continue;
                    }
                    // we can derive an automorphism!
                    bijection leaf;
                    leaf.read_from_coloring(&c);
                    *automorphism = leaf;
                    automorphism->inverse();
                    automorphism->compose(*canon_leaf);
                    assert(g->certify_automorphism(*automorphism));
                    return;
                } else {
                    I.push_level();
                    //R.complete_colorclass_invariant(g, &c, &I);
                    canon_leaf->read_from_coloring(&c);
                    *canon_I = I;
                    return;
                }
            }
        }

        if (last_op == OP_I) { // add new operations to trail...
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                REFINEMENT
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            changes.clear();
            I.push_level();
            R.refine_coloring(g, &c, &changes, &I);
            //R.complete_colorclass_invariant(g, &c, &I);
            last_op = OP_R;
            if (compare) {
                // compare invariant
                if (I.level_is_eq(canon_I, I.current_level())) {
                    continue;
                } else {
                    backtrack = true;
                    continue;
                }
            }
        } else if (last_op == OP_R) {
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                             INDIVIDUALIZATION
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // collect all elements of color s
            int rpos = s + ((*re)() % (c.ptn[s] + 1));
            int v = c.lab[rpos];
            // individualize random vertex of class
            R.individualize_vertex(g, &c, v);
            last_op = OP_I;
        }
    }
}

void auto_blaster::sample(sgraph* g, bool master, bool* done) {
    // find comparison leaf
    //std::cout << &g << ", " << this << std::endl;
    invariant canon_I;
    std::vector<std::thread> work_threads;
    bijection canon_leaf;
    bijection trash;
    bool trash_bool = false;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine re = std::default_random_engine(seed);

    refinement R;
    std::set<std::pair<int, int>> changes;
    if(master) {
        start_I.push_level();
        g->initialize_coloring(&start_c);
        R.refine_coloring(g, &start_c, &changes, &start_I);
        std::cout << "Launching workers..." << std::endl;
        for(int i = 0; i < CONFIG_THREADS - 1; i++)
            work_threads.emplace_back(std::thread(&auto_blaster::sample, this, g, false, done));
    }
    int trash_int = 0;
    find_automorphism_prob(g, false, &canon_I, &canon_leaf, &trash, &re, &trash_int, &trash_bool);
    std::cout << "Found canonical leaf." << std::endl;

    int abort_counter = 0;
    int sampled_paths = 0;
    int sampled_paths_all = 0;
    int restarts = 0;

    // sample for automorphisms
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                      MASTER THREAD
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(master) {
        // initialize automorphism group
        group G(g->v.size());
        // run algorithm
        while (abort_counter <= 4) {
            sampled_paths += 1;
            // sample myself
            bijection automorphism;
            find_automorphism_prob(g, true, &canon_I, &canon_leaf, &automorphism, &re, &restarts, &trash_bool);
            bool added = G.add_permutation(&automorphism);
            if (added) {
                abort_counter = 0;
            } else {
                abort_counter += 1;
            }

            // add samples of other threads
            bool d;
            do {
                d = Q.try_dequeue(automorphism);
                if (d) {
                    sampled_paths_all += 1;
                    added = G.add_permutation(&automorphism);
                    if (added) {
                        abort_counter = 0;
                    } else {
                        abort_counter += 1;
                    }
                }
            } while(d);
        }
        *done = true;
        std::cout << "Sampled leaves (local thread): \t" << sampled_paths << std::endl;
        std::cout << "Sampled leaves (all threads): \t"  << sampled_paths_all + sampled_paths << std::endl;
        std::cout << "Restart probing (local thread):\t" << restarts << std::endl;
        std::cout << "Group size: ";
        G.print_group_size();
        while(!work_threads.empty()) {
            work_threads[work_threads.size()-1].join();
            work_threads.pop_back();
        }
    } else {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                      SLAVE THREAD
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        while(!(*done)) {
            bijection automorphism;
            find_automorphism_prob(g, true, &canon_I, &canon_leaf, &automorphism, &re, &restarts, done);
            Q.enqueue(automorphism);
        }
        return;
    }
}