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
#include "sequential_group.h"
#include "concurrentqueue.h"
#include "invariant_acc.h"
#include "pipeline_group.h"
#include <pthread.h>

extern long mmultcount;
extern long mfiltercount;

// ToDo: recreate backtracking version
void auto_blaster::find_automorphism_prob(sgraph* g, bool compare, invariant* canon_I, bijection* canon_leaf,
        bijection* automorphism, std::default_random_engine* re, int *restarts, bool *done, int selector_seed) {
    bool backtrack = false;
    std::list<std::pair<int, int>> changes;
    refinement R;
    selector S;
    coloring c;
    invariant I;
    std::list<int> init_color_class;
    *restarts = 0;
    ir_operation last_op = OP_R;
    I = start_I;
    c = start_c;
    if(compare)
        I.set_compare_invariant(canon_I);
    int startlevel = I.current_level();
    

    while (true) {
        if (backtrack) {
            // initialize a search state
            if(*done) {
                return;
            }
            *restarts += 1;
            //c.rewrite_ptn(&start_c);
            c = start_c;

            while(I.current_level() != startlevel)
                I.pop_level();
            // invariant, hopefully becomes complete in leafs such that automorphisms can be found
            init_color_class.clear();
            last_op = OP_R;
            backtrack = false;
        }

        int s;
        if (!backtrack) {
            s = S.select_color(g, &c, selector_seed);
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
                    //std::cout << "Found automorphism." << *restarts << std::endl;
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
            //assert(init_color_class.size() == 2);
            bool comp = R.refine_coloring(g, &c, &changes, &I, &init_color_class, false);
            //R.complete_colorclass_invariant(g, &c, &I);
            last_op = OP_R;
            if (compare) {
                // compare invariant
                if(!comp || !I.compare_sizes()) {
                    backtrack = true;
                    continue;
                }
                //if (I.level_is_eq(canon_I, I.current_level())) {
                    continue;
                //} else {
                //    backtrack = true;
               //     continue;
               // }
            }
        } else if (last_op == OP_R) {
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                             INDIVIDUALIZATION
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // collect all elements of color s
            init_color_class.clear();
            int rpos = s + ((*re)() % (c.ptn[s] + 1));
            int v = c.lab[rpos];
            // individualize random vertex of class
            R.individualize_vertex(g, &c, v);
            last_op = OP_I;
            assert(init_color_class.empty());
            int labpos = c.vertex_to_lab[v];
            assert(labpos == c.vertex_to_col[v]);
            init_color_class.push_back(labpos);
            assert(c.vertex_to_col[v] > 0);
            //init_color_class.push_back(c.vertex_to_col[c.lab[labpos - 1]]);
            if (!compare) {
                automorphism->map.push_back(v);
            }
        }
    }
}

void auto_blaster::find_automorphism_bt(sgraph* g, bool compare, invariant* canon_I, bijection* canon_leaf,
        bijection* automorphism, std::default_random_engine* re, int *restarts, bool *done, int selector_seed) {
    bool backtrack;
    int backtrack_to_level = -1;
    std::list<std::pair<int, int>> changes;
    refinement R;
    selector S;
    coloring c;
    invariant I;
    std::list<int> init_color_class;
    *restarts -= 1;

    trail T(g->v.size());
    c = start_c;
    I = start_I; // invariant, hopefully becomes complete in leafs such that automorphisms can be found
    if(compare)
        I.set_compare_invariant(canon_I);
    init_color_class.clear();
    backtrack = false;
    T.push_op_r(&changes);

    while (true) {
        if(T.last_op() == OP_END) {
            assert(false);
        }

        if (backtrack) {
            //std::cout << I.current_level() << std::endl;
            // initialize a search state
            if(*done) {
                break;
            }
            if(*restarts >= 0) {
                //std::cout << "Restart." << *restarts << std::endl;
            }
        }

        int s;
        if (!backtrack) {
            s = S.select_color(g, &c, selector_seed);
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
                    //std::cout << "Found automorphism." << *restarts << std::endl;
                    assert(g->certify_automorphism(*automorphism));
                    break;
                } else {
                    I.push_level();
                    //R.complete_colorclass_invariant(g, &c, &I);
                    canon_leaf->read_from_coloring(&c);
                    *canon_I = I;
                    break;
                }
            }
        }

        if (!backtrack && T.last_op() == OP_I) { // add new operations to trail...
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                REFINEMENT
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            changes.clear();
            I.push_level();
            //assert(init_color_class.size() == 2);
            bool comp = R.refine_coloring(g, &c, &changes, &I, &init_color_class, true);
            //assert(R.assert_is_equitable(g, &c));
            T.push_op_r(&changes);
            //R.complete_colorclass_invariant(g, &c, &I);
            //last_op = OP_R;
            if (compare) {
                // compare invariant
                if(!comp || !I.compare_sizes()) {
                    if(config.CONFIG_IR_BACKTRACK_RANDOM) {
                        backtrack_to_level = ((*re)() % I.current_level());
                    }
                    backtrack = true;
                    continue;
                }
                continue;
            }
        } else if (!backtrack && T.last_op() == OP_R) {
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                             INDIVIDUALIZATION
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // collect all elements of color s
            init_color_class.clear();
            std::list<int> color_s;
            int rpos = s + ((*re)() % (c.ptn[s] + 1));
            int v = c.lab[rpos];
            int i = s;
            if(!config.CONFIG_IR_BACKTRACK_RANDOM) {
                while ((i == s) || (i == 0) || c.ptn[i - 1] != 0) {
                    if (i != rpos) {
                        color_s.push_front(c.lab[i]);
                    }
                    i += 1;
                }
                assert(color_s.size() > 0);
            }
            // individualize random vertex of class
            R.individualize_vertex(g, &c, v);
            // last_op = OP_I;
            assert(init_color_class.empty());
            int labpos = c.vertex_to_lab[v];
            assert(labpos == c.vertex_to_col[v]);
            init_color_class.push_back(labpos);
            assert(c.vertex_to_col[v] > 0);
            //init_color_class.push_back(c.vertex_to_col[c.lab[labpos - 1]]);
            if (!compare) {
                automorphism->map.push_back(v);
            }
            T.push_op_i(&color_s, v);
        } else if(T.last_op() == OP_I && backtrack) { // backtrack trail, undo operations...
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                       BACKTRACK INDIVIDUALIZATION
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            init_color_class.clear();
            // backtrack until we can do a new operation
            int v = T.top_op_i_v();
            R.undo_individualize_vertex(g, &c, v);
            T.pop_op_i_v();
            //assert(R.assert_is_equitable(g, &c));

            // undo individualization
            if((T.top_op_i_class().empty() && !config.CONFIG_IR_BACKTRACK_RANDOM) || backtrack_to_level > 0) {
                backtrack_to_level -= 1;
                // we tested the entire color class, need to backtrack further
                T.pop_op_i_class();
                continue;
            } else {
                // there is another vertex we have to try, so we are done backtracking
                if(!config.CONFIG_IR_BACKTRACK_RANDOM) {
                    T.shuffle_top_i_class(re);
                    v = T.top_op_i_class().front();
                    T.top_op_i_class().pop_front();
                } else {
                    s = c.vertex_to_col[v];
                    int rpos = s + ((*re)() % (c.ptn[s] + 1));
                    v = c.lab[rpos];
                }
                T.push_op_i_v(v);
                R.individualize_vertex(g, &c, v);
                int labpos = c.vertex_to_lab[v];
                assert(labpos == c.vertex_to_col[v]);
                init_color_class.push_back(labpos);
                backtrack = false;
            }
        }  else if(T.last_op() == OP_R && backtrack) {
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                             UNDO REFINEMENT
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // we are backtracking, so we have to undo refinements
            I.pop_level();
           // assert(R.assert_is_equitable(g, &c));
            R.undo_refine_color_class(g, &c, &T.top_op_r());
            //assert(T.top_op_r().empty()?R.assert_is_equitable(g, &c):!R.assert_is_equitable(g, &c));
            T.pop_op_r();
        }
    }
    T.free_path();
}


void auto_blaster::sample(sgraph* g, bool master, bool* done) {
    // find comparison leaf
    invariant canon_I;
    std::vector<std::thread> work_threads;
    bijection canon_leaf;
    bijection base_points;
    bool trash_bool = false;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine re = std::default_random_engine(seed);
    int selector_seed = re() % INT32_MAX;

    refinement R;
    std::list<std::pair<int, int>> changes;
    if(master) {
        start_I.push_level();
        g->initialize_coloring(&start_c);
        std::list<int> init_color_class;
        R.refine_coloring(g, &start_c, &changes, &start_I, &init_color_class, false);
        std::cout << "Launching workers..." << std::endl;
        for(int i = 0; i < config.CONFIG_THREADS_NO_PIPELINE - 1; i++)
            work_threads.emplace_back(std::thread(&auto_blaster::sample, this, g, false, done));
    }
    int trash_int = 0;
    find_automorphism_prob(g, false, &canon_I, &canon_leaf, &base_points, &re, &trash_int, &trash_bool, selector_seed);
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
        sequential_group G(g->v.size(), &base_points);
        // run algorithm
        while (abort_counter <= 4) {
            sampled_paths += 1;
            // sample myself
            bijection automorphism;
            bool added;
            find_automorphism_prob(g, true, &canon_I, &canon_leaf, &automorphism, &re, &restarts, &trash_bool, selector_seed);
            added = G.add_permutation(&automorphism);
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
            } while(d && abort_counter <= 3);
        }
        *done = true;
        std::cout << "Sampled leaves (local thread): \t" << sampled_paths << std::endl;
        std::cout << "Sampled leaves (all threads): \t"  << sampled_paths_all + sampled_paths << std::endl;
        std::cout << "Restart probing (local thread):\t" << restarts << std::endl;
        //std::cout << "Schreier Filtercount: \t" << mfiltercount << std::endl;
        //std::cout << "Schreier Multcount: \t"   << mmultcount << std::endl;
        std::cout << "Base size:  " << G.base_size << std::endl;
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
            find_automorphism_prob(g, true, &canon_I, &canon_leaf, &automorphism, &re, &restarts, done, selector_seed);
            Q.enqueue(automorphism);
            if(Q.size_approx() > 20) {
                if(Q.size_approx() < 100) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(50));
                } else {
                    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
                }
            }
        }
        return;
    }
}

void auto_blaster::sample_pipelined(sgraph* g, bool master, bool* done, pipeline_group* G) {
    // find comparison leaf
    invariant canon_I;
    std::vector<std::thread> work_threads;
    bijection canon_leaf;
    bijection base_points;
    bool trash_bool = false;
    int trash_int = 0;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine re = std::default_random_engine(seed);
    int selector_seed = re() % INT32_MAX;

    refinement R;
    std::list<std::pair<int, int>> changes;
    if(master) {
        start_I.push_level();
        g->initialize_coloring(&start_c);
        std::list<int> init_color_class;
        R.refine_coloring(g, &start_c, &changes, &start_I, &init_color_class, false);
        G = new pipeline_group();
        std::cout << "Launching refinement worker..." << std::endl;
        for(int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++)
            work_threads.emplace_back(std::thread(&auto_blaster::sample_pipelined, this, g, false, done, G));
    }
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
        find_automorphism_prob(g, false, &canon_I, &canon_leaf, &base_points, &re, &trash_int, &trash_bool, selector_seed);
        G->initialize(g->v.size(), &base_points, config.CONFIG_THREADS_PIPELINE_DEPTH);
        G->launch_pipeline_threads(done);
        G->pipeline_stage(0, done);
        // run algorithm
        std::cout << "Sampled leaves (local thread): \t" << sampled_paths << std::endl;
        std::cout << "Sampled leaves (all threads): \t"  << sampled_paths_all + sampled_paths << std::endl;
        std::cout << "Restart probing (local thread):\t" << restarts << std::endl;
        //std::cout << "Schreier Filtercount: \t" << mfiltercount << std::endl;
        //std::cout << "Schreier Multcount: \t"   << mmultcount << std::endl;
        std::cout << "Base size:  " << G->base_size << std::endl;
        std::cout << "Group size: ";
        G->print_group_size();
        G->join_threads();
        while(!work_threads.empty()) {
            work_threads[work_threads.size()-1].join();
            work_threads.pop_back();
        }
        delete G;
    } else {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                      SLAVE THREAD
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int idle_ms = 0;
        find_automorphism_prob(g, false, &canon_I, &canon_leaf, &base_points, &re, &trash_int, &trash_bool, selector_seed);
        while(!(*done)) {
            bijection automorphism;
            if(config.CONFIG_IR_BACKTRACK) {
                find_automorphism_bt(g, true, &canon_I, &canon_leaf, &automorphism, &re, &restarts, done, selector_seed);
            } else {
                find_automorphism_prob(g, true, &canon_I, &canon_leaf, &automorphism, &re, &restarts, done, selector_seed);
            }
            G->add_permutation(&automorphism, &idle_ms, done);
        }
        std::cout << "Refinement worker idle: " << idle_ms << "ms" << std::endl;
        return;
    }
}