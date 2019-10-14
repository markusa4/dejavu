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
#include "refinement_bucket.h"
#include <pthread.h>
#include <chrono>

extern long mmultcount;
extern long mfiltercount;

/*void auto_blaster::find_automorphism_prob_bucket(sgraph* g, bool compare, invariant* canon_I, bijection* canon_leaf,
                                          bijection* automorphism, std::default_random_engine* re, int *restarts, bool *done, int selector_seed, refinement_bucket* R) {
    bool backtrack = false;
    std::list<std::pair<int, int>> changes;
    selector S;
    coloring_bucket c;
    invariant I;
    std::list<int> init_color_class;
    *restarts = 0;
    ir_operation last_op = OP_R;
    I = start_I;
    c.copy(&start_cb);
    int level = 2;

    if(compare)
        I.set_compare_invariant(canon_I);
    int startlevel = I.current_level();

    while (true) {
        if (backtrack) {
            //std::cout << "backtrack" << std::endl;
            // initialize a search state
            //I.print();
            //assert(false);
            if(*done) {
                return;
            }
            *restarts += 1;
            c.copy(&start_cb);
            while(I.current_level() != startlevel)
                I.pop_level();
            // invariant, hopefully becomes complete in leafs such that automorphisms can be found
            init_color_class.clear();
            last_op = OP_R;
            backtrack = false;
            level = 2;
        }

        std::pair<int, int> s;
        if (!backtrack) {
            s = S.select_color_bucket(g, &c, selector_seed, level);
            if (s.first == -1) {
                if (compare) {
                    // we can derive an automorphism!
                    bijection leaf;
                    leaf.read_from_coloring_bucket(&c);
                    //I.print();
                    //c.print();
                    *automorphism = leaf;
                    automorphism->inverse();
                    automorphism->compose(*canon_leaf);
                    if(!g->certify_automorphism(*automorphism)) {
                        std::cout << "Restart (leaf)." << *restarts << std::endl;
                        backtrack = true;
                        continue;
                    }
                    std::cout << "Found automorphism." << *restarts << std::endl;
                    return;
                } else {
                    I.push_level();//I.print();

                    canon_leaf->read_from_coloring_bucket(&c);
                    *canon_I = I;
                    std::cout << "Finished first" << std::endl;
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
            bool comp = R->refine_coloring(g, &c, &I, level);
            last_op = OP_R;
            if (compare) {
                // compare invariant
                if(!comp || !I.compare_sizes()) {
                    backtrack = true;
                    continue;
                }
                //if (I.level_is_eq(canon_I, I.current_level())) {
                continue;
            }
        } else if (last_op == OP_R) {
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                             INDIVIDUALIZATION
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // collect all elements of color s //
            level += 1;
            int rpos = s.first + ((*re)() % (s.second));
            //std::cout << s.first << ", " << s.second << std::endl;
            assert(rpos < c.lab_sz);
            // cell size?
            int v = c.lab[rpos];
            // individualize random vertex of class
            R->individualize_vertex(g, &c, s.first, v, level);
            last_op = OP_I;
            if (!compare) { // base points
                automorphism->map.push_back(v);
            }
        }
    }
}*/

void auto_blaster::find_automorphism_prob(sgraph* g, bool compare, invariant* canon_I, bijection* canon_leaf,
        bijection* automorphism, std::default_random_engine* re, int *restarts, bool *done, int selector_seed, auto_workspace* w) {
    bool backtrack = false;
    std::list<std::pair<int, int>> changes;

    refinement* R = &w->R;
    selector* S = &w->S;
    coloring* c = &w->c;
    invariant* I = &w->I;
    work_set* first_level_fail = &w->first_level_fail;
    coloring* start_c = w->start_c;
    invariant start_I;

    S->empty_cache();

    std::list<int> init_color_class;
    *restarts = 0;
    int level = 1;
    ir_operation last_op = OP_R;
    *I = start_I;
    c->copy(start_c);
    if(compare)
        I->set_compare_invariant(canon_I);

    int base;


    while (true) {
        if(*done) return;
        if (backtrack) {
            S->empty_cache();
            // initialize a search state
            *restarts += 1;
            //c.rewrite_ptn(&start_c);
            c->copy(start_c);
            //while(I->current_level() != startlevel)
            //    I->pop_level();
            *I = start_I;
            if(compare)
                I->set_compare_invariant(canon_I);

            // invariant, hopefully becomes complete in leafs such that automorphisms can be found
            init_color_class.clear();
            last_op = OP_R;
            if(level == 2) first_level_fail->set(base);
            //std::cout << "level " << level << "base " << base << std::endl;
            backtrack = false;
            level = 1;
        }

        int s;
        if (!backtrack) {
            s = S->select_color(g, c, selector_seed);
            if (s == -1) {
                if (compare) {
                    //assert(I.level_is_eq(canon_I, I.current_level()));
                    //I->push_level();
                    //I->write_top_and_compare(INT32_MAX);
                    // we can derive an automorphism!
                    bijection leaf;
                    leaf.read_from_coloring(c);
                    *automorphism = leaf;
                    automorphism->inverse();
                    automorphism->compose(canon_leaf);
                    //std::cout << "Found automorphism." << *restarts << std::endl;
                    assert(g->certify_automorphism(*automorphism));
                    return;
                } else {
                    //I->push_level();
                    canon_leaf->read_from_coloring(c);
                    *canon_I = *I;
                    return;
                }
            }
        }

        if (last_op == OP_I) { // add new operations to trail...
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                REFINEMENT
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            changes.clear();
            bool comp = I->write_top_and_compare(INT32_MAX);
            //assert(init_color_class.size() == 2);
            comp = comp && R->refine_coloring(g, c, &changes, I, &init_color_class, false);
            comp = comp && I->write_top_and_compare(INT32_MAX);
            comp = comp && I->write_top_and_compare(INT32_MIN);
            //R.complete_colorclass_invariant(g, &c, &I);
            last_op = OP_R;
            if (compare) {
                // compare invariant
                if(!comp) { // || !I->compare_sizes()
                    backtrack = true;
                    continue;
                }
                continue;
            }
        } else if (last_op == OP_R) {
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                             INDIVIDUALIZATION
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // collect all elements of color s
            init_color_class.clear();
            int rpos = s + ((*re)() % (c->ptn[s] + 1));
            int v = c->lab[rpos];
            //assert(rpos == c->vertex_to_lab[v]);
            // first level fail prevention
            if(level == 1) {
                while(first_level_fail->get(v)) {
                    rpos = s + ((*re)() % (c->ptn[s] + 1));
                    v = c->lab[rpos];
                }
            }

            // individualize random vertex of class
            int newpos = R->individualize_vertex(g, c, v);
            last_op = OP_I;
            assert(init_color_class.empty());
            init_color_class.push_back(newpos);
            assert(c->vertex_to_col[v] > 0);
            //init_color_class.push_back(c.vertex_to_col[c.lab[labpos - 1]]);
            if (!compare) { // base points
                automorphism->map.push_back(v);
            }
            base = v;
            level += 1;

            bool comp = I->write_top_and_compare(INT32_MIN);
            comp && I->write_top_and_compare(INT32_MIN);
            if(!comp) {
                backtrack = true;
                continue;
            }
        }
    }
}

void auto_blaster::find_automorphism_bt(sgraph* g, bool compare, invariant* canon_I, bijection* canon_leaf,
       bijection* automorphism, std::default_random_engine* re, int *restarts, bool *done, int selector_seed) {
    /*bool backtrack;
    int backtrack_to_level = -1;
    std::list<std::pair<int, int>> changes;
    refinement R;
    selector S;
    coloring c;
    invariant I;
    std::list<int> init_color_class;
    *restarts -= 1;

    trail T(g->v_size);
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
            init_color_class.push_back(rpos);
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
    T.free_path();*/
}


/*void auto_blaster::sample(sgraph* g, bool master, bool* done) {
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
        //std::cout << "Launching workers..." << std::endl;
        for(int i = 0; i < config.CONFIG_THREADS_NO_PIPELINE - 1; i++)
            work_threads.emplace_back(std::thread(&auto_blaster::sample, this, g, false, done));
    }
    int trash_int = 0;
    auto_workspace W;
    W.first_level_fail.initialize(g->v_size);
    find_automorphism_prob(g, false, &canon_I, &canon_leaf, &base_points, &re, &trash_int, &trash_bool, selector_seed, &W);
    //std::cout << "Found canonical leaf." << std::endl;

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
        sequential_group G(g->v_size, &base_points);
        // run algorithm
        while (abort_counter <= 4) {
            sampled_paths += 1;
            // sample myself
            bijection automorphism;
            bool added;
            find_automorphism_prob(g, true, &canon_I, &canon_leaf, &automorphism, &re, &restarts, &trash_bool, selector_seed, &W);
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
            find_automorphism_prob(g, true, &canon_I, &canon_leaf, &automorphism, &re, &restarts, done, selector_seed, &W);
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
}*/

void auto_blaster::sample_pipelined(sgraph* g_, bool master, bool* done, pipeline_group* G, coloring* start_c, bijection* canon_leaf, invariant* canon_I) {
    // find comparison leaf
    sgraph* g = g_;

    if(config.CONFIG_THREADS_COPYG && !master) {
        g = new sgraph;
        g->copy_graph(g_);
    }

    std::vector<std::thread> work_threads;
    bijection base_points;
    bool trash_bool = false;
    int trash_int = 0;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine re = std::default_random_engine(seed);
    int selector_seed = re() % INT32_MAX;

    std::list<std::pair<int, int>> changes;

    invariant start_I;
   // start_I.push_level();

    auto_workspace W;
    W.first_level_fail.initialize(g->v_size);

    if(master) {
        canon_I    = new invariant;
        canon_leaf = new bijection;
        start_c = new coloring;
        std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
        g->initialize_coloring(start_c);
        std::list<int> init_color_class;
        W.start_c = start_c;
        W.R.refine_coloring(g, start_c, &changes, &start_I, &init_color_class, false);
        double cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "Color ref: " << cref / 1000000.0 << "ms" << std::endl;
        find_automorphism_prob(g, false, canon_I, canon_leaf, &base_points, &re, &trash_int, &trash_bool, selector_seed, &W);

        G = new pipeline_group();
        //std::cout << "Launching refinement worker..." << std::endl;
        for(int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++)
            work_threads.emplace_back(std::thread(&auto_blaster::sample_pipelined, auto_blaster(), g, false, done, G, start_c, canon_leaf, canon_I));
        std::cout << "Refinement workers (" << config.CONFIG_THREADS_REFINEMENT_WORKERS << ")" << std::endl;
    }
    //std::cout << "Found canonical leaf." << std::endl;

    int sampled_paths = 0;
    int restarts = 0;

    // sample for automorphisms
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                      MASTER THREAD
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(master) {
        // initialize automorphism group
        G->initialize(g->v_size, &base_points, config.CONFIG_THREADS_PIPELINE_DEPTH);
        G->launch_pipeline_threads(done);
        G->pipeline_stage(0, done);
        // run algorithm
        //std::cout << "Sampled leaves (local thread): \t" << sampled_paths << std::endl;
        //std::cout << "Sampled leaves (all threads): \t"  << sampled_paths_all + sampled_paths << std::endl;
        //std::cout << "Restart probing (local thread):\t" << restarts << std::endl;
        //std::cout << "Schreier Filtercount: \t" << mfiltercount << std::endl;
        //std::cout << "Schreier Multcount: \t"   << mmultcount << std::endl;
        std::cout << "Base size:  " << G->base_size << std::endl;
        std::cout << "Group size: ";
        G->print_group_size();
        std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
        G->join_threads();
        while(!work_threads.empty()) {
            work_threads[work_threads.size()-1].join();
            work_threads.pop_back();
        }
        double cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "Join: " << cref / 1000000.0 << "ms" << std::endl;
        delete G;
    } else {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                      SLAVE THREAD
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int idle_ms = 0;
        W.start_c = start_c;
       // std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();

       // make my own canonical leaf...
        canon_I    = new invariant;
        canon_leaf = new bijection;
        find_automorphism_prob(g, false, canon_I, canon_leaf, &base_points, &re, &trash_int, &trash_bool, selector_seed, &W);

        //double inner_t = 0;
        while(!(*done)) {
            bijection automorphism;
            //std::chrono::high_resolution_clock::time_point inner_timer = std::chrono::high_resolution_clock::now();
            if(config.CONFIG_IR_BACKTRACK) {
                //find_automorphism_bt(g, true, canon_I, canon_leaf, &automorphism, &re, &restarts, done, selector_seed);
            } else {
                find_automorphism_prob(g, true, canon_I, canon_leaf, &automorphism, &re, &restarts, done, selector_seed, &W);
            }
            //inner_t += (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - inner_timer).count());
            G->add_permutation(&automorphism, &idle_ms, done);
            sampled_paths += 1;
        }
        //double cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
        //std::cout << "Refinement speed: " << sampled_paths / (cref / 1000000.0) << "l/s" << std::endl;
        //std::cout << "Refinement speed (raw): " << sampled_paths / (inner_t / 1000000.0) << "l/s" << std::endl;
        //std::cout << "Dif: " << sampled_paths / (inner_t / 1000000.0)  - sampled_paths / (cref / 1000000.0) << "d" << std::endl;
        //std::cout << "Refinement worker idle: " << idle_ms << "ms" << std::endl;
        return;
    }
}


/*void auto_blaster::sample_pipelined_bucket(sgraph* g, bool master, bool* done, pipeline_group* G) {
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

    refinement_bucket R;
    std::list<std::pair<int, int>> changes;
    if(master) {
        start_I.push_level();

        std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
        g->initialize_coloring_bucket(&start_cb);
        std::list<int> init_color_class;
        R.initialize_active(g);
        R.refine_coloring(g, &start_cb, &start_I, 1);
       // c.print();
        double cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "Color ref: " << cref / 1000000.0 << "ms" << std::endl;

        G = new pipeline_group();
        for(int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++)
            work_threads.emplace_back(std::thread(&auto_blaster::sample_pipelined_bucket, this, g, false, done, G));
        std::cout << "Refinement workers (" << config.CONFIG_THREADS_REFINEMENT_WORKERS << ")" << std::endl;
    }

    int restarts = 0;

    // sample for automorphisms
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                      MASTER THREAD
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(master) {
        // initialize automorphism group
        find_automorphism_prob_bucket(g, false, &canon_I, &canon_leaf, &base_points, &re, &trash_int, &trash_bool, selector_seed, &R);
        //sleep(10);
        G->initialize(g->v_size, &base_points, config.CONFIG_THREADS_PIPELINE_DEPTH);
        G->launch_pipeline_threads(done);
        G->pipeline_stage(0, done);
        // run algorithm
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
        find_automorphism_prob_bucket(g, false, &canon_I, &canon_leaf, &base_points, &re, &trash_int, &trash_bool, selector_seed, &R);
        while(!(*done)) {
            bijection automorphism;
            if(config.CONFIG_IR_BACKTRACK) {
                //find_automorphism_bt(g, true, &canon_I, &canon_leaf, &automorphism, &re, &restarts, done, selector_seed);
            } else {
                find_automorphism_prob_bucket(g, true, &canon_I, &canon_leaf, &automorphism, &re, &restarts, done, selector_seed, &R);
            }
            G->add_permutation(&automorphism, &idle_ms, done);
        }
        //std::cout << "Refinement worker idle: " << idle_ms << "ms" << std::endl;
        return;
    }
}*/