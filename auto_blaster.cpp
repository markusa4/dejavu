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
#include "diy_group.h"
#include <pthread.h>
#include <chrono>
#include <tuple>

extern long mmultcount;
extern long mfiltercount;

int intRand(const int & min, const int & max, int seed) {
    static thread_local std::mt19937 generator(seed);
    std::uniform_int_distribution<int> distribution(min,max);
    return distribution(generator);
}

bool auto_blaster::proceed_state(auto_workspace* w, sgraph* g, coloring* c, invariant* I, int v, change_tracker* changes) {
    if(changes != nullptr)
        changes->track(c->vertex_to_col[v]);
    int init_color_class = w->R.individualize_vertex(c, v);
    bool comp = I->write_top_and_compare(INT32_MIN);
    comp && I->write_top_and_compare(INT32_MIN);
    //assert(comp);
    //w->first_level += 1;
    //assert(init_color_class.size() == 2);
    comp = comp && I->write_top_and_compare(INT32_MAX);
    if(!comp) return comp;
    //assert(comp);
    comp = comp && w->R.refine_coloring(g, c, changes, I, init_color_class, changes != nullptr);
    comp = comp && I->write_top_and_compare(INT32_MAX);
    comp = comp && I->write_top_and_compare(INT32_MIN);
    return comp;
}

void auto_blaster::find_automorphism_from_bfs(auto_workspace *w, sgraph *g, bool compare, invariant *canon_I,
                                          bijection *canon_leaf, bijection *automorphism, int *restarts,
                                          shared_switches *switches, int selector_seed) {
    bool backtrack = false;
    bool* done = &switches->done;

    refinement *R = &w->R;
    selector *S = &w->S;
    coloring *c = &w->c;
    invariant *I = &w->I;

    S->empty_cache();

    int init_color_class;
    *restarts = 0;
    int reguess = 0;
    int level = w->first_level;
    int first_level_point = -1;
    int second_level_point = -1;

    int enqueue_fail_point_sz = 0;

    if (compare) {
        automorphism->map    = new int[g->v_size];
        automorphism->map_sz = 0;
    }
    ir_operation last_op;

    // Pick start from BFS level instead!

    int bfs_level    = w->BW->BW.current_level - 1;
    int bfs_level_sz = w->BW->BW.level_sizes[bfs_level];
    int rand_pos     = intRand(0, bfs_level_sz - 1, selector_seed);
    bfs_element* picked_elem = w->BW->BW.level_states[bfs_level][rand_pos];
    //std::cout << "picking level " << bfs_level << " element " << rand_pos << std::endl;
    *I = *picked_elem->I;
    c->copy_force(picked_elem->c);
    I->set_compare_invariant(canon_I);
    last_op = OP_R;
    level = bfs_level + 1;
    S->empty_cache();
    int base;

    while (true) {
        if(*done) return;
        if (backtrack) {
            bfs_level    = w->BW->BW.current_level - 1;
            bfs_level_sz = w->BW->BW.level_sizes[bfs_level];
            rand_pos     = intRand(0, bfs_level_sz - 1, selector_seed);
            picked_elem = w->BW->BW.level_states[bfs_level][rand_pos];
            *I = *picked_elem->I;
            c->copy_force(picked_elem->c);
            I->set_compare_invariant(canon_I);
            backtrack = false;
            last_op = OP_R;
            level = bfs_level + 1;
            S->empty_cache();

            *restarts += 1;
            //if(*restarts % 1000 == 0)
            //    std::cout << *restarts << ", " << reguess << ", " << level << std::endl;
        }

        int s;
        if (!backtrack) {
            s = S->select_color(g, c, selector_seed);
            if (s == -1) {
                if (compare) {
                    // we can derive an automorphism!
                    w->measure2 += 1;
                    bijection leaf;
                    leaf.read_from_coloring(c);
                    leaf.not_deletable();
                    *automorphism = leaf;
                    automorphism->inverse();
                    automorphism->compose(canon_leaf);//enqueue_fail_point_sz
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
            //changes.clear();
            bool comp = I->write_top_and_compare(INT32_MAX);
            //assert(init_color_class.size() == 2);
            comp = comp && R->refine_coloring(g, c, nullptr, I, init_color_class, false);
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
            init_color_class = -1;
            int rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
            int v = c->lab[rpos];
            // if(switches->done_shared_group && config.CONFIG_THREADS_COLLABORATE) // ToDo why is this problematic?
            //    v = (*w->shared_orbit)[v];

            //assert(rpos == c->vertex_to_lab[v]);

            if(level == w->first_level)
                first_level_point = v;

            if(level == w->first_level + 1)
                second_level_point = v;

            // individualize random vertex of class
            int newpos = R->individualize_vertex(c, v);
            last_op = OP_I;
//            assert(init_color_class >= 0);
            init_color_class = newpos;
            assert(c->vertex_to_col[v] > 0);
            //init_color_class.push_back(c.vertex_to_col[c.lab[labpos - 1]]);
            if (!compare) { // base points
                //automorphism->map.push_back(v);
                automorphism->map[automorphism->map_sz] = v;
                automorphism->map_sz += 1;
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


void auto_blaster::find_automorphism_prob(auto_workspace *w, sgraph *g, bool compare, invariant *canon_I,
                                          bijection *canon_leaf, bijection *automorphism, int *restarts,
                                          shared_switches *switches, int selector_seed) {
    bool backtrack = false;
    bool* done = &switches->done;
    //std::list<std::pair<int, int>> changes;

    refinement *R = &w->R;
    selector *S = &w->S;
    coloring *c = &w->c;
    invariant *I = &w->I;

    work_set *first_level_fail = &w->first_level_fail;
    work_set *first_level_succ = &w->first_level_succ;

    coloring *start_c = w->start_c;
    invariant *start_I = &w->start_I;

    S->empty_cache();

    int init_color_class;
    *restarts = 0;
    int reguess = 0;
    int level = w->first_level;
    int first_level_point = -1;
    int second_level_point = -1;

    int enqueue_fail_point_sz = 0;

    if (compare) {
        start_I->set_compare_invariant(canon_I);
    } else {
        start_I->create_vector();
        automorphism->map = new int[g->v_size];
        automorphism->map_sz = 0;
    }
    ir_operation last_op = OP_R;
    *I = *start_I;
    c->copy(start_c);
    int base;


    while (true) {
        if(*done) return;
        if (backtrack) {
            w->measure1 += 1;
            // check for new messages
            bool receive_advance = false;
            size_t num = 1;
            size_t max_repeat = 0;
            if(config.CONFIG_THREADS_COLLABORATE) {
                while(num > 0) {
                    num = (*w->communicator_pad)[w->communicator_id].try_dequeue_bulk(*w->ctok, w->dequeue_space, w->dequeue_space_sz);
                    for (int j = 0; j < num; ++j) {
                        if (abs(std::get<0>(w->dequeue_space[j])) == w->first_level) {
                            if (std::get<0>(w->dequeue_space[j]) > 0) {
                                if(!first_level_succ->get(std::get<1>(w->dequeue_space[j]))) {
                                    int vertex = std::get<1>(w->dequeue_space[j]);
                                    first_level_succ->set(vertex);
                                    w->first_level_succ_point = vertex;
                                    w->first_level_sz += 1;
                                }
                            } else {
                                if(!first_level_fail->get(std::get<1>(w->dequeue_space[j]))) {
                                    first_level_fail->set(std::get<1>(w->dequeue_space[j]));
                                    w->first_level_sz += 1;
                                }
                            }
                        } else if (abs(std::get<0>(w->dequeue_space[j])) == 0) {
                            if(w->first_level == 1) {
                                receive_advance = true;
                            }
                        }
                    }
                }
            }
            if(*done) return;
            S->empty_cache();
            // initialize a search state
            *restarts += 1;
            c->copy(start_c);
            *I = *start_I;

            // invariant, hopefully becomes complete in leafs such that automorphisms can be found
            if (level == w->first_level + 1) {
                if (!first_level_fail->get(base)) {
                    if(config.CONFIG_THREADS_COLLABORATE) {
                        if(enqueue_fail_point_sz < w->enqueue_space_sz) {
                            w->enqueue_space[enqueue_fail_point_sz] = std::pair<int, int>(-w->first_level, base);
                            enqueue_fail_point_sz += 1;
                        } else {
                            //for (int i = 0; i < w->communicator_pad->size(); ++i) {
                                //if (i != w->communicator_id) {
                                    //(*w->communicator_pad)[i].try_enqueue_bulk(*w->ptoks[i], w->enqueue_space, enqueue_fail_point_sz);
                                    w->G->first_level_points.try_enqueue_bulk(*w->ptok, w->enqueue_space, enqueue_fail_point_sz);
                                    enqueue_fail_point_sz = 0;
                                //}
                            //}
                        }
                    }

                    //std::cout << "failed" << first_level_point << std::endl;
                    first_level_fail->set(base);
                    w->first_level_sz += 1;
                }
            }

            backtrack = false;
            level = w->first_level;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                             INDIVIDUALIZATION (1)
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //S->empty_cache();
            int s = S->select_color(g, c, selector_seed);
            if((receive_advance || w->first_level_sz == c->ptn[s] + 1) && w->first_level < w->base_size) {
                /*assert(s != -1);
                // proceed first level
                assert(c->vertex_to_col[w->first_level_succ_point] == s);
                init_color_class = R->individualize_vertex(start_c, w->first_level_succ_point);
                bool comp = start_I->write_top_and_compare(INT32_MIN);
                comp && start_I->write_top_and_compare(INT32_MIN);
                assert(comp);
                w->first_level += 1;
                //assert(init_color_class.size() == 2);
                comp = comp && start_I->write_top_and_compare(INT32_MAX);
                assert(comp);
                comp = comp && R->refine_coloring(g, start_c, nullptr, start_I, init_color_class, false);
                comp = comp && start_I->write_top_and_compare(INT32_MAX);
                comp = comp && start_I->write_top_and_compare(INT32_MIN);*/
                //std::cout << "proceed receive" << std::endl;
                proceed_state(w, g, start_c, start_I, w->G->b[w->first_level - 1], nullptr);
                w->first_level += 1;

                first_level_fail->reset();
                first_level_succ->reset();
                w->first_level_succ_point = -1;
                w->first_level_sz = 0;

                level = w->first_level;
                c->copy_force(start_c);
                *I = *start_I;
                //S->empty_cache();
                s  = S->select_color(g, c, selector_seed);
                assert(s != -1);
                assert(start_I->cur_pos == I->cur_pos);
            } else if (w->first_level_sz == c->ptn[s] + 1) {
                // I added all automorphisms!
                std::cout << "Levels explored..." << std::endl;
                *done = true;
            }

            // collect all elements of color s
            int rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
            int v = c->lab[rpos];
            //assert(rpos == c->vertex_to_lab[v]);
            // first level fail prevention done_shared_group
            if(!switches->done_shared_group) {
                while (first_level_fail->get(v)) {
                    if (*done) return;
                    reguess += 1;
                    //rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
                    rpos = s + ((rpos - s + 1) % (c->ptn[s] + 1));
                    v = c->lab[rpos];
                }
            } else {
               // std::cout << "using shared orbit " << std::endl;
                v = (*w->shared_orbit)[v];
                while (first_level_fail->get(v)) {
                    if (*done) return;
                    reguess += 1;
                    //rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
                    rpos = s + ((rpos - s + 1) % (c->ptn[s] + 1));
                    v = (*w->shared_orbit)[c->lab[rpos]];
                }
            }

            first_level_point = v;
            // individualize random vertex of class
            int newpos = R->individualize_vertex(c, v);
            last_op = OP_I;
            assert(init_color_class >= 0);
            init_color_class = newpos;
            assert(c->vertex_to_col[v] > 0);
            //init_color_class.push_back(c.vertex_to_col[c.lab[labpos - 1]]);
            if (!compare) { // base points
                automorphism->map[automorphism->map_sz] = v;
                automorphism->map_sz += 1;
                //automorphism->map.push_back(v);
            }
            base = v;
            level += 1;

            bool comp = I->write_top_and_compare(INT32_MIN);
            comp && I->write_top_and_compare(INT32_MIN);
            assert(comp);
        }

        int s;
        if (!backtrack) {
            s = S->select_color(g, c, selector_seed);
            if (s == -1) {
                if (compare) {
                    // we can derive an automorphism!
                    //std::cout << *restarts << ", " << reguess << ", " << w->first_level << std::endl;
                    w->measure2 += 1;
                    bijection leaf;
                    leaf.read_from_coloring(c);
                    leaf.not_deletable();
                    *automorphism = leaf;
                    automorphism->inverse();
                    automorphism->compose(canon_leaf);//enqueue_fail_point_sz
                    if(enqueue_fail_point_sz > 0) {
                        if(config.CONFIG_THREADS_COLLABORATE) {
                                    w->G->first_level_points.try_enqueue_bulk(*w->ptok, w->enqueue_space, enqueue_fail_point_sz);
                        }
                    }
                    if(!first_level_succ->get(first_level_point)) {
                        if(config.CONFIG_THREADS_COLLABORATE) {
                                w->G->first_level_points.try_enqueue(*w->ptok, std::tuple<int, int>(w->first_level,first_level_point));
                        }
                        w->first_level_succ_point = first_level_point;
                        w->first_level_sz += 1;
                        first_level_succ->set(first_level_point);
                    }
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
            //changes.clear();
            bool comp = I->write_top_and_compare(INT32_MAX);
            //assert(init_color_class.size() == 2);
            comp = comp && R->refine_coloring(g, c, nullptr, I, init_color_class, false);
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
            init_color_class = -1;
            int rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
            int v = c->lab[rpos];
           // if(switches->done_shared_group && config.CONFIG_THREADS_COLLABORATE) // ToDo why is this problematic?
            //    v = (*w->shared_orbit)[v];

            //assert(rpos == c->vertex_to_lab[v]);

            if(level == w->first_level)
                first_level_point = v;

            if(level == w->first_level + 1)
                second_level_point = v;

            // individualize random vertex of class
            int newpos = R->individualize_vertex(c, v);
            last_op = OP_I;
//            assert(init_color_class >= 0);
            init_color_class = newpos;
            assert(c->vertex_to_col[v] > 0);
            //init_color_class.push_back(c.vertex_to_col[c.lab[labpos - 1]]);
            if (!compare) { // base points
                //automorphism->map.push_back(v);
                automorphism->map[automorphism->map_sz] = v;
                automorphism->map_sz += 1;
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

void auto_blaster::fast_automorphism_non_uniform(sgraph* g, bool compare, invariant* canon_I, bijection* canon_leaf,
        bijection* automorphism, int *restarts, bool *done, int selector_seed, auto_workspace* w) {
    bool backtrack = false;
    bool skipped_level = false;
    //std::list<std::pair<int, int>> changes;

    refinement *R = &w->R;
    selector *S = &w->S;
    coloring *c = &w->c;
    invariant *I = &w->I;

    coloring  *start_c     = &w->skip_c;
    invariant *start_I     = &w->skip_I;
    mschreier *group_level = w->skip_schreier_level;

    automorphism->non_uniform = false;

    bool base_aligned = true;

    S->empty_cache();
    int init_color_class;

    *restarts = 0;
    int level = w->first_skiplevel;


    if (compare) {
        start_I->set_compare_invariant(canon_I);
    } else {
        start_I->create_vector();
        automorphism->map = new int[g->v_size];
        automorphism->map_sz = 0;
    }
    ir_operation last_op = OP_R;
    *I = *start_I;
    c->copy_force(start_c);


    {
        S->empty_cache();
        while(w->first_skiplevel <= w->skiplevels) {
            proceed_state(w, g, start_c, start_I, w->my_base_points[w->first_skiplevel - 1], nullptr);
            w->first_skiplevel += 1;
            if(!w->is_foreign_base) {
                w->skip_schreier_level = w->skip_schreier_level->next;
            }
        }

        // initialize a search state
        c->copy_force(start_c);
        *I = *start_I;

        backtrack    = false;
        base_aligned = true;
        level = w->first_skiplevel;
        if(!w->is_foreign_base)
            group_level = w->skip_schreier_level;
    }

    while (true) {
        if(*done) return;
        if (backtrack) {
            //if(w->skiplevels > 0)
             //   std::cout << "wrong guess" << w->skiplevels << std::endl;
            if(*done) return;
            if(*restarts % 3 == 0) {
                w->skiplevels += 1;
            }
            S->empty_cache();
            if(w->first_skiplevel <= w->skiplevels) {
                proceed_state(w, g, start_c, start_I, w->my_base_points[w->first_skiplevel - 1], nullptr);
                w->first_skiplevel += 1;
                if(!w->is_foreign_base) {
                    w->skip_schreier_level = w->skip_schreier_level->next;
                }
            }

            // initialize a search state
            *restarts += 1;
            c->copy_force(start_c);
            *I = *start_I;

            backtrack    = false;
            base_aligned = true;
            level = w->first_skiplevel;
            if(!w->is_foreign_base)
                group_level = w->skip_schreier_level;
        }

        int s;
        if (!backtrack) {
            s = S->select_color(g, c, selector_seed);
            if (s == -1) {
                if (compare) {
                    // we can derive an automorphism!
                    //std::cout << *restarts << ", " << reguess << ", " << w->first_level << std::endl;
                    w->measure2 += 1;
                    bijection leaf;
                    leaf.read_from_coloring(c);
                    leaf.not_deletable();
                    *automorphism = leaf;
                    automorphism->inverse();
                    automorphism->compose(canon_leaf);//enqueue_fail_point_sz
                    automorphism->non_uniform = (w->skiplevels > 0);
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
            //changes.clear();
            bool comp = I->write_top_and_compare(INT32_MAX);
            //assert(init_color_class.size() == 2);
            comp = comp && R->refine_coloring(g, c, nullptr, I, init_color_class, false);
            comp = comp && I->write_top_and_compare(INT32_MAX);
            comp = comp && I->write_top_and_compare(INT32_MIN);
            //R.complete_colorclass_invariant(g, &c, &I);
            last_op = OP_R;
            assert(skipped_level?comp:true);
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
            int rpos;
            int v;

            if(level <= w->skiplevels) {
                v = w->my_base_points[level - 1];
                assert(c->vertex_to_col[v] == s);
                skipped_level = true;
                base_aligned  = true;
            } else {
                rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
                v = c->lab[rpos];
                skipped_level = false;
            }

            // check if base point can be chosen instead
            if(!w->is_foreign_base) {
                if (group_level->vec[v] && base_aligned) {
                    v = group_level->fixed;// choose base point
                    if(level == w->skiplevels + 1) {
                        bool total_orbit = (c->ptn[s] + 1 == group_level->fixed_orbit_sz);
                        //bool total_orbit = true;
                        /*for(int i = s; i < s + c->ptn[s] + 1; ++i)
                            total_orbit = total_orbit && group_level->vec[c->lab[i]];*/
                        if(total_orbit)
                            w->skiplevels += 1;
                    }
                } else {
                    base_aligned = false;
                }
            }

            int newpos = R->individualize_vertex(c, v);
            last_op = OP_I;
//            assert(init_color_class >= 0);
            init_color_class = newpos;
            assert(c->vertex_to_col[v] > 0);
            //init_color_class.push_back(c.vertex_to_col[c.lab[labpos - 1]]);
            if (!compare) { // base points
                //automorphism->map.push_back(v);
                automorphism->map[automorphism->map_sz] = v;
                automorphism->map_sz += 1;
            }
            level += 1;

            bool comp = I->write_top_and_compare(INT32_MIN);
            comp && I->write_top_and_compare(INT32_MIN);
            assert(skipped_level?comp:true);
            if(!comp) {
                backtrack = true;
                continue;
            }
            if(!w->is_foreign_base)
                group_level = group_level->next;
        }
    }
}

// if canonical vertex is v, true is returned and orbit contains vertices of the orbit of v
// if the canonical vertex is not v, false is returned

bool auto_blaster::get_orbit(auto_workspace* w, int* base, int base_sz, int v, work_list* orbit, bool reuse_generators) {
    orbit->reset();
    // test if identity
    if(*w->shared_generators_size == 0) {
        orbit->push_back(v);
        return true;
    }

    // if level == 1 we can return shared orbit
    if(base_sz == 0) {
        int map_v = (*w->shared_orbit)[v];
        assert(v >= 0);
        assert(v <= w->G_->domain_size);
        assert(map_v >= 0);
        assert(map_v <= w->G_->domain_size);

        if(map_v != v)
            return false;
    } else { // if level > 1, we collect generators that fix base and perform orbit algorithm on v
        // collect generators
        if(!reuse_generators) {
            w->generator_fix_base_size = 0;
            mpermnode *it = *w->shared_generators;
            do {
                // does it fix base?
                bool it_fixes_base = true; // do not need this variable
                for (int i = 0; i < base_sz; ++i) {
                    int b = base[i];
                    if (it->p[b] != b) {
                        it_fixes_base = false;
                        break;
                    }
                }

                if (it_fixes_base) {
                    w->generator_fix_base[w->generator_fix_base_size] = it;
                    w->generator_fix_base_size += 1;
                }
                it = it->next;
            } while (it != *w->shared_generators);
        } //else {
          //  if(w->canonical_v != v)
        //        return false;
        //}

        // do orbit algorithm on v
        if(w->generator_fix_base_size > 0) {
            int min_v = v; // find canonical v
            w->orbit_vertex_worklist.reset();
            w->orbit_considered.reset();
            w->orbit_vertex_worklist.push_back(v);
            w->orbit_considered.set(v);
            orbit->push_back(v);

            while(!w->orbit_vertex_worklist.empty()) {
                int next_v = w->orbit_vertex_worklist.pop_back();
                if(next_v < min_v)
                    min_v = next_v;

                // apply all generators exhaustively on v
                for(int j = 0; j < w->generator_fix_base_size; ++j) {
                    int mapped_v = w->generator_fix_base[j]->p[next_v];
                    if(!w->orbit_considered.get(mapped_v)) {
                        w->orbit_considered.set(mapped_v);
                        w->orbit_vertex_worklist.push_back(mapped_v);
                        orbit->push_back(mapped_v);
                    }
                }
            }
            w->canonical_v = min_v;
            if(min_v == v)
                return true;
            //std::cout << "safed" << std::endl;
            return false;
        } // else is identity again (below)
    }

    // return identity
    orbit->push_back(v);
    return true;
}

bool auto_blaster::bfs_chunk(sgraph* g, invariant* canon_I, bijection* canon_leaf, bool *done, int selector_seed, auto_workspace* w) {
    bfs* BFS = w->BW;
    int level        = BFS->BW.current_level;
    int target_level = BFS->BW.target_level;
    if(level == target_level) return false; // we are done with BFS!

    // initialize structures
    if(!w->init_bfs) {
        int chunk_sz = w->BW->BW.chunk_size;
        w->todo_dequeue = new std::tuple<bfs_element*, int>[chunk_sz];
        w->todo_deque_sz = chunk_sz;
        w->todo_elements = new std::tuple<bfs_element *, int>[chunk_sz * 8];
        w->todo_elements_sz = chunk_sz * 8;
        w->finished_elements = new std::pair<bfs_element *, int>[chunk_sz];
        w->finished_elements_sz = chunk_sz;
        w->init_bfs = true;
        w->orbit.initialize(g->v_size);
        w->orbit_considered.initialize(g->v_size);
        w->orbit_vertex_worklist.initialize(g->v_size);
        w->generator_fix_base = new mpermnode*[*w->shared_generators_size];
        w->changes.initialize((int) floor(sqrt(g->v_size)) + 1, g->v_size);
    }

    // try to dequeue a chunk of work
    // ToDo: try out next levels as well! (but should check target_level, though)
    size_t num = BFS->BW.bfs_level_todo[level].try_dequeue_bulk(w->todo_dequeue, w->BW->BW.chunk_size);
    int finished_elements_sz = 0;
    int todo_elements_sz = 0;

    for(int i = 0; i < num; ++i) {
        bfs_element *elem = std::get<0>(w->todo_dequeue[i]);
        int v             = std::get<1>(w->todo_dequeue[i]);

        // check orbit
        bool comp = get_orbit(w, elem->base, elem->base_sz, v, &w->orbit, w->prev_bfs_element == elem); // todo: need to consider weights
        // I could prune the todos relatively cheaply now (just a further restriction of generators, and they will not
        // appear in queues at all)

        if(comp) {
            // copy to workspace
            if (w->prev_bfs_element != elem) { // <-> last computed base is the same!
                w->work_c->copy_force(elem->c);
                w->prev_bfs_element = elem;
                *w->work_I = *elem->I;
                w->work_I->set_compare_invariant(canon_I);
            } else {
                // ToDo: check if w->changes did not overflow
                *w->work_I = *elem->I;
                w->work_I->set_compare_invariant(canon_I);
                //if(w->changes.did_overflow()) {
                    w->work_c->copy(elem->c);
                //} else {
                  // w->R.undo_changes_with_reference(&w->changes, w->work_c, elem->c);
                //};
            }

            // compute next coloring
            // ToDo: track changes and undo if !comp
            w->changes.reset();
            comp = comp && proceed_state(w, g, w->work_c, w->work_I, v, nullptr); // &w->changes
        }
        // not equal to canonical invariant?
        if (!comp) {
            // throw this node away, but keep track of that we computed it!
            w->finished_elements[finished_elements_sz] = std::pair<bfs_element *, int>(nullptr, 0);
            finished_elements_sz += 1;
            continue;
        }

        // create node
        bfs_element *next_elem = new bfs_element;
        next_elem->c = w->work_c;
        next_elem->I = w->work_I;
        next_elem->init_c = true;
        next_elem->init_I = true;
        next_elem->level = level + 1;
        next_elem->base_sz = elem->base_sz + 1;
        next_elem->base = new int[next_elem->base_sz];
        for(int j = 0; j < elem->base_sz; ++j)
            next_elem->base[j] = elem->base[j];
        next_elem->base[next_elem->base_sz - 1] = v;

        // create todos for this node
        int sz = 0;
        //if(level != target_level - 1) {
            w->S.empty_cache();
            int c = w->S.select_color(g, w->work_c, selector_seed);
            //assert(c != -1);
            for (int i = c; i < c + w->work_c->ptn[c] + 1; ++i) {
                int next_v = w->work_c->lab[i];
                sz += 1;
                if (todo_elements_sz == w->todo_elements_sz) {
                    BFS->BW.bfs_level_todo[level + 1].enqueue_bulk(w->todo_elements, todo_elements_sz);
                    todo_elements_sz = 0;
                }
                w->todo_elements[todo_elements_sz] = std::tuple<bfs_element *, int>(next_elem, next_v);
                todo_elements_sz += 1;
            }
        //}

        w->finished_elements[finished_elements_sz] = std::pair<bfs_element *, int>(next_elem, sz);
        finished_elements_sz += 1;
        w->work_c = new coloring;
        w->work_I = new invariant;
    }

    BFS->BW.bfs_level_finished_elements[level].enqueue_bulk(w->finished_elements, finished_elements_sz);
    BFS->BW.bfs_level_todo[level + 1].enqueue_bulk(w->todo_elements, todo_elements_sz);

    return true;
}

void auto_blaster::sample_pipelined(sgraph* g_, bool master, shared_switches* switches, pipeline_group* G, coloring* start_c, bijection* canon_leaf, invariant* canon_I,
                                    com_pad* communicator_pad, int communicator_id, int** shared_orbit, bfs* bwork, mpermnode** gens, int* shared_group_size) {
    // find comparison leaf
    std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
    sgraph* g = g_;

    bool* done      = &switches->done;
    bool* done_fast = &switches->done_fast;
    int _shared_group_size = false;
    mpermnode* _gens = nullptr;

    if(config.CONFIG_THREADS_COPYG && !master) {
        g = new sgraph;
        g->copy_graph(g_);
    }

    std::vector<std::thread> work_threads;
    bijection base_points;
    bool trash_bool = false;
    int trash_int = 0;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() * ((communicator_id + 3)* 5135235);
    //std::cout << seed << std::endl;
    //std::default_random_engine re = std::default_random_engine(seed);
    int selector_seed = seed;
    com_pad pad;

   // start_I.push_level();

    auto_workspace W;
    invariant start_I;
   // bool* p = new bool[g->v_size * 2];
    //W.first_level_fail.initialize_from_array(p, g->v_size);
    //W.first_level_succ.initialize_from_array(p + g->v_size, g->v_size);

    W.first_level_fail.initialize(g->v_size);
    W.first_level_succ.initialize(g->v_size);

    W.dequeue_space    = new std::pair<int, int>[2048];
    W.dequeue_space_sz = 2048;
    W.enqueue_space    = new std::pair<int, int>[128];
    W.enqueue_space_sz = 128; // <- choose this dynamic?

    int* shrd_orbit;

    if(master) {
        shrd_orbit = new int[g->v_size];
        for(int i = 0; i < g->v_size; ++i) {
            shrd_orbit[i] = i;
        }

        canon_I    = new invariant;
        canon_leaf = new bijection;
        start_c = new coloring;
        g->initialize_coloring(start_c);
        W.start_c = start_c;
        start_I.create_vector();
        W.R.refine_coloring_first(g, start_c, -1);
        double cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "[T] Color ref: " << cref / 1000000.0 << "ms" << std::endl;


        find_automorphism_prob(&W, g, false, canon_I, canon_leaf, &base_points, &trash_int, switches,
                               selector_seed);
        cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "[T] Canonical leaf found: " << cref / 1000000.0 << "ms" << std::endl;
        G = new pipeline_group(g->v_size);
        base_points.not_deletable();
        G->base_size = base_points.map_sz;
        G->b = base_points.map;
        cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "[T] Group created: " << cref / 1000000.0 << "ms" << std::endl;

        // initialize BFS
        bfs_element* root_elem = new bfs_element;
        root_elem->id = 0;
        root_elem->c = new coloring;
        root_elem->I = new invariant;
        root_elem->c->copy_force(start_c);
        root_elem->base_sz = 0;
        *root_elem->I = start_I;
        W.S.empty_cache();
        int init_c = W.S.select_color(g, start_c, selector_seed);

        W.BW = new bfs();
        W.BW->initialize(root_elem, init_c, g->v_size, G->base_size);

        //std::cout << "Launching refinement worker..." << std::endl;
        for(int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++)
            pad.emplace_back(moodycamel::ConcurrentQueue<std::pair<int, int>>(g->v_size, config.CONFIG_THREADS_REFINEMENT_WORKERS, config.CONFIG_THREADS_REFINEMENT_WORKERS));
        for(int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++)
            work_threads.emplace_back(std::thread(&auto_blaster::sample_pipelined, auto_blaster(), g, false, switches, G, start_c, canon_leaf, canon_I, &pad, i, &shrd_orbit, W.BW, &_gens, &_shared_group_size));
        //std::cout << "Refinement workers (" << config.CONFIG_THREADS_REFINEMENT_WORKERS << ")" << std::endl;
        cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "[T] Refinement workers created: " << cref / 1000000.0 << "ms" << std::endl;

        W.communicator_pad = &pad;
        W.communicator_id  = -1;
        W.shared_orbit = &shrd_orbit;
        W.shared_generators = &_gens;
        W.shared_generators_size = &_shared_group_size;
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
        G->launch_pipeline_threads(switches, &W);
        double cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
       // std::cout << "Pipeline group created: " << cref / 1000000.0 << "ms" << std::endl;
        G->pipeline_stage(0, switches, &W);
        // run algorithm
        //std::cout << "Sampled leaves (local thread): \t" << sampled_paths << std::endl;
        //std::cout << "Sampled leaves (all threads): \t"  << sampled_paths_all + sampled_paths << std::endl;
        //std::cout << "Restart probing (local thread):\t" << restarts << std::endl;
        //std::cout << "Schreier Filtercount: \t" << mfiltercount << std::endl;
        //std::cout << "Schreier Multcount: \t"   << mmultcount << std::endl; // <- ToDo: pick elements to minimize this?
        malloc_lock.lock();
        std::cout << "Base size:  " << G->base_size << std::endl;
        std::cout << "Group size: ";
        G->print_group_size();
        std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
        malloc_lock.unlock();
        G->join_threads();
        while(!work_threads.empty()) {
            work_threads[work_threads.size()-1].join();
            work_threads.pop_back();
        }
        cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "Join: " << cref / 1000000.0 << "ms" << std::endl;
        delete G;
    } else {
        W.shared_orbit = shared_orbit;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                      SLAVE THREAD
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int idle_ms = 0;
        W.start_c = start_c;
       // std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();

       // make my own canonical leaf...
       // canon_I    = new invariant;
       // canon_leaf = new bijection;
        coloring* start_c = W.start_c;
        W.start_c = new coloring;
        W.start_c->copy_force(start_c);
        W.skip_c.copy_force(start_c);
        W.work_c = new coloring;
        W.work_I = new invariant;

        W.communicator_pad = communicator_pad;
        W.communicator_id  = communicator_id;
        W.ctok = new moodycamel::ConsumerToken((*communicator_pad)[communicator_id]);
        W.ptok = new moodycamel::ProducerToken(G->first_level_points);
        W.base_size = G->base_size;
        W.G = G;
        W.BW = bwork;
        W.shared_generators = gens;
        W.shared_generators_size = shared_group_size;
        /*for(int i = 0; i < W.communicator_pad->size(); ++i) {
            if (communicator_id != i) {
                W.ptoks.emplace_back(new moodycamel::ProducerToken((*communicator_pad)[i]));
            } else {
                W.ptoks.emplace_back(nullptr);
            }
        }*/

        double cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
        //std::cout << "[" << communicator_id << "]: Startup: " << cref / 1000000.0 << "ms" << std::endl;

        invariant* _canon_I    = new invariant;
        _canon_I->create_vector();
        bijection* _canon_leaf = new bijection;

        W.skiplevels = 0;
        bool switched1 = false;
        bool switched2 = false;

        if(W.skiplevels == G->base_size - 1) {
            if(communicator_id == 0) std::cout << "[N] Skipped non-uniform automorphism search" << std::endl;
            G->skip_shared_group();
            switched1 = true;
            *done_fast   = true;
        } else {
            find_automorphism_prob(&W, g, false, _canon_I, _canon_leaf, &base_points, &trash_int, switches,
                                   selector_seed);
            W.my_base_points    = base_points.map;
            W.my_base_points_sz = base_points.map_sz;
            std::cout << "Found canonical leaf." << std::endl;
        }

        int earliest_found = -1;
        int n_found = 0;
        while(!(*done)) {
            bijection automorphism;
            //std::chrono::high_resolution_clock::time_point inner_timer = std::chrono::high_resolution_clock::now();
            if(config.CONFIG_IR_FAST_AUTOPRE && !(*done_fast)) { // detect when done
                fast_automorphism_non_uniform(g, true, _canon_I, _canon_leaf, &automorphism, &restarts, done_fast, selector_seed, &W); // <- we should already safe unsuccessfull / succ first level stuff here
                //std::cout << "Found automorphism." << std::endl;
                n_found += 1;
                if(n_found % 3 == 0 && W.skiplevels < W.my_base_points_sz)
                    W.skiplevels += 1;
                if((*done_fast && !automorphism.non_uniform )) continue;
                // set target level to earliest found automorphism skiplevel
                if(earliest_found < 0 && communicator_id == 0) {
                    int proposed_level = W.skiplevels + 1;
                    if(proposed_level == G->base_size)
                        proposed_level += 1;
                    W.BW->BW.target_level.store(proposed_level);
                    earliest_found = W.skiplevels;
                }
            } else if(W.BW->BW.current_level != W.BW->BW.target_level) {
                if(communicator_id == 0 && W.BW->BW.target_level < 0) {
                    int proposed_level = W.skiplevels + 1;
                    if(proposed_level == G->base_size)
                        proposed_level += 1;
                    W.BW->BW.target_level.store(proposed_level);
                }
                if(switches->done_shared_group && W.BW->BW.target_level >= 0) {
                    if(W.communicator_id == 0 && !switched1) {
                        switched1 = true;
                        cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
                        std::cout << "[N] Finished non-uniform automorphism search (" << *W.shared_generators_size << " generators)" << std::endl;
                        std::cout << "[N] Ended in skiplevel " << W.skiplevels << std::endl;
                        std::cout << "[T] " << cref / 1000000.0 << "ms" << std::endl;
                        std::cout << "[B] Determined target level: " << W.BW->BW.target_level << "" << std::endl;
                    }
                    bfs_chunk(g, canon_I, canon_leaf, done, selector_seed, &W);
                }
                continue;
            } else {
                if(W.communicator_id == 0 && !switched2) {
                    switched2 = true;
                    cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() - timer).count());
                    std::cout << "[T] " << cref / 1000000.0 << "ms" << std::endl;
                }
                //std::cout << W.BW->BW.current_level << std::endl;
                //std::cout << "Searching automorphism." << std::endl;
                find_automorphism_from_bfs(&W, g, true, canon_I, canon_leaf, &automorphism, &restarts, switches, selector_seed);
                //std::cout << "Found automorphism." << std::endl;
            }
            /*else {
                if(!switched1 && communicator_id == 0) {
                    cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
                    std::cout << "[T] Fast automorphisms done: " << cref / 1000000.0 << "ms" << std::endl;
                }
                if(!switched2 && (communicator_id == 0) && switches->done_shared_group) {
                    cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
                    std::cout << "[T] Shared group done: " << cref / 1000000.0 << "ms" << std::endl;
                    switched2 = true;
                }
                switched1 = true;
                find_automorphism_prob(&W, g, true, canon_I, canon_leaf, &automorphism, &restarts, switches, selector_seed);
            }*/


            if(sampled_paths == 0) {
                //cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
              //  std::cout << "[" << communicator_id << "]: First automorphism (probed): " << cref / 1000000.0 << "ms, " << W.first_level << ":" << W.first_level_sz << ":" << W.skiplevels << std::endl;
            }
            automorphism.not_deletable();
            bool test = G->add_permutation(&automorphism, &idle_ms, done);

            automorphism.map = new int[g->v_size];
            sampled_paths += 1;
        }
        delete W.start_c;
        return;
    }
}

void auto_blaster::sample_shared(sgraph* g_, bool master, shared_switches* switches, diy_group* G, coloring* start_c, bijection** canon_leaf, invariant** canon_I,
                                    com_pad* communicator_pad, int communicator_id, int** shared_orbit, bfs* bwork, mpermnode** gens, int* shared_group_size) {
    // find comparison leaf
    std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
    sgraph *g = g_;
    double cref;

    bool *done      = &switches->done;
    bool *done_fast = &switches->done_fast;
    int _shared_group_size = false;
    mpermnode *_gens       = nullptr;

    int* shrd_orbit;
    std::vector<std::thread> work_threads;
    bijection base_points;
    int trash_int = 0;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() * ((communicator_id * 5) * 5135235);
    int selector_seed = seed;
    com_pad pad;

    auto_workspace W;
    invariant start_I;

    W.first_level_fail.initialize(g->v_size);
    W.first_level_succ.initialize(g->v_size);

    W.dequeue_space = new std::pair<int, int>[2048];
    W.dequeue_space_sz = 2048;
    W.enqueue_space = new std::pair<int, int>[128];
    W.enqueue_space_sz = 128; // <- choose this dynamic?

    if (master) {
        shrd_orbit = new int[g->v_size];
        for(int i = 0; i < g->v_size; ++i) {
            shrd_orbit[i] = i;
        }

        int** shrd_orbit_ = new (int*);
        *shrd_orbit_ = shrd_orbit;

        switches->current_mode = modes::MODE_LEAF_TOURNAMENT;

        canon_I    = new invariant*;
        canon_leaf = new bijection*;
        start_c    = new coloring;
        g->initialize_coloring(start_c);
        W.start_c = start_c;
        start_I.create_vector();
        W.R.refine_coloring_first(g, start_c, -1);
        cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "[T] Color ref: " << cref / 1000000.0 << "ms" << std::endl;

        // create some objects that are initialized after tournament
        G = new diy_group(g->v_size);
        W.BW = new bfs();
        bwork = W.BW;

        // launch worker threads
        for (int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++)
            work_threads.emplace_back(
                    std::thread(&auto_blaster::sample_shared, auto_blaster(), g, false, switches, G, start_c,
                                canon_leaf, canon_I, &pad, i, shrd_orbit_, W.BW, &_gens, &_shared_group_size));
        cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "[T] Refinement workers created: " << cref / 1000000.0 << "ms" << std::endl;

        // set some workspace variables
        W.start_c = new coloring;
        W.start_c->copy_force(start_c);
        W.communicator_pad       = &pad;
        W.communicator_id        = -1;
        W.shared_orbit           = shrd_orbit_;
        W.shared_generators      = &_gens;
        W.shared_generators_size = &_shared_group_size;
    }

    int sampled_paths = 0;
    int restarts = 0;
    int idle_ms = 0;
    invariant *_canon_I;
    bijection *_canon_leaf;
    bool switched1 = false;
    bool switched2 = false;

    W.skip_c.copy_force(start_c);
    W.work_c = new coloring;
    W.work_I = new invariant;
    W.G_ = G;
    W.BW = bwork;
    W.skiplevels = 0;
    W.skip_schreier_level = G->gp;

    if (!master) {
        W.shared_orbit = shared_orbit;
        W.start_c = new coloring;
        W.start_c->copy_force(start_c);
        W.communicator_pad = communicator_pad;
        W.communicator_id = communicator_id;
        W.shared_generators = gens;
        W.shared_generators_size = shared_group_size;
    }

    {
        //W.base_size = G->base_size;
        _canon_I = new invariant;
        _canon_I->create_vector();
        _canon_leaf = new bijection;
        {
            find_automorphism_prob(&W, g, false, _canon_I, _canon_leaf, &base_points, &trash_int, switches, selector_seed);
            W.my_base_points    = base_points.map;
            W.my_base_points_sz = base_points.map_sz;
            W.is_foreign_base   = true;
        }

        // skip immediately if base size is 1
        if(master && W.my_base_points_sz == 1) {
            *canon_I    = _canon_I;
            *canon_leaf = _canon_leaf;
            G->initialize(g->v_size, &base_points);
            bfs_element *root_elem = new bfs_element;
            root_elem->id = 0;
            root_elem->c = new coloring;
            root_elem->I = new invariant;
            root_elem->c->copy_force(start_c);
            root_elem->base_sz = 0;
            *root_elem->I = start_I;
            W.S.empty_cache();
            int init_c = W.S.select_color(g, start_c, selector_seed);
            W.BW->initialize(root_elem, init_c, g->v_size, G->base_size);
            switches->done_created_group = true;
            int proposed_level = W.skiplevels + 1;
            if(proposed_level == G->base_size)
                proposed_level += 1;
            W.BW->BW.target_level.store(proposed_level);
            W.is_foreign_base = false;
            W.skip_schreier_level = G->gp;
            for(int i = 0; i < W.skiplevels; ++i)
                W.skip_schreier_level = W.skip_schreier_level->next;
            W.base_size = G->base_size;
            *W.shared_generators_size = 0;
            switches->done_shared_group = true;
            switches->current_mode = modes::MODE_BFS;
            std::cout << "[T] No tournament, created group by " << communicator_id << " with restarts " << restarts << std::endl;
            // ToDo: create BFS and group here
        } else if(W.my_base_points_sz == 1) {
            // wait until shared group created
            while(!switches->done_shared_group) continue;
            while(!switches->current_mode == modes::MODE_BFS || !switches->current_mode == modes::MODE_UNIFORM_PROBE) continue;
        } else {
            // why is fast_automorphism called in ranreg?
            //std::cout << W.my_base_points_sz << std::endl;
        }
    }

    int earliest_found = -1;
    int n_found = 0;
    int n_restarts = 0;
    bool foreign_base_done = false;

    while(!(*done)) {
        if(switches->done_fast)
            G->ack_done_shared();
        bijection automorphism;
        //std::chrono::high_resolution_clock::time_point inner_timer = std::chrono::high_resolution_clock::now();


        switch(switches->current_mode) {
            case modes::MODE_LEAF_TOURNAMENT:
                fast_automorphism_non_uniform(g, true, _canon_I, _canon_leaf, &automorphism, &restarts, done_fast, selector_seed, &W); // <- we should already safe unsuccessfull / succ first level stuff here
                if(n_found == 0) { // check if I won
                    // wait until everyone checked
                    while(!switches->check_leaf_tournament(communicator_id, restarts) && !switches->done_created_group) continue;
                    // check if I won, if yes: create group
                    if(switches->done_created_group) continue;
                    if(switches->win_id == communicator_id) {
                        *canon_I    = _canon_I;
                        *canon_leaf = _canon_leaf;
                        G->initialize(g->v_size, &base_points);
                        bfs_element *root_elem = new bfs_element;
                        root_elem->id = 0;
                        root_elem->c = new coloring;
                        root_elem->I = new invariant;
                        root_elem->c->copy_force(start_c);
                        root_elem->base_sz = 0;
                        *root_elem->I = start_I;
                        W.S.empty_cache();
                        int init_c = W.S.select_color(g, start_c, selector_seed);
                        W.BW->initialize(root_elem, init_c, g->v_size, G->base_size);
                        switches->done_created_group = true;
                        int proposed_level = W.skiplevels + 1;
                        if(proposed_level == G->base_size)
                            proposed_level += 1;
                        W.BW->BW.target_level.store(proposed_level);
                        W.is_foreign_base = false;

                        W.skip_schreier_level = G->gp;
                        for(int i = 0; i < W.skiplevels; ++i)
                            W.skip_schreier_level = W.skip_schreier_level->next;

                        W.base_size = G->base_size;
                        foreign_base_done = true;
                        switches->current_mode = modes::MODE_NON_UNIFORM_PROBE;
                        std::cout << "[T] Created group by " << communicator_id << " with restarts " << restarts << std::endl;
                    }

                    while(!(switches->done_created_group)) continue;
                }
                automorphism.foreign_base = true;
                n_restarts += restarts;
                n_found += 1;
                if(n_found % 3 == 0 && W.skiplevels < W.my_base_points_sz)
                    W.skiplevels += 1;
                if((*done_fast && !automorphism.non_uniform )) continue;
                break;

            case modes::MODE_NON_UNIFORM_PROBE:
                if(!foreign_base_done) {
                    fast_automorphism_non_uniform(g, true, _canon_I, _canon_leaf, &automorphism, &restarts, done_fast, selector_seed, &W);
                    automorphism.foreign_base = true;
                    n_restarts += restarts;
                } else {
                    fast_automorphism_non_uniform(g, true, *canon_I, *canon_leaf, &automorphism, &restarts, done_fast, selector_seed, &W);
                    automorphism.foreign_base = false;
                    n_restarts += restarts;
                }
                n_found += 1;
                if(n_found % 3 == 0 && W.skiplevels < W.my_base_points_sz)
                    W.skiplevels += 1;
                if((*done_fast && !automorphism.non_uniform )) continue;
                    break;

            case modes::MODE_BFS:
                if(W.BW->BW.current_level != W.BW->BW.target_level) {
                    if (communicator_id == -1 && W.BW->BW.target_level < 0) {
                        int proposed_level = W.skiplevels + 1;
                        if (proposed_level == G->base_size)
                            proposed_level += 1;
                        W.BW->BW.target_level.store(proposed_level);
                    }
                    if(switches->done_shared_group && W.BW->BW.target_level >= 0) {
                        if(W.communicator_id == 0 && !switched1) {
                            switched1 = true;
                            cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
                            std::cout << "[N] Finished non-uniform automorphism search (" << *W.shared_generators_size << " generators, " << n_restarts << " restarts)" << std::endl;
                            std::cout << "[N] Ended in skiplevel " << W.skiplevels << ", found " << n_found << std::endl;
                            std::cout << "[T] " << cref / 1000000.0 << "ms" << std::endl;
                            std::cout << "[B] Determined target level: " << W.BW->BW.target_level << "" << std::endl;
                        }
                        bfs_chunk(g, *canon_I, *canon_leaf, done, selector_seed, &W);
                        if(master)
                            bwork->work_queues();
                    }
                } else {
                    if(master)
                        switches->current_mode = modes::MODE_UNIFORM_PROBE;
                }
                continue;
                break;

            case modes::MODE_UNIFORM_PROBE:
                if(W.communicator_id == 0 && !switched2) {
                    switched2 = true;
                    cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() - timer).count());
                    std::cout << "[T] " << cref / 1000000.0 << "ms" << std::endl;
                }
                find_automorphism_from_bfs(&W, g, true, *canon_I, *canon_leaf, &automorphism, &restarts, switches, selector_seed);
                break;

            case modes::MODE_WAIT:
                continue;
        }

        if(switches->done)
            return;

        automorphism.not_deletable();
        //std::cout << automorphism.foreign_base << std::endl;
        bool test = true;
        if(switches->done_created_group) {
            test = G->add_permutation(&automorphism, &idle_ms, done);
        }

        if(!test && !foreign_base_done) {
            //std::cout << "foreign base done" << std::endl;
            W.skip_c.copy_force(start_c);
            W.skip_I = W.start_I;
            W.skiplevels = 0;
            if(W.BW->BW.target_level > 0)
                W.skiplevels = W.BW->BW.target_level;
            W.skip_schreier_level = G->gp;
            W.first_skiplevel = 1;
            W.my_base_points    = G->b;
            W.my_base_points_sz = G->base_size;
            W.is_foreign_base   = false;
            foreign_base_done = true;
            //std::cout << "[N] Switching to canonical search (" << W.communicator_id << ", " << n_found << " generators)" << std::endl;
        }

        // ToDo: sift some random elements!

        delete[] automorphism.map;
        automorphism.map = new int[g->v_size];
        automorphism.foreign_base = false;
        sampled_paths += 1;

        if(master) {
            //std::cout << "checking" << std::endl;
            G->manage_results(switches);
            if(switches->done_fast && !switches->done_shared_group) {
                // wait for ack of done_fast
                std::cout << "[N] Waiting for ACK" << std::endl;
                switches->current_mode = modes::MODE_WAIT;
                G->ack_done_shared();
                G->wait_for_ack_done_shared(config.CONFIG_THREADS_REFINEMENT_WORKERS + 1);
                *W.shared_generators      = G->gens;
                *W.shared_orbit           = G->gp->orbits;
                *W.shared_generators_size = G->number_of_generators();
                switches->done_shared_group = true;
                switches->current_mode = modes::MODE_BFS;
            }
            if(switches->done) {
                std::cout << "Base size:  " << G->base_size << std::endl;
                while(!work_threads.empty()) {
                    work_threads[work_threads.size()-1].join();
                    work_threads.pop_back();
                }
                std::cout << "Group size: ";
                G->print_group_size();
                std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
                cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
                std::cout << "Join: " << cref / 1000000.0 << "ms" << std::endl;
                break;
            }
        }
    }

    if(master)
        delete G;
    //delete W.start_c;
    return;
}