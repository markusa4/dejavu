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
#include "diy_group.h"
#include "lowdeg.h"
#include <pthread.h>
#include <tuple>

extern long mmultcount;
extern long mfiltercount;

int intRand(const int & min, const int & max, int seed) {
    static thread_local std::mt19937 generator(seed);
    std::uniform_int_distribution<int> distribution(min,max);
    return distribution(generator);
}

double doubleRand(const double & min, const double & max, int seed) {
    static thread_local std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distribution(min,max);
    return floor(distribution(generator));
}

bool auto_blaster::proceed_state(auto_workspace* w, sgraph* g, coloring* c, invariant* I, int v, change_tracker* changes) {
    if(changes != nullptr)
        changes->track(c->vertex_to_col[v]);

    int init_color_class = w->R.individualize_vertex(c, v);
    bool comp = I->write_top_and_compare(INT32_MIN);
    comp && I->write_top_and_compare(INT32_MIN);
    comp = comp && I->write_top_and_compare(INT32_MAX);

    if(!comp) return comp;

    comp = comp && w->R.refine_coloring(g, c, changes, I, init_color_class, changes != nullptr);
    bool tempcomp = comp;
    comp = comp && I->write_top_and_compare(INT32_MAX);
    comp = comp && I->write_top_and_compare(INT32_MIN);
    assert(tempcomp == comp);
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
    int level = w->first_level;

    if (compare) {
        automorphism->map    = new int[g->v_size];
        automorphism->map_sz = 0;
    }
    ir_operation last_op;

    // Pick start from BFS level instead!
    int bfs_level    = w->BW->BW.current_level - 1;
    int bfs_level_sz = w->BW->BW.level_sizes[bfs_level];

    mschreier* start_group_level = w->G_->gp;
    for(int i = 0; i < bfs_level; i++)
        start_group_level = start_group_level->next;
    mschreier* group_level = start_group_level;

    int rand_pos     = intRand(0, bfs_level_sz - 1, selector_seed);
    bfs_element* picked_elem = w->BW->BW.level_states[bfs_level][rand_pos];
    bool base_aligned = picked_elem->is_identity;

    *I = *picked_elem->I;
    c->copy_force(picked_elem->c);
    I->set_compare_invariant(canon_I);
    last_op = OP_R;
    level = bfs_level + 1;
    S->empty_cache();
    int base;
    double picked_weight, max_weight, rand_weight;

    while (true) {
        if(*done) return;
        if (backtrack) {
            bfs_level    = w->BW->BW.current_level - 1;
            bfs_level_sz = w->BW->BW.level_sizes[bfs_level];
            rand_pos     = intRand(0, bfs_level_sz - 1, selector_seed);
            picked_elem = w->BW->BW.level_states[bfs_level][rand_pos];
            group_level = start_group_level;
            base_aligned = picked_elem->is_identity;

            // consider the weight by redrawing
            picked_weight = picked_elem->weight;
            max_weight    = w->BW->BW.level_maxweight[bfs_level];
            assert(max_weight > 0);
            rand_weight   = doubleRand(1, max_weight, selector_seed);
            if(rand_weight > picked_weight) continue; // need to redraw

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

            if (group_level->vec[v] && base_aligned) {
                v = group_level->fixed;// choose base point
                if(level == w->skiplevels + 1 &&  (w->skiplevels < w->my_base_points_sz - 1)) {
                    bool total_orbit = (c->ptn[s] + 1 == group_level->fixed_orbit_sz);
                    if(total_orbit)
                        w->skiplevels += 1;
                }
            } else {
                base_aligned = false;
            }


            // individualize random vertex of class
            int newpos = R->individualize_vertex(c, v);
            last_op = OP_I;
            init_color_class = newpos;
            assert(c->vertex_to_col[v] > 0);
            if (!compare) { // base points
                automorphism->map[automorphism->map_sz] = v;
                automorphism->map_sz += 1;
            }

            group_level = group_level->next;
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
            if(*done) return;
            S->empty_cache();
            // initialize a search state
            *restarts += 1;
            c->copy(start_c);
            *I = *start_I;

            // invariant, hopefully becomes complete in leafs such that automorphisms can be found
            if (level == w->first_level + 1) {
                if (!first_level_fail->get(base)) {
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
            if (s == -1 && last_op == OP_R) {
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
        bijection* automorphism, int *restarts, bool *done, int selector_seed, auto_workspace* w, int tolerance) {
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

    if(w->skiplevels >= w->my_base_points_sz - 1)
        w->skiplevels = w->my_base_points_sz - 2;

    {
        S->empty_cache();
        while(w->first_skiplevel <= w->skiplevels) {
            if(*done) return;
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
            if((*restarts % (5 * tolerance) == 0) && (w->skiplevels < w->my_base_points_sz - 1)) {
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
            if (s == -1 && last_op == OP_R) {
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
                base_aligned  = true;
            } else {
                rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
                v = c->lab[rpos];
            }

            // check if base point can be chosen instead
            if(!w->is_foreign_base) {
                if (group_level->vec[v] && base_aligned) {
                    v = group_level->fixed;// choose base point
                    if(level == w->skiplevels + 1 &&  (w->skiplevels < w->my_base_points_sz - 1)) {
                        bool total_orbit = (c->ptn[s] + 1 == group_level->fixed_orbit_sz);
                        if(total_orbit)
                            w->skiplevels += 1;
                    }
                } else {
                    base_aligned = false;
                }
            }

            int newpos = R->individualize_vertex(c, v);
            last_op = OP_I;
            init_color_class = newpos;
            assert(c->vertex_to_col[v] > 0);

            if (!compare) { // base points
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

bool auto_blaster::get_orbit(auto_workspace* w, int* base, int base_sz, int v, int v_base, work_list* orbit, bool reuse_generators) {
    orbit->reset();
    // test if identity
    if(*w->shared_generators_size == 0) {
        orbit->push_back(v);
        return true;
    }

    // if level == 1 we can return shared orbit
    if(base_sz == 0) {
        // ToDo: difficult with weights here, actually... need one pass of orbits to create orbit sizes
        // ToDo: and also map the orbit of base_v to base_v...
        int map_v = (*w->shared_orbit)[v];
        assert(v >= 0);
        assert(v <= w->G_->domain_size);
        assert(map_v >= 0);
        assert(map_v <= w->G_->domain_size);

        orbit->push_back(v);
        orbit->cur_pos = (*w->shared_orbit_weights)[map_v];
        assert(orbit->cur_pos > 0);
        return (map_v == v);
    } else { // if level > 1, we collect generators that fix base and perform orbit algorithm on v
        // collect generators
        if(!reuse_generators) {
            w->generator_fix_base_size = 0;
            mpermnode *it = *w->shared_generators;
            do {
                // does it fix base?
                // do not need this variable
                int i;
                for (i = 0; i < base_sz; ++i) {
                    int b = base[i];
                    assert(b < w->G_->domain_size && b >= 0);
                    if (it->p[b] != b) {
                        break;
                    }
                }

                if (i == base_sz) {
                    assert(w->generator_fix_base_size < w->generator_fix_base_alloc);
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
            int  min_v = v; // find canonical v
            w->orbit_vertex_worklist.reset();
            w->orbit_considered.reset();
            w->orbit_vertex_worklist.push_back(v);
            w->orbit_considered.set(v);
            orbit->push_back(v);

            while(!w->orbit_vertex_worklist.empty()) {
                int next_v = w->orbit_vertex_worklist.pop_back();
                if((next_v < min_v && min_v != v_base) || next_v == v_base)
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
            return (min_v == v);
        } // else is identity again (below)
    }

    // return identity
    orbit->push_back(v);
    return true;
}

bool auto_blaster::bfs_chunk(sgraph* g, invariant* canon_I, bijection* canon_leaf, bool *done, int selector_seed, auto_workspace* w) {
    thread_local bool done_test = false;

    bfs* BFS = w->BW;
    int level        = BFS->BW.current_level;
    int target_level = BFS->BW.target_level;
    if(level == target_level) return false; // we are done with BFS!

    // initialize structures
    if(!w->init_bfs) {
        int chunk_sz = w->BW->BW.chunk_size;
        w->todo_dequeue = new std::pair<bfs_element*, int>[chunk_sz];
        w->todo_deque_sz = chunk_sz;
        w->todo_elements = new std::pair<bfs_element *, int>[chunk_sz * 8];
        w->todo_elements_sz = chunk_sz * 8;
        w->finished_elements = new std::pair<bfs_element *, int>[chunk_sz + 1];
        w->finished_elements_sz = chunk_sz;
        w->init_bfs = true;
        w->orbit.initialize(g->v_size);
        w->orbit_considered.initialize(g->v_size);
        w->orbit_vertex_worklist.initialize(g->v_size);
        w->generator_fix_base = new mpermnode*[*w->shared_generators_size];
        w->generator_fix_base_alloc = *w->shared_generators_size;
    }

    if(w->generator_fix_base_alloc < *w->shared_generators_size) {
        delete[] w->generator_fix_base;
        w->generator_fix_base = new mpermnode*[*w->shared_generators_size];
        w->generator_fix_base_alloc = *w->shared_generators_size;
        w->prev_bfs_element = nullptr;
    }

    // try to dequeue a chunk of work
    // ToDo: try out next levels as well! (but should check target_level, though)
    size_t num = BFS->BW.bfs_level_todo[level].try_dequeue_bulk(w->todo_dequeue, w->BW->BW.chunk_size);
    int finished_elements_sz = 0;
    int todo_elements_sz = 0;

    int finished_elements_null_buffer = 0;

    for(int i = 0; i < num; ++i) {
        bfs_element *elem = w->todo_dequeue[i].first;
        int v             = w->todo_dequeue[i].second;
        bool is_identity  = elem->is_identity && (v == w->my_base_points[elem->base_sz]);

        // check orbit
        bool comp = get_orbit(w, elem->base, elem->base_sz, v, w->my_base_points[elem->base_sz], &w->orbit, w->prev_bfs_element == elem);
        assert(!comp?(!is_identity) : true);
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
                *w->work_I = *elem->I;
                w->work_I->set_compare_invariant(canon_I);
                w->work_c->copy(elem->c);
            }

            // compute next coloring
            comp = comp && proceed_state(w, g, w->work_c, w->work_I, v, nullptr); // &w->changes
        }

        assert(elem->base_sz < w->G_->base_size);
        assert(!comp?(!is_identity) : true);

        // not equal to canonical invariant?
        if (!comp) {
            // throw this node away, but keep track of that we computed it!
            finished_elements_null_buffer += 1;
            continue;
        }

        // create node
        bfs_element *next_elem = new bfs_element;
        next_elem->c = w->work_c;
        next_elem->I = w->work_I;
        next_elem->init_c = true;
        next_elem->init_I = true;
        next_elem->is_identity = is_identity;
        next_elem->level = level + 1;
        next_elem->base_sz = elem->base_sz + 1;
        next_elem->base = new int[next_elem->base_sz];
        next_elem->init_base = true;
        for(int j = 0; j < elem->base_sz; ++j) {
            assert(elem->base[j] >= 0 && elem->base[j] < g->v_size);
            next_elem->base[j] = elem->base[j];
        }
        assert(v >= 0 && v < g->v_size);
        next_elem->base[next_elem->base_sz - 1] = v;
        next_elem->weight = elem->weight * w->orbit.cur_pos;
        next_elem->parent_weight = elem->weight;
        next_elem->parent = elem;

        // create todos for this node
        int sz = 0;
        //if(level != target_level - 1) {
            w->S.empty_cache();
            int c = w->S.select_color(g, w->work_c, selector_seed);
            next_elem->target_color = c;
            // ToDo: sort lab?
            //std::sort(w->work_c->lab + c, w->work_c->lab + c + w->work_c->ptn[c]);
            //for(int ii = c; ii < c + w->work_c->ptn[c]; ++ii)
            //    w->work_c->vertex_to_lab[w->work_c->lab[ii]] = ii;

            //assert(c != -1);
            sz += w->work_c->ptn[c] + 1;
            /*for (int i = c; i < c + w->work_c->ptn[c] + 1; ++i) {
                int next_v = w->work_c->lab[i];
                sz += 1;
                if (todo_elements_sz == w->todo_elements_sz) {
                    BFS->BW.bfs_level_todo[level + 1].enqueue_bulk(w->todo_elements, todo_elements_sz);
                    todo_elements_sz = 0;
                }
                w->todo_elements[todo_elements_sz] = std::pair<bfs_element *, int>(next_elem, next_v);
                todo_elements_sz += 1;
            }*/
        //}

        w->finished_elements[finished_elements_sz] = std::pair<bfs_element *, int>(next_elem, sz);
        finished_elements_sz += 1;
        w->work_c = new coloring;
        w->work_I = new invariant;
    }

    if(finished_elements_null_buffer > 0) {
        w->finished_elements[finished_elements_sz] = std::pair<bfs_element *, int>(nullptr, finished_elements_null_buffer);
        finished_elements_sz += 1;
    }

    if(finished_elements_sz > 0)
        BFS->BW.bfs_level_finished_elements[level].enqueue_bulk(w->finished_elements, finished_elements_sz);
    //BFS->BW.bfs_level_todo[level + 1].enqueue_bulk(w->todo_elements, todo_elements_sz);

    return true;
}

void auto_blaster::bfs_reduce_tree(auto_workspace* w) {
    if(!w->init_bfs) {
        int chunk_sz = w->BW->BW.chunk_size;
        w->todo_dequeue = new std::pair<bfs_element*, int>[chunk_sz];
        w->todo_deque_sz = chunk_sz;
        w->todo_elements = new std::pair<bfs_element *, int>[chunk_sz * 8];
        w->todo_elements_sz = chunk_sz * 8;
        w->finished_elements = new std::pair<bfs_element *, int>[chunk_sz + 1];
        w->finished_elements_sz = chunk_sz;
        w->init_bfs = true;
        w->orbit.initialize(w->G_->domain_size);
        w->orbit_considered.initialize(w->G_->domain_size);
        w->orbit_vertex_worklist.initialize(w->G_->domain_size);
        w->generator_fix_base = new mpermnode*[*w->shared_generators_size];
        w->generator_fix_base_alloc = *w->shared_generators_size;
        std::cout << "late init" << std::endl;
    }

    if(w->generator_fix_base_alloc < *w->shared_generators_size) {
        delete[] w->generator_fix_base;
        w->generator_fix_base = new mpermnode*[*w->shared_generators_size];
        w->generator_fix_base_alloc = *w->shared_generators_size;
        w->prev_bfs_element = nullptr;
    }

    // step1: check on level i if orbits can be reduced
    int i, j, v;
    bool updated = false;
    bfs_workspace* BW = &w->BW->BW;
    int active = 0;
    bool found_canon = false;
    for (int j = 0; j < BW->level_sizes[BW->current_level - 1]; ++j) {
        bfs_element *elem = BW->level_states[BW->current_level - 1][j];
        active += (elem->weight > 0);
        found_canon = found_canon || (elem->base[elem->base_sz - 1] == w->G_->b[elem->base_sz - 1]);
    }
    assert(found_canon);
    std::cout << "top level active: " << active << std::endl;


    for(i = 1; i < BW->current_level; ++i) {
        //std::cout << "checking level " << i << ", sz " << BW->level_sizes[i] << std::endl;
        if(i == 1) {
            for(j = 0; j < BW->level_sizes[i]; ++j) {
                bfs_element* elem = BW->level_states[i][j];
                if(elem->weight > 0) {
                    v = elem->base[elem->base_sz - 1];
                    assert(v >= 0 && v < w->G_->domain_size);
                    if ((*w->shared_orbit)[v] != v) {
                        elem->weight = 0;
                        assert(v != w->G_->b[elem->base_sz - 1]);
                        //std::cout << "reducing " << elem->weight << " " << std::endl;
                    } else {
                        if (elem->weight != (*w->shared_orbit_weights)[v]) {
                            elem->weight = (*w->shared_orbit_weights)[v];
                        }
                    }
                }
            }
        } else {
            // update weights...
            for(j = 0; j < BW->level_sizes[i]; ++j) {
                bfs_element *elem = BW->level_states[i][j];
                assert(elem->parent != NULL);
                if(elem->parent->weight == 0) {
                    elem->weight = 0;
                    continue;
                }
                if(elem->parent_weight != elem->parent->weight) {
                    elem->weight        = elem->parent->weight * (elem->weight / elem->parent_weight);
                    elem->parent_weight = elem->parent->weight;
                }
            }

            // reduce using orbits
            for(j = 0; j < BW->level_sizes[i]; ++j) {
                bfs_element *elem = BW->level_states[i][j];
                assert(elem->parent != NULL);
                if(elem->weight != 0) {
                    //int* orbits = mgetorbits(elem->base, elem->base_sz - 1, w->G_->gp, &w->G_->gens, w->G_->domain_size); // ToDo: not the right thing on canonical base...
                    assert(elem->base_sz == elem->level);
                    assert(elem->level == i);
                    assert(i > 1);
                    v = elem->base[elem->base_sz - 1];
                    assert((v >= 0) && (v < w->G_->domain_size));
                    w->orbit.reset();
                    assert(elem->base_sz - 1 > 0);
                    bool comp = get_orbit(w, elem->base, elem->base_sz - 1, v, w->my_base_points[elem->base_sz - 1], &w->orbit, false);
                    if(!comp) {
                        assert(v != w->G_->b[elem->base_sz - 1]);
                        elem->weight = 0;
                    } else {
                        elem->weight = elem->parent_weight * w->orbit.cur_pos;
                    }
                    w->orbit.reset();
                    w->orbit_vertex_worklist.reset();
                    w->orbit_considered.reset();
                }
            }
        }
    }

    active = 0;
    for (j = 0; j < BW->level_sizes[BW->current_level - 1]; ++j) {
        bfs_element *elem = BW->level_states[BW->current_level - 1][j];
        active += (elem->weight > 0);
    }

    std::cout << "top level active (after): " << active << std::endl;
    std::cout << "[B] Reducing queue from sz " << BW->bfs_level_todo[BW->current_level].size_approx();
    moodycamel::ConcurrentQueue<std::pair<bfs_element *, int>> throwaway_queue;
    BW->bfs_level_todo[BW->current_level].swap(throwaway_queue);


    int expecting_finished = 0;
    for (int j = 0; j < BW->level_sizes[BW->current_level - 1]; ++j) {
        bfs_element *elem = BW->level_states[BW->current_level - 1][j];
        if (elem->weight > 0) {
            int c = elem->target_color;
            int c_size = elem->c->ptn[c] + 1;
            for (i = c; i < c + c_size; ++i) {
                expecting_finished += 1;
                BW->bfs_level_todo[BW->current_level].enqueue(std::pair<bfs_element *, int>(elem, elem->c->lab[i]));
            }
        }
    }
    BW->level_expecting_finished[BW->current_level] = expecting_finished;

    std::cout << " to " << BW->bfs_level_todo[BW->current_level].size_approx() << std::endl;
    w->prev_bfs_element = nullptr;
    w->orbit.reset();
    w->orbit_vertex_worklist.reset();
    w->orbit_considered.reset();

    //mgetorbits(w->G_->b, w->G_->base_size, w->G_->gp, &w->G_->gens, w->G_->domain_size);
}

void auto_blaster::reset_skiplevels(auto_workspace* w) {
    w->skip_c.copy_force(w->start_c);
    w->skip_I = w->start_I;
    w->skiplevels = 0;
    //if(w->BW->BW.target_level > 0)
    //    w->skiplevels = w->BW->BW.target_level;
    w->skip_schreier_level = w->G_->gp;
    w->first_skiplevel = 1;
    w->my_base_points    = w->G_->b;
    w->my_base_points_sz = w->G_->base_size;
    w->is_foreign_base   = false;
}

void auto_blaster::sample_shared(sgraph* g_, bool master, shared_switches* switches, diy_group* G, coloring* start_c, bijection** canon_leaf, invariant** canon_I,
                                    com_pad* communicator_pad, int communicator_id, int** shared_orbit, int** shared_orbit_weights, bfs* bwork, mpermnode** gens, int* shared_group_size) {
    // find comparison leaf
    std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
    sgraph *g = g_;
    lowdeg* L;

    if(master) {
        L = new lowdeg();
        std::pair<sgraph*, coloring*> preprocessed_graph = L->preprocess(g);
        g = preprocessed_graph.first;
    }

    double cref;

    bool *done      = &switches->done;
    bool *done_fast = &switches->done_fast;
    int _shared_group_size = false;
    mpermnode *_gens       = nullptr;

    int* shrd_orbit;
    int* shrd_orbit_weights;

    std::vector<std::thread> work_threads;
    bijection base_points;
    bijection actual_base;
    int trash_int = 0;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() * ((communicator_id * 5) * 5135235);
    int selector_seed = seed;
    com_pad pad;

    auto_workspace W;
    invariant start_I;

    //W.first_level_fail.initialize(g->v_size);
    //W.first_level_succ.initialize(g->v_size);

    W.dequeue_space = new std::pair<int, int>[2048];
    W.dequeue_space_sz = 2048;
    W.enqueue_space = new std::pair<int, int>[128];
    W.enqueue_space_sz = 128; // <- choose this dynamic?

    if (master) {
        config.CONFIG_IR_DENSE = !(g->e_size < g->v_size || g->e_size / g->v_size < g->v_size / (g->e_size / g->v_size));
        std::cout << "[R] Dense graph: " << (config.CONFIG_IR_DENSE?"true":"false") << std::endl;

        shrd_orbit = new int[g->v_size];
        for(int i = 0; i < g->v_size; ++i)
            shrd_orbit[i] = i;

        shrd_orbit_weights = new int[g->v_size];
        memset(shrd_orbit_weights, 0, g->v_size * sizeof(int));

        int** shrd_orbit_ = new (int*);
        *shrd_orbit_ = shrd_orbit;

        int** shrd_orbit_weights_ = new (int*);
        *shrd_orbit_weights_ = shrd_orbit_weights;

        switches->current_mode = modes::MODE_TOURNAMENT;

        // first color refinement
        canon_I    = new invariant*;
        canon_leaf = new bijection*;
        start_c    = new coloring;
        g->initialize_coloring(start_c);
        W.start_c = start_c;
        start_I.create_vector();
        //W.R.old_refine_coloring_first(g, start_c, -1);
        //invariant trashi;
        //trashi.create_vector();
        W.R.old_refine_coloring_first(g, start_c, -1); // ToDo: why does new color ref not work properly on rantree-2000?
        //W.R.refine_coloring(g, start_c, nullptr, &trashi, -1, false);
        cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "[T] Color ref: " << cref / 1000000.0 << "ms" << std::endl;
        // ToDo: check if start_c now discrete

        // create some objects that are initialized after tournament
        G = new diy_group(g->v_size);
        W.BW = new bfs();
        bwork = W.BW;

        int init_c = W.S.select_color(g, start_c, selector_seed);
        if(init_c == -1) {
            *done = true;
            std::cout << "First coloring discrete." << std::endl;
            std::cout << "Base size: 0" << std::endl;
            std::cout << "Group size: 1" << std::endl;
            return;
        }
        W.S.empty_cache();
        // launch worker threads
        for (int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++)
            work_threads.emplace_back(
                    std::thread(&auto_blaster::sample_shared, auto_blaster(), g, false, switches, G, start_c,
                                canon_leaf, canon_I, &pad, i, shrd_orbit_, shrd_orbit_weights_, W.BW, &_gens, &_shared_group_size));
        cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "[T] Refinement workers created: " << cref / 1000000.0 << "ms" << std::endl;

        // set some workspace variables
        W.start_c = new coloring;
        W.start_c->copy_force(start_c);
        W.communicator_pad       = &pad;
        W.communicator_id        = -1;
        W.shared_orbit           = shrd_orbit_;
        W.shared_orbit_weights   = shrd_orbit_weights_;
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
        W.shared_orbit_weights = shared_orbit_weights;
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
        _canon_I->has_compare = false;
        _canon_I->compare_vec = nullptr;
        _canon_I->compareI    = nullptr;
        _canon_leaf = new bijection;
        {
            W.S.empty_cache();
            find_automorphism_prob(&W, g, false, _canon_I, _canon_leaf, &base_points, &trash_int, switches, selector_seed);
            W.my_base_points    = base_points.map;
            W.my_base_points_sz = base_points.map_sz;
            W.is_foreign_base   = true;
        }

        // skip immediately if base size is 1
        if(master && W.my_base_points_sz == 1) {
            std::cout << "[N] Base size 1 skip" << std::endl;
            *canon_I    = _canon_I;
            *canon_leaf = _canon_leaf;
            actual_base = base_points;
            base_points.not_deletable();
            G->initialize(g->v_size, &actual_base);
            bfs_element *root_elem = new bfs_element;
            root_elem->id = 0;
            root_elem->c = new coloring;
            root_elem->I = new invariant;
            root_elem->c->copy_force(start_c);
            root_elem->base_sz = 0;
            root_elem->is_identity = true;
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
            W.my_base_points    = W.G_->b;
            W.my_base_points_sz = W.G_->base_size;
            W.is_foreign_base   = false;
        }
    }

    int earliest_found = -1;
    int n_found = 0;
    int n_restarts = 0;
    bool foreign_base_done = false;
    bool reset_non_uniform_switch = true;

    // main loop...
    while(true) {
        if(switches->done_fast) {
            G->ack_done_shared();
            reset_non_uniform_switch = true;
        }

        if(switches->done) {
            //while (!switches->ack_done()) continue;
            return;
        }

        bijection automorphism;

        // in what phase are we in?
        switch(switches->current_mode) {
            case modes::MODE_TOURNAMENT:
                fast_automorphism_non_uniform(g, true, _canon_I, _canon_leaf, &automorphism, &restarts, done_fast, selector_seed, &W, switches->tolerance); // <- we should already safe unsuccessfull / succ first level stuff here
                if(n_found == 0) { // check if I won
                    // wait until everyone checked
                    while(!switches->check_leaf_tournament(communicator_id, restarts) && !switches->done_created_group) continue;
                    // check if I won, if yes: create group
                    if(switches->done_created_group) continue;
                    if(switches->win_id == communicator_id) {
                        *canon_I    = _canon_I;
                        *canon_leaf = _canon_leaf;
                        actual_base = base_points;
                        base_points.not_deletable();
                        G->initialize(g->v_size, &actual_base);
                        bfs_element *root_elem = new bfs_element;
                        root_elem->id = 0;
                        root_elem->c = new coloring;
                        root_elem->I = new invariant;
                        root_elem->c->copy_force(start_c);
                        root_elem->base_sz = 0;
                        root_elem->is_identity = true;
                        *root_elem->I = start_I;
                        W.S.empty_cache();
                        int init_c = W.S.select_color(g, start_c, selector_seed);
                        W.BW->initialize(root_elem, init_c, g->v_size, G->base_size);
                        switches->done_created_group = true;
                        int proposed_level = W.skiplevels + 1;
                        if(proposed_level == G->base_size) // ToDo: do this better...
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
                if(n_found % 3 == 0 && (W.skiplevels < W.my_base_points_sz - 1))
                    W.skiplevels += 1;
                if((*done_fast && !automorphism.non_uniform )) continue;
                break;

            case modes::MODE_NON_UNIFORM_PROBE:
                if(!foreign_base_done) {
                    fast_automorphism_non_uniform(g, true, _canon_I, _canon_leaf, &automorphism, &restarts, done_fast, selector_seed, &W, switches->tolerance);
                    automorphism.foreign_base = true;
                    n_restarts += restarts;
                } else {
                    fast_automorphism_non_uniform(g, true, *canon_I, *canon_leaf, &automorphism, &restarts, done_fast, selector_seed, &W, switches->tolerance);
                    automorphism.foreign_base = false;
                    n_restarts += restarts;
                }
                n_found += 1;
                if(n_found % 3 == 0 && (W.skiplevels < W.my_base_points_sz - 1))
                    W.skiplevels += 1;
                if((*done_fast && !automorphism.non_uniform )) continue;
                    break;

            case modes::MODE_NON_UNIFORM_PROBE_IT:
                // fast automorphism search, but from initial bfs pieces
                if(!*done_fast) {
                    if (reset_non_uniform_switch) {
                        reset_skiplevels(&W);
                        if(!master) { // guess a new leaf
                            //delete _canon_I->compare_vec;
                            _canon_I = new invariant;
                            _canon_I->create_vector();
                            //delete _canon_leaf;
                            _canon_leaf = new bijection;
                            base_points = bijection();
                            find_automorphism_prob(&W, g, false, _canon_I, _canon_leaf, &base_points, &trash_int, switches, selector_seed);
                            W.my_base_points    = base_points.map;
                            W.my_base_points_sz = base_points.map_sz;
                            W.is_foreign_base   = true;
                            foreign_base_done = false;
                        }
                        reset_non_uniform_switch = false;
                    }

                    if (!foreign_base_done) {
                        fast_automorphism_non_uniform(g, true, _canon_I, _canon_leaf, &automorphism, &restarts,
                                                      done_fast, selector_seed, &W, switches->tolerance);
                        automorphism.foreign_base = true;
                        n_restarts += restarts;
                    } else {
                        fast_automorphism_non_uniform(g, true, *canon_I, *canon_leaf, &automorphism, &restarts,
                                                      done_fast, selector_seed, &W, switches->tolerance);
                        automorphism.foreign_base = false;
                        n_restarts += restarts;
                    }

                    if (master && n_found == 0) {
                        int proposed_level = W.skiplevels + 1;
                        if (proposed_level == G->base_size)
                            proposed_level += 1;
                        if (proposed_level > W.BW->BW.target_level)
                            W.BW->BW.target_level.store(proposed_level);
                    }

                    n_found += 1;
                    if (n_found % (3 * switches->tolerance) == 0 && (W.skiplevels < W.my_base_points_sz - 1))
                        W.skiplevels += 1;
                    if ((*done_fast && !automorphism.non_uniform)) continue;
                } else continue;
                break;

            case modes::MODE_BFS:
                reset_non_uniform_switch = true;
                if(W.BW->BW.current_level != W.BW->BW.target_level) {
                    if (communicator_id == -1 && W.BW->BW.target_level < 0) {
                        int proposed_level = W.skiplevels + 1;
                        if (proposed_level == G->base_size)
                            proposed_level += 1;
                        W.BW->BW.target_level.store(proposed_level);
                    }
                    if(switches->done_shared_group && W.BW->BW.target_level >= 0) {
                        if(master && !switched1) {
                            switched1 = true;
                            cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
                            std::cout << "[N] Finished non-uniform automorphism search (" << *W.shared_generators_size << " generators, " << n_restarts << " restarts)" << std::endl;
                            std::cout << "[N] Ended in skiplevel " << W.skiplevels << ", found " << n_found << std::endl;
                            std::cout << "[T] " << cref / 1000000.0 << "ms" << std::endl;
                            std::cout << "[B] Determined target level: " << W.BW->BW.target_level << "" << std::endl;
                        }
                        bfs_chunk(g, *canon_I, *canon_leaf, done, selector_seed, &W);
                        if(master)
                            bwork->work_queues(switches->tolerance);
                    }
                } else {
                    if(master) {
                        bwork->work_queues(switches->tolerance);
                        if(bwork->BW.reached_initial_target) {
                            // reached the desired target level? go to next phase!
                            std::cout << "[A] Starting uniform probe, tolerance: " << switches->tolerance << std::endl;
                            switches->current_mode = modes::MODE_UNIFORM_PROBE;
                        } else {
                            // did not reach the target level within tolerance? iterate!
                            switches->iterate_tolerance();
                            switches->done_fast = false;
                            switches->done_shared_group = false;
                            G->non_uniform_abort_counter = 0;
                            n_found = 0;
                            switched1 = false;
                            bwork->reset_initial_target();
                            W.skiplevels = 0;
                            std::cout << "[A] Iterating, tolerance: " << switches->tolerance << std::endl;
                            // switches->reset_leaf_tournament();
                            switches->current_mode = modes::MODE_NON_UNIFORM_PROBE_IT; // ToDo: actually should go to leaf tournament
                        }
                    }
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

        if(switches->done) {
            //while (!switches->ack_done()) continue;
            return;
        }
        automorphism.not_deletable();
        //std::cout << automorphism.foreign_base << std::endl;
        bool test = true;
        if(switches->done_created_group) {
            test = G->add_permutation(&automorphism, &idle_ms, done);
            if(test && foreign_base_done) {
                G->sift_random();
            }
        }

        if(!test && !foreign_base_done) {
            // switch this worker to canonical search
            reset_skiplevels(&W);
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
                // ToDo: if current_level != 0, re-filter elements according to orbits
                // ToDo: need one pass of orbits to create orbit sizes for bfs_element weights (array orbit weights)
                std::cout << "[N] Waiting for ACK" << std::endl;
                switches->current_mode = modes::MODE_WAIT;

                G->ack_done_shared();
                reset_non_uniform_switch = true;
                G->wait_for_ack_done_shared(config.CONFIG_THREADS_REFINEMENT_WORKERS + 1);

                std::cout << "[N] Creating shared orbit and generators" << std::endl;
                G->sift_random();

                *W.shared_generators        = G->gens;

                memset(shrd_orbit_weights, 0, g->v_size * sizeof(int));

                int base_v_orbit = G->gp->orbits[G->b[0]];
                for(int i = 0; i < g->v_size; ++i)
                    if(G->gp->orbits[i] != base_v_orbit) {
                        (*W.shared_orbit)[i] = G->gp->orbits[i];
                        (*W.shared_orbit_weights)[G->gp->orbits[i]]++;
                    } else {
                        (*W.shared_orbit)[i] = G->b[0];
                        (*W.shared_orbit_weights)[G->b[0]]++;
                    }
                *W.shared_generators_size   = G->number_of_generators();

                    if(W.BW->BW.current_level > 1) {
                        bool tree_reduce = false;
                        if(*W.shared_generators_size > 0) {
                            std::cout << "[B] Reducing tree (" << n_found << ")" << std::endl;
                            assert(master && communicator_id == -1);
                            bfs_reduce_tree(&W);
                            tree_reduce = true;
                        }
                        // check if expected size is still too large...
                        if(W.BW->BW.level_expecting_finished[W.BW->BW.current_level] >= config.CONFIG_IR_SIZE_FACTOR * g->v_size * switches->tolerance) {
                            std::cout << "[B] Expected size still too large, not going into BFS" << std::endl;
                            W.BW->BW.reached_initial_target = (W.BW->BW.target_level == W.BW->BW.current_level);
                            W.BW->BW.target_level.store(W.BW->BW.current_level);
                        } else if(!tree_reduce) {
                            std::cout << "[B] No tree reduction, but filling queue..." << std::endl;
                            moodycamel::ConcurrentQueue<std::pair<bfs_element *, int>> throwaway_queue;
                            W.BW->BW.bfs_level_todo[W.BW->BW.current_level].swap(throwaway_queue);
                            int expected = 0;
                            for (int j = 0; j < W.BW->BW.level_sizes[W.BW->BW.current_level - 1]; ++j) {
                                bfs_element *elem = W.BW->BW.level_states[W.BW->BW.current_level - 1][j];
                                if (elem->weight > 0) {
                                    int c = elem->target_color;
                                    int c_size = elem->c->ptn[c] + 1;
                                    for (int i = c; i < c + c_size; ++i) {
                                        expected += 1;
                                        W.BW->BW.bfs_level_todo[W.BW->BW.current_level].enqueue(std::pair<bfs_element *, int>(elem, elem->c->lab[i]));
                                    }
                                }
                            }
                            W.BW->BW.level_expecting_finished[W.BW->BW.current_level] = expected;
                        }
                    }

                switches->done_shared_group = true;
                switches->current_mode      = modes::MODE_BFS;
            }
            if(switches->done) {
                std::cout << "Base size:  " << G->base_size << std::endl;

                //while(!switches->ack_done()) continue;

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