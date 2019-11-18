//
// Created by markus on 23.09.19.
//
//#define NDEBUG

#include "dejavu.h"
#include <stack>
#include <iostream>
#include <assert.h>
#include <chrono>
#include "schreier_sequential.h"
#include "refinement.h"
#include "selector.h"
#include "invariant.h"
#include "configuration.h"
#include "group_sequential.h"
#include "concurrentqueue.h"
#include "group_diy.h"
#include "lowdeg.h"
#include <pthread.h>
#include <tuple>

bool dejavu::proceed_state(dejavu_workspace* w, sgraph* g, coloring* c, invariant* I, int v, change_tracker* changes, strategy_metrics* m) {
    //std::cout << "proc" << std::endl;
    if(changes != nullptr)
        changes->track(c->vertex_to_col[v]);

    int init_color_class = w->R.individualize_vertex(c, v);
    bool comp = I->write_top_and_compare(INT32_MIN);
    comp && I->write_top_and_compare(INT32_MIN);

    comp = comp && I->write_top_and_compare(INT32_MAX);
    //if(!comp) return comp;
    comp = comp && w->R.refine_coloring(g, c, changes, I, init_color_class, changes != nullptr, m);
    comp = comp && I->write_top_and_compare(INT32_MAX);
    comp = comp && I->write_top_and_compare(INT32_MIN);
    //assert(tempcomp == comp);
    return comp;
}


abort_code dejavu::find_automorphism_from_bfs(dejavu_workspace *w, sgraph *g, bool compare, strategy* canon_strategy, bijection *automorphism, int *restarts,
                                        shared_switches *switches, int selector_seed) {
    bool backtrack = false;
    bool* done = &switches->done;

    refinement *R = &w->R;
    selector *S = &w->S;
    coloring *c = &w->c;
    invariant *I = &w->I;

    invariant* canon_I    = canon_strategy->I;
    bijection* canon_leaf = canon_strategy->leaf;

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

    mschreier* start_group_level = w->G->gp;
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
        if(*done) return abort_code();
        if(switches->current_mode != modes::MODE_UNIFORM_PROBE) return abort_code(2);
        if (backtrack) {
            *restarts += 1;
            if(w->communicator_id == -1) {
                if(*restarts > (switches->tolerance * 10)) {
                    return abort_code(1);
                }
                w->G->manage_results(switches); // ToDo: make this work
                if(*done) {
                    //std::cout << "early cancel" << std::endl;
                    return abort_code(2);
                }
            }

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
            //if(*restarts % 1000 == 0)
            //    std::cout << *restarts << ", " << reguess << ", " << level << std::endl;
        }

        int s;
        if (!backtrack) {
            s = S->select_color_dynamic(g, c, canon_strategy);
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
                    return abort_code();
                } else {
                    //I->push_level();
                    canon_leaf->read_from_coloring(c);
                    *canon_I = *I;
                    return abort_code();
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
            comp = comp && R->refine_coloring(g, c, nullptr, I, init_color_class, false, nullptr);
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


void dejavu::find_automorphism_prob(dejavu_workspace *w, sgraph *g, bool compare, invariant *canon_I,
                                    bijection *canon_leaf, strategy* canon_strategy, bijection *automorphism, int *restarts,
                                    shared_switches *switches, int selector_seed) {
    bool backtrack = false;
    bool* done = &switches->done;

    // workspace
    refinement *R = &w->R;
    selector *S = &w->S;
    coloring *c = &w->c;
    invariant *I = &w->I;

    coloring *start_c  = w->start_c;
    invariant *start_I = &w->start_I;

    S->empty_cache();

    int init_color_class;
    *restarts = 0;
    int level = 1;

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
    int s;

    while (true) {
        if(*done) return;
        if (!backtrack) {
            s = S->select_color_dynamic(g, c, canon_strategy);
            if (s == -1 && last_op == OP_R) {
                if (compare) {
                    // we can derive an automorphism!
                    w->measure2 += 1;
                    bijection leaf;
                    leaf.read_from_coloring(c);
                    leaf.not_deletable();
                    *automorphism = leaf;
                    automorphism->inverse();
                    automorphism->compose(canon_leaf);
                    assert(g->certify_automorphism(*automorphism));
                    return;
                } else {
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
            bool comp = I->write_top_and_compare(INT32_MAX);
            comp = comp && R->refine_coloring(g, c, nullptr, I, init_color_class, false, nullptr);
            comp = comp && I->write_top_and_compare(INT32_MAX);
            comp = comp && I->write_top_and_compare(INT32_MIN);
            last_op = OP_R;
            if (compare) {
                // compare invariant
                if(!comp) {
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
            int rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
            int v = c->lab[rpos];

            // individualize random vertex of class
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
            if(!comp) {
                backtrack = true;
                continue;
            }
        }
    }
}

void dejavu::fast_automorphism_non_uniform(sgraph* g, bool compare, strategy* canon_s,
                                           bijection* automorphism, strategy_metrics *m, bool *done, shared_switches* switches, int selector_seed, dejavu_workspace* w, int tolerance) {
    bool backtrack = false;
    bool skipped_level = false;

    // workspace
    refinement *R = &w->R;
    selector *S   = &w->S;
    coloring *c   = &w->c;
    invariant *I  = &w->I;

    coloring  *start_c     = &w->skip_c;
    invariant *start_I     = &w->skip_I;
    mschreier *group_level = w->skip_schreier_level;

    invariant* canon_I    = canon_s->I;
    bijection* canon_leaf = canon_s->leaf;

    automorphism->non_uniform = false;

    bool base_aligned = true;

    S->empty_cache();
    int init_color_class;

    m->restarts = 0;
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
            proceed_state(w, g, start_c, start_I, w->my_base_points[w->first_skiplevel - 1], nullptr, m);
            w->first_skiplevel += 1;
            if(!w->is_foreign_base)
                w->skip_schreier_level = w->skip_schreier_level->next;
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

    m->expected_bfs_size = 1;
    m->expected_level = -1;

    skipped_level = w->first_skiplevel > 1;

    while (true) {
        if(*done) return;
        if (backtrack) {
            if(*done) return;
            if((m->restarts % (5 * tolerance) == ((5 * tolerance) - 1)) && (w->skiplevels < w->my_base_points_sz))
                w->skiplevels += 1;

            S->empty_cache();
            if(w->first_skiplevel <= w->skiplevels) {
                proceed_state(w, g, start_c, start_I, w->my_base_points[w->first_skiplevel - 1], nullptr, m);
                w->first_skiplevel += 1;
                if(!w->is_foreign_base) {
                    w->skip_schreier_level = w->skip_schreier_level->next;
                }
            }

            // initialize a search state
            m->restarts += 1;
            if(switches->current_mode == modes::MODE_TOURNAMENT)
                switches->check_strategy_tournament(w->communicator_id, m, true);
            // ToDo: master should check results here and put done fast
            if(w->communicator_id == -1) // but need to be able to reach proper state afterwads
                w->G->manage_results(switches);

            c->copy_force(start_c);
            *I = *start_I;

            backtrack    = false;
            base_aligned = true;
            level = w->first_skiplevel;
            if(!w->is_foreign_base)
                group_level = w->skip_schreier_level;
            last_op = OP_R;

            skipped_level = w->first_skiplevel > 1;
        }

        int s;
        if (!backtrack) {
            s = S->select_color_dynamic(g, c, canon_s);
            if (s == -1 && last_op == OP_R) {
                if (compare) {
                    // we can derive an automorphism!
                    w->measure2 += 1;
                    bijection leaf;
                    leaf.read_from_coloring(c);
                    leaf.not_deletable();
                    *automorphism = leaf;
                    automorphism->inverse();
                    automorphism->compose(canon_leaf);
                    automorphism->non_uniform = skipped_level;
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
            comp = comp && R->refine_coloring(g, c, nullptr, I, init_color_class, false, m);
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
            int rpos;
            int v;

            if(level <= w->skiplevels) {
                skipped_level = true;
                v = w->my_base_points[level - 1];
                assert(c->vertex_to_col[v] == s);
                base_aligned  = true;
            } else {
                rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
                v = c->lab[rpos];
            }

            if(level > m->expected_level) {
                m->expected_bfs_size *= c->ptn[s] + 1;
                m->expected_level = level;
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
bool dejavu::get_orbit(dejavu_workspace* w, int* base, int base_sz, int v, int v_base, work_list* orbit, bool reuse_generators) {
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
        assert(v <= w->G->domain_size);
        assert(map_v >= 0);
        assert(map_v <= w->G->domain_size);

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
                    assert(b < w->G->domain_size && b >= 0);
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
        }

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

void dejavu::bfs_assure_init(dejavu_workspace* w) {
    if(!w->init_bfs) {
        int chunk_sz = w->BW->BW.chunk_size;
        w->todo_dequeue = new std::tuple<bfs_element*, int, int>[chunk_sz];
        w->todo_deque_sz = chunk_sz;
        w->todo_elements = new std::pair<bfs_element *, int>[chunk_sz * 8];
        w->todo_elements_sz = chunk_sz * 8;
        w->finished_elements = new std::pair<bfs_element *, int>[chunk_sz + 1];
        w->finished_elements_sz = chunk_sz;
        w->init_bfs = true;
        w->orbit.initialize(w->G->domain_size);
        w->orbit_considered.initialize(w->G->domain_size);
        w->orbit_vertex_worklist.initialize(w->G->domain_size);
        w->generator_fix_base = new mpermnode*[*w->shared_generators_size];
        w->generator_fix_base_alloc = *w->shared_generators_size;
    }
}

bool dejavu::bfs_chunk(sgraph* g, strategy* canon_strategy, bool *done, int selector_seed, dejavu_workspace* w) {
    thread_local bool done_test = false;

    bfs* BFS = w->BW;
    int level        = BFS->BW.current_level;
    int target_level = BFS->BW.target_level;
    if(level == target_level) return false; // we are done with BFS!

    // initialize bfs structures
    bfs_assure_init(w);

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
        bfs_element *elem = std::get<0>(w->todo_dequeue[i]);
        int v             = std::get<1>(w->todo_dequeue[i]);
        int weight        = std::get<2>(w->todo_dequeue[i]);
        bool is_identity  = elem->is_identity && (v == w->my_base_points[elem->base_sz]);

        // check orbit
        bool comp = elem->weight > 0;
        if(weight == -1) {
            comp = comp && get_orbit(w, elem->base, elem->base_sz, v, w->my_base_points[elem->base_sz], &w->orbit,
                                  w->prev_bfs_element == elem);
            assert(!comp ? (!is_identity) : true);
            // I could prune the todos relatively cheaply now (just a further restriction of generators, and they will not
            // appear in queues at all)
        }
        if(comp) {
            // copy to workspace
            if (w->prev_bfs_element != elem) { // <-> last computed base is the same!
                w->work_c->copy_force(elem->c);
                w->prev_bfs_element = elem;
                *w->work_I = *elem->I;
                w->work_I->set_compare_invariant(canon_strategy->I);
            } else {
                *w->work_I = *elem->I;
                w->work_I->set_compare_invariant(canon_strategy->I);
                w->work_c->copy(elem->c);
            }

            // compute next coloring
            comp = comp && proceed_state(w, g, w->work_c, w->work_I, v, nullptr, nullptr); // &w->changes

            // ToDo: if !comp consider abort map
            if(comp && elem->is_identity && level > 1) {
                // decrease abort map done...
                BFS->BW.level_abort_map_mutex[level]->lock();
                BFS->BW.level_abort_map_done[level]--;
                BFS->BW.level_abort_map_mutex[level]->unlock();
            }
            if(!comp && level > 1) {
                //std::cout << w->work_I->comp_fail_pos << ", " << w->work_I->comp_fail_val << std::endl;
                if(elem->is_identity) { // save to abort map...
                    BFS->write_abort_map(level, w->work_I->comp_fail_pos, w->work_I->comp_fail_val);
                } else { // if abort map done, check abort map...
                    bool comp_ = BFS->read_abort_map(level, w->work_I->comp_fail_pos, w->work_I->comp_fail_val);
                    if(!comp_) {
                        elem->weight = 0;
                    }
                }
            }
        }

        assert(elem->base_sz < w->G->base_size);
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
        // ToDo: use memcpy
        for(int j = 0; j < elem->base_sz; ++j) {
            assert(elem->base[j] >= 0 && elem->base[j] < g->v_size);
            next_elem->base[j] = elem->base[j];
        }
        assert(v >= 0 && v < g->v_size);
        next_elem->base[next_elem->base_sz - 1] = v;
        if(weight == -1)
            next_elem->weight = elem->weight * w->orbit.cur_pos;
        else
            next_elem->weight = weight;
        next_elem->parent_weight = elem->weight;
        next_elem->parent = elem;

        // create todos for this node
        int sz = 0;
        w->S.empty_cache();
        int c = w->S.select_color_dynamic(g, w->work_c, canon_strategy);
        next_elem->target_color = c;
        sz += w->work_c->ptn[c] + 1;

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

void dejavu::sequential_init_copy(dejavu_workspace* w) {
    if(!w->sequential_init) {
        _newgroup(&w->sequential_gp, &w->sequential_gens, w->G->domain_size);
        w->sequential_init = true;
    }

    // copy generators that have not been copied yet
    mpermnode *it = w->G->gens;
    do {
        if(it->copied == 0) {
            _addpermutation(&w->sequential_gens, it->p, w->G->domain_size);
            it->copied = 1;
        }
        it = it->next;
    } while (it != w->G->gens);
}

void dejavu::bfs_reduce_tree(dejavu_workspace* w) {
    thread_local bool init_group = false;

    bfs_assure_init(w);

    int domain_size = w->G->domain_size;
    sequential_init_copy(w);

    if(w->generator_fix_base_alloc < *w->shared_generators_size) {
        delete[] w->generator_fix_base;
        w->generator_fix_base = new mpermnode*[*w->shared_generators_size];
        w->generator_fix_base_alloc = *w->shared_generators_size;
        w->prev_bfs_element = nullptr;
    }

    // step1: check on level i if orbits can be reduced
    int i, j, v;
    bfs_workspace* BW = &w->BW->BW;
    int active = 0;
    bool found_canon = false;
    for (int j = 0; j < BW->level_sizes[BW->current_level - 1]; ++j) {
        bfs_element *elem = BW->level_states[BW->current_level - 1][j];
        active += (elem->weight > 0);
        found_canon = found_canon || (elem->base[elem->base_sz - 1] == w->G->b[elem->base_sz - 1]);
    }
    assert(found_canon);
    std::cout << "[B] Top level active (before): " << active << std::endl;


    for(i = 1; i < BW->current_level; ++i) {
        //std::cout << "checking level " << i << ", sz " << BW->level_sizes[i] << std::endl;
        if(i == 1) {
            for(j = 0; j < BW->level_sizes[i]; ++j) {
                bfs_element* elem = BW->level_states[i][j];
                if(elem->weight > 0) {
                    v = elem->base[elem->base_sz - 1];
                    assert(v >= 0 && v < w->G->domain_size);
                    if ((*w->shared_orbit)[v] != v) {
                        elem->weight = 0;
                        assert(v != w->G->b[elem->base_sz - 1]);
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
                    int* orbits_sz = nullptr;
                    // ToDo: find all other elements with same base -1 and prune using exactly this orbit, then its fine
                    int* orbits = _getorbits(elem->base, elem->base_sz - 1, w->sequential_gp, &w->sequential_gens, domain_size, w->G->b, &orbits_sz); // ToDo: weights incorrect?

                    int calc_sz = 0;
                    for(int ii = 0; ii < domain_size; ++ii) {
                        assert(orbits[ii] >= 0 && orbits[ii] < domain_size);
                        if(ii == orbits[ii]) {
                            //std::cout << (orbits_sz)[ii] << ", " << ii << std::endl;
                            calc_sz += (orbits_sz)[ii];
                        }
                    }
                    assert(calc_sz == domain_size);
                    assert(elem->base_sz == elem->level);
                    assert(elem->level == i);
                    assert(i > 1);
                    v = elem->base[elem->base_sz - 1];
                    assert((v >= 0) && (v < w->G->domain_size));
                    assert(elem->base_sz - 1 > 0);
                    if(orbits[v] != v) {
                        assert(v != w->G->b[elem->base_sz - 1]);
                        elem->weight = 0;
                    } else {
                        elem->weight = elem->parent_weight * orbits_sz[v];
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

    std::cout << "[B] Top level active (after): " << active << std::endl;
    std::cout << "[B] Emptying queue..." << std::endl;
    moodycamel::ConcurrentQueue<std::tuple<bfs_element *, int, int>> throwaway_queue;
    BW->bfs_level_todo[BW->current_level].swap(throwaway_queue);


    // do not fill this if there is no hope...
    int expecting_finished = 0;
    for (int j = 0; j < BW->level_sizes[BW->current_level - 1]; ++j) {
        bfs_element *elem = BW->level_states[BW->current_level - 1][j];
        int added = 0;
        if (elem->weight > 0) {
            int c      = elem->target_color;
            int c_size = elem->c->ptn[c] + 1;
            expecting_finished += c_size;
        }
    }
    BW->level_expecting_finished[BW->current_level] = expecting_finished;

    //std::cout << " to " << BW->bfs_level_todo[BW->current_level].size_approx() << std::endl;
    w->prev_bfs_element = nullptr;
    w->orbit.reset();
    w->orbit_vertex_worklist.reset();
    w->orbit_considered.reset();
}

void dejavu::reset_skiplevels(dejavu_workspace* w) {
    w->skip_c.copy_force(w->start_c);
    w->skip_I = w->start_I;
    w->skiplevels = 0;
    w->skip_schreier_level = w->G->gp;
    w->first_skiplevel = 1;
    w->my_base_points    = w->G->b;
    w->my_base_points_sz = w->G->base_size;
    w->is_foreign_base   = false;
}

void dejavu::bfs_fill_queue(dejavu_workspace* w) {
    if(w->BW->BW.current_level == w->BW->BW.target_level)
        return;
    moodycamel::ConcurrentQueue<std::tuple<bfs_element *, int, int>> throwaway_queue;
    w->BW->BW.bfs_level_todo[w->BW->BW.current_level].swap(throwaway_queue);

    int expected = 0;

    if(*w->shared_generators_size > 0)
        sequential_init_copy(w);

    if(!w->sequential_init) {
        for (int j = 0; j < w->BW->BW.level_sizes[w->BW->BW.current_level - 1]; ++j) {
            bfs_element *elem = w->BW->BW.level_states[w->BW->BW.current_level - 1][j];
            if (elem->weight > 0) {
                int c = elem->target_color;
                int c_size = elem->c->ptn[c] + 1;
                for (int i = c; i < c + c_size; ++i) {
                    expected += 1;
                    w->BW->BW.bfs_level_todo[w->BW->BW.current_level].enqueue(
                            std::tuple<bfs_element *, int, int>(elem, elem->c->lab[i], -1)); // ToDo: prune this, too!
                }
                if (elem->is_identity) {
                    std::cout << "Abort map expecting: " << c_size << std::endl;
                    w->BW->BW.level_abort_map_done[w->BW->BW.current_level] = c_size;
                }
            }
        }
    } else {
        std::cout << "Filling with orbits..." << std::endl;
        int i;
        for (int j = 0; j < w->BW->BW.level_sizes[w->BW->BW.current_level - 1]; ++j) {
            bfs_element *elem = w->BW->BW.level_states[w->BW->BW.current_level - 1][j];
            int added = 0;
            if (elem->weight > 0) {
                int c = elem->target_color;
                int c_size = elem->c->ptn[c] + 1;
                int *orbits_sz = nullptr;
                int *orbits = _getorbits(elem->base, elem->base_sz, w->sequential_gp, &w->sequential_gens,
                                         w->G->domain_size, w->G->b, &orbits_sz);
                for (i = c; i < c + c_size; ++i) {
                    if (orbits[elem->c->lab[i]] == elem->c->lab[i]) {
                        w->BW->BW.bfs_level_todo[w->BW->BW.current_level].enqueue(
                                std::tuple<bfs_element *, int, int>(elem, elem->c->lab[i], orbits_sz[i]));
                        expected += 1;
                        added += 1;
                    }
                }
            }
            if (elem->is_identity) {
                w->BW->BW.level_abort_map_done[w->BW->BW.current_level] = added;
            }
        }
    }

    w->BW->BW.level_expecting_finished[w->BW->BW.current_level] = expected;
}

void dejavu::sample_shared(sgraph* g_, bool master, shared_switches* switches, group_diy* G, coloring* start_c, strategy* canon_strategy,
                           int communicator_id, int** shared_orbit, int** shared_orbit_weights, bfs* bwork, mpermnode** gens, int* shared_group_size) {
    // find comparison leaf
    std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
    sgraph *g = g_;
    lowdeg* L;
    dejavu_workspace W;

    if(master) {
        config.CONFIG_IR_DENSE = !(g->e_size < g->v_size || g->e_size / g->v_size < g->v_size / (g->e_size / g->v_size));
        L = new lowdeg();
        start_c = new coloring;
        g->initialize_coloring(start_c);
        std::pair<sgraph*, coloring*> preprocessed_graph = L->preprocess(start_c, g, &W.R);
        g       = preprocessed_graph.first;
        start_c = preprocessed_graph.second;
        assert(start_c->check());
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

    invariant start_I;

    if (master) {
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
        canon_strategy = new strategy;

        W.start_c = start_c;
        start_I.create_vector();
        W.R.refine_coloring_first(g, start_c, -1); // ToDo: why does new color ref not work properly on rantree-2000?
        cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "[T] Color ref: " << cref / 1000000.0 << "ms" << std::endl;
        // ToDo: check if start_c now discrete

        // create some objects that are initialized after tournament
        G = new group_diy(g->v_size);
        W.BW = new bfs();
        bwork = W.BW;

        int init_c = W.S.select_color(g, start_c, selector_seed);
        if(init_c == -1) {
            *done = true;
            std::cout << "First coloring discrete." << std::endl;
            std::cout << "Base size: 0" << std::endl;
            std::cout << "Group size: 1" << std::endl;
            L->postprocess(nullptr);
            return;
        }
        W.S.empty_cache();
        // launch worker threads
        for (int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++)
            work_threads.emplace_back(
                    std::thread(&dejavu::sample_shared, dejavu(), g, false, switches, G, start_c,
                                canon_strategy, i, shrd_orbit_, shrd_orbit_weights_, W.BW, &_gens, &_shared_group_size));
        cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "[T] Refinement workers created: " << cref / 1000000.0 << "ms" << std::endl;

        // set some workspace variables
        W.start_c = new coloring;
        W.start_c->copy_force(start_c);
        W.communicator_id        = -1;
        W.shared_orbit           = shrd_orbit_;
        W.shared_orbit_weights   = shrd_orbit_weights_;
        W.shared_generators      = &_gens;
        W.shared_generators_size = &_shared_group_size;
    }

    int sampled_paths = 0;
    int restarts = 0;
    int idle_ms = 0;
    invariant *my_canon_I;
    bijection *my_canon_leaf;
    bool switched1 = false;
    bool switched2 = false;

    W.skip_c.copy_force(start_c);
    W.work_c = new coloring;
    W.work_I = new invariant;
    W.G = G;
    W.BW = bwork;
    W.skiplevels = 0;
    W.skip_schreier_level = G->gp;

    if (!master) {
        W.shared_orbit = shared_orbit;
        W.shared_orbit_weights = shared_orbit_weights;
        W.start_c = new coloring;
        W.start_c->copy_force(start_c);
        W.communicator_id = communicator_id;
        W.shared_generators = gens;
        W.shared_generators_size = shared_group_size;
    }

    strategy* my_strategy;
    {
        //W.base_size = G->base_size;
        my_canon_I = new invariant;
        my_canon_I->create_vector();
        my_canon_I->has_compare = false;
        my_canon_I->compare_vec = nullptr;
        my_canon_I->compareI    = nullptr;
        my_canon_leaf = new bijection;
        selector_type rst = (selector_type) intRand(0, 2, selector_seed);
        my_strategy = new strategy(my_canon_leaf, my_canon_I, rst, -1);

        {
            W.S.empty_cache();
            find_automorphism_prob(&W, g, false, my_canon_I, my_canon_leaf, my_strategy, &base_points, &trash_int, switches, selector_seed);
            W.my_base_points    = base_points.map;
            W.my_base_points_sz = base_points.map_sz;
            W.is_foreign_base   = true;
        }

        // skip immediately if base size is 1
        if(master && W.my_base_points_sz == 1 && !config.CONFIG_IR_FULLBFS) { // ToDo: can end in deadlock if threads dont aggree on base size 1
            std::cout << "[N] Base size 1 skip" << std::endl;
            //*canon_I    = my_canon_I;
            //*canon_leaf = my_canon_leaf;
            canon_strategy->replace(my_strategy);
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
            int init_c = W.S.select_color_dynamic(g, start_c, my_strategy);
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
        } else if(W.my_base_points_sz == 1 && !config.CONFIG_IR_FULLBFS) {
            // wait until shared group created
            while(!switches->done_shared_group) continue;
            while(!switches->current_mode == modes::MODE_BFS || !switches->current_mode == modes::MODE_UNIFORM_PROBE) continue;
            W.my_base_points    = W.G->b;
            W.my_base_points_sz = W.G->base_size;
            W.is_foreign_base   = false;
        }
    }

    int earliest_found = -1;
    int n_found = 0;
    int n_restarts = 0;
    strategy_metrics m;
    bool foreign_base_done = false;
    bool reset_non_uniform_switch = true;
    int required_level = -1;

    // main loop...
    while(true) {
        if(master) {
            // sifting results
            G->manage_results(switches);

            // non-uniform search over, fix a group state for collaborative bfs
            if(switches->done_fast && !switches->done_shared_group) {
                // wait for ack of done_fast
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
                    } else {
                        std::cout << "[B] Filling queue..." << std::endl;
                        bfs_fill_queue(&W);
                    }
                } else {
                    /*if(*W.shared_generators_size > 0) { // ToDo: deletes identity?
                        std::cout << "[B] Pre-re-fill" << std::endl;
                        bfs_fill_queue(&W);
                    }*/
                }

                switches->done_shared_group = true;
                switches->current_mode      = modes::MODE_BFS;
            }

            // search is done
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
                L->postprocess(G);
                break;
            }
        }

        if(switches->done_fast) {
            G->ack_done_shared();
            reset_non_uniform_switch = true;
        }

        if(switches->done) {
            //while (!switches->ack_done()) continue;
            return;
        }

        bijection automorphism;
        abort_code A;
        automorphism.mark = false;

        // in what phase are we in?
        switch(switches->current_mode) {
            case modes::MODE_TOURNAMENT:
                m.restarts = 0;
                m.expected_bfs_size = 0;
                fast_automorphism_non_uniform(g, true, my_strategy, &automorphism, &m, done_fast, switches, selector_seed, &W, switches->tolerance); // <- we should already safe unsuccessfull / succ first level stuff here
                if(n_found == 0) { // check if I won
                    // wait until everyone checked
                    while(!switches->check_strategy_tournament(communicator_id, &m, false) && !switches->done_created_group) continue;
                    // check if I won, if yes: create group
                    if(switches->done_created_group) continue;
                    if(switches->win_id == communicator_id) {
                        //*canon_I    = my_canon_I;
                        //*canon_leaf = my_canon_leaf;
                        canon_strategy->replace(my_strategy);
                        actual_base = base_points;
                        base_points.not_deletable();
                        G->initialize(g->v_size, &actual_base);
                        std::cout << "Strategy: " << canon_strategy->cell_selector_type << std::endl;
                        std::cout << "Base size: " << actual_base.map_sz << std::endl;
                        bfs_element *root_elem = new bfs_element;
                        root_elem->id = 0;
                        root_elem->c = new coloring;
                        root_elem->I = new invariant;
                        root_elem->c->copy_force(start_c);
                        root_elem->base_sz = 0;
                        root_elem->is_identity = true;
                        *root_elem->I = start_I;
                        W.S.empty_cache();
                        int init_c = W.S.select_color_dynamic(g, start_c, my_strategy);
                        W.BW->initialize(root_elem, init_c, g->v_size, G->base_size);
                        switches->done_created_group = true;
                        int proposed_level = W.skiplevels + 1;
                        //if(proposed_level == G->base_size) // ToDo: do this better...
                        //    proposed_level += 1;
                        if(config.CONFIG_IR_FULLBFS)
                            proposed_level = G->base_size + 1;
                        W.BW->BW.target_level.store(proposed_level);
                        std::cout << "proposed level: " << proposed_level << std::endl;
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
                n_restarts += m.restarts;
                n_found += 1;
                if(n_found % 3 == 0 && (W.skiplevels < W.my_base_points_sz - 1))
                    W.skiplevels += 1;
                if((*done_fast && !automorphism.non_uniform )) continue;
                break;

            case modes::MODE_NON_UNIFORM_PROBE:
                if(!foreign_base_done) {
                    fast_automorphism_non_uniform(g, true, my_strategy, &automorphism, &m, done_fast, switches, selector_seed, &W, switches->tolerance);
                    automorphism.foreign_base = true;
                    n_restarts += m.restarts;
                    automorphism.mark = true;
                } else {
                    fast_automorphism_non_uniform(g, true, canon_strategy, &automorphism, &m, done_fast, switches, selector_seed, &W, switches->tolerance);
                    automorphism.foreign_base = false;
                    n_restarts += m.restarts;
                    automorphism.mark = true;
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
                            //delete my_canon_I->compare_vec;
                            my_canon_I = new invariant; // delete old if non canon
                            my_canon_I->create_vector();
                            //delete my_canon_leaf;
                            my_canon_leaf = new bijection;
                            base_points = bijection();
                            selector_type rst = (selector_type) intRand(0, 2, selector_seed);
                            my_strategy = new strategy(my_canon_leaf, my_canon_I, rst, -1);
                            find_automorphism_prob(&W, g, false, my_canon_I, my_canon_leaf, my_strategy, &base_points, &trash_int, switches, selector_seed);
                            W.my_base_points    = base_points.map;
                            W.my_base_points_sz = base_points.map_sz;
                            W.is_foreign_base   = true;
                            foreign_base_done = false;
                        }
                        reset_non_uniform_switch = false;
                    }

                    if (!foreign_base_done) {
                        fast_automorphism_non_uniform(g, true, my_strategy, &automorphism, &m,
                                                      done_fast, switches, selector_seed, &W, switches->tolerance);
                        automorphism.foreign_base = true;
                        automorphism.mark = true;
                        n_restarts += m.restarts;
                    } else {
                        fast_automorphism_non_uniform(g, true, canon_strategy, &automorphism, &m,
                                                      done_fast, switches, selector_seed, &W, switches->tolerance);
                        automorphism.foreign_base = false;
                        automorphism.mark = true;
                        n_restarts += m.restarts;
                    }

                    if (master && n_found == 0) {
                        int proposed_level = std::max(W.skiplevels + 1, required_level);
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
                if ((*done_fast && !automorphism.non_uniform)) continue;
                break;

            case modes::MODE_BFS:
                //std::cout << "b" << std::endl;
                reset_non_uniform_switch = true;
                if(W.is_foreign_base) {
                    reset_skiplevels(&W);
                    foreign_base_done = true;
                }
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
                        bfs_chunk(g, canon_strategy, done, selector_seed, &W);
                        if(master) {
                            bool fill = bwork->work_queues(switches->tolerance); // ToDo: maybe current level should be incremented later for sync...
                            if(fill)
                                bfs_fill_queue(&W);
                        }
                    }
                } else {
                    if(master) {
                        bool fill = bwork->work_queues(switches->tolerance);
                        if(fill)
                            bfs_fill_queue(&W);
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
                            std::cout << "[A] Iterating, tolerance: " << switches->tolerance << std::endl;
                            // switches->reset_leaf_tournament();
                            reset_skiplevels(&W);
                            foreign_base_done = true;
                            switches->current_mode = modes::MODE_NON_UNIFORM_PROBE_IT; // ToDo: actually should go to leaf tournament
                            continue;
                        }
                    }
                }
                continue;
                break;

            case modes::MODE_UNIFORM_PROBE:
                reset_non_uniform_switch = true;
                if(W.communicator_id == 0 && !switched2) {
                    switched2 = true;
                    cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() - timer).count());
                    std::cout << "[T] " << cref / 1000000.0 << "ms" << std::endl;
                }
                A = find_automorphism_from_bfs(&W, g, true, canon_strategy, &automorphism, &restarts, switches, selector_seed);
                if(A.reason == 2)
                    continue;
                if(A.reason == 1) {
                    // go back to bfs?
                    switches->current_mode = MODE_WAIT;
                    // manage?
                    switches->iterate_tolerance();
                    switches->done_fast = false;
                    switches->done_shared_group = false;
                    G->non_uniform_abort_counter = 0;
                    n_found = 0;
                    switched1 = false;
                    bwork->reset_initial_target();
                    std::cout << "[A] Iterating, tolerance: " << switches->tolerance << std::endl;
                    reset_skiplevels(&W);
                    foreign_base_done = true;
                    switches->current_mode = MODE_NON_UNIFORM_PROBE_IT;
                    required_level = W.BW->BW.current_level + 1;
                    std::cout << "[B] Requiring level " << required_level << std::endl;
                    continue;
                }
                automorphism.mark = true;
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
        if(switches->done_created_group && automorphism.mark) {
           // std::cout << "found auto" << std::endl;
            test = G->add_permutation(&automorphism, &idle_ms, done);
            if(test && foreign_base_done) {
                G->sift_random();
            }
        }

        automorphism.mark = false;

        if(!test && !foreign_base_done) {
            // switch this worker to canonical search
            reset_skiplevels(&W);
            foreign_base_done = true;
            std::cout << "[N] Switching to canonical search (" << W.communicator_id << ", " << n_found << " generators)" << std::endl;
        }

        delete[] automorphism.map;
        automorphism.map = new int[g->v_size];
        automorphism.foreign_base = false;
        sampled_paths += 1;

        // master thread managing sifting results, bfs, ...
    }

    if(master)
        delete G;
    //delete W.start_c;
    return;
}