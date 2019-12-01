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
#include "concurrentqueue.h"
#include "group_shared.h"
#include "lowdeg.h"
#include <pthread.h>
#include <tuple>

// individualize and refine
bool dejavu::proceed_state(dejavu_workspace* w, sgraph* g, coloring* c, invariant* I, int v, change_tracker* changes, strategy_metrics* m) {
    if(changes != nullptr)
        changes->track(c->vertex_to_col[v]);

    int init_color_class = w->R.individualize_vertex(c, v);
    bool comp = I->write_top_and_compare(INT32_MIN);
    comp && I->write_top_and_compare(INT32_MIN);
    comp = comp && I->write_top_and_compare(INT32_MAX);

    comp = comp && w->R.refine_coloring(g, c, changes, I, init_color_class, changes != nullptr, m);
    comp = comp && I->write_top_and_compare(INT32_MAX);
    comp = comp && I->write_top_and_compare(INT32_MIN);
    return comp;
}

// search uniformly from bfs structure, no additional leaf storage
abort_code dejavu::uniform_from_bfs_search(dejavu_workspace *w, sgraph *g, bool compare, strategy* canon_strategy, bijection *automorphism, int *restarts,
                                           shared_workspace *switches, int selector_seed) {
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


    automorphism->certified = false;

    // pick start from BFS level
    int bfs_level    = w->BW->current_level - 1;
    int bfs_level_sz = w->BW->level_sizes[bfs_level];

    mschreier* start_group_level = w->G->gp;
    for(int i = 0; i < bfs_level; i++)
        start_group_level = start_group_level->next;
    mschreier* group_level = start_group_level;

    int rand_pos     = intRand(0, bfs_level_sz - 1, selector_seed);
    bfs_element* picked_elem = w->BW->level_states[bfs_level][rand_pos];
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

            // make some global checks
            *restarts += 1;
            if(w->id == -1) {
                // too many restarts? abort and try bfs_workspace again...
                if(*restarts > (switches->tolerance * 10)) {
                    return abort_code(1);
                }

                // manage sifting results too detect if other threads finished the task
                w->G->manage_results(switches);
                if(*done) {
                    return abort_code(2);
                }
            }

            // do uniform search
            bfs_level    = w->BW->current_level - 1;
            bfs_level_sz = w->BW->level_sizes[bfs_level];
            rand_pos     = intRand(0, bfs_level_sz - 1, selector_seed);
            picked_elem = w->BW->level_states[bfs_level][rand_pos];
            group_level = start_group_level;
            base_aligned = picked_elem->is_identity;

            // consider the weight by redrawing
            picked_weight = picked_elem->weight;
            max_weight    = w->BW->level_maxweight[bfs_level];
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
        }

        int s;
        if (!backtrack) {
            s = S->select_color_dynamic(g, c, canon_strategy);
            if (s == -1) {
                if (compare) {
                    // we can derive an automorphism!
                    bijection leaf;
                    leaf.read_from_coloring(c);
                    leaf.not_deletable();
                    *automorphism = leaf;
                    automorphism->inverse();
                    automorphism->compose(canon_leaf);//enqueue_fail_point_sz
                    //std::cout << "Found automorphism." << *restarts << std::endl;
                    if(!config.CONFIG_IR_FULL_INVARIANT && !R->certify_automorphism(g, automorphism)) {
                        std::cout << "late backtrack2" << std::endl;
                        backtrack = true;
                        continue;
                    }
                    automorphism->certified = true;
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


void dejavu::find_first_leaf(dejavu_workspace *w, sgraph *g, bool compare, invariant *canon_I,
                             bijection *canon_leaf, strategy* canon_strategy, bijection *automorphism, int *restarts,
                             shared_workspace *switches, int selector_seed) {
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
        start_I->create_vector(g->v_size * 2);
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

void
dejavu::base_aligned_search(dejavu_workspace *w, sgraph *g, strategy *canon_strategy, bijection *automorphism,
                            strategy_metrics *m, bool *done, shared_workspace *switches, int selector_seed) {
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

    invariant* canon_I    = canon_strategy->I;
    bijection* canon_leaf = canon_strategy->leaf;

    automorphism->non_uniform = false;
    automorphism->certified   = false;

    bool base_aligned = true;

    S->empty_cache();

    m->restarts = 0;
    int level = w->first_skiplevel;

    start_I->set_compare_invariant(canon_I);
    *I = *start_I;
    c->copy_force(start_c);

    if(w->skiplevels >= w->my_base_points_sz - 1)
        w->skiplevels = w->my_base_points_sz - 2;

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

    m->expected_bfs_size = 1;
    m->expected_level = -1;

    skipped_level = w->first_skiplevel > 1;

    int it = 0;

    while (true) {
        if(*done) return;

        ++it;
        if(it % 3 == 0) {
            if(switches->current_mode == modes::MODE_TOURNAMENT)
                switches->check_strategy_tournament(w->id, m, true);
            if(w->id == -1) // but need to be able to reach proper state afterwads
                w->G->manage_results(switches);
        }

        if (backtrack) {
            if(*done) return;
            if((m->restarts % (5 * switches->tolerance) == ((5 * switches->tolerance) - 1)) && (w->skiplevels < w->my_base_points_sz))
                w->skiplevels += 1;

            S->empty_cache();
            if(w->first_skiplevel <= w->skiplevels) {
                m->expected_bfs_size *= start_c->ptn[start_c->vertex_to_col[w->my_base_points[w->first_skiplevel - 1]]] + 1;
                proceed_state(w, g, start_c, start_I, w->my_base_points[w->first_skiplevel - 1], nullptr, m);
                w->first_skiplevel += 1;
                if(!w->is_foreign_base) {
                    w->skip_schreier_level = w->skip_schreier_level->next;
                }
            }

            // initialize a search state
            m->restarts += 1;

            c->copy_force(start_c);
            *I = *start_I;

            backtrack    = false;
            base_aligned = true;
            level = w->first_skiplevel;
            if(!w->is_foreign_base)
                group_level = w->skip_schreier_level;

            skipped_level = w->first_skiplevel > 1;
        }

        int s = S->select_color_dynamic(g, c, canon_strategy);
        if (s == -1) {
            // we can derive an automorphism!
            bijection leaf;
            leaf.read_from_coloring(c);
            leaf.not_deletable();
            *automorphism = leaf;
            automorphism->inverse();
            automorphism->compose(canon_leaf);
            automorphism->non_uniform = skipped_level;
            if(!config.CONFIG_IR_FULL_INVARIANT && !R->certify_automorphism(g, automorphism)) {
                //leaf.deletable();
                std::cout << "late backtrack1" << std::endl;
                backtrack = true;
                continue;
            }
            automorphism->certified = true;
            assert(g->certify_automorphism(*automorphism));
            return;
        }

        // individualize and refine now
        int rpos, v;

        if(level <= w->skiplevels) {
            skipped_level = true;
            v = w->my_base_points[level - 1];
            assert(c->vertex_to_col[v] == s);
            base_aligned  = true;
        } else {
            rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
            v = c->lab[rpos];
        }

       // if(level > m->expected_level) {
     //       m->expected_bfs_size *= c->ptn[s] + 1;
      //      m->expected_level = level;
      //  }

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

        bool comp = proceed_state(w, g, c, I, v, nullptr, m);

        level += 1;

        if(!comp) {
            backtrack = true;
            continue;
        }

        if(!w->is_foreign_base)
            group_level = group_level->next;
    }
}

bool dejavu::uniform_from_bfs_search_with_storage(dejavu_workspace* w, sgraph* g, shared_workspace* switches, bfs_element* elem, int selector_seed, strategy* strat, bijection* automorphism, bool look_close) {
    coloring*  c = &w->c;
    invariant* I = &w->I;

    c->copy_force(elem->c);
    *I = *elem->I;

    I->set_compare_invariant(strat->I);

    int col, v, rpos;
    bool comp;

    // first individualization
    col  = elem->target_color;
    rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[col] + 1));
    v    = c->lab[rpos];

    I->reset_deviation();
    if(look_close)
        I->never_fail = true;
    //I->never_fail = true;
    comp = proceed_state(w, g, c, I, v, nullptr, nullptr);

    if(!comp) { // fail on first level, set abort_val and abort_pos in elem
        // std::cout << "saving..." << std::endl;
        // ToDo: need to sync this write?
        ++switches->experimental_deviation;
        elem->deviation_pos    = I->comp_fail_pos;
        elem->deviation_val    = I->comp_fail_val;
        elem->deviation_acc    = I->comp_fail_acc;
        elem->deviation_vertex = v;
        return false;
    }

    ++switches->experimental_paths;

    w->S.empty_cache();

    I->never_fail = true;
    do {
        if(switches->done_fast) return false;

        col = w->S.select_color_dynamic(g, c, strat);
        if(col == -1) break;
        rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[col] + 1));
        v    = c->lab[rpos];

        comp = proceed_state(w, g, c, I, v, nullptr, nullptr);
    } while(comp);

    if(comp && (strat->I->acc == I->acc)) { // automorphism computed
        bijection leaf;
        leaf.read_from_coloring(c);
        leaf.not_deletable();
        *automorphism = leaf;
        automorphism->inverse();
        automorphism->compose(strat->leaf);
        automorphism->non_uniform = false;
        if(!config.CONFIG_IR_FULL_INVARIANT && !w->R.certify_automorphism(g, automorphism)) {
            comp = false;
        } else {
            automorphism->certified = true;
        }
    } else {
        bijection leaf;
        leaf.read_from_coloring(c);
        leaf.not_deletable();
        // consider leaf store...
        std::vector<int*> pointers;

        switches->leaf_store_mutex.lock();
        auto range = switches->leaf_store.equal_range(I->acc);
        for (auto it = range.first; it != range.second; ++it)
            pointers.push_back(it->second);
        if(pointers.empty()) {
            switches->leaf_store.insert(std::pair<long, int *>(I->acc, leaf.map));
        }
        switches->leaf_store_mutex.unlock();

        comp = false;

        for(int i = 0; i < pointers.size(); ++i) {
            *automorphism = leaf;
            automorphism->inverse();
            bijection fake_leaf;
            fake_leaf.map = pointers[i];
            fake_leaf.map_sz = g->v_size;
            fake_leaf.not_deletable();
            automorphism->compose(&fake_leaf);

            int j;
            for(j = 0; j < automorphism->map_sz; ++j)
                if(automorphism->map[j] != j) break;

            if(j == automorphism->map_sz) {comp = false; break;}

            if(w->R.certify_automorphism(g, automorphism)) {
                // std::cout << "found auto" << std::endl;
                automorphism->certified = true;
                automorphism->non_uniform = false;
                if(i > 0) {
                    std::cout << "found on second" << std::endl;
                }
                comp = true;
                break;
            } else {
                std::cout << "should add / check more" << std::endl;
                comp = false;
            }
        }

    }

    return comp;
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
        int chunk_sz = w->BW->chunk_size;
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

bool dejavu::bfs_chunk(dejavu_workspace *w, sgraph *g, strategy *canon_strategy, bool *done, int selector_seed) {
    bfs_workspace *BFS = w->BW;
    int level = BFS->current_level;
    int target_level = BFS->target_level;
    if (level == target_level) return false; // we are done with BFS!

    // initialize bfs_workspace structures
    bfs_assure_init(w);

    if (w->generator_fix_base_alloc < *w->shared_generators_size) {
        delete[] w->generator_fix_base;
        w->generator_fix_base = new mpermnode *[*w->shared_generators_size];
        w->generator_fix_base_alloc = *w->shared_generators_size;
        w->prev_bfs_element = nullptr;
    }

    // try to dequeue a chunk of work
    // ToDo: try out next levels as well! (but should check target_level, though)
    size_t num = BFS->bfs_level_todo[level].try_dequeue_bulk(w->todo_dequeue, w->BW->chunk_size);
    int finished_elements_sz = 0;
    int todo_elements_sz = 0;

    int finished_elements_null_buffer = 0;

    for (int i = 0; i < num; ++i) {
        bfs_element *elem = std::get<0>(w->todo_dequeue[i]);
        int v = std::get<1>(w->todo_dequeue[i]);
        int weight = std::get<2>(w->todo_dequeue[i]);
        bool is_identity = elem->is_identity && (v == w->my_base_points[elem->base_sz]);

        // check orbit
        bool comp = elem->weight != 0;
        comp = comp && (elem->deviation_vertex != v);
        if(elem->deviation_pos > 0) {
            // check in abort map
            if (!elem->is_identity) {
                bool comp_ = BFS->read_abort_map(level, elem->deviation_pos, elem->deviation_acc);
                if(!comp_) {
                    //std::cout << "synergy" << std::endl;
                    elem->weight = 0;
                }
                comp = comp && comp_;
            }
        }

        if(elem->is_identity) {
            comp = true;
        }

        if(!comp) {
            //BFS->BW.abort_map_prune++;
        }

        if (weight == -1 && comp) {
            comp = comp && get_orbit(w, elem->base, elem->base_sz, v, w->my_base_points[elem->base_sz], &w->orbit,
                                     w->prev_bfs_element == elem);
            assert(!comp ? (!is_identity) : true);
        }

        if (comp) {
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
            w->work_I->reset_deviation();
            comp = comp && proceed_state(w, g, w->work_c, w->work_I, v, nullptr, nullptr); // &w->changes

            // manage abort map counter
            if (comp && elem->is_identity && level > 1) {
                // decrease abort map done...
                BFS->level_abort_map_mutex[level]->lock();
                BFS->level_abort_map_done[level]--;
                BFS->level_abort_map_mutex[level]->unlock();
            }

            // if !comp consider abort map
            if (!comp && level > 1) {
                if (elem->is_identity) { // save to abort map...
                    BFS->write_abort_map(level, w->work_I->comp_fail_pos, w->work_I->comp_fail_acc);
                } else { // if abort map done, check abort map...
                    bool comp_ = BFS->read_abort_map(level, w->work_I->comp_fail_pos, w->work_I->comp_fail_acc);
                    if (!comp_) {
                        elem->weight = 0; // element is pruned
                    }
                }
            }
        }

        assert(elem->base_sz < w->G->base_size);
        assert(!comp ? (!is_identity) : true);

        if (!comp) {
            // throw this node away, but keep track of that we computed it
            finished_elements_null_buffer += 1;
            continue;
        }

        // still looks equal to canonical base, so create a node
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

        for (int j = 0; j < elem->base_sz; ++j) {
            assert(elem->base[j] >= 0 && elem->base[j] < g->v_size);
            next_elem->base[j] = elem->base[j];
        }
        assert(v >= 0 && v < g->v_size);
        next_elem->base[next_elem->base_sz - 1] = v;
        if (weight == -1)
            next_elem->weight = elem->weight * w->orbit.cur_pos;
        else
            next_elem->weight = weight;
        next_elem->parent_weight = elem->weight;
        next_elem->parent = elem;

        // compute target color for this node
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

    if (finished_elements_null_buffer > 0) {
        w->finished_elements[finished_elements_sz] = std::pair<bfs_element *, int>(nullptr,
                                                                                   finished_elements_null_buffer);
        finished_elements_sz += 1;
    }

    if (finished_elements_sz > 0)
        BFS->bfs_level_finished_elements[level].enqueue_bulk(w->finished_elements, finished_elements_sz);

    return true;
}


bool dejavu::sequential_init_copy(dejavu_workspace* w) {
    // update sequential group in workspace
    bool new_gen = false;
    if(!w->sequential_init) {
        _newgroup(&w->sequential_gp, &w->sequential_gens, w->G->domain_size);
        w->sequential_init = true;
    }

    // copy generators that have not been copied yet
    mpermnode *it = w->G->gens;
    do {
        if(it->copied == 0) {
            new_gen = true;
            _addpermutation(&w->sequential_gens, it->p, w->G->domain_size);
            it->copied = 1;
        }
        it = it->next;
    } while (it != w->G->gens);

    return new_gen;
}

bool bfs_element_parent_sorter(bfs_element* const& lhs, bfs_element* const& rhs) {
    if(lhs->parent < rhs->parent)
        return true;
    if(lhs->parent == rhs->parent) {
        return(lhs->parent->parent < rhs->parent->parent);
    }
    return false;
}

void dejavu::bfs_reduce_tree(dejavu_workspace* w) {
    thread_local bool init_group = false;
    thread_local bool first_call = true;

    _schreier_fails(2);
    bfs_assure_init(w);

    int domain_size = w->G->domain_size;
    bool new_gens = sequential_init_copy(w);

    if(!new_gens && !first_call) {
        std::cout << "[B] Skipping reduce tree (no new generators)." << std::endl;
        return;
    }

    first_call = false;

    if(w->generator_fix_base_alloc < *w->shared_generators_size) {
        delete[] w->generator_fix_base;
        w->generator_fix_base = new mpermnode*[*w->shared_generators_size];
        w->generator_fix_base_alloc = *w->shared_generators_size;
        w->prev_bfs_element = nullptr;
    }

    // step1: check on level i if orbits can be reduced
    int i, j, v;
    bfs_workspace* BFS = w->BW;
    int active = 0;
    bool found_canon = false;
    for (int j = 0; j < BFS->level_sizes[BFS->current_level - 1]; ++j) {
        bfs_element *elem = BFS->level_states[BFS->current_level - 1][j];
        active += (elem->weight > 0);
        found_canon = found_canon || (elem->base[elem->base_sz - 1] == w->G->b[elem->base_sz - 1]);
    }
    assert(found_canon);
    std::cout << "[B] Top level active (before): " << active << std::endl;


    for(i = 1; i < BFS->current_level; ++i) {
        //std::cout << "Sorting with respect to parents..." << std::endl;
        bfs_element** arr = w->BW->level_states[i];
        int arr_sz        = w->BW->level_sizes[i];
        std::sort(arr, arr + arr_sz, &bfs_element_parent_sorter);

        //std::cout << "checking level " << i << ", sz " << BW->level_sizes[i] << std::endl;
        if(i == 1) {
            for(j = 0; j < BFS->level_sizes[i]; ++j) {
                bfs_element* elem = BFS->level_states[i][j];
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
            for(j = 0; j < BFS->level_sizes[i]; ++j) {
                bfs_element *elem = BFS->level_states[i][j];
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
            for(j = 0; j < BFS->level_sizes[i]; ++j) {
                bfs_element *elem = BFS->level_states[i][j];
                if(elem == nullptr)
                    continue;
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
    for (j = 0; j < BFS->level_sizes[BFS->current_level - 1]; ++j) {
        bfs_element *elem = BFS->level_states[BFS->current_level - 1][j];
        active += (elem->weight > 0);
    }

    std::cout << "[B] Top level active (after): " << active << std::endl;
    std::cout << "[B] Emptying queue..." << std::endl;
    moodycamel::ConcurrentQueue<std::tuple<bfs_element *, int, int>> throwaway_queue;
    BFS->bfs_level_todo[BFS->current_level].swap(throwaway_queue);


    // do not fill this if there is no hope...
    int expecting_finished = 0;
    for (int j = 0; j < BFS->level_sizes[BFS->current_level - 1]; ++j) {
        bfs_element *elem = BFS->level_states[BFS->current_level - 1][j];
        int added = 0;
        if (elem->weight > 0) {
            int c      = elem->target_color;
            int c_size = elem->c->ptn[c] + 1;
            expecting_finished += c_size;
        }
    }
    BFS->level_expecting_finished[BFS->current_level] = expecting_finished;

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
    if(w->BW->current_level == w->BW->target_level)
        return;
    moodycamel::ConcurrentQueue<std::tuple<bfs_element *, int, int>> throwaway_queue;
    w->BW->bfs_level_todo[w->BW->current_level].swap(throwaway_queue);

    int expected = 0;

    if(*w->shared_generators_size > 0)
        sequential_init_copy(w);

    // swap identity to first position...
    for (int j = 0; j < w->BW->level_sizes[w->BW->current_level - 1]; ++j) {
        bfs_element *elem = w->BW->level_states[w->BW->current_level - 1][j];
        if(elem->is_identity) {
            bfs_element *first_elem = w->BW->level_states[w->BW->current_level - 1][0];
            w->BW->level_states[w->BW->current_level - 1][j] = first_elem;
            w->BW->level_states[w->BW->current_level - 1][0] = elem;
            break;
        }
    }

    if(!w->sequential_init) {
            // then rest...
        for (int j = 0; j < w->BW->level_sizes[w->BW->current_level - 1]; ++j) {
            bfs_element *elem = w->BW->level_states[w->BW->current_level - 1][j];
            if (elem->weight > 0) {
                int c = elem->target_color;
                int c_size = elem->c->ptn[c] + 1;
                for (int i = c; i < c + c_size; ++i) {
                    expected += 1;
                    w->BW->bfs_level_todo[w->BW->current_level].enqueue(
                            std::tuple<bfs_element *, int, int>(elem, elem->c->lab[i], -1)); // ToDo: prune this, too!
                }
                if (elem->is_identity) {
                    std::cout << "Abort map expecting: " << c_size << std::endl;
                    w->BW->level_abort_map_done[w->BW->current_level] = c_size;
                }
            }
        }
    } else {
        std::cout << "Filling with orbits..." << std::endl;

        // swap identity to first position...
        for (int j = 0; j < w->BW->level_sizes[w->BW->current_level - 1]; ++j) {
            bfs_element *elem = w->BW->level_states[w->BW->current_level - 1][j];
            if(elem->is_identity) {
                bfs_element *first_elem = w->BW->level_states[w->BW->current_level - 1][0];
                w->BW->level_states[w->BW->current_level - 1][j] = first_elem;
                w->BW->level_states[w->BW->current_level - 1][0] = elem;
                break;
            }
        }

        int i;
        // could parallelize this easily?
        for (int j = 0; j < w->BW->level_sizes[w->BW->current_level - 1]; ++j) {
            bfs_element *elem = w->BW->level_states[w->BW->current_level - 1][j];
            int added = 0;
            if (elem->weight > 0) {
                int c = elem->target_color;
                int c_size = elem->c->ptn[c] + 1;
                int * orbits_sz;
                int * orbits;
                if(w->BW->current_level > 1) {
                     orbits_sz = nullptr;
                     orbits = _getorbits(elem->base, elem->base_sz, w->sequential_gp, &w->sequential_gens,
                                             w->G->domain_size, w->G->b, &orbits_sz);
                } else {
                     orbits_sz = *w->shared_orbit_weights; // <- something wrong here
                     orbits    = *w->shared_orbit;
                }
                for (i = c; i < c + c_size; ++i) {
                    if (orbits[elem->c->lab[i]] == elem->c->lab[i]) {
                        w->BW->bfs_level_todo[w->BW->current_level].enqueue(
                                std::tuple<bfs_element *, int, int>(elem, elem->c->lab[i], orbits_sz[elem->c->lab[i]])); // should be elem->c->lab[i], i no error?
                        expected += 1;
                        added += 1;
                    }
                }
                if (elem->is_identity) {
                    std::cout << "Abort map expecting: " << added << std::endl;
                    w->BW->level_abort_map_done[w->BW->current_level] = added;
                }
            }
        }
    }

    w->BW->level_expecting_finished[w->BW->current_level] = expected;
}

void dejavu::worker_thread(sgraph* g_, bool master, shared_workspace* switches, group_shared* G, coloring* start_c, strategy* canon_strategy,
                           int communicator_id, int** shared_orbit, int** shared_orbit_weights, bfs_workspace* bwork, mpermnode** gens, int* shared_group_size) {
    // find comparison leaf
    std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
    sgraph *g = g_;
    lowdeg* L;
    dejavu_workspace W;

    if(master) {
        config.CONFIG_IR_DENSE = !(g->e_size < g->v_size || g->e_size / g->v_size < g->v_size / (g->e_size / g->v_size));
        L       = new lowdeg();
        start_c = new coloring;
        g->initialize_coloring(start_c);
        if(config.CONFIG_PREPROCESS) {
          //  std::pair<sgraph *, coloring *> preprocessed_graph = L->preprocess(start_c, g, &W.R);
          //  g = preprocessed_graph.first;
          //  start_c = preprocessed_graph.second;
        }
        assert(start_c->check());
    }

    double cref;

    bool *done      = &switches->done;
    bool *done_fast = &switches->done_fast;
    int _shared_group_size = false;
    mpermnode *_gens       = nullptr;

    int*  shrd_orbit;
    int*  shrd_orbit_weights;
    int** shrd_orbit_;
    int** shrd_orbit_weights_;

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

        shrd_orbit_ = new (int*);
        *shrd_orbit_      = shrd_orbit;

        shrd_orbit_weights_ = new (int*);
        *shrd_orbit_weights_      = shrd_orbit_weights;

        switches->current_mode = modes::MODE_TOURNAMENT;

        // first color refinement
        canon_strategy = new strategy;

        W.start_c = start_c;
        //start_I.create_vector();
        W.R.refine_coloring_first(g, start_c, -1);
        cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "[T] Color ref: " << cref / 1000000.0 << "ms" << std::endl;

        // create some objects that are initialized after tournament
        G    = new group_shared(g->v_size);
        W.BW = new bfs_workspace();
        bwork = W.BW;

        int init_c = W.S.select_color(g, start_c, selector_seed);
        if(init_c == -1) {
            *done = true;
            std::cout << "First coloring discrete." << std::endl;
            std::cout << "Base size: 0" << std::endl;
            std::cout << "Group size: 1" << std::endl;
            //if(config.CONFIG_PREPROCESS)
            //    L->postprocess(nullptr);
            return;
        }
        W.S.empty_cache();
        // launch worker threads
        for (int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++)
            work_threads.emplace_back(
                    std::thread(&dejavu::worker_thread, dejavu(), g, false, switches, G, start_c,
                                canon_strategy, i, shrd_orbit_, shrd_orbit_weights_, W.BW, &_gens, &_shared_group_size));
        cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - timer).count());
        std::cout << "[T] Refinement workers created: " << cref / 1000000.0 << "ms" << std::endl;

        // set some workspace variables
        W.start_c = new coloring;
        W.start_c->copy_force(start_c);
        W.id        = -1;
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
        W.id = communicator_id;
        W.shared_generators = gens;
        W.shared_generators_size = shared_group_size;
    }

    strategy* my_strategy;
    {
        //W.base_size = G->base_size;
        my_canon_I = new invariant;
        //my_canon_I->create_vector();
        my_canon_I->has_compare = false;
        my_canon_I->compare_vec = nullptr;
        my_canon_I->compareI    = nullptr;
        my_canon_leaf = new bijection;
        selector_type rst = (selector_type) intRand(0, 3, selector_seed);
        my_strategy = new strategy(my_canon_leaf, my_canon_I, rst, -1);

        {
            W.S.empty_cache();
            find_first_leaf(&W, g, false, my_canon_I, my_canon_leaf, my_strategy, &base_points, &trash_int, switches,
                            selector_seed);
            std::cout << "I_sz: " << my_canon_I->vec_invariant->size() << std::endl;
            W.my_base_points    = base_points.map;
            W.my_base_points_sz = base_points.map_sz;
            W.is_foreign_base   = true;
        }

        // skip immediately if base size is 1
        if(master && W.my_base_points_sz == 1 && !config.CONFIG_IR_FULLBFS) { // ToDo: can end in deadlock if threads dont aggree on base size 1
            std::cout << "[N] Base size 1 skip" << std::endl;
            switches->base1_skip.store(2);
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
            W.BW->target_level.store(proposed_level);
            W.is_foreign_base = false;
            W.skip_schreier_level = G->gp;
            for(int i = 0; i < W.skiplevels; ++i)
                W.skip_schreier_level = W.skip_schreier_level->next;
            W.base_size = G->base_size;
            *W.shared_generators_size = 0;
            switches->done_shared_group = true;
            switches->current_mode = modes::MODE_BFS;
            std::cout << "[T] No tournament, created group by " << communicator_id << " with restarts " << restarts << std::endl;
        } else if(W.my_base_points_sz == 1 && !config.CONFIG_IR_FULLBFS) {
            while(switches->base1_skip.load() == 0) continue;
            if(switches->base1_skip == 2) {
                // wait until shared group created
                while (!switches->done_shared_group)
                    continue;
                while (!switches->current_mode == modes::MODE_BFS ||
                       !switches->current_mode == modes::MODE_UNIFORM_PROBE)
                    continue;
                W.my_base_points    = W.G->b;
                W.my_base_points_sz = W.G->base_size;
                W.is_foreign_base = false;
            }
        } else if(master) {
            switches->base1_skip.store(1);
        }
    }

    int n_found = 0;
    int n_restarts = 0;
    int rotate_i = 0;
    strategy_metrics m;
    bool foreign_base_done = false;
    bool reset_non_uniform_switch = true;
    bool increase_budget = true;
    int required_level = -1;

    // main loop...
    while(true) {
        if(master) {
            // sifting results
            G->manage_results(switches);
            if(*done) std::cout << "abort in head" << std::endl;

            // non-uniform search over, fix a group state for collaborative bfs_workspace
            if(switches->done_fast && !switches->done_shared_group && !switches->done) {
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
                    } else {
                        (*W.shared_orbit)[i] = G->b[0];
                    }
                for(int i = 0; i < g->v_size; ++i)
                    (*W.shared_orbit_weights)[(*W.shared_orbit)[i]]++;

                *W.shared_generators_size   = G->number_of_generators();

                if(W.BW->current_level > 1) { // > 1
                    bool tree_reduce = false;
                    if(*W.shared_generators_size > 0) {
                        std::cout << "[B] Reducing tree (" << n_found << ")" << std::endl;
                        assert(master && communicator_id == -1);
                        bfs_reduce_tree(&W);
                        tree_reduce = true;
                    }
                    // check if expected size is still too large...
                    if(W.BW->level_expecting_finished[W.BW->current_level] >= config.CONFIG_IR_SIZE_FACTOR * g->v_size * switches->tolerance) {
                        std::cout << "[B] Expected size still too large, not going into BFS" << std::endl;
                        W.BW->reached_initial_target = (W.BW->target_level == W.BW->current_level);
                        W.BW->target_level.store(W.BW->current_level);
                    } else {
                        switches->reset_tolerance(W.BW->level_expecting_finished[W.BW->current_level], g->v_size);
                        std::cout << "[T] tolerance " << switches->tolerance << std::endl;
                        std::cout << "[B] Filling queue..." << W.BW->current_level << " -> " << W.BW->target_level << std::endl;
                        bfs_fill_queue(&W);
                    }
                } else {
                    if(*W.shared_generators_size > 0) {
                        std::cout << "[B] Pre-re-fill" << std::endl;
                    //    bfs_fill_queue(&W);
                    }
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
                std::cout << "Abort map prune: " << W.BW->abort_map_prune << std::endl;
                std::cout << "Group size: ";
                /*G->sift_random(); // should sift random again here... or add generators from sequential group
                G->sift_random();
                G->sift_random();
                G->sift_random();*/
                G->print_group_size();
                std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
                cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
                std::cout << "Join: " << cref / 1000000.0 << "ms" << std::endl;
                //if(config.CONFIG_PREPROCESS)
                //    L->postprocess(G);
                break;
            }
        }

        if(switches->done_fast) {
            G->ack_done_shared();
            reset_non_uniform_switch = true;
        }

        if(switches->done) {
            if(!master)
                return;
            else
                continue;
        }

        bijection automorphism;
        abort_code A;
        automorphism.mark = false;

        // in what phase are we in?
        switch(switches->current_mode) {
            case modes::MODE_TOURNAMENT:
                m.restarts = 0;
                m.expected_bfs_size = 0;
                base_aligned_search(&W, g, my_strategy, &automorphism, &m, done_fast, switches,
                                    selector_seed); // <- we should already safe unsuccessfull / succ first level stuff here
                if(n_found == 0) { // check if I won
                    // wait until everyone checked
                    while(!switches->check_strategy_tournament(communicator_id, &m, false) && !switches->done_created_group) continue;
                    // check if I won, if yes: create group
                    if(switches->done_created_group) continue;
                    if(switches->win_id == communicator_id) {
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
                        W.BW->target_level.store(proposed_level);
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
                    base_aligned_search(&W, g, my_strategy, &automorphism, &m, done_fast, switches,
                                        selector_seed);
                    automorphism.foreign_base = true;
                    n_restarts += m.restarts;
                    automorphism.mark = true;
                } else {
                    base_aligned_search(&W, g, canon_strategy, &automorphism, &m, done_fast, switches,
                                        selector_seed);
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
                // fast automorphism search, but from initial bfs_workspace pieces
                if(!*done_fast) {
                    if (reset_non_uniform_switch) {
                        reset_skiplevels(&W);
                        if(!master) { // guess a new leaf
                            //delete my_canon_I->compare_vec;
                            my_canon_I = new invariant; // delete old if non canon
                            //my_canon_I->create_vector();
                            //delete my_canon_leaf;
                            my_canon_leaf = new bijection;
                            base_points = bijection();
                            selector_type rst = (selector_type) intRand(0, 2, selector_seed);
                            my_strategy = new strategy(my_canon_leaf, my_canon_I, rst, -1);
                            find_first_leaf(&W, g, false, my_canon_I, my_canon_leaf, my_strategy, &base_points,
                                            &trash_int, switches, selector_seed);
                            W.my_base_points    = base_points.map;
                            W.my_base_points_sz = base_points.map_sz;
                            W.is_foreign_base   = true;
                            foreign_base_done = false;
                        }
                        reset_non_uniform_switch = false;
                    }

                    if (!foreign_base_done) {
                        base_aligned_search(&W, g, my_strategy, &automorphism, &m,
                                            done_fast, switches, selector_seed);
                        automorphism.foreign_base = true;
                        automorphism.mark = true;
                        n_restarts += m.restarts;
                    } else {
                        base_aligned_search(&W, g, canon_strategy, &automorphism, &m,
                                            done_fast, switches, selector_seed);
                        automorphism.foreign_base = false;
                        automorphism.mark = true;
                        n_restarts += m.restarts;
                    }

                    if (master && n_found == 0) {
                        int proposed_level = std::max(W.skiplevels + 1, required_level);
                        if (proposed_level == G->base_size)
                            proposed_level += 1;
                        if (proposed_level > W.BW->target_level)
                            W.BW->target_level.store(proposed_level);
                    }

                    n_found += 1;
                    if (n_found % (3 * switches->tolerance) == 0 && (W.skiplevels < W.my_base_points_sz - 1))
                        W.skiplevels += 1;
                    if ((*done_fast && !automorphism.non_uniform)) continue;
                } else continue;
                if(*done && master) {
                    std::cout << "abort within non-uniform" << std::endl;
                    continue;
                }
                if ((*done_fast && !automorphism.non_uniform)) continue;
                break;

            case modes::MODE_NON_UNIFORM_FROM_BFS:
                {
                    // pick initial path from BFS level that is allocated to me
                    --switches->experimental_budget;
                    if(master && (switches->experimental_budget <= 0 || switches->done_fast)) {
                        if(!switches->done_fast) {
                            if(switches->experimental_paths > switches->experimental_deviation) {
                                if(!switches->experimental_look_close) {
                                    switches->experimental_look_close = true;
                                    switches->experimental_budget += W.BW->level_sizes[W.BW->current_level - 1];
                                    std::cout << "Switching to close look..." << std::endl;
                                    continue;
                                }
                            }
                        }

                        switches->experimental_budget = -1;
                        switches->current_mode = modes::MODE_WAIT;
                        switches->done_fast = true;
                        switches->done_shared_group = false;
                        continue;
                    }

                    bfs_element *elem;
                    int bfs_level    = W.BW->current_level - 1;
                    int max_weight   = W.BW->level_maxweight[bfs_level];
                    int bfs_level_sz = W.BW->level_sizes[bfs_level];
                    if(reset_non_uniform_switch) {
                        rotate_i = bfs_level_sz / (config.CONFIG_THREADS_REFINEMENT_WORKERS + 1);
                        rotate_i = rotate_i * (communicator_id + 1);
                        increase_budget = true;
                        reset_non_uniform_switch = false;
                        if(master) {
                            int proposed_level = std::max(bfs_level + 1, required_level);
                            if (proposed_level == G->base_size)
                                proposed_level += 1;
                            if (proposed_level > W.BW->target_level) {
                                W.BW->target_level.store(proposed_level);
                            }
                        }
                    }
                    int picked_weight, rand_weight;
                    do {
                        int pick_elem = intRand(0, bfs_level_sz - 1, selector_seed);
                        elem = W.BW->level_states[bfs_level][pick_elem];
                        picked_weight = elem->weight;
                        assert(max_weight > 0);
                        rand_weight   = doubleRand(1, max_weight, selector_seed);
                        if(rand_weight > picked_weight) continue;
                    } while (elem->weight <= 0 && !switches->done_fast && !switches->done); // && elem->deviation_vertex == -1
                    // compute one experimental path
                    bool comp = uniform_from_bfs_search_with_storage(&W, g, switches, elem, selector_seed,
                                                                     canon_strategy, &automorphism,
                                                                     switches->experimental_look_close);

                    if (!comp) {
                        // if failed, deduct experimental_budget and continue
                        continue;
                    } else {
                        if(increase_budget) {
                            increase_budget = false;
                            int budget_fac = switches->experimental_look_close?std::max(switches->tolerance, 10):1;
                            switches->experimental_budget += ((bfs_level_sz * 2 * budget_fac) / (config.CONFIG_THREADS_REFINEMENT_WORKERS + 1));
                            //std::cout << "Increasing experimental_budget..." << std::endl;
                        }
                        // otherwise add automorphism, if it exists...
                        automorphism.mark = true;
                    }
                }
                break;

            case modes::MODE_BFS:
                //std::cout << "b" << std::endl;
                reset_non_uniform_switch = true;
                if(W.is_foreign_base) {
                    reset_skiplevels(&W);
                    foreign_base_done = true;
                }
                if(W.BW->current_level != W.BW->target_level) {
                    if (communicator_id == -1 && W.BW->target_level < 0) {
                        int proposed_level = W.skiplevels + 1;
                        if (proposed_level == G->base_size)
                            proposed_level += 1;
                        W.BW->target_level.store(proposed_level);
                    }
                    if(switches->done_shared_group && W.BW->target_level >= 0) {
                        if(master && !switched1) {
                            switched1 = true;
                            cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
                            std::cout << "[N] Finished non-uniform automorphism search (" << *W.shared_generators_size << " generators, " << n_restarts << " restarts)" << std::endl;
                            std::cout << "[N] Ended in skiplevel " << W.skiplevels << ", found " << n_found << std::endl;
                            std::cout << "[T] " << cref / 1000000.0 << "ms" << std::endl;
                            std::cout << "[B] Determined target level: " << W.BW->target_level << "" << std::endl;
                        }
                        bfs_chunk(&W, g, canon_strategy, done, selector_seed);
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
                        if(bwork->reached_initial_target) {
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
                            //switches->current_mode = modes::MODE_NON_UNIFORM_PROBE_IT; // ToDo: actually should go to leaf tournament
                            int budget_fac = switches->experimental_look_close?std::max(switches->tolerance, 10):1;

                            std::cout << "[B] Uniform extensions, experimental_budget " << bwork->level_sizes[bwork->current_level - 1] * budget_fac << std::endl;
                            switches->experimental_budget.store(bwork->level_sizes[bwork->current_level - 1] * budget_fac);
                            switches->experimental_paths.store(0);
                            switches->experimental_deviation.store(0);
                            switches->current_mode = modes::MODE_NON_UNIFORM_FROM_BFS;
                            continue;
                        }
                    }
                }
                continue;
                break;

            case modes::MODE_UNIFORM_PROBE:
                reset_non_uniform_switch = true;
                if(W.id == 0 && !switched2) {
                    switched2 = true;
                    cref = (std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() - timer).count());
                    std::cout << "[T] " << cref / 1000000.0 << "ms" << std::endl;
                }
                A = uniform_from_bfs_search(&W, g, true, canon_strategy, &automorphism, &restarts, switches,
                                            selector_seed);
                if(A.reason == 2) {
                    //std::cout << "early abort" << std::endl;
                    continue;
                }
                if(A.reason == 1) {
                    // go back to bfs_workspace?
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
                    //std::cout << "[B] Non-uniform extensions, experimental_budget " << bwork->BW.level_sizes[bwork->BW.current_level - 1] << std::endl;
                    //switches->experimental_budget.store(bwork->BW.level_sizes[bwork->BW.current_level - 1]);
                    //switches->current_mode = MODE_NON_UNIFORM_FROM_BFS;
                    required_level = W.BW->current_level + 1;
                    std::cout << "[B] Requiring level " << required_level << std::endl;
                    continue;
                }
                automorphism.mark = true;
                break;

            case modes::MODE_WAIT:
                continue;
        }

        if(switches->done) {
            if(!master)
                return;
            else
                continue;
        }

        automorphism.not_deletable();
        //std::cout << automorphism.foreign_base << std::endl;
        bool test = true;
        if(switches->done_created_group && automorphism.mark && automorphism.certified) {
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
            //std::cout << "[N] Switching to canonical search (" << W.id << ", " << n_found << " generators)" << std::endl;
        }

        delete[] automorphism.map;
        automorphism.map = new int[g->v_size];
        automorphism.foreign_base = false;
        sampled_paths += 1;

        // master thread managing sifting results, bfs_workspace, ...
    }

    if(master) {
        delete[] shrd_orbit;
        delete[] shrd_orbit_weights;

        delete shrd_orbit_;
        delete shrd_orbit_weights_;

        delete G;
        delete bwork;
    }
    return;
}

void dejavu::automorphisms(sgraph *g) {
    shared_workspace switches;
    worker_thread(g, true, &switches, nullptr, nullptr, nullptr, -1,
                  nullptr, nullptr, nullptr, nullptr, nullptr);
}
