#ifndef DEJAVU_VUJADE_H
#define DEJAVU_VUJADE_H


#include <random>
#include "sgraph.h"
#include "invariant.h"
#include "concurrentqueue.h"
#include "selector.h"
#include "bfs.h"

struct abort_code {
    abort_code()=default;
    abort_code(int reason):reason(reason){};
    int reason = 0;
};

enum vujade_modes {VU_MODE_BIDIRECTIONAL, VU_MODE_BFS};

template <class vertex_t, class degree_t, class edge_t>
struct alignas(64) dejavu_workspace {
    // workspace for normal search
    refinement<vertex_t, degree_t, edge_t> R;
    selector<vertex_t, degree_t, edge_t>   S;
    coloring<vertex_t> c;
    invariant  I;

    // workspace for bfs_workspace
    coloring<vertex_t>*  work_c;
    invariant* work_I;

    // workspace for base aligned search
    int first_level = 1;
    int base_size = 0;
    int skiplevels = 1;
    int first_skiplevel = 1;
    coloring<vertex_t> skip_c;
    invariant  skip_I;
    bool       skiplevel_is_uniform = false;

    int*         my_base_points;
    int          my_base_points_sz;
    bool is_foreign_base;

    coloring<vertex_t>* start_c1;
    coloring<vertex_t>* start_c2;
    invariant start_I;

    // indicates which thread this is
    int id;

    // deprecated workspace for simple orbit method
    work_set  orbit_considered;
    work_list orbit_vertex_worklist;
    work_list orbit;
    int canonical_v;
    int         generator_fix_base_size;

    // bfs_workspace workspace
    bfs_workspace<vertex_t>* BW1;
    bfs_workspace<vertex_t>* BW2;
    std::tuple<bfs_element<vertex_t>*, int, int>* todo_dequeue;
    int todo_deque_sz        = -1;
    std::pair<bfs_element<vertex_t> *, int>* finished_elements;
    int finished_elements_sz = -1;
    std::pair<bfs_element<vertex_t> *, int>* todo_elements;
    int todo_elements_sz     = -1;
    bfs_element<vertex_t>* prev_bfs_element = nullptr;
    bool init_bfs = false;

    ~dejavu_workspace() {
        if(init_bfs) {
            delete[] todo_dequeue;
            delete[] todo_elements,
            delete[] finished_elements;
        }

        delete work_c;
        delete work_I;
        delete start_c1;
        delete start_c2;
    };
};

template<class vertex_t>
bool bfs_element_parent_sorter(bfs_element<vertex_t>* const& lhs, bfs_element<vertex_t>* const& rhs) {
    if(lhs->parent < rhs->parent)
        return true;
    if(lhs->parent == rhs->parent) {
        return(lhs->parent->parent < rhs->parent->parent);
    }
    return false;
}

template<class vertex_t, class degree_t, class edge_t>
class vujade_t {
public:
    bool iso(sgraph_t<vertex_t, degree_t, edge_t> *g1, sgraph_t<vertex_t, degree_t, edge_t> *g2) {
        bool simple_check = (g1->v_size == g2->v_size) && (g1->e_size == g2->e_size) && (g1->d_size == g2->d_size);
        if(!simple_check)
            return false;

        if(config.CONFIG_THREADS_REFINEMENT_WORKERS == -1) {
            const int max_threads = std::thread::hardware_concurrency();
            if (g1->v_size <= 100) {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = std::min(0, max_threads - 1);
            } else if(g1->v_size <= 150) {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = std::min(1, max_threads - 1);
            } else if(g1->v_size <= 200) {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = std::min(3, max_threads - 1);
            } else {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = max_threads - 1;
            }
        }

        shared_iso_workspace<vertex_t> switches;
        return worker_thread(g1, g2, true, &switches, nullptr, nullptr, nullptr,
                -1,nullptr, nullptr);
    }

private:
    bool worker_thread(sgraph_t<vertex_t, degree_t, edge_t> *g1_, sgraph_t<vertex_t, degree_t, edge_t> *g2_, bool master,
                       shared_iso_workspace<vertex_t> *switches, coloring<vertex_t> *start_c1,
                       coloring<vertex_t> *start_c2, strategy<vertex_t>* canon_strategy,
                       int communicator_id, bfs_workspace<vertex_t> *bwork1, bfs_workspace<vertex_t> *bwork2) {
        sgraph_t<vertex_t, degree_t, edge_t> *g1 = g1_;
        sgraph_t<vertex_t, degree_t, edge_t> *g2 = g2_;
        dejavu_workspace<vertex_t, degree_t, edge_t> W;

        numnodes  = 0;
        colorcost = 0;

        // preprocessing
        if(master) {
            config.CONFIG_IR_DENSE = !(g1->e_size<g1->v_size||g1->e_size/g1->v_size<g1->v_size/(g1->e_size/g1->v_size));
            start_c1 = new coloring<vertex_t>;
            g1->initialize_coloring(start_c1);
            start_c2 = new coloring<vertex_t>;
            g2->initialize_coloring(start_c2);
            if(config.CONFIG_PREPROCESS) {
                //  add preprocessing here
            }
            assert(start_c->check());
        }

        double cref;

        std::vector<std::thread> work_threads;
        bijection<vertex_t> base_points;
        bijection<vertex_t> actual_base;
        int trash_int = 0;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() * ((communicator_id * 5) * 5135235);
        int selector_seed = seed;

        invariant start_I;
        invariant *my_canon_I;
        bijection<vertex_t> *my_canon_leaf;
        my_canon_I = new invariant;
        my_canon_I->has_compare = false;
        my_canon_I->compare_vec = nullptr;
        my_canon_I->compareI    = nullptr;
        my_canon_leaf = new bijection<vertex_t>;

        // first color refinement, initialize some more shared structures, launch threads
        if (master) {
            PRINT("[Vuj] Dense graph: " << (config.CONFIG_IR_DENSE?"true":"false"));
            switches->current_mode = modes::MODE_TOURNAMENT;

            // first color refinement
            canon_strategy = new strategy<vertex_t>;
            my_canon_I->create_vector(g1->v_size);
            W.start_c1 = start_c1;
            strategy_metrics m;
            bool comp;
            W.R.refine_coloring(g1, start_c1, my_canon_I, -1, &m, -1);
            W.start_I.set_compare_invariant(my_canon_I);
            W.start_c2 = start_c2;
            comp = W.R.refine_coloring(g2, start_c2, &W.start_I, -1, &m, -1);
            delete my_canon_I;
            my_canon_I = new invariant;
            if(!comp)
                return comp;
            PRINT("[Vuj] First refinement: " << cref / 1000000.0 << "ms");

            int init_c = W.S.select_color(g1, start_c1, selector_seed);
            if(init_c == -1) {
                std::cout << "First coloring discrete." << std::endl;
                std::cout << "Base size: 0" << std::endl;
                std::cout << "Group size: 1" << std::endl;
                W.work_c = new coloring<vertex_t>;
                W.work_I = new invariant;
                delete canon_strategy;
                return false;
            }

            if(config.CONFIG_PREPROCESS_EDGELIST_SORT) {
                if (start_c1->cells == 1) {
                    for (int i = 0; i < start_c1->lab_sz; ++i) {
                        start_c1->lab[i] = i;
                        start_c1->vertex_to_lab[i] = i;
                    }
                    for (int i = 0; i < start_c2->lab_sz; ++i) {
                        start_c2->lab[i] = i;
                        start_c2->vertex_to_lab[i] = i;
                    }
                }

                g1->sort_edgelist();
                g2->sort_edgelist();
            }

            // create some objects that are initialized after tournament
            W.BW1 = new bfs_workspace<vertex_t>();
            bwork1 = W.BW1;

            W.BW2 = new bfs_workspace<vertex_t>();
            bwork2 = W.BW2;

            W.S.empty_cache();
            // launch worker threads
            for (int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++)
                work_threads.emplace_back(
                        std::thread(&vujade_t<vertex_t, degree_t, edge_t>::worker_thread,
                                    vujade_t<vertex_t, degree_t, edge_t>(), g1, g2, false, switches, start_c1, start_c2,
                                    canon_strategy, i, W.BW1, W.BW2));
            PRINT("[Vuj] Refinement workers created (" << config.CONFIG_THREADS_REFINEMENT_WORKERS << " threads)");

            // set some workspace variables
            W.start_c1 = new coloring<vertex_t>;
            W.start_c1->copy_force(start_c1);
            W.start_c2 = new coloring<vertex_t>;
            W.start_c2->copy_force(start_c2);
            W.id        = -1;
        }

        int sampled_paths = 0;
        int restarts = 0;
        int idle_ms = 0;
        int base_sz = 0;
        bool switched1 = false;
        bool switched2 = false;

        W.skip_c.copy_force(start_c1);
        W.work_c = new coloring<vertex_t>;
        W.work_I = new invariant;
        W.BW1 = bwork1;
        W.BW2 = bwork2;
        W.skiplevels = 0;

        if (!master) {
            W.start_c1 = new coloring<vertex_t>;
            W.start_c1->copy_force(start_c1);
            W.start_c2 = new coloring<vertex_t>;
            W.start_c2->copy_force(start_c2);
            W.id = communicator_id;
        }

        strategy<vertex_t>* my_strategy;

        auto rst = (selector_type) ((communicator_id + 2) % 3);
        if(config.CONFIG_IR_FORCE_SELECTOR)
            rst = (selector_type) config.CONFIG_IR_CELL_SELECTOR;

        my_strategy = new strategy<vertex_t>(my_canon_leaf, my_canon_I, rst, -1);

        W.S.empty_cache();
        find_first_leaf(&W, g1, my_canon_I, my_canon_leaf, my_strategy, &base_points, switches, selector_seed);
        base_sz = base_points.map_sz;
        if(std::is_same<vertex_t, int>::value) {
            W.my_base_points = (int*) base_points.map;
        } else {
            W.my_base_points = new int[base_points.map_sz];
            for(int i = 0; i < base_points.map_sz; ++i) {
                W.my_base_points[i] = static_cast<int>(base_points.map[i]);
            }
            base_points.deletable();
        }
        W.my_base_points_sz = base_points.map_sz;
        W.is_foreign_base   = true;

        if(master) {
            canon_strategy->replace(my_strategy);
            actual_base = base_points;
            base_points.not_deletable();
            PRINT("[Strat] Chosen strategy: " << canon_strategy->cell_selector_type);
            PRINT("[Strat] Combinatorial base size:" << W.my_base_points_sz);
            bfs_element<vertex_t> *root_elem1 = new bfs_element<vertex_t>;
            root_elem1->id = 0;
            root_elem1->c = new coloring<vertex_t>;
            root_elem1->I = new invariant;
            root_elem1->c->copy_force(start_c1);
            root_elem1->base_sz = 0;
            root_elem1->is_identity = true;
            *root_elem1->I = start_I;

            bfs_element<vertex_t> *root_elem2 = new bfs_element<vertex_t>;
            root_elem2->id = 0;
            root_elem2->c = new coloring<vertex_t>;
            root_elem2->I = new invariant;
            root_elem2->c->copy_force(start_c2);
            root_elem2->base_sz = 0;
            root_elem2->is_identity = true;
            *root_elem2->I = start_I;
            W.S.empty_cache();
            int init_c = W.S.select_color_dynamic(g1, start_c1, my_strategy);
            W.BW1->initialize(root_elem1, init_c, g1->v_size, base_sz);
            W.BW2->initialize(root_elem2, init_c, g2->v_size, base_sz);
            W.base_size = base_sz;
            // suppress extended deviation on last level, if there is only one level...
            if(W.base_size == 1)
                config.CONFIG_IR_EXPAND_DEVIATION = 0;
            switches->current_mode = modes::MODE_NON_UNIFORM_PROBE;
            switches->done_created_group = true;
        }

        while(!switches->done_created_group) {
            continue;
        }

        int n_found    = 0;
        int n_restarts = 0;
        int rotate_i   = 0;
        strategy_metrics m;
        bool foreign_base_done        = false;
        bool reset_non_uniform_switch = true;
        bool increase_budget   = true;
        bool is_canon_strategy = false;
        int required_level     = -1;

        vujade_modes mode = VU_MODE_BIDIRECTIONAL;
        int g_id = communicator_id % 2;

        // main loop
        while(!switches->done && (switches->noniso_counter <= config.CONFIG_RAND_ABORT)) {
            switch(mode) {
                case VU_MODE_BIDIRECTIONAL:
                {
                    g_id = (g_id + 1) % 2;
                    bfs_workspace<vertex_t>* bwork = (g_id == 0)?W.BW1:W.BW2;

                    // pick initial path from BFS level that is allocated to me
                    --switches->experimental_budget;
                    bijection<vertex_t> automorphism;
                    bfs_element<vertex_t> *elem;
                    int bfs_level    = bwork->current_level - 1;
                    int max_weight   = bwork->level_maxweight[bfs_level];
                    int bfs_level_sz = bwork->level_sizes[bfs_level];
                    double picked_weight, rand_weight;
                    do {
                        int pick_elem = intRand(0, bfs_level_sz - 1, selector_seed);
                        elem = bwork->level_states[bfs_level][pick_elem];
                        picked_weight = elem->weight;
                        assert(max_weight > 0);
                        rand_weight   = doubleRand(1, max_weight, selector_seed);
                        if(rand_weight > picked_weight) continue;
                    } while (elem->weight <= 0 && !switches->done); // && elem->deviation_vertex == -1
                    // compute one experimental path
                    bool comp = uniform_from_bfs_search_with_storage(&W, g1, g2, g_id, switches, elem, selector_seed,
                                                                     canon_strategy, &automorphism,
                                                                     true);
                    if(comp) {
                        switches->done = true;
                        break;
                    }
                }
                    break;

                case VU_MODE_BFS:

                    break;
            }
        }

        if(master && !dejavu_kill_request) {
            while (!work_threads.empty()) {
                work_threads[work_threads.size() - 1].join();
                work_threads.pop_back();
            }

            PRINT("[Vuj] Store " << switches->leaf_store[0].size() << ":" << switches->leaf_store[1].size());
            PRINT("[Vuj] Cleanup...");
        }
        return switches->found_iso;
    }

    void find_first_leaf(dejavu_workspace<vertex_t, degree_t, edge_t> *w,
                         sgraph_t<vertex_t, degree_t, edge_t> *g, invariant *canon_I,
                         bijection<vertex_t> *canon_leaf, strategy<vertex_t>* canon_strategy,
                         bijection<vertex_t> *automorphism, shared_iso_workspace<vertex_t> *switches,
                         int selector_seed) {
        const bool* done = &switches->done;

        // workspace
        refinement<vertex_t, degree_t, edge_t> *R = &w->R;
        selector<vertex_t, degree_t, edge_t> *S = &w->S;
        coloring<vertex_t> *c = &w->c;
        invariant *I = &w->I;
        coloring<vertex_t> *start_c  = w->start_c1;
        invariant *start_I = &w->start_I;

        S->empty_cache();
        start_I->reset_compare_invariant();
        start_I->create_vector(g->v_size * 2);
        automorphism->map    = new vertex_t[g->v_size];
        automorphism->map_sz = 0;

        *I = *start_I;
        c->copy(start_c);

        while (true) {
            if(*done) return;
            const int s = S->select_color_dynamic(g, c, canon_strategy);
            if (s == -1) {
                canon_leaf->read_from_coloring(c);
                *canon_I = *I;
                return;
            }

            // choose random vertex of class
            const int rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
            const int v = c->lab[rpos];

            // individualize and refine
            proceed_state(w, g, c, I, v, nullptr, nullptr, -1);
            assert(c->vertex_to_col[v] > 0);

            // base point
            automorphism->map[automorphism->map_sz] = v;
            automorphism->map_sz += 1;
        }
    }

    void reset_skiplevels(dejavu_workspace<vertex_t, degree_t, edge_t> *w) {
        w->skip_c.copy_force(w->start_c);
        w->skip_I = w->start_I;
        w->skiplevels = 0;
        w->skip_schreier_level = w->G->gp;
        w->first_skiplevel = 1;
        w->skiplevel_is_uniform = true;
        w->my_base_points    = w->G->b;
        w->my_base_points_sz = w->G->base_size;
        w->is_foreign_base   = false;
    }

    int extract_selector(invariant* I, int base_point) {
        if((I->compareI->vec_selections)->size() <= base_point)
            return - 1;
        return (*I->compareI->vec_selections)[base_point];
    }

    bool proceed_state(dejavu_workspace<vertex_t, degree_t, edge_t>* w,
                       sgraph_t<vertex_t, degree_t, edge_t> * g, coloring<vertex_t>* c,
                       invariant* I, int v, change_tracker* changes, strategy_metrics* m, int cell_early) {
        if(!config.CONFIG_IR_IDLE_SKIP)
            cell_early = -1;

        // protocol of selector choice
        I->selection_write(c->vertex_to_col[v]);
        const int init_color_class = w->R.individualize_vertex(c, v);
        bool comp = true;
        comp = comp && w->R.refine_coloring(g, c, I, init_color_class, m, cell_early);
        return comp;
    }

    bool bfs_chunk(dejavu_workspace<vertex_t, degree_t, edge_t> *w,
                   sgraph_t<vertex_t, degree_t, edge_t>  *g, strategy<vertex_t> *canon_strategy,
                   bool *done, int selector_seed) {
        bfs_workspace<vertex_t> *BFS = w->BW;
        int level = BFS->current_level;
        int target_level = BFS->target_level;
        if (level == target_level) return false; // we are done with BFS!

        // initialize bfs_workspace structures
        bfs_assure_init(w);

        // try to dequeue a chunk of work
        size_t num = BFS->bfs_level_todo[level].try_dequeue_bulk(w->todo_dequeue, w->BW->chunk_size);
        int finished_elements_sz = 0;
        int finished_elements_null_buffer = 0;

        for (size_t i = 0; i < num; ++i) {
            bfs_element<vertex_t> *elem = std::get<0>(w->todo_dequeue[i]);
            int v      = std::get<1>(w->todo_dequeue[i]);
            int weight = std::get<2>(w->todo_dequeue[i]);
            bool is_identity  = elem->is_identity && (v == w->my_base_points[elem->base_sz]);

            // check orbit
            bool comp = elem->weight != 0;
            comp = comp && (elem->deviation_vertex != v);
            if(elem->deviation_pos > 0) {
                // check in abort map
                if (!elem->is_identity) {
                    bool comp_ = BFS->read_abort_map(level, elem->deviation_pos, elem->deviation_acc);

                    if(!comp_) {
                        if(elem->weight != 0 && elem->deviation_write.try_lock()) {
                            elem->weight = 0;
                            elem->deviation_write.unlock();
                        }
                    }
                    comp = comp && comp_;
                }
            }

            if(elem->is_identity) {
                comp = true;
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

                numnodes++;
                // compute next coloring
                w->work_I->reset_deviation();
                const int cells_early = (*w->work_I->compareI->vec_cells)[elem->base_sz];
                comp = comp && proceed_state(w, g, w->work_c, w->work_I, v, nullptr, nullptr, cells_early); // &w->changes

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
            bfs_element<vertex_t> *next_elem = new bfs_element<vertex_t>;
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

            w->finished_elements[finished_elements_sz] = std::pair<bfs_element<vertex_t> *, int>(next_elem, sz);
            finished_elements_sz += 1;
            w->work_c = new coloring<vertex_t>;
            w->work_I = new invariant;
        }

        if (finished_elements_null_buffer > 0) {
            w->finished_elements[finished_elements_sz] = std::pair<bfs_element<vertex_t> *, int>(nullptr,
                                                                                                 finished_elements_null_buffer);
            finished_elements_sz += 1;
        }

        if (finished_elements_sz > 0)
            BFS->bfs_level_finished_elements[level].enqueue_bulk(w->finished_elements, finished_elements_sz);

        return true;
    }

    void bfs_fill_queue(dejavu_workspace<vertex_t, degree_t, edge_t> *w) {
        if(w->BW->current_level == w->BW->target_level)
            return;
        moodycamel::ConcurrentQueue<std::tuple<bfs_element<vertex_t> *, int, int>> throwaway_queue;
        w->BW->bfs_level_todo[w->BW->current_level].swap(throwaway_queue);

        int expected = 0;

        if(*w->shared_generators_size > 0)
            sequential_init_copy(w);

        // swap identity to first position...
        for (int j = 0; j < w->BW->level_sizes[w->BW->current_level - 1]; ++j) {
            bfs_element<vertex_t> *elem = w->BW->level_states[w->BW->current_level - 1][j];
            if(elem->is_identity) {
                bfs_element<vertex_t> *first_elem = w->BW->level_states[w->BW->current_level - 1][0];
                w->BW->level_states[w->BW->current_level - 1][j] = first_elem;
                w->BW->level_states[w->BW->current_level - 1][0] = elem;
                break;
            }
        }

        if(!w->sequential_init) {
            // then rest...
            for (int j = 0; j < w->BW->level_sizes[w->BW->current_level - 1]; ++j) {
                bfs_element<vertex_t> *elem = w->BW->level_states[w->BW->current_level - 1][j];
                if (elem->weight > 0) {
                    int c = elem->target_color;
                    int c_size = elem->c->ptn[c] + 1;
                    for (int i = c; i < c + c_size; ++i) {
                        expected += 1;
                        w->BW->bfs_level_todo[w->BW->current_level].enqueue(
                                std::tuple<bfs_element<vertex_t> *, int, int>(elem, elem->c->lab[i], -1));
                    }
                    if (elem->is_identity) {
                        PRINT("[BFS] Abort map expecting: " << c_size);
                        w->BW->level_abort_map_done[w->BW->current_level] = c_size;
                    }
                }
            }
        } else {
            PRINT("[BFS] Filling with orbits...");

            // swap identity to first position...
            for (int j = 0; j < w->BW->level_sizes[w->BW->current_level - 1]; ++j) {
                bfs_element<vertex_t> *elem = w->BW->level_states[w->BW->current_level - 1][j];
                if(elem->is_identity) {
                    bfs_element<vertex_t> *first_elem = w->BW->level_states[w->BW->current_level - 1][0];
                    w->BW->level_states[w->BW->current_level - 1][j] = first_elem;
                    w->BW->level_states[w->BW->current_level - 1][0] = elem;
                    break;
                }
            }

            int i;
            // could parallelize this easily?
            for (int j = 0; j < w->BW->level_sizes[w->BW->current_level - 1]; ++j) {
                bfs_element<vertex_t> *elem = w->BW->level_states[w->BW->current_level - 1][j];
                int added = 0;
                if (elem->weight > 0) {
                    int c = elem->target_color;
                    int c_size = elem->c->ptn[c] + 1;
                    int * orbits_sz;
                    int * orbits;
                    for (i = c; i < c + c_size; ++i) {
                        w->BW->bfs_level_todo[w->BW->current_level].enqueue(
                                std::tuple<bfs_element<vertex_t> *, int, int>(
                                        elem, elem->c->lab[i], orbits_sz[elem->c->lab[i]]));
                        expected += 1;
                        added += 1;
                    }
                    if (elem->is_identity) {
                        PRINT("[BFS] Abort map expecting: " << added);
                        w->BW->level_abort_map_done[w->BW->current_level] = added;
                    }
                }
            }
        }

        w->BW->level_expecting_finished[w->BW->current_level] = expected;
    }

    void bfs_assure_init(dejavu_workspace<vertex_t, degree_t, edge_t> *w) {
        if(!w->init_bfs) {
            int chunk_sz = w->BW->chunk_size;
            w->todo_dequeue = new std::tuple<bfs_element<vertex_t>*, int, int>[chunk_sz];
            w->todo_deque_sz = chunk_sz;
            w->todo_elements = new std::pair<bfs_element<vertex_t> *, int>[chunk_sz * 8];
            w->todo_elements_sz = chunk_sz * 8;
            w->finished_elements = new std::pair<bfs_element<vertex_t> *, int>[chunk_sz + 1];
            w->finished_elements_sz = chunk_sz;
            w->init_bfs = true;
        }
    }

    bool uniform_from_bfs_search_with_storage(dejavu_workspace<vertex_t, degree_t, edge_t> *w,
                                              sgraph_t<vertex_t, degree_t, edge_t>  *g1, sgraph_t<vertex_t, degree_t, edge_t>  *g2,
                                              int g_id, shared_iso_workspace<vertex_t>* switches, bfs_element<vertex_t> *elem,
                                              int selector_seed, strategy<vertex_t> *strat,
                                              bijection<vertex_t> *automorphism, bool look_close) {
        sgraph_t<vertex_t, degree_t, edge_t>* g;
        if(g_id == 0) {
            g = g1;
        } else {
            g = g2;
        }

        coloring<vertex_t>* c = &w->c;
        invariant* I = &w->I;
        c->copy_force(elem->c);
        *I = *elem->I;
        I->set_compare_invariant(strat->I);

        bool comp;

        // first individualization
        {
            const int col = elem->target_color;
            const int rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[col] + 1));
            const int v = c->lab[rpos];
            int cell_early = -1;
            if (look_close) {
                I->never_fail = true;
            } else {
                I->never_fail = false;
                cell_early = (*I->compareI->vec_cells)[elem->base_sz];
            }
            I->reset_deviation();
            if(v != elem->deviation_vertex) {
                comp = proceed_state(w, g, c, I, v, nullptr, nullptr, cell_early);
                if (!comp) { // fail on first level, set deviation acc, pos and vertex in elem
                    ++switches->experimental_deviation;
                    if (elem->deviation_write.try_lock() && elem->deviation_pos == -1) {
                        elem->deviation_pos = I->comp_fail_pos;
                        elem->deviation_acc = I->comp_fail_acc;
                        elem->deviation_vertex = v;
                        elem->deviation_write.unlock();
                    }
                    return false;
                }
            } else {
                return false;
            }
        }

        ++switches->experimental_paths;
        w->S.empty_cache();

        I->never_fail = true;
        do {
            if(switches->done_fast) return false;
            const int col = w->S.select_color_dynamic(g, c, strat);
            if(col == -1) break;
            const int rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[col] + 1));
            const int v    = c->lab[rpos];

            comp = proceed_state(w, g, c, I, v, nullptr, nullptr, -1);
        } while(comp);

        /*if(comp && (strat->I->acc == I->acc)) { // automorphism computed
            bijection<vertex_t> leaf;
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
        } else {*/
            bijection<vertex_t> leaf;
            leaf.read_from_coloring(c);
            leaf.not_deletable();
            // consider leaf store...

            for (int gi = 0; gi <= 1; ++gi) {
                std::vector<vertex_t*> pointers;
                switches->leaf_store_mutex[gi].lock();
                auto range = switches->leaf_store[gi].equal_range(I->acc);
                for (auto it = range.first; it != range.second; ++it)
                    pointers.push_back(it->second);
                if (pointers.empty()) {
                    if(gi == g_id)
                        switches->leaf_store[gi].insert(std::pair<long, vertex_t *>(I->acc, leaf.map));
                }
                switches->leaf_store_mutex[gi].unlock();

                comp = false;

                for (size_t i = 0; i < pointers.size(); ++i) {
                    automorphism->copy(&leaf);
                    automorphism->inverse();
                    bijection<vertex_t> fake_leaf;
                    fake_leaf.map = pointers[i];
                    fake_leaf.map_sz = g->v_size;
                    fake_leaf.not_deletable();
                    automorphism->compose(&fake_leaf);

                    if(gi == g_id) {
                        if (w->R.certify_automorphism_iso(g, automorphism)) {
                            switches->noniso_counter++;
                            PRINT("[Bid] Found uniform automorphism. (" << g_id << ")");
                            automorphism->certified = true;
                            automorphism->non_uniform = false;
                            break;
                        }
                    } else {
                        if (w->R.certify_isomorphism(g_id?g2:g1, (1-g_id)?g2:g1, automorphism)) {
                            PRINT("[Bid] Found isomorphism. (" << g_id << ")");
                            switches->found_iso.store(true);
                            switches->noniso_counter.store(-10);
                            return true;
                        }
                    }
                }

                if(!pointers.empty() && gi == g_id) {
                    switches->leaf_store_mutex[gi].lock();
                    switches->leaf_store[gi].insert(std::pair<long, vertex_t *>(I->acc, leaf.map));
                    switches->leaf_store_mutex[gi].unlock();
                }
            }
        return false;
    }
};

typedef vujade_t<int, int, int> vujade;

/*void vujade_automorphisms_dispatch(dynamic_sgraph *sgraph, shared_permnode **gens) {
    switch(sgraph->type) {
        case sgraph_type::DSG_INT_INT_INT: {
            PRINT("[Dispatch] <int32, int32, int32>");
            dejavu_t<int, int, int> d;
            d.automorphisms(sgraph->sgraph_0, gens);
        }
            break;
        case sgraph_type::DSG_SHORT_SHORT_INT: {
            PRINT("[Dispatch] <int16, int16, int>");
            dejavu_t<int16_t, int16_t, int> d;
            d.automorphisms(sgraph->sgraph_1, gens);
        }
            break;
        case sgraph_type::DSG_SHORT_SHORT_SHORT: {
            PRINT("[Dispatch] <int16, int16, int16>");
            dejavu_t<int16_t, int16_t, int16_t> d;
            d.automorphisms(sgraph->sgraph_2, gens);
        }
            break;
        case sgraph_type::DSG_CHAR_CHAR_SHORT:{
            PRINT("[Dispatch] <int8, int8, int16>");
            dejavu_t<int8_t, int8_t, int16_t> d;
            d.automorphisms(sgraph->sgraph_3, gens);
        }
            break;
        case sgraph_type::DSG_CHAR_CHAR_CHAR: {
            PRINT("[Dispatch] <int8, int8, int8>");
            dejavu_t<int8_t, int8_t, int8_t> d;
            d.automorphisms(sgraph->sgraph_4, gens);
        }
            break;
    }
}*/

bool vujade_iso(sgraph_t<int, int, int> *g1, sgraph_t<int, int, int> *g2) {
    vujade v;
    bool res = v.iso(g1, g2);
    if(res) {
        PRINT("ISOMORPHIC");
    } else {
        PRINT("NON_ISOMORPHIC");
    }
    return res;
}

#endif //DEJAVU_VUJADE_H
