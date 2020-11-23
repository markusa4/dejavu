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

enum uniform_outcome {OUT_NONE, OUT_AUTO, OUT_ISO, OUT_AUTO_DEV, OUT_ISO_DEV};

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

template<class vertex_t, class degree_t, class edge_t>
class vujade_t {
public:
    bool iso(sgraph_t<vertex_t, degree_t, edge_t> *g1, sgraph_t<vertex_t, degree_t, edge_t> *g2) {
        if(config.CONFIG_THREADS_REFINEMENT_WORKERS == -1) {
            const int max_threads = std::thread::hardware_concurrency();
            if (g1->v_size <= 150) {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = std::min(0, max_threads - 1);
            } else if(g1->v_size <= 200) {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = std::min(1, max_threads - 1);
            } else if(g1->v_size <= 250) {
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

        config.CONFIG_VUJADE = true; // some things ought to work differently now...

        numnodes  = 0;
        colorcost = 0;

        // preprocessing
        if(master) {
            config.CONFIG_IR_FORCE_EXPAND_DEVIATION = true;
            config.CONFIG_IR_DENSE = !(g1->e_size<g1->v_size||g1->e_size/g1->v_size<g1->v_size/(g1->e_size/g1->v_size));
            start_c1 = new coloring<vertex_t>;
            g1->initialize_coloring_raw(start_c1);
            start_c2 = new coloring<vertex_t>;
            g2->initialize_coloring_raw(start_c2);
            if(config.CONFIG_PREPROCESS) {
                //  add preprocessing here
            }
            assert(start_c->check());
        }

        double cref;

        std::vector<std::thread> work_threads;
        bijection<vertex_t> base_points;
        bijection<vertex_t> actual_base;
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
            PRINT("[vuj] Dense graph: " << (config.CONFIG_IR_DENSE?"true":"false"));
            switches->current_mode = vujade_modes ::VU_MODE_BIDIRECTIONAL_DEVIATION;

            // first color refinement
            canon_strategy = new strategy<vertex_t>;
            my_canon_I->create_vector(g1->v_size);
            W.start_c1 = start_c1;
            strategy_metrics m;
            bool comp;
            W.R.refine_coloring(g1, start_c1, my_canon_I, -1, &m, -1);
            W.start_I.set_compare_invariant(my_canon_I);
            W.start_c2 = start_c2;
            comp = W.R.refine_coloring(g2, start_c2, &W.start_I, -1, &m, start_c2->cells);
            delete my_canon_I;
            my_canon_I = new invariant;
            if(!comp) {
                W.work_c = new coloring<vertex_t>;
                W.work_I = new invariant;
                delete canon_strategy;
                return comp;
            }
            PRINT("[vuj] First refinement: " << cref / 1000000.0 << "ms");

            int init_c = W.S.select_color(g1, start_c1, selector_seed);
            if(init_c == -1) {
                std::cout << "First coloring discrete, checking isomorphism." << std::endl;
                bijection<vertex_t> automorphism;
                automorphism.read_from_coloring(start_c1);
                bijection<vertex_t> leaf2;
                leaf2.read_from_coloring(start_c2);
                automorphism.inverse();
                automorphism.compose(&leaf2);
                W.work_c = new coloring<vertex_t>;
                W.work_I = new invariant;
                delete canon_strategy;

                return (W.R.certify_isomorphism(g1, g2, &automorphism));
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

            const int master_sched = sched_getcpu();
            
            W.S.empty_cache();
            {
                cpu_set_t cpuset;
                CPU_ZERO(&cpuset);
                CPU_SET(master_sched, &cpuset);
                int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
            }
            // launch worker threads
            for (int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++) {
                work_threads.emplace_back(
                        std::thread(&vujade_t<vertex_t, degree_t, edge_t>::worker_thread,
                                    vujade_t<vertex_t, degree_t, edge_t>(), g1, g2, false, switches, start_c1, start_c2,
                                    canon_strategy, i, W.BW1, W.BW2));
                cpu_set_t cpuset;
                CPU_ZERO(&cpuset);
                CPU_SET(i + (i >= master_sched), &cpuset);
                int rc = pthread_setaffinity_np(work_threads[i].native_handle(),
                                                sizeof(cpu_set_t), &cpuset);
            }
            PRINT("[vuj] Refinement workers created (" << config.CONFIG_THREADS_REFINEMENT_WORKERS << " threads)");

            // set some workspace variables
            W.start_c1 = new coloring<vertex_t>;
            W.start_c1->copy_force(start_c1);
            W.start_c2 = new coloring<vertex_t>;
            W.start_c2->copy_force(start_c2);
            W.id        = -1;
        }

        int base_sz = 0;
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
        const int gid_first = communicator_id % 2 == 0;
        if(gid_first == 0) {
            find_first_leaf(&W, g1, my_canon_I, my_canon_leaf, my_strategy, &base_points, switches, selector_seed);
        } else {
            find_first_leaf(&W, g2, my_canon_I, my_canon_leaf, my_strategy, &base_points, switches, selector_seed);
        }
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

        // add first leaf to leaf store
        switches->leaf_store_mutex[gid_first].lock();
        my_canon_leaf->not_deletable();
        switches->leaf_store[gid_first].insert(std::pair<long, stored_leaf<vertex_t>>(my_canon_I->acc,
                stored_leaf<vertex_t>(my_canon_leaf->map, g1->v_size, true)));
        switches->leaf_store_mutex[gid_first].unlock();

        if(master) {
            canon_strategy->replace(my_strategy);
            actual_base = base_points;
            base_points.not_deletable();
            PRINT("[strat] Chosen strategy: " << canon_strategy->cell_selector_type);
            PRINT("[strat] Combinatorial base size:" << W.my_base_points_sz);
            W.S.empty_cache();
            int init_c = W.S.select_color_dynamic(g1, start_c1, my_strategy);
            W.BW1->initialize(bfs_element<vertex_t>::root_element(start_c1, &start_I), init_c, g1->v_size, 2);
            W.BW2->initialize(bfs_element<vertex_t>::root_element(start_c2, &start_I), init_c, g2->v_size, 2);
            W.base_size = base_sz;
            // suppress extended deviation on last level, if there is only one level...
            if(W.base_size == 0) {
                switches->current_mode = vujade_modes::VU_MODE_BIDIRECTIONAL;
                config.CONFIG_IR_EXPAND_DEVIATION = 0;
            } else {
                config.CONFIG_IR_EXPAND_DEVIATION = floor(1.5*sqrt(g1->v_size));
                PRINT("[dev] Expansion: " << config.CONFIG_IR_EXPAND_DEVIATION);
            }
            switches->done_created_group = true;
        }

        while(!switches->done_created_group) {
            continue;
        }

        strategy_metrics m;
        int  g_id = (communicator_id + 1) % 2;
        int  end_test = (g_id + 1) % 2;
        bool found_auto = false;

        // main loop
        while(!switches->done && (switches->noniso_counter <= config.CONFIG_RAND_ABORT)) {
            switch(switches->current_mode ) {
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
                    int pick_elem = intRand(0, bfs_level_sz - 1, selector_seed);
                    elem = bwork->level_states[bfs_level][pick_elem];

                    // compute one experimental path
                    uniform_outcome res = uniform_from_bfs_search_with_storage(&W, g1, g2, g_id, switches, elem, selector_seed,
                                                                     canon_strategy, &automorphism,true);

                    switch(res) {
                        case OUT_AUTO:
                            found_auto = true;
                            if(g_id == end_test && found_auto) {
                                switches->noniso_counter++;
                                found_auto = false;
                            }
                            break;
                        case OUT_ISO:
                            switches->done = true;
                            switches->found_iso.store(true);
                            switches->noniso_counter.store(-10);
                            break;
                        default:
                            break;
                    }
                }
                    break;

                case VU_MODE_BIDIRECTIONAL_DEVIATION:
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
                    int pick_elem = intRand(0, bfs_level_sz - 1, selector_seed);
                    elem = bwork->level_states[bfs_level][pick_elem];

                    // compute one experimental path
                    uniform_outcome res = uniform_deviation_from_bfs_search_with_storage(&W, g1, g2, g_id, switches, elem, selector_seed,
                                                                               canon_strategy, &automorphism,false);
                    switch(res) {
                        case OUT_AUTO:
                            found_auto = true;
                            if(g_id == end_test && found_auto) {
                                switches->noniso_counter++;
                                found_auto = false;
                            }
                            break;
                        case OUT_ISO:
                            switches->done = true;
                            switches->found_iso.store(true);
                            switches->noniso_counter.store(-10);
                            break;
                        case OUT_AUTO_DEV:
                            found_auto = true;
                            if(g_id == end_test && found_auto) {
                                switches->noniso_counter++;
                                found_auto = false;
                            }
                            break;
                        case OUT_ISO_DEV:
                            switches->switch_mode_mutex.lock();
                            // config.CONFIG_RAND_ABORT = 1;
                            //config.CONFIG_IR_EXPAND_DEVIATION = floor(sqrt(g1->v_size));
                            if(switches->current_mode != VU_MODE_BIDIRECTIONAL) {
                                switches->current_mode = VU_MODE_BIDIRECTIONAL;
                                switches->noniso_counter.store(0);
                                switches->switch_mode_mutex.unlock();
                                continue;
                            }
                            //} else {
                            switches->switch_mode_mutex.unlock();
                            //}
                            break;
                        default:
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

            PRINT("[vuj] Store: leafs(" << switches->leaf_store[0].size() << ":" << switches->leaf_store[1].size()
            << "), " << "devs(" << switches->deviation_store[0].size() << ":" << switches->deviation_store[1].size() << ")");
            PRINT("[vuj] Cleanup...");
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

    void reconstruct_leaf(dejavu_workspace<vertex_t, degree_t, edge_t> *w,
                          sgraph_t<vertex_t, degree_t, edge_t>  *g, coloring<vertex_t>* start_c, vertex_t* base,
                          int base_sz, bijection<vertex_t> *leaf) {
        coloring<vertex_t>* c = &w->c;
        invariant* I = &w->I;
        c->copy_force(start_c);
        I->never_fail = true;
        for(int pos = 0; pos < base_sz; ++pos) {
            const int v    = base[pos];
            proceed_state(w, g, c, I, v, nullptr, nullptr, -1);
        };
        leaf->read_from_coloring(c);
    }

    uniform_outcome uniform_from_bfs_search_with_storage(dejavu_workspace<vertex_t, degree_t, edge_t> *w,
                                              sgraph_t<vertex_t, degree_t, edge_t>  *g1, sgraph_t<vertex_t, degree_t, edge_t>  *g2,
                                              int g_id, shared_iso_workspace<vertex_t>* switches, bfs_element<vertex_t> *elem,
                                              int selector_seed, strategy<vertex_t> *strat,
                                              bijection<vertex_t> *automorphism, bool look_close) {
        sgraph_t<vertex_t, degree_t, edge_t>* g;
        thread_local work_list_t<vertex_t> collect_base;
        thread_local int collect_base_sz;

        if(!collect_base.init)
            collect_base.initialize(g1->v_size);
        collect_base.reset();
        collect_base_sz = 0;

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
        bool found_auto = false;

        // first individualization
        {
            const int col = elem->target_color;
            const int rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[col] + 1));
            const int v = c->lab[rpos];
            collect_base.push_back(v);
            collect_base_sz += 1;
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
                    return OUT_NONE;
                }
            } else {
                return OUT_NONE;
            }
        }

        ++switches->experimental_paths;
        w->S.empty_cache();

        I->never_fail = true;
        do {
            if(switches->done) return OUT_NONE;
            const int col = w->S.select_color_dynamic(g, c, strat);
            if(col == -1) break;
            const int rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[col] + 1));
            const int v    = c->lab[rpos];
            collect_base.push_back(v);
            collect_base_sz += 1;
            comp = proceed_state(w, g, c, I, v, nullptr, nullptr, -1);
        } while(comp);

        bijection<vertex_t> leaf;
        leaf.read_from_coloring(c);
        leaf.not_deletable();

        // consider leaf store...
        for (int gi = 0; gi <= 1; ++gi) {
            std::vector<stored_leaf<vertex_t>> pointers;
            switches->leaf_store_mutex[gi].lock();
            auto range = switches->leaf_store[gi].equal_range(I->acc);
            for (auto it = range.first; it != range.second; ++it)
                pointers.push_back(it->second);
            if (pointers.empty()) {
                if(gi == g_id) {
                    if (switches->leaf_store_explicit <= config.CONFIG_IR_LEAF_STORE_LIMIT) {
                        switches->leaf_store_explicit++;
                        switches->leaf_store[gi].insert(std::pair<long,
                                stored_leaf<vertex_t>>(I->acc, stored_leaf<vertex_t>(leaf.map, g1->v_size, true)));
                    } else {
                        vertex_t* base = new vertex_t[collect_base_sz];
                        memcpy(base, collect_base.arr, sizeof(vertex_t) * collect_base_sz);
                        switches->leaf_store[gi].insert(std::pair<long,
                                stored_leaf<vertex_t>>(I->acc, stored_leaf<vertex_t>(base, collect_base_sz, false)));
                        leaf.deletable();
                    }
                }
            }
            switches->leaf_store_mutex[gi].unlock();

            for (size_t i = 0; i < pointers.size(); ++i) {
                automorphism->copy(&leaf);
                automorphism->inverse();
                bijection<vertex_t> fake_leaf;
                if(pointers[i].explicit_leaf) {
                    fake_leaf.map = pointers[i].map;
                    fake_leaf.not_deletable();
                } else {
                    reconstruct_leaf(w, (gi==0)?g1:g2, (gi==0)?(w->start_c1):(w->start_c2),
                            pointers[i].map,pointers[i].map_sz, &fake_leaf);
                    fake_leaf.deletable();
                }
                fake_leaf.map_sz = g->v_size;
                automorphism->compose(&fake_leaf);

                if(gi == g_id) {
                    if (w->R.certify_automorphism_iso(g, automorphism)) {

                        int j;
                        for(j = 0; j < automorphism->map_sz; ++j)
                            if(automorphism->map[j] != j) break;

                        PRINT("[bid] Found uniform automorphism. (" << g_id << "), (" << (j == automorphism->map_sz) << ")");
                        automorphism->certified = true;
                        automorphism->non_uniform = false;
                        found_auto = true;
                    }
                } else {
                    if (w->R.certify_isomorphism(g_id?g2:g1, (1-g_id)?g2:g1, automorphism)) {
                        PRINT("[bid] Found isomorphism. (" << g_id << ")");
                        return OUT_ISO;
                    }
                }
            }

            if(!pointers.empty() && gi == g_id && !found_auto) {
                if (switches->leaf_store_explicit <= config.CONFIG_IR_LEAF_STORE_LIMIT) {
                    switches->leaf_store_explicit++;
                    switches->leaf_store[gi].insert(std::pair<long,
                            stored_leaf<vertex_t>>(I->acc, stored_leaf<vertex_t>(leaf.map, g1->v_size, true)));
                } else {
                    vertex_t* base = new vertex_t[collect_base_sz];
                    memcpy(base, collect_base.arr, sizeof(vertex_t) * collect_base_sz);
                    switches->leaf_store[gi].insert(std::pair<long,
                            stored_leaf<vertex_t>>(I->acc, stored_leaf<vertex_t>(base, collect_base_sz, false)));
                    leaf.deletable();
                }
            }
        }
        return (found_auto?OUT_AUTO:OUT_NONE);
    }

    uniform_outcome uniform_deviation_from_bfs_search_with_storage(dejavu_workspace<vertex_t, degree_t, edge_t> *w,
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
        I->reset_deviation();
        I->never_fail = false;

        bool comp = true;
        bool found_auto = false;
        w->S.empty_cache();

        int level = 0;

        do {
            if(switches->done) return OUT_NONE;
            if(c->cells == g->v_size) break;
            //const int col = w->S.select_color_dynamic(g, c, strat);
            const int col = strat->I->selection_read(level);
            ++level;
            const int rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[col] + 1));
            const int v    = c->lab[rpos];
            const int cell_early = (*I->compareI->vec_cells)[elem->base_sz];
            comp = proceed_state(w, g, c, I, v, nullptr, nullptr, cell_early);
        } while(comp);

        if(comp) {
            bijection<vertex_t> leaf;
            leaf.read_from_coloring(c);
            leaf.not_deletable();

            // consider leaf store...
            for (int gi = 0; gi <= 1; ++gi) {
                std::vector<stored_leaf<vertex_t>> pointers;
                switches->leaf_store_mutex[gi].lock();
                auto range = switches->leaf_store[gi].equal_range(I->acc);
                for (auto it = range.first; it != range.second; ++it)
                    pointers.push_back(it->second);
                if (pointers.empty()) {
                    if (gi == g_id) {
                        if(switches->leaf_store_explicit <= config.CONFIG_IR_LEAF_STORE_LIMIT) {
                            switches->leaf_store_explicit++;
                            switches->leaf_store[gi].insert(std::pair<long,
                                    stored_leaf<vertex_t>>(I->acc, stored_leaf<vertex_t>(leaf.map, g1->v_size, true)));
                        } else {

                        }
                    }
                }
                switches->leaf_store_mutex[gi].unlock();

                for (size_t i = 0; i < pointers.size(); ++i) {
                    automorphism->copy(&leaf);
                    automorphism->inverse();
                    bijection<vertex_t> fake_leaf;
                    if(pointers[i].explicit_leaf)
                        fake_leaf.map = pointers[i].map;
                    fake_leaf.map_sz = g->v_size;
                    fake_leaf.not_deletable();
                    automorphism->compose(&fake_leaf);

                    if (gi == g_id) {
                        if (w->R.certify_automorphism_iso(g, automorphism)) {

                            int j;
                            for (j = 0; j < automorphism->map_sz; ++j)
                                if (automorphism->map[j] != j) break;

                            PRINT("[bid] Found uniform automorphism. (" << g_id << "), (" << (j == automorphism->map_sz)
                                                                        << ")");
                            automorphism->certified = true;
                            automorphism->non_uniform = false;
                            found_auto = true;
                        }
                    } else {
                        if (w->R.certify_isomorphism(g_id ? g2 : g1, (1 - g_id) ? g2 : g1, automorphism)) {
                            PRINT("[bid] Found isomorphism. (" << g_id << ")");
                            return OUT_ISO;
                        }
                    }
                }

                if (!pointers.empty() && gi == g_id && !found_auto) {
                    switches->leaf_store_mutex[gi].lock();
                    switches->leaf_store[gi].insert(std::pair<long, stored_leaf<vertex_t>>(I->acc,
                            stored_leaf<vertex_t>(leaf.map, g1->v_size, true)));
                    switches->leaf_store_mutex[gi].unlock();
                }
            }
            return (found_auto ? OUT_AUTO : OUT_NONE);
        } else {
            long acc = I->comp_fail_acc;
            // consider deviation stores...
            for (int gi = 0; gi <= 1; ++gi) {
                switches->deviation_store_mutex[gi].lock();
                auto find = switches->deviation_store[gi].find(acc);
                bool found = !(find == switches->deviation_store[gi].end());
                if (!found) {
                    if (gi == g_id)
                        switches->deviation_store[gi].insert(acc);
                }
                switches->deviation_store_mutex[gi].unlock();
                if(!found) {
                    continue;
                } else {
                    if (gi == g_id) {
                        found_auto = true;
                        PRINT("[bid] Found uniform auto deviation. (" << g_id << ")" << "(" << level << ")");
                    } else {
                        PRINT("[bid] Found iso deviation. (" << g_id << ")" << "(" << level << ")");
                        return OUT_ISO_DEV;
                    }
                }
            }
            return (found_auto ? OUT_AUTO_DEV : OUT_NONE);
        }
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
    bool res = (g1->v_size == g2->v_size) && (g1->e_size == g2->e_size) && (g1->d_size == g2->d_size);
    if(res) {
        vujade v;
        res = v.iso(g1, g2);
    }
    return res;
}

#endif //DEJAVU_VUJADE_H
