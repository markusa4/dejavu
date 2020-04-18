#ifndef DEJAVU_DEJAVU_H
#define DEJAVU_DEJAVU_H


#include <random>
#include "sgraph.h"
#include "invariant.h"
#include "concurrentqueue.h"
#include "selector.h"
#include "bfs.h"
#include "schreier_shared.h"
#include "group_shared.h"
#include "schreier_sequential.h"

struct abort_code {
    abort_code()=default;
    abort_code(int reason):reason(reason){};
    int reason = 0;
};

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
    shared_schreier* skip_schreier_level;
    bool       skiplevel_is_uniform = false;

    int*         my_base_points;
    int          my_base_points_sz;
    bool is_foreign_base;

    group_shared<vertex_t>* G;

    coloring<vertex_t>* start_c;
    invariant start_I;

    // indicates which thread this is
    int id;

    // shared orbit and generators
    int** shared_orbit;
    int** shared_orbit_weights;
    shared_permnode** shared_generators;
    int*  shared_generators_size;
    int   generator_fix_base_alloc = -1;

    // sequential, local group
    sequential_permnode*      sequential_gens;
    sequential_schreierlevel* sequential_gp;
    bool            sequential_init = false;

    // deprecated workspace for simple orbit method
    work_set  orbit_considered;
    work_list orbit_vertex_worklist;
    work_list orbit;
    int canonical_v;
    shared_permnode** generator_fix_base;
    int         generator_fix_base_size;

    // bfs_workspace workspace
    bfs_workspace<vertex_t>* BW;
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
            delete[] generator_fix_base;
        }
        if(sequential_init) {
            _freeschreier(&sequential_gp, &sequential_gens);
        }

        delete work_c;
        delete work_I;
        delete start_c;
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
class dejavu_t {
public:
    void automorphisms(sgraph_t<vertex_t, degree_t, edge_t> *g, vertex_t* colmap, shared_permnode **gens) {
        if(config.CONFIG_THREADS_REFINEMENT_WORKERS == -1) {
            const int max_threads = std::thread::hardware_concurrency();
            if (g->v_size <= 100) {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = std::min(0, max_threads - 1);
            } else if(g->v_size <= 150) {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = std::min(1, max_threads - 1);
            } else if(g->v_size <= 200) {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = std::min(3, max_threads - 1);
            } else {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = max_threads - 1;
            }
        }

        coloring<vertex_t> start_c;
        start_c.vertex_to_col = colmap;
        start_c.init = false;
        shared_workspace<vertex_t> switches;
        worker_thread(g, true, &switches, nullptr, &start_c, nullptr, -1,
                      nullptr, nullptr, nullptr, gens, nullptr);
    }

private:
    void worker_thread(sgraph_t<vertex_t, degree_t, edge_t> *g_, bool master,
                       shared_workspace<vertex_t> *switches, group_shared<vertex_t> *G,
                       coloring<vertex_t> *start_c, strategy<vertex_t>* canon_strategy,
                       int communicator_id, int **shared_orbit, int** shared_orbit_weights,
                       bfs_workspace<vertex_t> *bwork, shared_permnode **gens, int *shared_group_size) {
        sgraph_t<vertex_t, degree_t, edge_t> *g = g_;
        dejavu_workspace<vertex_t, degree_t, edge_t> W;

        numnodes  = 0;
        colorcost = 0;

        // preprocessing
        if(master) {
            config.CONFIG_IR_DENSE = !(g->e_size<g->v_size||g->e_size/g->v_size<g->v_size/(g->e_size/g->v_size));
            coloring<vertex_t>* init_c  = start_c;
            start_c = new coloring<vertex_t>;
            g->initialize_coloring(start_c, init_c->vertex_to_col);
            if(config.CONFIG_PREPROCESS) {
                //  add preprocessing here
            }
            assert(start_c->check());
        }

        double cref;

        bool *done             = &switches->done;
        bool *done_fast        = &switches->done_fast;

        int _shared_group_size   = false;
        shared_permnode *_gens   = nullptr;
        int*  shrd_orbit         = nullptr;
        int*  shrd_orbit_weights = nullptr;
        int** shrd_orbit_        = nullptr;
        int** shrd_orbit_weights_= nullptr;

        std::vector<std::thread> work_threads;
        bijection<vertex_t> base_points;
        bijection<vertex_t> actual_base;
        int trash_int = 0;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() * ((communicator_id * 5) * 5135235);
        int selector_seed = seed;

        invariant start_I;

        // first color refinement, initialize some more shared structures, launch threads
        if (master) {
            PRINT("[Dej] Dense graph: " << (config.CONFIG_IR_DENSE?"true":"false"));
            switches->current_mode = modes::MODE_TOURNAMENT;

            // first color refinement
            canon_strategy = new strategy<vertex_t>;
            W.start_c = start_c;
            W.R.refine_coloring_first(g, start_c, -1);
            PRINT("[Dej] First refinement: " << cref / 1000000.0 << "ms");

            int init_c = W.S.select_color(g, start_c, selector_seed);
            if(init_c == -1) {
                *done = true;
                std::cout << "First coloring discrete." << std::endl;
                std::cout << "Base size: 0" << std::endl;
                std::cout << "Group size: 1" << std::endl;
                W.work_c = new coloring<vertex_t>;
                W.work_I = new invariant;
                delete canon_strategy;
                return;
            }

            if(config.CONFIG_PREPROCESS_EDGELIST_SORT) {
                if (start_c->cells == 1) {
                    for (int i = 0; i < start_c->lab_sz; ++i) {
                        start_c->lab[i] = i;
                        start_c->vertex_to_lab[i] = i;
                    }
                }

                g->sort_edgelist();
            }

            shrd_orbit = new int[g->v_size];
            for(int i = 0; i < g->v_size; ++i)
                shrd_orbit[i] = i;

            shrd_orbit_weights = new int[g->v_size];
            memset(shrd_orbit_weights, 0, g->v_size * sizeof(int));

            shrd_orbit_ = new (int*);
            *shrd_orbit_= shrd_orbit;

            shrd_orbit_weights_ = new (int*);
            *shrd_orbit_weights_= shrd_orbit_weights;

            // create some objects that are initialized after tournament
            G    = new group_shared<vertex_t>(g->v_size);
            W.BW = new bfs_workspace<vertex_t>();
            bwork = W.BW;

            W.S.empty_cache();
            // launch worker threads
            for (int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++)
                work_threads.emplace_back(
                        std::thread(&dejavu_t<vertex_t, degree_t, edge_t>::worker_thread,
                                    dejavu_t<vertex_t, degree_t, edge_t>(), g, false, switches, G, start_c,
                                    canon_strategy, i, shrd_orbit_, shrd_orbit_weights_, W.BW, &_gens,
                                    &_shared_group_size));
            PRINT("[Dej] Refinement workers created (" << config.CONFIG_THREADS_REFINEMENT_WORKERS << " threads)");

            // set some workspace variables
            W.start_c = new coloring<vertex_t>;
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
        bijection<vertex_t> *my_canon_leaf;
        bool switched1 = false;
        bool switched2 = false;

        W.skip_c.copy_force(start_c);
        W.work_c = new coloring<vertex_t>;
        W.work_I = new invariant;
        W.G = G;
        W.BW = bwork;
        W.skiplevels = 0;
        W.skip_schreier_level = G->gp;

        if (!master) {
            W.shared_orbit = shared_orbit;
            W.shared_orbit_weights = shared_orbit_weights;
            W.start_c = new coloring<vertex_t>;
            W.start_c->copy_force(start_c);
            W.id = communicator_id;
            W.shared_generators = gens;
            W.shared_generators_size = shared_group_size;
        }

        strategy<vertex_t>* my_strategy;

        my_canon_I = new invariant;
        my_canon_I->has_compare = false;
        my_canon_I->compare_vec = nullptr;
        my_canon_I->compareI    = nullptr;
        my_canon_leaf = new bijection<vertex_t>;
        auto rst = (selector_type) ((communicator_id + 2) % 3);
        if(config.CONFIG_IR_FORCE_SELECTOR)
            rst = (selector_type) config.CONFIG_IR_CELL_SELECTOR;

        my_strategy = new strategy<vertex_t>(my_canon_leaf, my_canon_I, rst, -1);

        W.S.empty_cache();
        find_first_leaf(&W, g, my_canon_I, my_canon_leaf, my_strategy, &base_points, switches, selector_seed);
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

        int n_found    = 0;
        int n_restarts = 0;
        int rotate_i   = 0;
        strategy_metrics m;
        bool foreign_base_done        = false;
        bool reset_non_uniform_switch = true;
        bool increase_budget   = true;
        bool is_canon_strategy = false;
        int required_level     = -1;

        // main loop
        while(true) {
            // master thread management
            if(master) {
                if(dejavu_kill_request) *done = true;

                // manage sifting results
                if(!switches->done)
                    G->manage_results(switches);

                // base-aligned search is over, we now fix a group state for collaborative bfs_workspace
                if(switches->done_fast && !switches->done_shared_group && !switches->done) {
                    // wait for ack of done_fast
                    PRINT("[N] Waiting for ACK");
                    switches->current_mode = modes::MODE_WAIT;
                    G->ack_done_shared();
                    reset_non_uniform_switch = true;
                    G->wait_for_ack_done_shared(config.CONFIG_THREADS_REFINEMENT_WORKERS + 1, &switches->done);

                    if(switches->done)
                        continue;

                    PRINT("[N] Creating shared orbit and generators");
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
                        if(*W.shared_generators_size > 0) {
                            PRINT("[BFS] Reducing tree (" << n_found << ")");
                            assert(master && communicator_id == -1);
                            bfs_reduce_tree(&W);
                        }
                        // check if expected size is still too large...
                        if(W.BW->level_expecting_finished[W.BW->current_level] >=
                           config.CONFIG_IR_SIZE_FACTOR * g->v_size * switches->tolerance) {
                            PRINT("[BFS] Expected size still too large, not going into BFS")
                            W.BW->reached_initial_target = (W.BW->target_level == W.BW->current_level);
                            W.BW->target_level.store(W.BW->current_level);
                        } else {
                            switches->reset_tolerance(W.BW->level_expecting_finished[W.BW->current_level], g->v_size);
                            PRINT("[Dej] Tolerance: " << switches->tolerance);
                            PRINT("[BFS] Filling queue..." << W.BW->current_level << " -> " << W.BW->target_level)
                            bfs_fill_queue(&W);
                        }
                    } else {
                        if(*W.shared_generators_size > 0) {
                            // PRINT("[BFS] Reducing queue using orbits")
                            //    bfs_fill_queue(&W);
                        }
                    }

                    switches->done_shared_group = true;
                    switches->current_mode      = modes::MODE_BFS;
                }

                // we are done
                if(switches->done) {
                    if(!dejavu_kill_request) {
                        std::cout << "Base size:  " << G->base_size << std::endl;
                        while (!work_threads.empty()) {
                            work_threads[work_threads.size() - 1].join();
                            work_threads.pop_back();
                        }
                        std::cout << "Group size: ";
                        G->print_group_size();

                        if(config.CONFIG_WRITE_AUTOMORPHISMS) {
                            std::cout << "Generators: " << std::endl;
                            shared_permnode *it = G->gens;
                            if(it != nullptr) {
                                do {
                                    for (int i = 0; i < g->v_size; ++i) {
                                        std::cout << it->p[i] << " ";
                                    }
                                    std::cout << std::endl;
                                    it = it->next;
                                } while (it != G->gens);
                            }
                            if(gens != nullptr) {
                                *gens = G->gens;
                            }
                        }

                        std::cout << "Join: " << cref / 1000000.0 << "ms" << std::endl;
                        std::cout << "Numnodes (master): " << numnodes << std::endl;
                        std::cout << "Colorcost (master): " << colorcost << std::endl;
                    } else {
                        while (!work_threads.empty()) {
                            work_threads[work_threads.size() - 1].join();
                            work_threads.pop_back();
                        }
                        std::cout << "Killed" << std::endl;
                    }
                    break;
                }
            }

            if(switches->done_fast) {
                G->ack_done_shared();
                reset_non_uniform_switch = true;
            }

            if(switches->done) {
                if(!master)
                    break;
                else
                    continue;
            }

            bijection<vertex_t> automorphism;
            abort_code A;
            automorphism.mark = false;

            // high-level solver mode
            switch(switches->current_mode) {
                case modes::MODE_TOURNAMENT:
                    m.restarts = 0;
                    m.expected_bfs_size = 0;
                    W.skiplevel_is_uniform = false;
                    base_aligned_search(&W, g, my_strategy, &automorphism, &m, done_fast, switches, selector_seed);
                    // check if this thread won
                    if(n_found == 0) {
                        // wait until everyone checked
                        while(!switches->check_strategy_tournament(communicator_id, &m, false)
                              && !switches->done_created_group) continue;
                        // if this thread won, it will now create the group
                        if(switches->win_id == communicator_id) {
                            canon_strategy->replace(my_strategy);
                            is_canon_strategy = true;
                            actual_base = base_points;
                            base_points.not_deletable();
                            G->initialize(g->v_size, &actual_base);
                            PRINT("[Strat] Chosen strategy: " << canon_strategy->cell_selector_type);
                            W.S.empty_cache();
                            int init_c = W.S.select_color_dynamic(g, start_c, my_strategy);
                            W.BW->initialize(bfs_element<vertex_t>::root_element(start_c, &start_I), init_c,
                                             g->v_size, G->base_size);
                            int proposed_level = W.skiplevels + 1;
                            if(config.CONFIG_IR_FULL_BFS)
                                proposed_level = G->base_size + 1;
                            W.BW->target_level.store(proposed_level);
                            PRINT("[Strat] Proposed level for BFS: " << proposed_level);
                            W.is_foreign_base = false;
                            W.skiplevel_is_uniform = W.skiplevels == 0;
                            W.skip_schreier_level = G->gp;
                            for(int i = 0; i < W.skiplevels; ++i)
                                W.skip_schreier_level = W.skip_schreier_level->next;

                            W.base_size = G->base_size;
                            // suppress extended deviation on last level, if there is only one level...
                            if(W.base_size == 1)
                                config.CONFIG_IR_EXPAND_DEVIATION = 0;

                            foreign_base_done = true;
                            switches->current_mode = modes::MODE_NON_UNIFORM_PROBE;
                            switches->done_created_group = true;
                            PRINT("[Strat] Created shared group by " << communicator_id << " with restarts " << restarts);
                        }

                        while(!(switches->done_created_group)) continue;
                        if(switches->all_no_restart && W.is_foreign_base) {
                            reset_skiplevels(&W);
                            foreign_base_done = true;
                            W.skiplevel_is_uniform = (W.skiplevels == 0);
                        }
                    }
                    automorphism.foreign_base = true;
                    automorphism.mark = true;
                    n_restarts += m.restarts;
                    n_found += 1;
                    if((*done_fast && !automorphism.certified)) continue;
                    break;

                case modes::MODE_NON_UNIFORM_PROBE:
                    if(!foreign_base_done) {
                        base_aligned_search(&W, g, my_strategy, &automorphism, &m, done_fast, switches,
                                            selector_seed);
                        automorphism.foreign_base = true;
                        n_restarts += m.restarts;
                        automorphism.mark = true;
                    } else {
                        abort_code a = base_aligned_search(&W, g, canon_strategy, &automorphism, &m, done_fast, switches,
                                                           selector_seed);

                        if(a.reason == 1) {
                            *done      = true;
                            *done_fast = true;
                        }

                        automorphism.foreign_base = false;
                        n_restarts += m.restarts;
                        automorphism.mark = true;
                    }
                    n_found += 1;
                    if((*done_fast && !automorphism.non_uniform )) continue;
                    break;

                case modes::MODE_NON_UNIFORM_PROBE_IT:
                    // automorphism search from initial bfs_workspace pieces
                    if(!*done_fast) {
                        if (reset_non_uniform_switch) {
                            reset_skiplevels(&W);
                            if(!master) { // guess a new leaf
                                if(!is_canon_strategy) {
                                    delete my_canon_I;
                                    delete my_canon_leaf;
                                    delete my_strategy;
                                }

                                is_canon_strategy = false;
                                my_canon_I    = new invariant; // delete old if non canon
                                my_canon_leaf = new bijection<vertex_t>;
                                base_points   = bijection<vertex_t>();
                                auto rst      = (selector_type) intRand(0, 2, selector_seed);
                                my_strategy   = new strategy<vertex_t>(my_canon_leaf, my_canon_I, rst, -1);
                                find_first_leaf(&W, g, my_canon_I, my_canon_leaf, my_strategy, &base_points, switches,
                                                selector_seed);
                                if(std::is_same<vertex_t, int>::value) {
                                    W.my_base_points = (int*) base_points.map;
                                } else {
                                    W.my_base_points = new int[base_points.map_sz];
                                    for(int i = 0; i < base_points.map_sz; ++i) {
                                        W.my_base_points[i] = static_cast<int>(base_points.map[i]);
                                    }
                                    base_points.deletable();
                                }
                                W.my_base_points_sz    = base_points.map_sz;
                                W.skiplevel_is_uniform = false;
                                W.is_foreign_base      = true;
                                foreign_base_done      = false;
                            }
                            reset_non_uniform_switch = false;
                        }

                        if(*done && !master)
                            break;

                        if (!foreign_base_done) {
                            base_aligned_search(&W, g, my_strategy, &automorphism, &m,
                                                done_fast, switches, selector_seed);
                            automorphism.foreign_base = true;
                            automorphism.mark = true;
                            n_restarts += m.restarts;
                        } else {
                            abort_code a = base_aligned_search(&W, g, canon_strategy, &automorphism, &m,
                                                               done_fast, switches, selector_seed);
                            if(a.reason == 1) {
                                *done      = true;
                                *done_fast = true;
                            }
                            automorphism.foreign_base = false;
                            automorphism.mark = true;
                            n_restarts += m.restarts;
                        }

                        if (master && n_found == 0) {
                            int proposed_level = required_level; // consider skiplevel (or previous "proposed level") here?
                            if (proposed_level == G->base_size)
                                proposed_level += 1;
                            if (proposed_level > G->base_size + 1)
                                proposed_level = G->base_size + 1;
                            if (proposed_level > W.BW->target_level)
                                W.BW->target_level.store(proposed_level);
                        }

                        n_found += 1;
                        if ((*done_fast && !automorphism.non_uniform)) continue;
                    } else continue;
                    if(*done && master) {
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
                                    PRINT("[UStore] Switching to close look...");
                                    continue;
                                }
                            }
                        }

                        if(config.CONFIG_IR_FAST_TOLERANCE_INC) {
                            if (switches->experimental_paths * 16 < switches->experimental_deviation) {
                                switches->iterate_tolerance();
                                switches->iterate_tolerance();
                                switches->iterate_tolerance();
                            }
                        }

                        switches->experimental_budget = -1;
                        switches->current_mode = modes::MODE_WAIT;
                        switches->done_fast = true;
                        switches->done_shared_group = false;
                        continue;
                    }

                    bfs_element<vertex_t> *elem;
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
                            if (proposed_level > G->base_size + 1)
                                proposed_level = G->base_size + 1;
                            if (proposed_level > W.BW->target_level) {
                                W.BW->target_level.store(proposed_level);
                            }
                        }
                    }
                    double picked_weight, rand_weight;
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
                            const int budget_fac = switches->experimental_look_close?std::max(switches->tolerance, 10):1;
                            switches->experimental_budget += ((bfs_level_sz * 2 * budget_fac) /
                                                              (config.CONFIG_THREADS_REFINEMENT_WORKERS + 1));
                        }
                        // otherwise add automorphism, if it exists...
                        automorphism.mark = true;
                    }
                }
                    break;

                case modes::MODE_BFS:
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
                            if (proposed_level > G->base_size + 1)
                                proposed_level = G->base_size + 1;
                            W.BW->target_level.store(proposed_level);
                        }
                        if(switches->done_shared_group && W.BW->target_level >= 0) {
                            if(master && !switched1) {
                                switched1 = true;
                                PRINT("[BA] Finished non-uniform automorphism search (" << *W.shared_generators_size
                                                                                          << " generators, " << n_restarts << " restarts)")
                                PRINT("[BA] Ended in skiplevel " << W.skiplevels << ", found " << n_found)
                                PRINT("[BA] " << cref / 1000000.0 << "ms")
                                PRINT("[BFS] Determined target level: " << W.BW->target_level << "")
                            }
                            bfs_chunk(&W, g, canon_strategy, done, selector_seed);
                            if(master) {
                                bool fill = bwork->work_queues(switches->tolerance);
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
                                if(bwork->current_level - 1 >= 0 && bwork->level_sizes[bwork->current_level - 1] == 1
                                                                 && bwork->current_level == bwork->base_size + 1) {
                                    PRINT("[BFS] Early-out");
                                    *done = true;
                                    continue;
                                }
                                PRINT("[UTarget] Starting uniform probe, tolerance: " << switches->tolerance)
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
                                PRINT("[UTarget] Iterating, tolerance: " << switches->tolerance)
                                reset_skiplevels(&W);
                                foreign_base_done = true;
                                const int budget_fac = switches->experimental_look_close?std::max(switches->tolerance, 10):1;

                                PRINT("[UStore] Switching to uniform with leaf storage, budget "
                                              << bwork->level_sizes[bwork->current_level - 1] * budget_fac)
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
                        PRINT("[UTarget] " << cref / 1000000.0 << "ms")
                    }
                    A = uniform_from_bfs_search(&W, g, canon_strategy, &automorphism, &restarts, switches, selector_seed);
                    if(A.reason == 2) // abort
                        continue;

                    if(A.reason == 1) { // too many restarts
                        switches->current_mode = MODE_WAIT;
                        // manage
                        switches->iterate_tolerance();
                        switches->done_fast = false;
                        switches->done_shared_group = false;
                        G->non_uniform_abort_counter = 0;
                        n_found = 0;
                        switched1 = false;
                        bwork->reset_initial_target();
                        PRINT("[Dej] Tolerance: " << switches->tolerance)
                        reset_skiplevels(&W);
                        foreign_base_done = true;
                        switches->current_mode = MODE_NON_UNIFORM_PROBE_IT;
                        required_level = W.BW->current_level + 1;
                        PRINT("[BFS] Requiring level " << required_level);
                        continue;
                    }
                    automorphism.mark = true;
                    break;

                case modes::MODE_WAIT:
                    continue;
            }

            if(switches->done) {
                if(!master)
                    break;
                else
                    continue;
            }

            bool test = true;
            if(switches->done_created_group && automorphism.mark && automorphism.certified) {
                test = G->add_permutation(&automorphism, &idle_ms, done);
                automorphism.not_deletable();
                if(test && foreign_base_done) { //
                    G->sift_random();
                }
            }

            if(!test && !foreign_base_done) {
                // switch this worker to canonical search
                reset_skiplevels(&W);
                foreign_base_done = true;
            }

            sampled_paths += 1;
        }

        //if(switch_map_init) {
         //   switch_map_init = false;
         //   delete[] switch_map;
        //}

        if(master && !dejavu_kill_request) {
            PRINT("[Dej] Cleanup...")
            delete[] shrd_orbit;
            delete[] shrd_orbit_weights;

            delete shrd_orbit_;
            delete shrd_orbit_weights_;

            G->generators_persistent = (gens == nullptr);
            delete G;
            delete bwork;

            delete canon_strategy->I;
            delete canon_strategy->leaf;
            delete canon_strategy;
        }

        if(!is_canon_strategy) {
            delete my_canon_I;
            delete my_canon_leaf;
            delete my_strategy;
        }

        FreeAll();

        return;
    }

    void find_first_leaf(dejavu_workspace<vertex_t, degree_t, edge_t> *w,
                         sgraph_t<vertex_t, degree_t, edge_t> *g, invariant *canon_I,
                         bijection<vertex_t> *canon_leaf, strategy<vertex_t>* canon_strategy,
                         bijection<vertex_t> *automorphism, shared_workspace<vertex_t> *switches,
                         int selector_seed) {
        const bool* done = &switches->done;

        // workspace
        refinement<vertex_t, degree_t, edge_t> *R = &w->R;
        selector<vertex_t, degree_t, edge_t> *S = &w->S;
        coloring<vertex_t> *c = &w->c;
        invariant *I = &w->I;
        coloring<vertex_t> *start_c  = w->start_c;
        invariant *start_I = &w->start_I;

        S->empty_cache();

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

    abort_code base_aligned_search(dejavu_workspace<vertex_t, degree_t, edge_t> *w,
                                   sgraph_t<vertex_t, degree_t, edge_t> *g,
                                   strategy<vertex_t> *canon_strategy, bijection<vertex_t> *automorphism,
                                   strategy_metrics *m, bool *done, shared_workspace<vertex_t> *switches,
                                   int selector_seed) {
        bool backtrack = false;
        bool skipped_level = false;
        bool full_orbit_check = false;

        // workspace
        refinement<vertex_t, degree_t, edge_t> *R = &w->R;
        selector<vertex_t, degree_t, edge_t>   *S = &w->S;
        coloring<vertex_t> *c       = &w->c;
        invariant *I  = &w->I;

        coloring<vertex_t> *start_c = &w->skip_c;
        invariant *start_I     = &w->skip_I;
        shared_schreier *group_level = w->skip_schreier_level;

        invariant* canon_I    = canon_strategy->I;
        bijection<vertex_t>* canon_leaf = canon_strategy->leaf;

        automorphism->non_uniform = false;
        automorphism->certified   = false;
        bool base_aligned = true;
        S->empty_cache();
        m->restarts = 0;
        int level = w->first_skiplevel;
        start_I->set_compare_invariant(canon_I);
        *I = *start_I;
        c->copy_force(start_c);

        if(w->skiplevels > w->my_base_points_sz)
            w->skiplevels = w->my_base_points_sz;

        m->expected_bfs_size = 1;
        m->expected_level = -1;
        S->empty_cache();

        while(w->first_skiplevel <= w->skiplevels) {
            m->expected_bfs_size *= start_c->ptn[start_c->vertex_to_col[w->my_base_points[w->first_skiplevel - 1]]] + 1;
            if(*done) return abort_code(0);
            proceed_state(w, g, start_c, start_I, w->my_base_points[w->first_skiplevel - 1], nullptr, m,
                          (*start_I->compareI->vec_cells)[w->first_skiplevel - 1]);
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
        skipped_level = w->first_skiplevel > 1;
        full_orbit_check = w->skiplevel_is_uniform;

        int it = 0;
        while (true) {
            if(*done) return abort_code(0);

            ++it;
            if(it % 3 == 0) {
                if(switches->current_mode == modes::MODE_TOURNAMENT)
                    switches->check_strategy_tournament(w->id, m, true);
                if(w->id == -1) // but need to be able to reach proper state afterwads
                    w->G->manage_results(switches);
            }

            if (backtrack) {
                if(*done) return abort_code(0);
                if((m->restarts % (5 * switches->tolerance) == ((5 * switches->tolerance) - 1))
                   && (w->skiplevels < w->my_base_points_sz)) {
                    w->skiplevel_is_uniform = false;
                    w->skiplevels += 1;
                }

                S->empty_cache();
                if(w->skiplevels > w->my_base_points_sz)
                    w->skiplevels = w->my_base_points_sz;

                if(w->first_skiplevel <= w->skiplevels) {
                    m->expected_bfs_size *= start_c->ptn[start_c->vertex_to_col[w->my_base_points[w->first_skiplevel - 1]]] + 1;
                    proceed_state(w, g, start_c, start_I, w->my_base_points[w->first_skiplevel - 1],
                            nullptr, m, -1);
                    w->first_skiplevel += 1;
                    if(!w->is_foreign_base) {
                        w->skip_schreier_level = w->skip_schreier_level->next;
                    }
                }

                S->empty_cache();

                // initialize a search state
                m->restarts += 1;
                c->copy_force(start_c);
                *I = *start_I;
                backtrack    = false;

                full_orbit_check = w->skiplevel_is_uniform;
                base_aligned = true;
                level = w->first_skiplevel;
                if(!w->is_foreign_base)
                    group_level = w->skip_schreier_level;
                skipped_level = w->first_skiplevel > 1;
            }

            // extract cell selection from compare invariant, no need to use my own selector
            int s = extract_selector(I, level - 1);
            // int s = S->select_color_dynamic(g, c, canon_strategy);
            if (s == -1) {
                // we can derive an automorphism!
                bijection<vertex_t> leaf;
                leaf.read_from_coloring(c);
                leaf.not_deletable();
                *automorphism = leaf;
                automorphism->inverse();
                automorphism->compose(canon_leaf);
                automorphism->non_uniform = (skipped_level && !w->skiplevel_is_uniform);

                if(full_orbit_check && base_aligned && w->skiplevel_is_uniform && !w->is_foreign_base
                   && switches->current_mode != MODE_TOURNAMENT) {
                    PRINT("[BA] Orbit equals cell abort");
                    return abort_code(1);
                }

                if(!config.CONFIG_IR_FULL_INVARIANT && !R->certify_automorphism(g, automorphism)) {
                    backtrack = true;
                    continue;
                }

                automorphism->certified = true;
                assert(g->certify_automorphism(*automorphism));
                return abort_code(0);
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

            // check if base point can be chosen instead
            if(!w->is_foreign_base) {
                if (group_level->vec[v] && base_aligned) {
                    v = group_level->fixed;// choose base point
                    if(level == w->skiplevels + 1 &&  (w->skiplevels < w->my_base_points_sz)) {
                        bool total_orbit = (c->ptn[s] + 1 == group_level->fixed_orbit_sz);
                        if(total_orbit) {
                            w->skiplevels += 1;
                        } else {
                            full_orbit_check = false;
                        }
                    }
                } else {
                    base_aligned = false;
                }
            }

            const int cell_early = (*I->compareI->vec_cells)[level - 1];
            bool comp = proceed_state(w, g, c, I, v, nullptr, m, cell_early);
            level += 1;

            if(!comp) {
                backtrack = true;
                continue;
            }

            if(!w->is_foreign_base)
                group_level = group_level->next;
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

    abort_code uniform_from_bfs_search(dejavu_workspace<vertex_t, degree_t, edge_t> *w,
                                       sgraph_t<vertex_t, degree_t, edge_t> *g,
                                       strategy<vertex_t>* canon_strategy, bijection<vertex_t> *automorphism,
                                       int *restarts, shared_workspace<vertex_t> *switches, int selector_seed) {
        bool backtrack = false;
        bool* done = &switches->done;

        refinement<vertex_t, degree_t, edge_t>  *R = &w->R;
        selector<vertex_t,   degree_t, edge_t>  *S = &w->S;
        coloring<vertex_t> *c = &w->c;
        invariant *I = &w->I;
        invariant* canon_I    = canon_strategy->I;
        bijection<vertex_t>* canon_leaf = canon_strategy->leaf;

        S->empty_cache();

        *restarts = 0;
        int level;

        automorphism->map    = new vertex_t[g->v_size];
        automorphism->map_sz = 0;

        automorphism->certified = false;

        // pick start from BFS level
        const int bfs_level    = w->BW->current_level - 1;
        const int bfs_level_sz = w->BW->level_sizes[bfs_level];

        shared_schreier* start_group_level = w->G->gp;
        for(int i = 0; i < bfs_level; i++)
            start_group_level = start_group_level->next;
        shared_schreier* group_level = start_group_level;

        int rand_pos       = intRand(0, bfs_level_sz - 1, selector_seed);
        bfs_element<vertex_t>* picked_elem = w->BW->level_states[bfs_level][rand_pos];
        bool base_aligned        = picked_elem->is_identity;

        *I = *picked_elem->I;
        c->copy_force(picked_elem->c);
        I->set_compare_invariant(canon_I);
        level = bfs_level + 1;
        S->empty_cache();
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
                rand_pos    = intRand(0, bfs_level_sz - 1, selector_seed);
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
                level = bfs_level + 1;
                S->empty_cache();
            }

            const int s = extract_selector(I, level - 1);
            //const int s = S->select_color_dynamic(g, c, canon_strategy);
            if (s == -1) {
                // we can derive an automorphism!
                bijection<vertex_t> leaf;
                leaf.read_from_coloring(c);
                leaf.not_deletable();
                *automorphism = leaf;
                automorphism->inverse();
                automorphism->compose(canon_leaf);//enqueue_fail_point_sz
                if(!config.CONFIG_IR_FULL_INVARIANT && !R->certify_automorphism(g, automorphism)) {
                    backtrack = true;
                    continue;
                }
                automorphism->certified = true;
                // assert(g->certify_automorphism(*automorphism));
                return abort_code();
            }

            const int rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
            int v = c->lab[rpos];

            if (group_level->vec[v] && base_aligned) {
                v = group_level->fixed;// choose base point
                if(level == w->skiplevels + 1 &&  (w->skiplevels < w->my_base_points_sz - 1)) {
                    const bool total_orbit = (c->ptn[s] + 1 == group_level->fixed_orbit_sz);
                    if(total_orbit)
                        w->skiplevels += 1;
                }
            } else {
                base_aligned = false;
            }

            group_level = group_level->next;
            level += 1;
            const bool comp = proceed_state(w, g, c, I, v, nullptr, nullptr, -1);

            if(!comp) {
                backtrack = true;
                continue;
            }
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

    bool get_orbit(dejavu_workspace<vertex_t, degree_t, edge_t> *w, int *base, int base_sz, int v,
                   int v_base, work_list *orbit, bool reuse_generators) {
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
                shared_permnode *it = *w->shared_generators;
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

    bool bfs_chunk(dejavu_workspace<vertex_t, degree_t, edge_t> *w,
                   sgraph_t<vertex_t, degree_t, edge_t>  *g, strategy<vertex_t> *canon_strategy,
                   bool *done, int selector_seed) {
        bfs_workspace<vertex_t> *BFS = w->BW;
        int level = BFS->current_level;
        int target_level = BFS->target_level;
        if (level == target_level) return false; // we are done with BFS!

        // initialize bfs_workspace structures
        bfs_assure_init(w);

        if (w->generator_fix_base_alloc < *w->shared_generators_size) {
            delete[] w->generator_fix_base;
            w->generator_fix_base = new shared_permnode *[*w->shared_generators_size];
            w->generator_fix_base_alloc = *w->shared_generators_size;
            w->prev_bfs_element = nullptr;
        }

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

    void bfs_reduce_tree(dejavu_workspace<vertex_t, degree_t, edge_t> *w) {
        thread_local bool first_call = true;

        _schreier_fails(2);
        bfs_assure_init(w);

        int domain_size = w->G->domain_size;
        bool new_gens = sequential_init_copy(w);

        if(!new_gens && !first_call) {
            PRINT("[BFS] Skipping reduce tree (no new generators)");
            return;
        }
        first_call = false;

        if(w->generator_fix_base_alloc < *w->shared_generators_size) {
            delete[] w->generator_fix_base;
            w->generator_fix_base = new shared_permnode*[*w->shared_generators_size];
            w->generator_fix_base_alloc = *w->shared_generators_size;
            w->prev_bfs_element = nullptr;
        }

        // step1: check on level i if orbits can be reduced
        int i, j, v;
        bfs_workspace<vertex_t>* BFS = w->BW;
        int active = 0;
        bool found_canon = false;
        for (j = 0; j < BFS->level_sizes[BFS->current_level - 1]; ++j) {
            bfs_element<vertex_t> *elem = BFS->level_states[BFS->current_level - 1][j];
            active += (elem->weight > 0);
            found_canon = found_canon || (elem->base[elem->base_sz - 1] == w->G->b[elem->base_sz - 1]);
        }
        assert(found_canon);
        PRINT("[BFS] Top level active (before): " << active);


        for(i = 1; i < BFS->current_level; ++i) {
            bfs_element<vertex_t>** arr = w->BW->level_states[i];
            int arr_sz        = w->BW->level_sizes[i];
            std::sort(arr, arr + arr_sz, &bfs_element_parent_sorter<vertex_t>);

            if(i == 1) {
                for(j = 0; j < BFS->level_sizes[i]; ++j) {
                    bfs_element<vertex_t>* elem = BFS->level_states[i][j];
                    if(elem->weight > 0) {
                        v = elem->base[elem->base_sz - 1];
                        assert(v >= 0 && v < w->G->domain_size);
                        if ((*w->shared_orbit)[v] != v) {
                            elem->weight = 0;
                            assert(v != w->G->b[elem->base_sz - 1]);
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
                    bfs_element<vertex_t> *elem = BFS->level_states[i][j];
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
                    bfs_element<vertex_t> *elem = BFS->level_states[i][j];
                    if(elem == nullptr)
                        continue;
                    assert(elem->parent != NULL);
                    if(elem->weight != 0) {
                        int* orbits_sz = nullptr;
                        // depends on earlier sorting
                        int* orbits = _getorbits(elem->base, elem->base_sz - 1, w->sequential_gp, &w->sequential_gens,
                                                 domain_size, w->G->b, &orbits_sz);
#ifdef NDEBUG
                        int calc_sz = 0;
                        for(int ii = 0; ii < domain_size; ++ii) {
                            assert(orbits[ii] >= 0 && orbits[ii] < domain_size);
                            if(ii == orbits[ii]) {
                                calc_sz += (orbits_sz)[ii];
                            }
                        }
                        assert(calc_sz == domain_size);
#endif
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
            bfs_element<vertex_t> *elem = BFS->level_states[BFS->current_level - 1][j];
            active += (elem->weight > 0);
        }

        PRINT("[BFS] Top level active (after): " << active);
        PRINT("[BFS] Emptying queue...");
        moodycamel::ConcurrentQueue<std::tuple<bfs_element<vertex_t> *, int, int>> throwaway_queue;
        BFS->bfs_level_todo[BFS->current_level].swap(throwaway_queue);


        // do not fill this if there is no hope...
        int expecting_finished = 0;
        for (j = 0; j < BFS->level_sizes[BFS->current_level - 1]; ++j) {
            bfs_element<vertex_t> *elem = BFS->level_states[BFS->current_level - 1][j];
            if (elem->weight > 0) {
                int c      = elem->target_color;
                int c_size = elem->c->ptn[c] + 1;
                expecting_finished += c_size;
            }
        }
        BFS->level_expecting_finished[BFS->current_level] = expecting_finished;
        w->prev_bfs_element = nullptr;
        w->orbit.reset();
        w->orbit_vertex_worklist.reset();
        w->orbit_considered.reset();
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
                    if(w->BW->current_level > 1) {
                        orbits_sz = nullptr;
                        orbits = _getorbits(elem->base, elem->base_sz, w->sequential_gp, &w->sequential_gens,
                                            w->G->domain_size, w->G->b, &orbits_sz);
                    } else {
                        orbits_sz = *w->shared_orbit_weights;
                        orbits    = *w->shared_orbit;
                    }
                    for (i = c; i < c + c_size; ++i) {
                        if (orbits[elem->c->lab[i]] == elem->c->lab[i]) {
                            w->BW->bfs_level_todo[w->BW->current_level].enqueue(
                                    std::tuple<bfs_element<vertex_t> *, int, int>(
                                            elem, elem->c->lab[i], orbits_sz[elem->c->lab[i]]));
                            expected += 1;
                            added += 1;
                        }
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
            w->orbit.initialize(w->G->domain_size);
            w->orbit_considered.initialize(w->G->domain_size);
            w->orbit_vertex_worklist.initialize(w->G->domain_size);
            w->generator_fix_base = new shared_permnode*[*w->shared_generators_size];
            w->generator_fix_base_alloc = *w->shared_generators_size;
        }
    }

    bool sequential_init_copy(dejavu_workspace<vertex_t, degree_t, edge_t> *w) {
        // update sequential group in workspace
        bool new_gen = false;
        if(!w->sequential_init) {
            _newgroup(&w->sequential_gp, &w->sequential_gens, w->G->domain_size);
            w->sequential_init = true;
        }

        // copy generators that have not been copied yet
        shared_permnode *it = w->G->gens;
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

    void reconstruct_leaf(dejavu_workspace<vertex_t, degree_t, edge_t> *w,
                          sgraph_t<vertex_t, degree_t, edge_t>  *g, coloring<vertex_t>* start_c, vertex_t* base,
                          int base_sz, bijection<vertex_t> *leaf) {
        coloring<vertex_t>* c = &w->c;
        invariant* I = &w->I;
        c->copy_force(start_c);
        I->never_fail = true;
        for(int pos = 0; pos < base_sz; ++pos) {
            const int v = base[pos];
            proceed_state(w, g, c, I, v, nullptr, nullptr, -1);
        };
        leaf->read_from_coloring(c);
    }

    bool uniform_from_bfs_search_with_storage(dejavu_workspace<vertex_t, degree_t, edge_t> *w,
                                              sgraph_t<vertex_t, degree_t, edge_t>  *g,
                                              shared_workspace<vertex_t>* switches, bfs_element<vertex_t> *elem,
                                              int selector_seed, strategy<vertex_t> *strat,
                                              bijection<vertex_t> *automorphism, bool look_close) {
        thread_local work_list_t<vertex_t> collect_base;
        thread_local int collect_base_sz;

        if(!collect_base.init)
            collect_base.initialize(g->v_size);
        collect_base.reset();
        collect_base_sz = 0;

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
            collect_base.push_back(v);
            collect_base_sz += 1;
            comp = proceed_state(w, g, c, I, v, nullptr, nullptr, -1);
        } while(comp);

        if(comp && (strat->I->acc == I->acc)) { // automorphism computed
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
        } else {
            bijection<vertex_t> leaf;
            leaf.read_from_coloring(c);
            leaf.not_deletable();
            // consider leaf store...
            std::vector<stored_leaf<vertex_t>> pointers;

            switches->leaf_store_mutex.lock();
            auto range = switches->leaf_store.equal_range(I->acc);
            for (auto it = range.first; it != range.second; ++it)
                pointers.push_back(it->second);
            if(pointers.empty()) {
                if (switches->leaf_store_explicit <= config.CONFIG_IR_LEAF_STORE_LIMIT) {
                    switches->leaf_store_explicit++;
                    switches->leaf_store.insert(std::pair<long,
                            stored_leaf<vertex_t>>(I->acc, stored_leaf<vertex_t>(leaf.map, g->v_size, true)));
                } else {
                    vertex_t* base = new vertex_t[collect_base_sz];
                    memcpy(base, collect_base.arr, sizeof(vertex_t) * collect_base_sz);
                    switches->leaf_store.insert(std::pair<long,
                            stored_leaf<vertex_t>>(I->acc, stored_leaf<vertex_t>(base, collect_base_sz, false, elem)));
                    leaf.deletable();
                }
            }
            switches->leaf_store_mutex.unlock();

            comp = false;

            for(size_t i = 0; i < pointers.size(); ++i) {
                automorphism->copy(&leaf);
                automorphism->inverse();
                bijection<vertex_t> fake_leaf;

                if(pointers[i].explicit_leaf) {
                    fake_leaf.map = pointers[i].map;
                    fake_leaf.not_deletable();
                } else {
                    reconstruct_leaf(w, g, pointers[i].start_elem->c, pointers[i].map,pointers[i].map_sz, &fake_leaf);
                    fake_leaf.deletable();
                }

                fake_leaf.map_sz = g->v_size;
                automorphism->compose(&fake_leaf);

                if(w->R.certify_automorphism(g, automorphism)) {
                    automorphism->certified = true;
                    automorphism->non_uniform = false;
                    comp = true;
                } else {
                    comp = false;
                }
            }
        }

        return comp;
    }
};

typedef dejavu_t<int, int, int> dejavu;

void dejavu_automorphisms_dispatch(dynamic_sgraph *sgraph, int* colmap, shared_permnode **gens) {
    switch(sgraph->type) {
        case sgraph_type::DSG_INT_INT_INT: {
            PRINT("[Dispatch] <int32, int32, int32>");
            dejavu_t<int, int, int> d;

            d.automorphisms(sgraph->sgraph_0, colmap, gens);
        }
            break;
        case sgraph_type::DSG_SHORT_SHORT_INT: {
            PRINT("[Dispatch] <int16, int16, int>");
            dejavu_t<int16_t, int16_t, int> d;
            const int v_size = sgraph->sgraph_1->v_size;
            int16_t* colmap_16 = nullptr;
            if(colmap != nullptr) {
                colmap_16 = new int16_t[v_size];
                for (int i = 0; i < v_size; ++i)
                    colmap_16[i] = static_cast<int16_t>(colmap[i]);
            }
            d.automorphisms(sgraph->sgraph_1, colmap_16, gens);
        }
            break;
        case sgraph_type::DSG_SHORT_SHORT_SHORT: {
            PRINT("[Dispatch] <int16, int16, int16>");
            dejavu_t<int16_t, int16_t, int16_t> d;
            const int v_size = sgraph->sgraph_2->v_size;
            int16_t* colmap_16 = nullptr;
            if(colmap != nullptr) {
                colmap_16 = new int16_t[v_size];
                for (int i = 0; i < v_size; ++i)
                    colmap_16[i] = static_cast<int16_t>(colmap[i]);
            }
            d.automorphisms(sgraph->sgraph_2, colmap_16, gens);
        }
            break;
        case sgraph_type::DSG_CHAR_CHAR_SHORT:{
            PRINT("[Dispatch] <int8, int8, int16>");
            dejavu_t<int8_t, int8_t, int16_t> d;
            const int v_size = sgraph->sgraph_3->v_size;
            int8_t* colmap_8 = nullptr;
            if(colmap != nullptr) {
                colmap_8 = new int8_t[v_size];
                for (int i = 0; i < v_size; ++i)
                    colmap_8[i] = static_cast<int8_t>(colmap[i]);
            }
            d.automorphisms(sgraph->sgraph_3, colmap_8, gens);
        }
            break;
        case sgraph_type::DSG_CHAR_CHAR_CHAR: {
            PRINT("[Dispatch] <int8, int8, int8>");
            dejavu_t<int8_t, int8_t, int8_t> d;
            const int v_size = sgraph->sgraph_4->v_size;
            int8_t* colmap_8 = nullptr;
            if(colmap != nullptr) {
                colmap_8 = new int8_t[v_size];
                for (int i = 0; i < v_size; ++i)
                    colmap_8[i] = static_cast<int8_t>(colmap[i]);
            }
            d.automorphisms(sgraph->sgraph_4, colmap_8, gens);
        }
            break;
    }
}

void dejavu_automorphisms(sgraph_t<int, int, int> *g, int* colmap, shared_permnode **gens) {
    if(config.CONFIG_PREPROCESS_COMPRESS) {
        dynamic_sgraph sg;
        dynamic_sgraph::read(g, &sg);
        dejavu_automorphisms_dispatch(&sg, colmap, gens);
    } else {
        dejavu d;
        d.automorphisms(g, colmap, gens);
    }
}

#endif //DEJAVU_DEJAVU_H
