#ifndef DEJAVU_DEJAVU_AUTO_H
#define DEJAVU_DEJAVU_AUTO_H

#include <random>
#include <chrono>
#include "sgraph.h"
#include "invariant.h"
#include "concurrentqueue.h"
#include "sassy/preprocessor.h"
#include "selector.h"
#include "bfs.h"
#include "schreier_shared.h"
#include "group_shared.h"
#include "schreier_sequential.h"
#include "dfs.h"

typedef std::chrono::high_resolution_clock Clock;

struct abort_code {
    abort_code()=default;
    abort_code(int reason):reason(reason){};
    int reason = 0;
};

struct dejavu_workspace {
    // workspace for normal search
    refinement R;
    selector   S;
    coloring c;
    invariant  I;

    // workspace for bfs_workspace
    coloring*  work_c;
    invariant* work_I;

    // workspace for base aligned search
    int base_size = 0;
    int skiplevels = 1;
    int first_skiplevel = 1;
    coloring           skip_c;
    invariant          skip_I;
    shared_schreier*   skip_schreier_level;
    bool               skiplevel_is_uniform = false;

    bijection my_base_points; // this should be bijection
    bool      is_foreign_base;

    group_shared* G;

    coloring start_c1;
    coloring start_c2;
    invariant start_I;

    // indicates which thread this is
    int id;

    // shared orbit and generators
    int**             shared_orbit;
    int**             shared_orbit_weights;
    shared_permnode** shared_generators;
    int*              shared_generators_size;
    int               generator_fix_base_alloc = -1;

    // sequential, local group
    sequential_permnode*      sequential_gens;
    sequential_schreierlevel* sequential_gp;
    bool                      sequential_init = false;

    // deprecated workspace for simple orbit method
    work_set  orbit_considered;
    work_list orbit_vertex_worklist;
    work_list orbit;
    int canonical_v;
    shared_permnode** generator_fix_base;
    int         generator_fix_base_size;

    bool checked_tournament = false;
    bool seen_done_shared   = false;

    // bfs workspace
    std::tuple<bfs_element*, int, int>* todo_dequeue;
    std::pair<bfs_element *, int>* finished_elements;
    std::pair<bfs_element *, int>* todo_elements;
    bfs_element* prev_bfs_element = nullptr;
    bool init_bfs = false;

    work_list_t<int> _collect_base;
    std::vector<int> _collect_early_individualized;
    int _collect_base_sz;

    dejavu_workspace() {
        work_c = new coloring;
        work_I = new invariant;
    }
    ~dejavu_workspace() {
        if(init_bfs) {
            delete[] todo_dequeue;
            delete[] todo_elements,
            delete[] finished_elements;
            delete[] generator_fix_base;
        }
        if(sequential_init) {
            _freeschreier(&sequential_gp, &sequential_gens);
            _schreier_freedyn();
        }

        delete work_c;
        delete work_I;
        //delete start_c1;
    };
};

bool bfs_element_parent_sorter(bfs_element* const& lhs, bfs_element* const& rhs) {
    if(lhs->parent < rhs->parent)
        return true;
    if(lhs->parent == rhs->parent) {
        return(lhs->parent->parent < rhs->parent->parent);
    }
    return false;
}

class dejavu_auto_t {
public:
    dejavu_stats automorphisms(sgraph *g, int* colmap, shared_permnode **gens) {
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

        coloring start_c;
        start_c.vertex_to_col = colmap;
        start_c.init = false;
        shared_workspace_auto switches;
        bfs_workspace BW;

        return worker_thread(g, true, &switches, nullptr, &start_c, nullptr, -1,
                      nullptr, nullptr, &BW, gens, nullptr);
    }

private:
    dejavu_stats worker_thread(sgraph *g_, bool master,
                               shared_workspace_auto *switches, group_shared *G,
                               coloring *start_c, strategy* canon_strategy,
                               int communicator_id, int **shared_orbit, int** shared_orbit_weights,
                               bfs_workspace *bwork, shared_permnode **gens, int *shared_group_size) {
        // first order of business: try to glue this thread to currently assigned core
        #ifndef OS_WINDOWS
        int master_sched;
        if(master) {
            master_sched = sched_getcpu();
            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(master_sched, &cpuset);
            pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
        }
        #endif

        dejavu_stats a = dejavu_stats();
        sgraph *g = g_;
        dejavu_workspace W;
        strategy _canon_strategy;


        numnodes  = 0;
        colorcost = 0;
	    coloring* init_c;
        // preprocessing
        if(master) {
            init_c  = start_c;
            start_c = &W.start_c1;
            g->initialize_coloring(&W.start_c1, init_c->vertex_to_col);
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
        bijection base_points;
        bijection actual_base;
        int trash_int = 0;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() * ((communicator_id * 5) * 5135235);
        int selector_seed = seed;

        invariant start_I;

        // first color refinement, initialize some more shared structures, launch threads
        if (master) {
            PRINT("[Dej] Dense graph: " << (g->dense?"true":"false"));
            switches->current_mode = modes_auto::MODE_AUTO_TOURNAMENT;

            // first color refinement
            canon_strategy = &_canon_strategy;
            if(!config.CONFIG_IR_SKIP_FIRST_REFINEMENT)
                W.R.refine_coloring_first(g, start_c, -1);
            PRINT("[Dej] First refinement: " << cref / 1000000.0 << "ms");

            if(start_c->cells == g->v_size) {
                *done = true;
                std::cout << "First coloring discrete." << std::endl;
                std::cout << "Base size: 0" << std::endl;
                std::cout << "Group size: 1" << std::endl;
                a.gens = nullptr;
                a.grp_sz_exp = 0;
                a.grp_sz_man = 1;
                return a;
            }

            dejavu::dfs_ir D;
            D.setup(0, &W.R, nullptr);
            const int dfs_reached_level = D.do_dfs(g, start_c);
            if(dfs_reached_level == 0) {
                long double add_grp_sz = D.grp_sz_man;
                while(add_grp_sz > 10) {
                    a.grp_sz_exp += 1;
                    add_grp_sz = add_grp_sz / 10;
                }
                a.grp_sz_exp += D.grp_sz_exp;
                a.grp_sz_man *= add_grp_sz;
                std::cout << "DFS finished graph" << std::endl;
                return a;
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
            G    = new group_shared(g->v_size);

            W.S.empty_cache();

            // launch worker threads
            for (int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++) {
                work_threads.emplace_back(
                        std::thread(&dejavu_auto_t::worker_thread,
                                    dejavu_auto_t(), g, false, switches, G, start_c,
                                    canon_strategy, i, shrd_orbit_, shrd_orbit_weights_, bwork, &_gens,
                                    &_shared_group_size));
                #ifndef OS_WINDOWS
                cpu_set_t cpuset;
                CPU_ZERO(&cpuset);
                CPU_SET(i + (i >= master_sched), &cpuset);
                pthread_setaffinity_np(work_threads[i].native_handle(), sizeof(cpu_set_t), &cpuset);
                #endif
            }

            PRINT("[Dej] Refinement workers created (" << config.CONFIG_THREADS_REFINEMENT_WORKERS << " threads)");

            // set some workspace variables
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
        W.G = G;

        W.skiplevels = 0;
        W.skip_schreier_level = G->gp;

        if (!master) {
            W.shared_orbit = shared_orbit;
            W.shared_orbit_weights = shared_orbit_weights;
            W.start_c1.copy_force(start_c);
            W.id = communicator_id;
            W.shared_generators = gens;
            W.shared_generators_size = shared_group_size;
        }

        strategy* my_strategy;
        bijection _my_canon_leaf;

        invariant _my_canon_I;
        my_canon_I = &_my_canon_I;
        my_canon_I->has_compare = false;
        my_canon_I->compare_vec = nullptr;
        my_canon_I->compareI    = nullptr;
        my_canon_leaf = &_my_canon_leaf;
        auto rst = (selector_type) ((communicator_id + 2) % 3);
        if(config.CONFIG_IR_FORCE_SELECTOR)
            rst = (selector_type) config.CONFIG_IR_CELL_SELECTOR;

        strategy _my_strategy(my_canon_leaf, my_canon_I, rst, -1);
        my_strategy = &_my_strategy;

        W.S.empty_cache();
        find_first_leaf(&W, g, my_canon_I, my_canon_leaf, my_strategy, &base_points, switches, selector_seed);
        W.my_base_points.swap(&base_points);

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
                    switches->current_mode = modes_auto::MODE_AUTO_WAIT;
                    G->ack_done_shared(&W.seen_done_shared);
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

                    if(bwork->current_level > 1) { // > 1
                        if(*W.shared_generators_size > 0) {
                            //PRINT("[BFS] Reducing tree (" << n_found << ")");
                            assert(master && communicator_id == -1);
                            //bfs_reduce_tree(&W, bwork); // removed feature
                        }
                        // check if expected size is still too large...
                        if(bwork->level_expecting_finished[bwork->current_level] >=
                           config.CONFIG_IR_SIZE_FACTOR * g->v_size * switches->tolerance) {
                            PRINT("[BFS] Expected size still too large, not going into BFS")
                            bwork->reached_initial_target = (bwork->target_level == bwork->current_level);
                            bwork->target_level.store(bwork->current_level);
                        } else {
                            switches->reset_tolerance(bwork->level_expecting_finished[bwork->current_level], g->v_size);
                            PRINT("[Dej] Tolerance: " << switches->tolerance);
                            PRINT("[BFS] Filling queue..." << bwork->current_level << " -> " << bwork->target_level)
                            bfs_fill_queue(&W, bwork);
                        }
                    } else {
                        if(*W.shared_generators_size > 0) {
                            // PRINT("[BFS] Reducing queue using orbits")
                            //    bfs_fill_queue(&W);
                        }
                    }

                    switches->done_shared_group = true;
                    switches->current_mode      = modes_auto::MODE_AUTO_BFS;
                }

                // we are done
                if(switches->done) {
                    if(!dejavu_kill_request) {
                        PRINT("Numnodes (master): " << numnodes);
                        PRINT("Colorcost (master): " << colorcost);

                        PRINT("Base size:  " << G->base_size);
                        while (!work_threads.empty()) {
                            work_threads[work_threads.size() - 1].join();
                            work_threads.pop_back();
                        }
                        if(config.CONFIG_IR_WRITE_GROUPORDER) {
                            std::cout << "Group size: ";
                            G->print_group_size_stdout();
                        }

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
                        }
                        if(config.CONFIG_WRITE_AUTOMORPHISMS_GAP) {
                            std::cout << "g := Group(";
                            shared_permnode *it = G->gens;
                            mark_set touched;
                            touched.initialize(g->v_size);

                            if(it != nullptr) {
                                do {
                                    touched.reset();
                                    for (int i = 0; i < g->v_size; ++i) {
                                        if(it->p[i] == i)
                                            continue;
                                        if(touched.get(i))
                                            continue;
                                        std::cout << "(";
                                        int j = i;
                                        int length = 0;
                                        do {
                                            touched.set(j);
                                            std::cout << j+1;
                                            j = it->p[j];
                                            ++length;
                                            if(j != i)
                                                std::cout << ",";
                                        } while(j != i);
                                        std::cout << ")";
                                        if(length < 2)
                                            std::cout << "error" << std::endl;
                                    }

                                    it = it->next;
                                    if(it != G->gens) {
                                        std::cout << ", " << std::endl;
                                    }
                                } while (it != G->gens);
                            }
                            std::cout << ");" << std::endl;
                        }
                        a.gens = gens;
                        a.base = std::vector<int>(G->b, G->b+G->base_size);
                        G->compute_group_size(&a.grp_sz_man, &a.grp_sz_exp);
                        if(gens != nullptr) {
                            *gens = G->gens;
                        }
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
                G->ack_done_shared(&W.seen_done_shared);
                reset_non_uniform_switch = true;
            }

            if(switches->done) {
                if(!master)
                    break;
                else
                    continue;
            }

            bijection automorphism;
            abort_code A;
            automorphism.mark = false;

            // high-level solver mode
            switch(switches->current_mode) {
                // In the tournament mode, threads compete for the best strategy (cell selector, target leaf).
                // After the tournament, the best strategy will be chosen as the canonical strategy (canon_strategy).
                case modes_auto::MODE_AUTO_TOURNAMENT:
                    m.restarts = 0;
                    m.expected_bfs_size = 0;
                    W.skiplevel_is_uniform = false;
                    base_aligned_search(&W, g, my_strategy, &automorphism, &m, done_fast, switches, selector_seed);
                    // check if this thread won
                    if(n_found == 0) {
                        // wait until everyone checked
                        while(!switches->check_strategy_tournament(communicator_id, &m, false,
                                                                   &W.checked_tournament)
                              && !switches->done_created_group) continue;
                        // if this thread won, it will now create the group
                        if(switches->win_id == communicator_id) {
                            canon_strategy->replace(my_strategy);
                            is_canon_strategy = true;
                            actual_base.copy(&W.my_base_points);
                            G->initialize(g->v_size, &actual_base);
                            PRINT("[Strat] Chosen strategy: " << canon_strategy->cell_selector_type);
                            W.S.empty_cache();
                            int init_c = W.S.select_color_dynamic(g, start_c, my_strategy);
                            bwork->initialize(bfs_element::root_element(start_c, &start_I), init_c,
                                             g->v_size, G->base_size);
                            int proposed_level = W.skiplevels + 1;
                            if(config.CONFIG_IR_FULL_BFS)
                                proposed_level = G->base_size + 1;
                            bwork->target_level.store(proposed_level);
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
                            switches->current_mode = modes_auto::MODE_AUTO_NON_UNIFORM_PROBE;
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

                // In this mode, base-aligned (non-uniform) search is performed.
                case modes_auto::MODE_AUTO_NON_UNIFORM_PROBE:
                    // This heuristic balances finding automorphisms aligned to the strategy of this thread (automorphisms
                    // might be more distinct but might cause higher sifting cost) as opposed to the canonical strategy
                    // (automorphisms of threads might be quite similar, but sifting cost is reduced).
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

                // Same as "MODE_AUTO_NON_UNIFORM_PROBE", but after breadth-first search and uniform probing has been performed.
                case modes_auto::MODE_AUTO_NON_UNIFORM_PROBE_IT:
                    // automorphism search from initial bfs_workspace pieces
                    if(!*done_fast) {
                        if (reset_non_uniform_switch) {
                            reset_skiplevels(&W);
                            if(!master) { // guess a new leaf
                                _my_canon_I    = invariant();
                                _my_canon_leaf = bijection();
                                auto rst      = (selector_type) intRand(0, 2, selector_seed);
                                _my_strategy   = strategy(my_canon_leaf, my_canon_I, rst, -1);
                                is_canon_strategy = false;
                                base_points   = bijection();
                                find_first_leaf(&W, g, my_canon_I, my_canon_leaf, my_strategy, &base_points, switches,
                                                selector_seed);
                                W.my_base_points.swap(&base_points);
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
                            if (proposed_level > bwork->target_level)
                                bwork->target_level.store(proposed_level);
                        }

                        n_found += 1;
                        if ((*done_fast && !automorphism.non_uniform)) continue;
                    } else continue;
                    if(*done && master) {
                        continue;
                    }
                    if ((*done_fast && !automorphism.non_uniform)) continue;
                    break;
                // Probes paths uniformly from the current breadth-first level and stores additional leaves.
                case modes_auto::MODE_AUTO_UNIFORM_WITH_LEAF_STORAGE:
                {
                    // pick initial path from BFS level that is allocated to me
                    --switches->experimental_budget;
                    if(master && (switches->experimental_budget <= 0 || switches->done_fast)) {
                        if(!switches->done_fast) {
                            if(switches->experimental_paths > switches->experimental_deviation) {
                                if(!switches->experimental_look_close) {
                                    switches->experimental_look_close = true;
                                    switches->experimental_budget += bwork->level_sizes[bwork->current_level - 1];
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
                        switches->current_mode = modes_auto::MODE_AUTO_WAIT;
                        switches->done_fast = true;
                        switches->done_shared_group = false;
                        continue;
                    }

                    bfs_element *elem;
                    int bfs_level    = bwork->current_level - 1;
                    int max_weight   = bwork->level_maxweight[bfs_level];
                    int bfs_level_sz = bwork->level_sizes[bfs_level];
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
                            if (proposed_level > bwork->target_level) {
                                bwork->target_level.store(proposed_level);
                            }
                        }
                    }
                    double picked_weight, rand_weight;
                    do {
                        int pick_elem = intRand(0, bfs_level_sz - 1, selector_seed);
                        elem = bwork->level_states[bfs_level][pick_elem];
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

                // Performs breadth-first search.
                case modes_auto::MODE_AUTO_BFS:
                    reset_non_uniform_switch = true;
                    if(W.is_foreign_base) {
                        reset_skiplevels(&W);
                        foreign_base_done = true;
                    }
                    if(bwork->current_level != bwork->target_level) {
                        if (communicator_id == -1 && bwork->target_level < 0) {
                            int proposed_level = W.skiplevels + 1;
                            if (proposed_level == G->base_size)
                                proposed_level += 1;
                            if (proposed_level > G->base_size + 1)
                                proposed_level = G->base_size + 1;
                            bwork->target_level.store(proposed_level);
                        }
                        if(switches->done_shared_group && bwork->target_level >= 0) {
                            if(master && !switched1) {
                                switched1 = true;
                                PRINT("[BA] Finished non-uniform automorphism search (" << *W.shared_generators_size
                                                                                          << " generators, " << n_restarts << " restarts)")
                                PRINT("[BA] Ended in skiplevel " << W.skiplevels << ", found " << n_found)
                                PRINT("[BA] " << cref / 1000000.0 << "ms")
                                PRINT("[BFS] Determined target level: " << bwork->target_level << "")
                            }
                            bfs_chunk(&W, g, canon_strategy, bwork, done, selector_seed);
                            if(master) {
                                bool fill = bwork->work_queues(switches->tolerance);
                                if(fill)
                                    bfs_fill_queue(&W, bwork);
                            }
                        }
                    } else {
                        if(master) {
                            bool fill = bwork->work_queues(switches->tolerance);
                            if(fill)
                                bfs_fill_queue(&W, bwork);
                            if(bwork->reached_initial_target && !config.CONFIG_IR_ALWAYS_STORE) { // should never go without leaf storage?
                                // reached the desired target level? go to next phase!
                                if(bwork->current_level - 1 >= 0 && bwork->level_sizes[bwork->current_level - 1] == 1
                                                                 && bwork->current_level == bwork->base_size + 1) {
                                    PRINT("[BFS] Early-out");
                                    *done = true;
                                    continue;
                                }
                                PRINT("[UTarget] Starting uniform probe, tolerance: " << switches->tolerance)
                                switches->current_mode = modes_auto::MODE_AUTO_UNIFORM_PROBE;
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
                                if(config.CONFIG_IR_ALWAYS_STORE) {
                                    required_level = bwork->current_level + 1;
                                    PRINT("[BFS] Requiring level " << required_level);
                                }
                                int budget_fac = switches->experimental_look_close?std::max(switches->tolerance, 10):1;

                                PRINT("[UStore] Switching to uniform with leaf storage, budget "
                                              << bwork->level_sizes[bwork->current_level - 1] * budget_fac)
                                switches->experimental_budget.store(bwork->level_sizes[bwork->current_level - 1] * budget_fac);
                                switches->experimental_paths.store(0);
                                switches->experimental_deviation.store(0);
                                switches->current_mode = modes_auto::MODE_AUTO_UNIFORM_WITH_LEAF_STORAGE;
                                continue;
                            }
                        }
                    }
                    continue;
                    break;

                // Performs uniform probing without additional leaf storage.
                case modes_auto::MODE_AUTO_UNIFORM_PROBE:
                    reset_non_uniform_switch = true;
                    if(W.id == 0 && !switched2) {
                        switched2 = true;
                        PRINT("[UTarget] " << cref / 1000000.0 << "ms")
                    }
                    A = uniform_from_bfs_search(&W, g, canon_strategy, &automorphism, &restarts, switches, bwork, selector_seed);
                    if(A.reason == 2) // abort
                        continue;

                    if(A.reason == 1) { // too many restarts
                        switches->current_mode = MODE_AUTO_WAIT;
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
                        switches->current_mode = MODE_AUTO_NON_UNIFORM_PROBE_IT;
                        required_level = bwork->current_level + 1;
                        PRINT("[BFS] Requiring level " << required_level);
                        continue;
                    }
                    automorphism.mark = true;
                    break;

                // Do nothing, for now.
                case modes_auto::MODE_AUTO_WAIT:
                    continue;
            }

            if(switches->done) {
                if(!master)
                    break;
                else
                    continue;
            }

            // If the current iteration of the main mode produced an automorphism, we want to sift it into the shared
            // Schreier structure.
            bool test = true;
            if(switches->done_created_group && automorphism.mark && automorphism.certified) {
                test = G->add_permutation(&automorphism, &idle_ms, done);
                if(test && foreign_base_done && !(switches->current_mode == MODE_AUTO_NON_UNIFORM_PROBE)
                        && !(switches->current_mode == MODE_AUTO_NON_UNIFORM_PROBE_IT)) { //
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

        if(master && !dejavu_kill_request) {
            PRINT("[Dej] Cleanup...")
            delete[] shrd_orbit;
            delete[] shrd_orbit_weights;

            delete shrd_orbit_;
            delete shrd_orbit_weights_;

            G->generators_persistent = (gens != nullptr);
            delete G;

            garbage_collector<int>::free_trash();
        }

        return a;
    }

    // Probes a random leaf of the tree and writes down an invariant. The invariant will be utilized for blueprints and
    // other comparisons.
    void find_first_leaf(dejavu_workspace *w,
                         sgraph *g, invariant *canon_I,
                         bijection *canon_leaf, strategy* canon_strategy,
                         bijection *automorphism, shared_workspace_auto *switches,
                         int selector_seed) {
        const bool* done = &switches->done;

        // workspace
        refinement *R = &w->R;
        selector *S   = &w->S;
        coloring *c        = &w->c;
        invariant *I                 = &w->I;
        coloring *start_c  = &w->start_c1;
        invariant *start_I           = &w->start_I;

        S->empty_cache();

        start_I->create_vector(g->v_size * 2);
        //if(automorphism->init)
        //    delete[] automorphism->map;
        automorphism->initialize_empty(g->v_size);
        //automorphism->map    = new vertex_t[g->v_size];
        //automorphism->map_sz = 0;
        //automorphism->deletable();

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
            assert(rpos >= 0);
            assert(rpos < g->v_size);
            const int v = c->lab[rpos];
            assert(v >= 0);
            assert(v < g->v_size);

            // individualize and refine
            assert(c->vertex_to_col[c->lab[s]] == s);
            assert(c->ptn[s] > 0 && c->ptn[s] < c->ptn_sz);
            assert(c->vertex_to_col[v] == s);
            proceed_state(w, g, c, I, v, nullptr, -1, -1, nullptr);
            assert(c->vertex_to_col[v] > 0);

            // base point
            automorphism->append(v);
            //automorphism->map[automorphism->map_sz] = v;
            //automorphism->map_sz += 1;
        }
    }

    abort_code base_aligned_search(dejavu_workspace *w,
                                   sgraph *g,
                                   strategy *canon_strategy, bijection *automorphism,
                                   strategy_metrics *m, bool *done, shared_workspace_auto *switches,
                                   int selector_seed) {
        bool backtrack = false;
        bool skipped_level = false;
        bool full_orbit_check = false;

        // workspace
        refinement *R = &w->R;
        selector   *S = &w->S;
        coloring   *c = &w->c;
        invariant  *I = &w->I;

        coloring *start_c  = &w->skip_c;
        invariant *start_I           = &w->skip_I;
        shared_schreier *group_level = w->skip_schreier_level;

        invariant* canon_I              = canon_strategy->I;
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

        if(w->skiplevels > w->my_base_points.map_sz)
            w->skiplevels = w->my_base_points.map_sz;

        m->expected_bfs_size = 1;
        m->expected_level = -1;
        S->empty_cache();

        while(w->first_skiplevel <= w->skiplevels) {
            m->expected_bfs_size *= start_c->ptn[start_c->vertex_to_col[w->my_base_points.map_vertex(w->first_skiplevel - 1)]] + 1;
            if(*done) return abort_code(0);
            proceed_state(w, g, start_c, start_I, w->my_base_points.map_vertex(w->first_skiplevel - 1), m,
                          (start_I->compareI->vec_cells)[w->first_skiplevel - 1], -1, nullptr);
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
        skipped_level    = w->first_skiplevel > 1;
        full_orbit_check = w->skiplevel_is_uniform;

        int it = 0;
        while (true) {
            if(*done) return abort_code(0);

            ++it;
            if(it % 3 == 0) {
                if(switches->current_mode == modes_auto::MODE_AUTO_TOURNAMENT)
                    switches->check_strategy_tournament(w->id, m, true, &w->checked_tournament);
                if(w->id == -1) // but need to be able to reach proper state afterwads
                    w->G->manage_results(switches);
            }

            if (backtrack) {
                if(*done) return abort_code(0);
                if((m->restarts % (5 * switches->tolerance) == ((5 * switches->tolerance) - 1))
                   && (w->skiplevels < w->my_base_points.map_sz)) {
                    w->skiplevel_is_uniform = false;
                    w->skiplevels += 1;
                }

                S->empty_cache();
                if(w->skiplevels > w->my_base_points.map_sz)
                    w->skiplevels = w->my_base_points.map_sz;

                if(w->first_skiplevel <= w->skiplevels) {
                    m->expected_bfs_size *= start_c->ptn[start_c->vertex_to_col[w->my_base_points.map_vertex(w->first_skiplevel - 1)]] + 1;
                    proceed_state(w, g, start_c, start_I, w->my_base_points.map_vertex(w->first_skiplevel - 1), m, -1,
                                  -1, nullptr);
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
                automorphism->read_from_coloring(c);
                automorphism->inverse();
                automorphism->compose(canon_leaf);
                automorphism->non_uniform = (skipped_level && !w->skiplevel_is_uniform);
                automorphism->non_uniform = true;
                if(full_orbit_check && base_aligned && w->skiplevel_is_uniform && !w->is_foreign_base
                   && switches->current_mode != MODE_AUTO_TOURNAMENT) {
                    PRINT("[BA] Orbit equals cell abort");
                    return abort_code(1);
                }

                if(!config.CONFIG_IR_FULL_INVARIANT && !R->certify_automorphism(g, automorphism)) {
                    backtrack = true;
                    continue;
                }

                automorphism->certified = true;
                //assert(g->certify_automorphism(automorphism)); // warning this assert is broken! (copy semantics!)
                return abort_code(0);
            }

            // individualize and refine now
            int rpos, v;

            if(level <= w->skiplevels) {
                skipped_level = true;
                v = w->my_base_points.map_vertex(level - 1);
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
                    if(level == w->skiplevels + 1 &&  (w->skiplevels < w->my_base_points.map_sz)) {
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

            const int cell_early = (I->compareI->vec_cells)[level - 1];
            bool comp = proceed_state(w, g, c, I, v, m, cell_early, -1, nullptr);
            level += 1;

            if(!comp) {
                backtrack = true;
                continue;
            }

            if(!w->is_foreign_base)
                group_level = group_level->next;
        }
    }

    void reset_skiplevels(dejavu_workspace *w) {
        w->skip_c.copy_force(&w->start_c1);
        w->skip_I = w->start_I;
        w->skiplevels = 0;
        w->skip_schreier_level = w->G->gp;
        w->first_skiplevel = 1;
        w->skiplevel_is_uniform = true;
        //w->my_base_points    = w->G->b;
        //w->my_base_points_sz = w->G->base_size;
        w->my_base_points.read_from_array(w->G->b, w->G->base_size);
        w->is_foreign_base   = false;
    }

    abort_code uniform_from_bfs_search(dejavu_workspace *w,
                                       sgraph *g,
                                       strategy* canon_strategy, bijection *automorphism,
                                       int *restarts, shared_workspace_auto *switches,
                                       bfs_workspace *bwork, int selector_seed) {
        bool backtrack = false;
        bool* done = &switches->done;

        refinement *R = &w->R;
        selector   *S = &w->S;
        coloring *c   = &w->c;
        invariant *I  = &w->I;
        invariant* canon_I                         = canon_strategy->I;
        bijection* canon_leaf            = canon_strategy->leaf;

        S->empty_cache();

        *restarts = 0;
        int level;

        automorphism->certified = false;

        // pick start from BFS level
        const int bfs_level    = bwork->current_level - 1;
        const int bfs_level_sz = bwork->level_sizes[bfs_level];

        shared_schreier* start_group_level = w->G->gp;
        for(int i = 0; i < bfs_level; i++)
            start_group_level = start_group_level->next;
        shared_schreier* group_level = start_group_level;

        int rand_pos       = intRand(0, bfs_level_sz - 1, selector_seed);
        bfs_element* picked_elem = bwork->level_states[bfs_level][rand_pos];
        bool base_aligned        = picked_elem->is_identity;

        *I = *picked_elem->I;
        c->copy_force(picked_elem->c);
        I->set_compare_invariant(canon_I);
        level = bfs_level + 1;
        S->empty_cache();
        double picked_weight, max_weight, rand_weight;

        while (true) {
            if(*done) return abort_code();
            if(switches->current_mode != modes_auto::MODE_AUTO_UNIFORM_PROBE) return abort_code(2);
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
                picked_elem = bwork->level_states[bfs_level][rand_pos];
                group_level = start_group_level;
                base_aligned = picked_elem->is_identity;

                // consider the weight by redrawing
                picked_weight = picked_elem->weight;
                max_weight    = bwork->level_maxweight[bfs_level];
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
                automorphism->read_from_coloring(c);
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
                if(level == w->skiplevels + 1 &&  (w->skiplevels < w->my_base_points.map_sz - 1)) {
                    const bool total_orbit = (c->ptn[s] + 1 == group_level->fixed_orbit_sz);
                    if(total_orbit)
                        w->skiplevels += 1;
                }
            } else {
                base_aligned = false;
            }

            group_level = group_level->next;
            level += 1;
            const bool comp = proceed_state(w, g, c, I, v, nullptr, -1, -1, nullptr);

            if(!comp) {
                backtrack = true;
                continue;
            }
        }
    }

    //  Extracts which color was chosen in the path of the given invariant at a given base point.
    int extract_selector(invariant* I, int base_point) {
        if((I->compareI->vec_selections).size() <= base_point)
            return - 1;
        return (I->compareI->vec_selections)[base_point].first;
    }

    // Performs one individualization followed by one refinement on the given coloring with the given settings.
    bool proceed_state(dejavu_workspace* w,
                       sgraph* g, coloring* c,
                       invariant* I, int v, strategy_metrics* m, int cell_early, int individualize_early,
                       std::vector<int>* early_individualized) {
        if(!config.CONFIG_IR_IDLE_SKIP)
            cell_early = -1;

        // protocol of selector choice
        assert(v >= 0);
        assert(v < g->v_size);
        assert(c->vertex_to_col[v] >= 0);
        assert(c->vertex_to_col[v] < g->v_size);
        I->selection_write(c->vertex_to_col[v], c->ptn[c->vertex_to_col[v]]);
        const int init_color_class = w->R.individualize_vertex(c, v); // TODO: should allow vector of vertices to be individualized
        bool comp = true;
        comp = comp && w->R.refine_coloring(g, c, I, init_color_class, m, cell_early, individualize_early,
                                            early_individualized, nullptr, nullptr);
        return comp;
    }

    bool get_orbit(dejavu_workspace *w, int *base, int base_sz, int v,
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

    bool bfs_chunk(dejavu_workspace *w,
                   sgraph *g, strategy *canon_strategy,
                   bfs_workspace* bwork, bool *done, int selector_seed) {
        int level = bwork->current_level;
        int target_level = bwork->target_level;
        if (level == target_level) return false; // we are done with BFS!

        // initialize bfs_workspace structures
        bfs_assure_init(w, bwork);

        if (w->generator_fix_base_alloc < *w->shared_generators_size) {
            delete[] w->generator_fix_base;
            w->generator_fix_base = new shared_permnode *[*w->shared_generators_size];
            w->generator_fix_base_alloc = *w->shared_generators_size;
            w->prev_bfs_element = nullptr;
        }

        // try to dequeue a chunk of work
        size_t num = bwork->bfs_level_todo[level].try_dequeue_bulk(w->todo_dequeue, bwork->chunk_size);
        int finished_elements_sz = 0;
        int finished_elements_null_buffer = 0;

        for (size_t i = 0; i < num; ++i) {
            bfs_element *elem = std::get<0>(w->todo_dequeue[i]);
            int v      = std::get<1>(w->todo_dequeue[i]);
            int weight = std::get<2>(w->todo_dequeue[i]);
            bool is_identity  = elem->is_identity && (v == w->my_base_points.map_vertex(elem->base_sz));

            // check orbit
            bool comp = elem->weight != 0;
            comp = comp && (elem->deviation_vertex != v);
            if(elem->deviation_pos > 0) {
                // check in abort map
                if (!elem->is_identity) {
                    bool comp_ = bwork->read_abort_map(level, elem->deviation_pos, elem->deviation_acc);

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
                //comp = comp && get_orbit(w, elem->base, elem->base_sz, v, w->my_base_points.map_vertex(elem->base_sz), &w->orbit,
                //                         w->prev_bfs_element == elem);
                weight = 1;
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
                const int cells_early = (w->work_I->compareI->vec_cells)[elem->base_sz];
                comp = comp && proceed_state(w, g, w->work_c, w->work_I, v, nullptr, cells_early, -1,
                                             nullptr);

                // manage abort map counter
                if (comp && elem->is_identity && level > 1) {
                    // decrease abort map done...
                    bwork->level_abort_map_mutex[level]->lock();
                    bwork->level_abort_map_done[level]--;
                    bwork->level_abort_map_mutex[level]->unlock();
                }

                // if !comp consider abort map
                if (!comp && level > 1) {
                    if (elem->is_identity) { // save to abort map...
                        bwork->write_abort_map(level, w->work_I->comp_fail_pos, w->work_I->comp_fail_acc);
                    } else { // if abort map done, check abort map...
                        bool comp_ = bwork->read_abort_map(level, w->work_I->comp_fail_pos, w->work_I->comp_fail_acc);
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

            // still looks equal to canonical base, so create a node in bfs tree
            bfs_element *next_elem = new bfs_element;
            next_elem->c = w->work_c;
            next_elem->I = w->work_I;
            next_elem->init_c = true;
            next_elem->init_I = true;
            next_elem->is_identity = is_identity;
            next_elem->level = level + 1;
            next_elem->base_sz = elem->base_sz + 1;
            assert(next_elem->base == nullptr);
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
            bwork->bfs_level_finished_elements[level].enqueue_bulk(w->finished_elements, finished_elements_sz);

        return true;
    }

    void bfs_fill_queue(dejavu_workspace *w, bfs_workspace* bwork) {
        if(bwork->current_level == bwork->target_level)
            return;
        // throw away old content of queue
        bwork->reserve_current_level();
        //moodycamel::ConcurrentQueue<std::tuple<bfs_element *, int, int>> throwaway_queue;
        //bwork->bfs_level_todo[bwork->current_level].swap(throwaway_queue);
        bwork->bfs_level_todo[bwork->current_level].clear();

        // swap identity to first position...
        for (int j = 0; j < bwork->level_sizes[bwork->current_level - 1]; ++j) {
            bfs_element *elem = bwork->level_states[bwork->current_level - 1][j];
            if(elem->is_identity) {
                bfs_element *first_elem = bwork->level_states[bwork->current_level - 1][0];
                bwork->level_states[bwork->current_level - 1][j] = first_elem;
                bwork->level_states[bwork->current_level - 1][0] = elem;
                break;
            }
        }

        // fill queue
        int expected = 0;
        for (int j = 0; j < bwork->level_sizes[bwork->current_level - 1]; ++j) {
            bfs_element *elem = bwork->level_states[bwork->current_level - 1][j];
            if (elem->weight > 0) {
                int c = elem->target_color;
                int c_size = elem->c->ptn[c] + 1;
                for (int i = c; i < c + c_size; ++i) {
                    expected += 1;
                    bwork->bfs_level_todo[bwork->current_level].enqueue(
                            std::tuple<bfs_element *, int, int>(elem, elem->c->lab[i], -1));
                }
                if (elem->is_identity) {
                    PRINT("[BFS] Abort map expecting: " << c_size);
                    bwork->level_abort_map_done[bwork->current_level] = c_size;
                }
            }
        }
        bwork->level_expecting_finished[bwork->current_level] = expected;
    }

    void bfs_assure_init(dejavu_workspace *w, bfs_workspace* bwork) {
        if(!w->init_bfs) {
            int chunk_sz = bwork->chunk_size;
            w->todo_dequeue = new std::tuple<bfs_element*, int, int>[chunk_sz];
            w->todo_elements = new std::pair<bfs_element *, int>[chunk_sz * 8];
            w->finished_elements = new std::pair<bfs_element *, int>[chunk_sz + 1];
            w->init_bfs = true;
            w->orbit.initialize(w->G->domain_size);
            w->orbit_considered.initialize(w->G->domain_size);
            w->orbit_vertex_worklist.initialize(w->G->domain_size);
            w->generator_fix_base = new shared_permnode*[*w->shared_generators_size];
            w->generator_fix_base_alloc = *w->shared_generators_size;
        }
    }

    // Computes a leaf of the tree by following the path given in base.
    void reconstruct_leaf(dejavu_workspace *w,
                          sgraph *g, coloring* start_c, int* base,
                          int base_sz, bijection *leaf) {
        coloring* c = &w->c;
        invariant* I = &w->I;
        c->copy_force(start_c);
        I->never_fail = true;
        //PRINT("Reconstruction base_sz " << base_sz);
        for(int pos = 0; pos < base_sz; ++pos) {
            const int v = base[pos];
            proceed_state(w, g, c, I, v, nullptr, -1, -1, nullptr);
        };
        leaf->read_from_coloring(c);
    }

    // Performs uniform probing with additional leaf storage.
    bool uniform_from_bfs_search_with_storage(dejavu_workspace *w,
                                              sgraph *g,
                                              shared_workspace_auto* switches, bfs_element *elem,
                                              int selector_seed, strategy *strat,
                                              bijection *automorphism, bool look_close) {
        if(!w->_collect_base.init)
            w->_collect_base.initialize(g->v_size);
        w->_collect_base.reset();
        w->_collect_early_individualized.clear();
        w->_collect_base_sz = 0;

        coloring* c = &w->c;
        invariant* I = &w->I;
        c->copy_force(elem->c);
        *I = *elem->I;
        I->set_compare_invariant(strat->I);

        bool comp;

        // first* individualization (*except for random blueprint individualization during refinement)
        if(elem->c->cells != g->v_size) {
            const int col = elem->target_color;
            const int rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[col] + 1));
            const int v = c->lab[rpos];
            w->_collect_base.push_back(v);
            w->_collect_base_sz += 1;
            int cell_early = -1;
            int individualize_early = -1;
            if (look_close) {
                I->never_fail = true;
                // use random blueprint individualization during refinement
                individualize_early = elem->base_sz + 1;
            } else {
                I->never_fail = false;
                cell_early = (I->compareI->vec_cells)[elem->base_sz];
            }
            I->reset_deviation();
            if(v != elem->deviation_vertex || look_close) { // added || look_close
                comp = proceed_state(w, g, c, I, v, nullptr, cell_early, individualize_early,
                                     &w->_collect_early_individualized);
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

        for(int i = 0; i < w->_collect_early_individualized.size(); ++i) {
            w->_collect_base.push_back(w->_collect_early_individualized[i]);
            w->_collect_base_sz += 1;
        }

        ++switches->experimental_paths;
        w->S.empty_cache();

        // continue individualizing until discrete
        I->never_fail = true;
        do {
            if(switches->done) return false;
            const int col = w->S.select_color_dynamic(g, c, strat);
            if(col == -1) break;
            const int rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[col] + 1));
            const int v    = c->lab[rpos];
            w->_collect_base.push_back(v);
            w->_collect_base_sz += 1;
            comp = proceed_state(w, g, c, I, v, nullptr, -1, -1, nullptr);
        } while(comp);

        if(comp && (strat->I->acc == I->acc)) { // automorphism computed
            automorphism->read_from_coloring(c);
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
            // consider leaf store...
            std::vector<stored_leaf> pointers;

            switches->leaf_store_mutex.lock();
            auto range = switches->leaf_store.equal_range(I->acc);
            for (auto it = range.first; it != range.second; ++it)
                pointers.push_back(it->second);
            if(pointers.empty()) {
                if (switches->leaf_store_explicit <= config.CONFIG_IR_LEAF_STORE_LIMIT) {
                    switches->leaf_store_explicit++;
                    switches->leaf_store.insert(std::pair<long,
                            stored_leaf>(I->acc, stored_leaf(leaf.extract_map(), g->v_size, true, look_close)));
                } else {
                    int* base = new int[w->_collect_base_sz];
                    memcpy(base, w->_collect_base.get_array(), sizeof(int) * w->_collect_base_sz); // ToDo: better solution to this, use bijection?
                    switches->leaf_store.insert(std::pair<long,
                            stored_leaf>(I->acc, stored_leaf(base, w->_collect_base_sz, false, look_close, elem)));
                }
            }
            switches->leaf_store_mutex.unlock();

            comp = false;

            for(size_t i = 0; i < pointers.size(); ++i) {
                // use reconstruction method on actual leaf to make it isomorphism-invariant
                if(w->_collect_early_individualized.size() > 0) {
                    reconstruct_leaf(w, g, elem->c, w->_collect_base.get_array(), w->_collect_base_sz, &leaf);
                }

                // compute automorphism
                automorphism->swap(&leaf);
                automorphism->inverse();
                bijection fake_leaf;

                if(pointers[i].explicit_leaf) {
                    fake_leaf.read_from_array(pointers[i].map, g->v_size);
                } else {
                    // if other leaf was not stored explicitly, reconstruct!
                    if(switches->done) return false;
                    reconstruct_leaf(w, g, pointers[i].start_elem->c, pointers[i].map,pointers[i].map_sz, &fake_leaf);
                }

                automorphism->compose(&fake_leaf);

                if(w->R.certify_automorphism(g, automorphism)) {
                    automorphism->certified = true;
                    automorphism->non_uniform = false;
                    comp = true;
                } else {
                    // PRINT("Reconstruction failed.")
                    comp = false;
                }
            }
        }

        return comp;
    }
};

typedef dejavu_auto_t dejavu_auto;

// TODO: add automorphism consumer

dejavu_stats dejavu_automorphisms(sgraph *g, int* colmap, sassy::sassy_hook* hook) {
    g->dense = !(g->e_size<g->v_size||g->e_size/g->v_size<g->v_size/(g->e_size/g->v_size));
    Clock::time_point timer1 = Clock::now();
    sassy::preprocessor p;
    int* gen_colmap = nullptr;
    if(colmap == nullptr) {
        gen_colmap = new int[g->v_size]; // TODO: memory leak
        for(int i = 0; i < g->v_size; ++i)
            gen_colmap[i] = 0;
        colmap = gen_colmap;
    }
    config.CONFIG_BULK_ALLOCATOR = false;
    if(config.CONFIG_PREPROCESS) {
        //config.CONFIG_PREP_DEACT_PROBE = true;
        //config.CONFIG_PREP_DEACT_DEG2  = true;
        //config.CONFIG_PREP_DEACT_DEG01  = true;
        sassy::sgraph gg;
        gg.v = g->v;
        gg.d = g->d;
        gg.e = g->e;
        gg.v_size = g->v_size;
        gg.d_size = g->d_size;
        gg.e_size = g->e_size;
        p.reduce(&gg, colmap, hook);
        g->v_size = gg.v_size;
        g->d_size = gg.d_size;
        g->e_size = gg.e_size;
        config.CONFIG_IR_SKIP_FIRST_REFINEMENT = true;
    } else {
        config.CONFIG_IR_SKIP_FIRST_REFINEMENT = false;
    }
    config.CONFIG_BULK_ALLOCATOR = true;
    double prep_red_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer1).count());
    Clock::time_point timer2 = Clock::now();
    dejavu_stats a;
    if(g->v_size > 1) {
        dejavu_auto_t d;
        shared_permnode *gens = nullptr;
        a = d.automorphisms(g, colmap, &gens);

        // call consumer for retrieved automorphisms
        work_list_t<int> automorphism_supp;
        automorphism_supp.initialize(g->v_size);
        if (gens != nullptr) {
            shared_permnode *gen_next = gens;
            do {
                if(config.CONFIG_PREPROCESS) {
                    automorphism_supp.reset();
                    for (int i = 0; i < g->v_size; ++i) {
                        if (gen_next->p[i] != i)
                            automorphism_supp.push_back(i);
                    }
                    //std::cout << automorphism_supp.cur_pos << std::endl;
                    p.pre_hook_buffered(g->v_size, gen_next->p, automorphism_supp.cur_pos, automorphism_supp.get_array(), hook);
                } else {
                    (*hook)(g->v_size, gen_next->p, g->v_size, gen_next->p);
                }
                gen_next = gen_next->next;
            } while (gen_next != gens);
        }
    }

    double solve_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer2).count());
    Clock::time_point timer3 = Clock::now();
    double prep_res_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timer3).count());
    std::cout << "(prep-red)  " << prep_red_time / 1000000.0 << "ms" << std::endl;
    std::cout << "(solve) " << solve_time / 1000000.0 << "ms" << std::endl;
    std::cout << "(prep-res)  " << prep_res_time / 1000000.0 << "ms" << std::endl;

    // calculate automorphism group size
    p.exp += a.grp_sz_exp;
    p.multiply_to_group_size(a.grp_sz_man);
    std::cout << "Group size: " << p.base << "*10^" << p.exp << std::endl;

    if(gen_colmap != nullptr)
        delete[] gen_colmap;

    return a;
}

#endif //DEJAVU_DEJAVU_AUTO_H
