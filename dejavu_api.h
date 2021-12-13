#ifndef DEJAVU_DEJAVU_API_H
#define DEJAVU_DEJAVU_API_H

#include "dejavu_iso.h"
#include "dejavu_auto.h"

class dejavu_api {
public:
    bool random_paths(sgraph_t<int, int, int> *g, int* vertex_to_col, int max_length, int num, std::set<std::tuple<int*, int, int*, long>>* paths) {
        if(config.CONFIG_THREADS_REFINEMENT_WORKERS == -1) {
            const int max_threads = std::thread::hardware_concurrency();
            if (g->v_size <= 150) {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = std::min(0, max_threads - 1);
            } else if(g->v_size <= 200) {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = std::min(1, max_threads - 1);
            } else if(g->v_size <= 250) {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = std::min(3, max_threads - 1);
            } else {
                config.CONFIG_THREADS_REFINEMENT_WORKERS = max_threads - 1;
            }
        }

        bulk_domain_reset = true;

        bfs_workspace<int> BW1;
        bfs_workspace<int> BW2;

        shared_workspace_iso<int> switches;
        switches.node_store = paths;
        return worker_thread(g, vertex_to_col, true, &switches, nullptr, nullptr,
                             -1,&BW1, &BW2, max_length, num);
    }

private:
    bool worker_thread(sgraph_t<int, int, int> *g_, int* vertex_to_col, bool master,
                       shared_workspace_iso<int> *switches, coloring<int> *_start_c, strategy<int>* canon_strategy,
                       int communicator_id, bfs_workspace<int> *bwork1, bfs_workspace<int> *bwork2, int max_length,
                       int num) {
        sgraph_t<int, int, int> *g = g_;
        dejavu_workspace<int, int, int> W;

        config.CONFIG_SOLVE_ISO = true; // some things ought to work differently now...

        numnodes = 0;
        colorcost = 0;

        // preprocessing
        if (master) {
            config.CONFIG_IR_FORCE_EXPAND_DEVIATION = true;
            config.CONFIG_IR_DENSE = !(g->e_size < g->v_size ||
                                       g->e_size / g->v_size < g->v_size / (g->e_size / g->v_size));
            //start_c1 = new coloring<int>;
            g->initialize_coloring(&W.start_c1, vertex_to_col);
            //g->initialize_coloring_raw(start_c1);
            if (config.CONFIG_PREPROCESS) {
                //  add preprocessing here
            }
            // assert(start_c->check());
        }

        double cref;

        std::vector<std::thread> work_threads;
        bijection<int> base_points;
        bijection<int> actual_base;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() * ((communicator_id * 5) * 5135235);
        int selector_seed = seed;

        invariant start_I;
        invariant *my_canon_I;
        bijection<int> *my_canon_leaf;
        my_canon_I = new invariant;
        my_canon_I->has_compare = false;
        my_canon_I->compare_vec = nullptr;
        my_canon_I->compareI = nullptr;
        my_canon_leaf = new bijection<int>;

        // first color refinement, initialize some more shared structures, launch threads
        if (master) {
            PRINT("[api] Dense graph: " << (config.CONFIG_IR_DENSE ? "true" : "false"));
            switches->current_mode = modes_iso::MODE_ISO_BIDIRECTIONAL_DEVIATION;

            // first color refinement
            canon_strategy = new strategy<int>;
            my_canon_I->create_vector(g->v_size);
            //W.start_c1 = start_c1;
            //W.start_c2 = new coloring<int>;
            strategy_metrics m;
            bool comp;
            W.R.refine_coloring(g, &W.start_c1, my_canon_I, -1, &m, -1, -1, nullptr);
            //PRINT("[api] First refinement inv: " << my_canon_I->acc << "ms");
            //my_canon_I->purge();
           // delete my_canon_I;
            //my_canon_I = new invariant;
            PRINT("[api] First refinement: " << cref / 1000000.0 << "ms");

            int init_c = W.S.select_color(g, &W.start_c1, selector_seed);
            if (init_c == -1) {
                //std::cout << "First coloring discrete, special case." << std::endl;
                //return 0;
            }

            #ifndef OS_WINDOWS
            const int master_sched = sched_getcpu();
            #endif

            W.S.empty_cache();
            {
                #ifndef OS_WINDOWS
                cpu_set_t cpuset;
                CPU_ZERO(&cpuset);
                CPU_SET(master_sched, &cpuset);
                int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
                #endif
            }
            // launch worker threads
            for (int i = 0; i < config.CONFIG_THREADS_REFINEMENT_WORKERS; i++) {
                work_threads.emplace_back(
                        std::thread(&dejavu_api::worker_thread,
                                    dejavu_api(), g, nullptr, false, switches, &W.start_c1,
                                    canon_strategy, i, bwork1, bwork2, max_length, num));
                #ifndef OS_WINDOWS
                cpu_set_t cpuset;
                CPU_ZERO(&cpuset);
                CPU_SET(i + (i >= master_sched), &cpuset);
                int rc = pthread_setaffinity_np(work_threads[i].native_handle(),
                                                sizeof(cpu_set_t), &cpuset);
                #endif
            }
            PRINT("[api] Refinement workers created (" << config.CONFIG_THREADS_REFINEMENT_WORKERS << " threads)");

            // set some workspace variables
            //W.start_c1 = new coloring<int>;
            //W.start_c1->copy_force(start_c1);
            _start_c = &W.start_c1;
            W.id = -1;
        }

        //PRINT("[api] Initial management...");
        int base_sz = 0;
        W.skip_c.copy_force(_start_c);
        W.skiplevels = 0;

        if (!master) {
            //W.start_c1 = new coloring<int>;
            //W.start_c2 = new coloring<int>;
            W.start_c1.copy_force(_start_c);
            W.id = communicator_id;
        }

        strategy<int> *my_strategy;

        auto rst = SELECTOR_LARGEST;
        //auto rst = SELECTOR_SMALLEST;
        if (config.CONFIG_IR_FORCE_SELECTOR)
            rst = (selector_type) config.CONFIG_IR_CELL_SELECTOR;

        my_strategy = new strategy<int>(my_canon_leaf, my_canon_I, rst, -1);

        //PRINT("[api] Entering main loop...");
        W.S.empty_cache();
        while(!switches->done) {
            //PRINT("[api] Computing path...");
            bool res = random_path_bounded(&W, g, my_strategy, &base_points, switches, selector_seed, max_length);
            if(res) {
                ++switches->experimental_paths;
                switches->leaf_store_mutex->lock();
                int* save_c = new int[g->v_size];
                //int* save_b = new int[base_points.map_sz];
                int* save_b = base_points.extract_map();
                memcpy(save_c, W.c.vertex_to_col, g->v_size * sizeof(int));
                //memcpy(save_b, base_points.map, base_points.map_sz * sizeof(int));

                switches->node_store->insert(std::tuple<int*, int, int*, long>(save_b, base_points.map_sz, save_c, W.I.acc));
                switches->leaf_store_mutex->unlock();
            }
            if(switches->experimental_paths >= num) {
                switches->done = true;
            }
            //delete[] base_points.map;
        }
        // TODO: check cleanup
        if (master && !dejavu_kill_request) {
            //PRINT("[api] Joining threads...");
            while (!work_threads.empty()) {
                work_threads[work_threads.size() - 1].join();
                work_threads.pop_back();
            }
            PRINT("[api] Found " << switches->experimental_paths << " paths of maximum length " << max_length);
            //PRINT("[api] Cleanup...");

            delete my_strategy;
            delete canon_strategy;
            delete my_canon_leaf;
            delete my_canon_I;
        }
        return false;
    }

    bool random_path_bounded(dejavu_workspace<int, int, int> *w, sgraph *g, strategy<int>* canon_strategy,
                             bijection<int> *automorphism, shared_workspace_iso<int> *switches,
                             int selector_seed, int max_length) {
        const bool* done = &switches->done;
        //PRINT("[api] Entering path routine...");

        // workspace
        selector<int, int, int> *S = &w->S;
        coloring<int> *c = &w->c;
        invariant *I = &w->I;
        coloring<int> *start_c  = &w->start_c1;
        invariant *start_I = &w->start_I;

        S->empty_cache();
        start_I->reset_compare_invariant();
        start_I->create_vector(1);

        automorphism->initialize_empty(g->v_size);

        int length = 0;

        *I = *start_I;
        c->copy(start_c);

        while (true) {
            if(*done) {
                //PRINT("[api] Aborting path...");
                start_I->purge();
                return false;
            }
            const int s = S->select_color_dynamic(g, c, canon_strategy);
            if (s == -1 || length == max_length) {
                //canon_leaf->read_from_coloring(c);
                //*canon_I = *I;
                PRINT("[api] Found path of length " << length  << ", invariant " << I->acc);
                start_I->purge();
                return true;
            }

            // choose random vertex of class
            const int rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
            const int v = c->lab[rpos];

            // individualize and refine
            proceed_state(w, g, c, I, v, nullptr, -1);
            assert(c->vertex_to_col[v] > 0);

            // check if max_length reached
            automorphism->append(v);
            length += 1;
        }
    }

    bool proceed_state(dejavu_workspace<int, int, int>* w,
                       sgraph_t<int, int, int> * g, coloring<int>* c,
                       invariant* I, int v, strategy_metrics* m, int cell_early) {
        if(!config.CONFIG_IR_IDLE_SKIP)
            cell_early = -1;
        //PRINT("[api] Individualizing " << v << " and refining...");
        // protocol of selector choice
        I->selection_write(c->vertex_to_col[v], c->ptn[c->vertex_to_col[v]]);
        const int init_color_class = w->R.individualize_vertex(c, v);
        bool comp = true;
        comp = comp && w->R.refine_coloring(g, c, I, init_color_class, m, cell_early, -1, nullptr);
        return comp;
    }
};

extern "C"  {
    extern void initialize();

    extern int graph_create(int size);

    extern void graph_delete(int graph_handle);

    extern void graph_add_edge(int graph_handle, int v1, int v2);

    extern void graph_add_edge_labelled(int graph_handle, int v1, int v2, int l);

    extern void graph_label(int graph_handle, int v, int l);

    extern int path_get_num(int path_handle);

    extern int path_get_size(int path_handle, int path_id);

    extern int path_get_inv(int path_handle, int path_id);

    extern int path_get_point(int path_handle, int path_id, int path_pos);

    extern int path_get_vertex_color(int path_handle, int path_id, int v);

    extern int random_paths(int graph_handle, int max_length, int num, bool fill_paths);

    extern bijection<int> are_isomorphic(int graph_handle1, int graph_handle2);

    extern std::vector<bijection<int>> automorphisms(int graph_handle);
}

#endif //DEJAVU_DEJAVU_API_H
