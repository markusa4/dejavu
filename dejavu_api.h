#ifndef DEJAVU_DEJAVU_API_H
#define DEJAVU_DEJAVU_API_H

#include "dejavu_iso.h"
#include <boost/python.hpp>

class dejavu_api {
public:
public:
    bool random_paths(sgraph_t<int, int, int> *g, int max_length, int num) {
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

        shared_iso_workspace<int> switches;
        return worker_thread(g, true, &switches, nullptr, nullptr,
                             -1,nullptr, nullptr, max_length, num);
    }

private:
    bool worker_thread(sgraph_t<int, int, int> *g_, bool master,
                       shared_iso_workspace<int> *switches, coloring<int> *start_c, strategy<int>* canon_strategy,
                       int communicator_id, bfs_workspace<int> *bwork1, bfs_workspace<int> *bwork2, int max_length,
                       int num) {
        sgraph_t<int, int, int> *g = g_;
        dejavu_workspace<int, int, int> W;

        config.CONFIG_VUJADE = true; // some things ought to work differently now...

        numnodes = 0;
        colorcost = 0;

        // preprocessing
        if (master) {
            config.CONFIG_IR_FORCE_EXPAND_DEVIATION = true;
            config.CONFIG_IR_DENSE = !(g->e_size < g->v_size ||
                                       g->e_size / g->v_size < g->v_size / (g->e_size / g->v_size));
            start_c = new coloring<int>;
            g->initialize_coloring_raw(start_c);
            if (config.CONFIG_PREPROCESS) {
                //  add preprocessing here
            }
            assert(start_c->check());
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
            switches->current_mode = vujade_modes::VU_MODE_BIDIRECTIONAL_DEVIATION;

            // first color refinement
            canon_strategy = new strategy<int>;
            my_canon_I->create_vector(g->v_size);
            W.start_c1 = start_c;
            W.start_c2 = new coloring<int>;
            strategy_metrics m;
            bool comp;
            W.R.refine_coloring(g, start_c, my_canon_I, -1, &m, -1);
            W.start_I.set_compare_invariant(my_canon_I);
            delete my_canon_I;
            my_canon_I = new invariant;
            PRINT("[api] First refinement: " << cref / 1000000.0 << "ms");

            int init_c = W.S.select_color(g, start_c, selector_seed);
            if (init_c == -1) {
                std::cout << "First coloring discrete, special case." << std::endl;
                return 0;
            }

            // create some objects that are initialized after tournament
            W.BW1 = new bfs_workspace<int>();
            bwork1 = W.BW1;

            W.BW2 = new bfs_workspace<int>();
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
                        std::thread(&dejavu_api::worker_thread,
                                    dejavu_api(), g, false, switches, start_c,
                                    canon_strategy, i, W.BW1, W.BW2, max_length, num));
                cpu_set_t cpuset;
                CPU_ZERO(&cpuset);
                CPU_SET(i + (i >= master_sched), &cpuset);
                int rc = pthread_setaffinity_np(work_threads[i].native_handle(),
                                                sizeof(cpu_set_t), &cpuset);
            }
            PRINT("[api] Refinement workers created (" << config.CONFIG_THREADS_REFINEMENT_WORKERS << " threads)");

            // set some workspace variables
            W.start_c1 = new coloring<int>;
            W.start_c1->copy_force(start_c);
            W.id = -1;
        }

        int base_sz = 0;
        W.skip_c.copy_force(start_c);
        W.work_c = new coloring<int>;
        W.work_I = new invariant;
        W.BW1 = bwork1;
        W.BW2 = bwork2;
        W.skiplevels = 0;

        if (!master) {
            W.start_c1 = new coloring<int>;
            W.start_c2 = new coloring<int>;
            W.start_c1->copy_force(start_c);
            W.id = communicator_id;
        }

        strategy<int> *my_strategy;

        auto rst = SELECTOR_LARGEST;
        if (config.CONFIG_IR_FORCE_SELECTOR)
            rst = (selector_type) config.CONFIG_IR_CELL_SELECTOR;

        my_strategy = new strategy<int>(my_canon_leaf, my_canon_I, rst, -1);

        W.S.empty_cache();
        while(!switches->done) {
            bool res= random_path_bounded(&W, g, my_strategy, &base_points, switches, selector_seed, max_length);
            if(res) {
                ++switches->experimental_paths;
                switches->leaf_store_mutex->lock();
                int* save_c = new int[g->v_size];
                memcpy(save_c, W.c.vertex_to_col, g->v_size * sizeof(int));
                switches->node_store.insert(std::pair<int*, long>(save_c, W.I.acc));
                switches->leaf_store_mutex->unlock();
            }
            if(switches->experimental_paths >= num) {
                switches->done = true;
            }
        }
        // TODO: check cleanup
        if (master && !dejavu_kill_request) {
            while (!work_threads.empty()) {
                work_threads[work_threads.size() - 1].join();
                work_threads.pop_back();
            }
            PRINT("[api] Found " << switches->experimental_paths << " paths of maximum length " << max_length);
            PRINT("[api] Cleanup...");
        }
        return false;
    }

    bool random_path_bounded(dejavu_workspace<int, int, int> *w, sgraph *g, strategy<int>* canon_strategy,
                             bijection<int> *automorphism, shared_iso_workspace<int> *switches,
                             int selector_seed, int max_length) {
        const bool* done = &switches->done;

        // workspace
        selector<int, int, int> *S = &w->S;
        coloring<int> *c = &w->c;
        invariant *I = &w->I;
        coloring<int> *start_c  = w->start_c1;
        invariant *start_I = &w->start_I;

        S->empty_cache();
        start_I->reset_compare_invariant();
        start_I->create_vector(g->v_size * 2);
        automorphism->map    = new int[g->v_size];
        automorphism->map_sz = 0;

        int length = 0;

        *I = *start_I;
        c->copy(start_c);

        while (true) {
            if(*done) return false;
            const int s = S->select_color_dynamic(g, c, canon_strategy);
            if (s == -1 || length == max_length) {
                //canon_leaf->read_from_coloring(c);
                //*canon_I = *I;
                PRINT("[api] Found path of length " << length  << ", invariant " << I->acc);
                return true;
            }

            // choose random vertex of class
            const int rpos = s + (intRand(0, INT32_MAX, selector_seed) % (c->ptn[s] + 1));
            const int v = c->lab[rpos];

            // individualize and refine
            proceed_state(w, g, c, I, v, nullptr, -1);
            assert(c->vertex_to_col[v] > 0);

            // check if max_length reached
            length += 1;

            // base point
            //automorphism->map[automorphism->map_sz] = v;
            //automorphism->map_sz += 1;
        }
    }

    bool proceed_state(dejavu_workspace<int, int, int>* w,
                       sgraph_t<int, int, int> * g, coloring<int>* c,
                       invariant* I, int v, strategy_metrics* m, int cell_early) {
        if(!config.CONFIG_IR_IDLE_SKIP)
            cell_early = -1;

        // protocol of selector choice
        I->selection_write(c->vertex_to_col[v]);
        const int init_color_class = w->R.individualize_vertex(c, v);
        bool comp = true;
        comp = comp && w->R.refine_coloring(g, c, I, init_color_class, m, cell_early);
        return comp;
    }
};

void random_paths(sgraph* g, int max_length, int num) {
    dejavu_api v;
    v.random_paths(g, max_length, num);

    // ToDo how to return stuff?
}

#endif //DEJAVU_DEJAVU_API_H
