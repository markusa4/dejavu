#include <atomic>
#include <mutex>
#include <algorithm>
#include <random>
#include <unordered_map>
#include "configuration.h"

#ifndef DEJAVU_UTILITY_H
#define DEJAVU_UTILITY_H

int intRand(const int & min, const int & max, int seed);
double doubleRand(const double & min, const double & max, int seed);

#define INV_MARK_ENDREF    (INT32_MAX - 5)
#define INV_MARK_STARTCELL (INT32_MAX - 6)
#define INV_MARK_ENDCELL   (INT32_MAX - 7)

#define PRINT(str) std::cout << str << std::endl;
//#define PRINT(str) (void)0;

// modes of the solver
enum modes {MODE_TOURNAMENT, MODE_NON_UNIFORM_PROBE, MODE_NON_UNIFORM_FROM_BFS, MODE_NON_UNIFORM_PROBE_IT,
            MODE_UNIFORM_PROBE, MODE_BFS, MODE_WAIT};

// metrics used to compare strategies
struct strategy_metrics {
    int    restarts              = 0;
    double expected_bfs_size     = 0;
    int    expected_level        = 0;
    int    color_refinement_cost = 0;
};

template<class vertex_type>
class shared_workspace_temp {
public:
    shared_workspace_temp() {
        done_shared_group.store(false);
        done_created_group.store(false);
        experimental_look_close.store(false);
        base1_skip.store(0);
        _ack_done.store(0);
        win_id.store(-2);
        checked.store(0);
        exit_counter.store(0);
        experimental_paths.store(0);
        experimental_deviation.store(0);
    };

    bool done = false;
    bool done_fast = false;
    std::atomic_bool done_shared_group;
    std::atomic_bool done_created_group;

    // solver mode
    std::atomic<modes> current_mode;
    std::atomic_int    exit_counter;

    // tournament variables
    std::mutex         tournament_mutex;
    std::atomic_int    checked;
    strategy_metrics   win_metrics;
    std::atomic_int    win_id;
    std::atomic_int    _ack_done;
    bool               all_no_restart = true;

    // used for leaf storage decisions and the leaf storage
    std::atomic_int    experimental_budget;
    std::atomic_int    experimental_paths;
    std::atomic_int    experimental_deviation;
    std::atomic_bool   experimental_look_close;
    std::unordered_multimap<long, vertex_type*> leaf_store;

    std::mutex leaf_store_mutex;

    // variable used to synchronize base 1 tournament skip heuristic
    std::atomic_int    base1_skip;

    int tolerance = 1;

    void iterate_tolerance() {
        tolerance *= 2;
    }

    void reset_tolerance(int size, int domain_size) {
        tolerance = std::max(size / (config.CONFIG_IR_SIZE_FACTOR * domain_size), 1);
    }

    bool check_strategy_tournament(int id, strategy_metrics* m, bool early_check) {
        thread_local bool ichecked = false;

        if(!early_check) {
            if (!ichecked) {
                tournament_mutex.lock();
                //std::cout << "late check" << m->color_refinement_cost << std::endl;
                if(m->restarts > 0)
                    all_no_restart = false;

                if ((m->restarts < win_metrics.restarts) ||
                    (m->restarts == win_metrics.restarts && m->expected_bfs_size < win_metrics.expected_bfs_size) ||
                    (m->restarts == win_metrics.restarts && m->expected_bfs_size == win_metrics.expected_bfs_size &&
                     m->color_refinement_cost < win_metrics.color_refinement_cost) ||
                    win_id == -2) {
                    PRINT("[Strat] Best: " << m->restarts << ", " << m->expected_bfs_size << ", " << m->color_refinement_cost);
                    win_metrics = *m;
                    win_id = id;
                }

                checked++;
                tournament_mutex.unlock();
            }
            ichecked = true;
        } else {
            if (!ichecked) {
                if (win_id != -2 && (m->restarts > win_metrics.restarts || win_metrics.restarts == 0)) {
                    tournament_mutex.lock();
                    if(m->restarts > 0)
                        all_no_restart = false;
                    //std::cout << "early concede" << m->color_refinement_cost << std::endl;
                    checked++;
                    tournament_mutex.unlock();
                    ichecked = true;
                }
            }
        }
        return (checked == config.CONFIG_THREADS_REFINEMENT_WORKERS + 1);
    }

    bool ack_done() {
        thread_local bool ichecked = false;

        if(!ichecked) {
            _ack_done++;
        }

        ichecked = true;
        return (checked == config.CONFIG_THREADS_REFINEMENT_WORKERS + 1);
    }
};

#endif //DEJAVU_UTILITY_H
