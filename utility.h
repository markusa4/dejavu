#include <atomic>
#include <mutex>
#include <algorithm>
#include <random>
#include <unordered_map>

#ifndef DEJAVU_UTILITY_H
#define DEJAVU_UTILITY_H

int intRand(const int & min, const int & max, int seed);
double doubleRand(const double & min, const double & max, int seed);

// modes of the solver
enum modes {MODE_TOURNAMENT, MODE_NON_UNIFORM_PROBE, MODE_NON_UNIFORM_FROM_BFS, MODE_NON_UNIFORM_PROBE_IT,  MODE_UNIFORM_PROBE, MODE_BFS, MODE_WAIT};

// metrics used to compare strategies
struct strategy_metrics {
    int    restarts              = 0;
    double expected_bfs_size     = 0;
    int    expected_level        = 0;
    int    color_refinement_cost = 0;
};

enum ir_operation {
    OP_I, OP_R, OP_END
};

class shared_workspace {
public:
    shared_workspace();
    bool done = false;
    bool done_fast = false;
    std::atomic_bool done_shared_group;
    std::atomic_bool done_created_group;

    // solver mode
    std::atomic<modes> current_mode;

    // tournament variables
    std::mutex         tournament_mutex;
    std::atomic_int    checked;
    strategy_metrics   win_metrics;
    std::atomic_int    win_id;
    std::atomic_int    _ack_done;

    // used for leaf storage decisions and the leaf storage
    std::atomic_int    experimental_budget;
    std::atomic_int    experimental_paths;
    std::atomic_int    experimental_deviation;
    std::atomic_bool   experimental_look_close;
    std::unordered_multimap<int, int*> leaf_store;
    std::mutex leaf_store_mutex;

    // variable used to synchronize base 1 tournament skip heuristic
    std::atomic_int    base1_skip;

    int tolerance = 1;

    void iterate_tolerance();
    void reset_tolerance(int size, int domain_size);
    bool check_strategy_tournament(int id, strategy_metrics* m, bool early_check);
    bool ack_done();
};

#endif //DEJAVU_UTILITY_H
