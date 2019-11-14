//
// Created by markus on 14.10.19.
//

#include <atomic>
#include <mutex>
#include <algorithm>
#include <random>

#ifndef DEJAVU_UTILITY_H
#define DEJAVU_UTILITY_H

int intRand(const int & min, const int & max, int seed);
double doubleRand(const double & min, const double & max, int seed);

enum modes {MODE_TOURNAMENT, MODE_NON_UNIFORM_PROBE, MODE_TOURNAMENT_IT, MODE_NON_UNIFORM_PROBE_IT,  MODE_UNIFORM_PROBE, MODE_BFS, MODE_WAIT};

struct strategy_metrics {
    int    restarts              = 0;
    double expected_bfs_size     = 0;
    int    expected_level        = 0;
    int    color_refinement_cost = 0;
};

enum ir_operation {
    OP_I, OP_R, OP_END
};

class shared_switches {
public:
    shared_switches();
    bool done = false;
    bool done_fast = false;
    std::atomic_bool done_shared_group;
    std::atomic_bool done_created_group;

    std::atomic<modes> current_mode;
    std::atomic_int    checked;
    strategy_metrics   win_metrics;
    std::atomic_int    win_id;
    std::atomic_int    _ack_done;

    std::mutex         tournament_mutex;

    int tolerance = 1;

    void iterate_tolerance();
    void reset_leaf_tournament();
    bool check_leaf_tournament(int id, strategy_metrics* m);
    bool ack_done();
};

#endif //DEJAVU_UTILITY_H
