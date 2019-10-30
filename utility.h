//
// Created by markus on 14.10.19.
//

#include <atomic>
#include <mutex>

#ifndef DEJAVU_UTILITY_H
#define DEJAVU_UTILITY_H

enum modes {MODE_TOURNAMENT, MODE_NON_UNIFORM_PROBE, MODE_UNIFORM_PROBE, MODE_BFS, MODE_WAIT};

class shared_switches {
public:
    shared_switches();
    bool done = false;
    bool done_fast = false;
    std::atomic_bool done_shared_group;
    std::atomic_bool done_created_group;

    std::atomic<modes> current_mode;
    std::atomic_int    checked;
    std::atomic_int    win_num;
    std::atomic_int    win_id;
    std::mutex         tournament_mutex;

    bool check_leaf_tournament(int id, int restarts);
};

#endif //DEJAVU_UTILITY_H
