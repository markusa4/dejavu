//
// Created by markus on 14.10.19.
//

#include <atomic>

#ifndef DEJAVU_UTILITY_H
#define DEJAVU_UTILITY_H

struct shared_switches {
    shared_switches();
    bool done = false;
    bool done_fast = false;
    std::atomic_bool done_shared_group;
};

#endif //DEJAVU_UTILITY_H
