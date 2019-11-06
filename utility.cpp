//
// Created by markus on 14.10.19.
//

#include "utility.h"
#include "configuration.h"

shared_switches::shared_switches() {
    done_shared_group.store(false);
    done_created_group.store(false);
    win_id.store(-2);
    checked.store(0);
}

bool shared_switches::check_leaf_tournament(int id, int restarts) {
    thread_local bool ichecked = false;

    if(!ichecked) {
        tournament_mutex.lock();
        if(restarts < win_num || win_id == -2) {
            win_num = restarts;
            win_id  = id;
        }
        checked++;
        tournament_mutex.unlock();
    }

    ichecked = true;
    return (checked == config.CONFIG_THREADS_REFINEMENT_WORKERS + 1);
}

void shared_switches::reset_leaf_tournament() {
    done_shared_group.store(false);
    done_created_group.store(false);
    win_id.store(-2);
    checked.store(0);
}

void shared_switches::iterate_tolerance() {
    tolerance *= 2;
}