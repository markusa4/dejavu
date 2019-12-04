#include <iostream>
#include "utility.h"
#include "configuration.h"

int intRand(const int & min, const int & max, int seed) {
    static thread_local std::mt19937 generator(seed);
    std::uniform_int_distribution<int> distribution(min,max);
    return distribution(generator);
}

double doubleRand(const double & min, const double & max, int seed) {
    static thread_local std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distribution(min,max);
    return floor(distribution(generator));
}


shared_workspace::shared_workspace() {
    done_shared_group.store(false);
    done_created_group.store(false);
    experimental_look_close.store(false);
    base1_skip.store(0);
    _ack_done.store(0);
    win_id.store(-2);
    checked.store(0);
    experimental_paths.store(0);
    experimental_deviation.store(0);
}

bool shared_workspace::check_strategy_tournament(int id, strategy_metrics* m, bool early_check) {
    thread_local bool ichecked = false;

    if(!early_check) {
        if (!ichecked) {
            tournament_mutex.lock();
            //std::cout << "late check" << m->color_refinement_cost << std::endl;
            if(m->restarts > 0)
                all_no_restart = false;

            if ((m->restarts < win_metrics.restarts) ||
                (m->restarts == win_metrics.restarts && m->expected_bfs_size < win_metrics.expected_bfs_size) ||
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
            if (win_id != -2 && m->restarts > win_metrics.restarts) {
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

bool shared_workspace::ack_done() {
    thread_local bool ichecked = false;

    if(!ichecked) {
        _ack_done++;
    }

    ichecked = true;
    return (checked == config.CONFIG_THREADS_REFINEMENT_WORKERS + 1);
}

void shared_workspace::iterate_tolerance() {
    //tolerance = std::max(tolerance + (tolerance / 2), 2);
    tolerance *= 2;
}

void shared_workspace::reset_tolerance(int size, int domain_size) {
    tolerance = std::max(size / (config.CONFIG_IR_SIZE_FACTOR * domain_size), 1);
}
