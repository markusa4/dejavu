#include <atomic>
#include <iostream>
#include <mutex>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include "configuration.h"
#include <fstream>
#include <set>
#include <cstring>
#include <queue>
#include <memory>

#ifndef DEJAVU_UTILITY_H
#define DEJAVU_UTILITY_H

#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__))
    #define OS_WINDOWS
#endif

int intRand(const int & min, const int & max, int seed);
double doubleRand(const double & min, const double & max, int seed);

#define INV_MARK_ENDREF    (INT32_MAX - 5)
#define INV_MARK_STARTCELL (INT32_MAX - 6)
#define INV_MARK_ENDCELL   (INT32_MAX - 7)

#define MASH0(i) (i * (35235237 - i * 5))
#define MASH1(i) (i * (352355 - i * 3))
#define MASH2(i) ((i + 1) * (423733 - (i + 1)))
#define MASH3(i) ((i + 1) * (423233 - (i + 1)))
#define MASH4(i) ((i + 1) * (23524361 - i * 3))
#define MASH5(i) ((i + 1) * (23524361 - i * 3))

#define PRINT(str) std::cout << str << std::endl;
//#define PRINT(str) (void)0;

// modes_auto of the solver
enum modes_auto {MODE_AUTO_TOURNAMENT, MODE_AUTO_NON_UNIFORM_PROBE, MODE_AUTO_UNIFORM_WITH_LEAF_STORAGE,
                 MODE_AUTO_NON_UNIFORM_PROBE_IT, MODE_AUTO_UNIFORM_PROBE, MODE_AUTO_BFS, MODE_AUTO_WAIT};
enum modes_iso {MODE_ISO_BIDIRECTIONAL, MODE_ISO_BIDIRECTIONAL_DEVIATION, MODE_ISO_BFS};

// metrics used to compare strategies
struct strategy_metrics {
    int    restarts              = 0;
    double expected_bfs_size     = 0;
    int    expected_level        = 0;
    int    color_refinement_cost = 0;
};

template<class vertex_t>
struct bfs_element;

template<class vertex_t>
struct stored_leaf {
    stored_leaf(vertex_t* map, int map_sz, bool explicit_leaf, bool look_close) :
                map(map), map_sz(map_sz), explicit_leaf(explicit_leaf), look_close(look_close) {};
    stored_leaf(vertex_t* map, int map_sz, bool explicit_leaf, bool look_close, bfs_element<vertex_t>* start_elem) :
                map(map), map_sz(map_sz), explicit_leaf(explicit_leaf), look_close(look_close), start_elem(start_elem) {};
    bool      explicit_leaf;
    bool      look_close;
    int       map_sz;
    vertex_t* map;
    bfs_element<vertex_t>* start_elem;
};

template<class vertex_t>
class shared_workspace_auto {
public:
    shared_workspace_auto() {
        done_shared_group.store(false);
        done_created_group.store(false);
        experimental_look_close.store(false);
        delete_canon_strategy.store(false);
        _ack_done.store(0);
        win_id.store(-2);
        checked.store(0);
        exit_counter.store(0);
        experimental_paths.store(0);
        experimental_deviation.store(0);
        leaf_store_explicit.store(0);
    };

    ~shared_workspace_auto() {
        // clean up leaf store
        for (auto it = leaf_store.begin(); it != leaf_store.end(); ++it) {
            delete[] it->second.map;
        }
    }

    bool done = false;
    bool done_fast = false;
    std::atomic_bool done_shared_group;
    std::atomic_bool done_created_group;

    // solver mode
    std::atomic<modes_auto> current_mode;
    std::atomic_int    exit_counter;

    //std::vector<unsigned char*>* buffer_buffer;

    // tournament variables
    std::mutex         tournament_mutex;
    std::atomic_int    checked;
    strategy_metrics   win_metrics;
    std::atomic_int    win_id;
    std::atomic_int    _ack_done;
    std::atomic_bool   delete_canon_strategy;
    bool               all_no_restart = true;

    // used for leaf storage decisions and the leaf storage
    std::atomic_int    experimental_budget;
    std::atomic_int    experimental_paths;
    std::atomic_int    experimental_deviation;
    std::atomic_bool   experimental_look_close;
    std::unordered_multimap<long, stored_leaf<vertex_t>> leaf_store;
    std::atomic_int    leaf_store_explicit;
    std::mutex leaf_store_mutex;

    int tolerance = 1;

    void iterate_tolerance() {
        tolerance *= 2;
    }

    void reset_tolerance(int size, int domain_size) {
        tolerance = std::max(size / (config.CONFIG_IR_SIZE_FACTOR * domain_size), 1);
    }

    bool check_strategy_tournament(int id, strategy_metrics* m, bool early_check, bool* this_checked) {
        if(!early_check) {
            if (!*this_checked) {
                tournament_mutex.lock();
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
            *this_checked = true;
        } else {
            if (!*this_checked) {
                if (win_id != -2 && (m->restarts > win_metrics.restarts || win_metrics.restarts == 0)) {
                    // we can already concede to currently best strategy
                    tournament_mutex.lock();
                    if(m->restarts > 0)
                        all_no_restart = false;
                    checked++;
                    tournament_mutex.unlock();
                    *this_checked = true;
                }
            }
        }
        return (checked == config.CONFIG_THREADS_REFINEMENT_WORKERS + 1);
    }
};

template<class vertex_t>
class shared_workspace_iso {
public:
    ~shared_workspace_iso() {
    }

    shared_workspace_iso() {
        done_created_group.store(false);
        experimental_look_close.store(false);
        _ack_done.store(0);
        win_id.store(-2);
        checked.store(0);
        noniso_counter.store(0);
        found_iso.store(false);
        exit_counter.store(0);
        experimental_paths.store(0);
        experimental_deviation.store(0);
        leaf_store[0] = std::unordered_multimap<long, stored_leaf<vertex_t>>();
        leaf_store[1] = std::unordered_multimap<long, stored_leaf<vertex_t>>();
        leaf_store_explicit.store(0);
        deviation_store[0] = std::unordered_set<long>();
        deviation_store[1] = std::unordered_set<long>();
    };

    bool done = false;
    std::atomic_bool done_created_group;

    // solver mode
    std::atomic<modes_iso> current_mode;
    std::mutex switch_mode_mutex;
    std::atomic_int    exit_counter;
    std::atomic_int    noniso_counter;
    std::atomic_bool   found_iso;

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
    std::unordered_multimap<long, stored_leaf<vertex_t>> leaf_store[2];
    std::mutex leaf_store_mutex[2];
    std::atomic_int    leaf_store_explicit;

    std::unordered_set<long> deviation_store[2];
    std::mutex deviation_store_mutex[2];

    // used for API
    std::set<std::tuple<int*, int, int*, long>>* node_store;

    int tolerance = 1;
};

inline bool file_exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

// set specialized for quick resets
class mark_set {
    int mark = 0;
    int *s;
    int sz;
    bool init = false;
public:
    void initialize(int size) {
        s = new int[size];
        sz = size;
        init = true;
        memset(s, mark, sz * sizeof(int));
        reset();
    }
    void initialize_from_array(int* arr, int size) {
        s  = arr;
        sz = size;
        init = false;
        memset(s, mark, sz * sizeof(int));
        reset();
    }
    bool get(int pos) {
        return s[pos] == mark;
    }
    void set(int pos) {
        s[pos] = mark;
    }
    void unset(int pos) {
        s[pos] = mark - 1;
    }
    void reset() {
        if(mark == -1) {
            memset(s, mark, sz * sizeof(int));
        }
        ++mark;
    }
    ~mark_set() {
        if(init)
            delete[] s;
    }
};

template<class T>
class concurrent_queue {
    std::unique_ptr<std::mutex> lock;
    std::deque<T> content;
public:
    concurrent_queue() {
        lock =  std::make_unique<std::mutex>();
    }

    void enqueue(T item) {
        lock->lock();
        content.emplace_back(item);
        lock->unlock();
    }

    void enqueue_bulk(T* items, int size) {
        lock->lock();
        for(int i = 0; i < size; ++i) {
            content.emplace_back(items[i]);
        }
        lock->unlock();
    }

    T dequeue() {
        lock->lock();
        T item = content.front();
        content.pop_front();
        lock->unlock();
        return item;
    }

    bool try_dequeue(T& item)  {
        int i = 0;
        lock->lock();
        if(content.size() > 0) {
            item = content.front();
            content.pop_front();
            ++i;
        }
        lock->unlock();
        return (i > 0);
    }

    int try_dequeue_bulk(T* arr, int chunk_size)  {
        int i = 0;
        lock->lock();
        while(chunk_size - i > 0 && content.size() > 0) {
            arr[i] = content.front();
            content.pop_front();
            ++i;
        }
        lock->unlock();
        return i;
    }

    void clear() {
        lock->lock();
        content.clear();
        lock->unlock();
    }

    int size() {
        int sz = 0;
        lock->lock();
        sz = content.size();
        lock->unlock();
        return sz;
    }
};

#endif //DEJAVU_UTILITY_H
