#ifndef DEJAVU_BFS_H
#define DEJAVU_BFS_H


#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include "concurrentqueue.h"
#include "coloring.h"
#include "invariant.h"

struct pair_hash {
    inline std::size_t operator()(const std::pair<int,int> & v) const {
        return v.first*31+v.second;
    }
};

class bfs_element {
public:
    // coloring and invariant for the specified path / base
    bfs_element* parent = nullptr;
    coloring*  c = nullptr;
    invariant* I = nullptr;
    int* base = nullptr;
    int  base_sz = -1;

    // position of element in level_states of workspace
    int level        = -1;
    int id           = -1;
    int in_orbit_id  = -1;
    int target_color = -1;
    bool is_identity = false;

    // probability weight for probing in level
    double weight        = -1;
    double parent_weight = -1;

    // synergy information for fast extension and deviation maps
    int deviation_pos    = -1; // save more than one?
    int deviation_val    = -1;
    int deviation_acc    = -1;
    int deviation_vertex = -1;

    // memory management
    bool init_c = false;
    bool init_I = false;
    bool init_base = false;
    ~bfs_element();
};

class bfs_workspace {
public:
    // level array that keeps a queue with tasks
    // bfs_element* is a pointer to the state where int has to be individualized
    moodycamel::ConcurrentQueue<std::tuple<bfs_element*, int, int>>* bfs_level_todo;

    // commit finished elements (or nullptr if element was deleted)
    // integer determines how many todos were added for this element
    moodycamel::ConcurrentQueue<std::pair<bfs_element*, int>>* bfs_level_finished_elements;

    // elements of all levels
    // delete elements of level once work of level + 1 is fully commited, then it is safe
    bfs_element*** level_states;
    int*           level_sizes;
    int*           level_reserved_sizes;
    int*           level_expecting_finished;
    std::unordered_set<std::pair<int, long>, pair_hash>* level_abort_map;
    int*                           level_abort_map_done;
    std::mutex**                   level_abort_map_mutex;

    double*           level_maxweight;
    double*           level_minweight;
    bool done = false;

    // current level where work has to be done, and where experimential paths should probe
    std::atomic_int current_level;

    // the level which was determined to be reached via BFS
    std::atomic_int target_level;

    std::atomic_int abort_map_prune;

    int domain_size;
    int base_size;
    int chunk_size = 16; // ToDo: dynamically adapt this
    bool reached_initial_target = true;
    int initial_target_level;

    std::pair<bfs_element*, int>* finished_elems;
    int finished_elems_sz = -1;

    bfs_workspace();
    ~bfs_workspace();
    void initialize(bfs_element* root_node, int init_c, int domain_size, int base_size);
    bool work_queues(int tolerance);
    void reset_initial_target();
    void write_abort_map(int level, int pos, long val);
    bool read_abort_map(int level, int pos, long val);
};


#endif //DEJAVU_BFS_H