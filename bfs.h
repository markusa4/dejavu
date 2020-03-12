#ifndef DEJAVU_BFS_H
#define DEJAVU_BFS_H


#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include "concurrentqueue.h"
#include "coloring.h"
#include "invariant.h"
#include "configuration.h"

struct pair_hash {
    inline std::size_t operator()(const std::pair<int,int> & v) const {
        return v.first*31+v.second;
    }
};

template<class vertex_t>
class bfs_element {
public:
    // coloring and invariant for the specified path / base
    bfs_element<vertex_t>* parent = nullptr;
    coloring<vertex_t>*    c = nullptr;
    invariant* I = nullptr;
    int* base    = nullptr;
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
    std::mutex deviation_write;
    int deviation_pos    = -1; // save more than one?
    int deviation_val    = -1;
    long deviation_acc    = -1;
    int deviation_vertex = -1;

    // memory management
    bool init_c = false;
    bool init_I = false;
    bool init_base = false;
    ~bfs_element() {
        if(init_c)
            delete c;
        if(init_I)
            delete I;
        if(init_base)
            delete base;
    }
};

template<class vertex_t>
class bfs_workspace {
public:
    // level array that keeps a queue with tasks
    // bfs_element* is a pointer to the state where int has to be individualized
    moodycamel::ConcurrentQueue<std::tuple<bfs_element<vertex_t>*, int, int>>* bfs_level_todo;

    // commit finished elements (or nullptr if element was deleted)
    // integer determines how many todos were added for this element
    moodycamel::ConcurrentQueue<std::pair<bfs_element<vertex_t>*, int>>* bfs_level_finished_elements;

    // elements of all levels
    // delete elements of level once work of level + 1 is fully commited, then it is safe
    bfs_element<vertex_t>***         level_states;
    int*                             level_sizes;
    int*                             level_reserved_sizes;
    int*                             level_expecting_finished;
    std::unordered_set<std::pair<int, long>, pair_hash>* level_abort_map;
    int*                             level_abort_map_done;
    std::mutex**                     level_abort_map_mutex;

    double* level_maxweight;
    double* level_minweight;
    bool    done = false;

    // current level where work has to be done, and where experimential paths should probe
    std::atomic_int current_level;

    // the level which was determined to be reached via BFS
    std::atomic_int target_level;

    std::atomic_int abort_map_prune;

    int domain_size;
    int base_size;
    int chunk_size = 32; // ToDo: dynamically adapt this
    bool reached_initial_target = true;
    int initial_target_level;

    std::pair<bfs_element<vertex_t>*, int>* finished_elems;
    int finished_elems_sz = -1;

    bfs_workspace() {}
    ~bfs_workspace() {
        for(int i = 0; i < base_size + 2; ++i)
            delete level_abort_map_mutex[i];

        for(int i = 0; i < current_level; ++i) {
            for(int j = 0; j < level_sizes[i]; ++j)
                delete level_states[i][j];
            delete[] level_states[i];
        }

        delete[] level_states;
        delete[] level_sizes;
        //delete[] level_reserved_sizes;
        delete[] level_maxweight;
        // delete[] level_minweight;
        // delete[] level_abort_map_done;
        delete[] level_abort_map_mutex;
        delete[] level_abort_map;
        // delete[] level_expecting_finished;
    }

    void initialize(bfs_element<vertex_t>* root_elem, int init_c, int domain_size, int base_size) {
        bfs_level_todo              = new moodycamel::ConcurrentQueue<std::tuple<bfs_element<vertex_t>*, int, int>>[base_size + 2];
        bfs_level_finished_elements = new moodycamel::ConcurrentQueue<std::pair<bfs_element<vertex_t>*, int>>[base_size + 2];

        this->domain_size = domain_size;
        this->base_size   = base_size;
        current_level = 1;
        target_level  = -1;
        level_states  = new bfs_element<vertex_t>**[base_size + 2];
        level_sizes          = new int[(base_size + 2) * 4];
        level_reserved_sizes = level_sizes + base_size + 2;
        level_maxweight = new double[(base_size + 2) * 2];
        level_minweight = level_maxweight + base_size + 2;
        level_abort_map_done  = level_sizes + (base_size + 2) * 2;
        level_abort_map_mutex = new std::mutex*[base_size + 2];
        level_abort_map = new std::unordered_set<std::pair<int, long>, pair_hash>[base_size + 2];

        abort_map_prune.store(0);

        level_expecting_finished = level_sizes + (base_size + 2) * 3;
        for(int i = 0; i < base_size + 2; ++i) {
            bfs_level_todo[i] =
                    moodycamel::ConcurrentQueue<std::tuple<bfs_element<vertex_t>*, int, int>>(chunk_size, 0,
                            config.CONFIG_THREADS_REFINEMENT_WORKERS + 1);
            bfs_level_finished_elements[i] =
                    moodycamel::ConcurrentQueue<std::pair<bfs_element<vertex_t>*, int>>(chunk_size, 0,
                            config.CONFIG_THREADS_REFINEMENT_WORKERS + 1);
            level_expecting_finished[i] = 0;
            level_maxweight[i] = 1;
            level_minweight[i] = INT32_MAX;
            level_abort_map[i] = std::unordered_set<std::pair<int, long>, pair_hash>();
            level_abort_map_done[i] = -1;
            level_abort_map_mutex[i] = new std::mutex();
        }

        level_states[0]    = new bfs_element<vertex_t>*[1];
        level_states[0][0] = root_elem;
        root_elem->weight = 1;
        root_elem->target_color = init_c;

        int sz = 0;
        for(int i = init_c; i < init_c + root_elem->c->ptn[init_c] + 1; ++i) {
            int next_v = root_elem->c->lab[i];
            sz += 1;
            bfs_level_todo[current_level].enqueue(
                    std::tuple<bfs_element<vertex_t>*, int, int>(root_elem, next_v, -1));
        }

        // PRINT("[BFS] Abort map expecting: " <<  root_elem->c->ptn[init_c] + 1);
        level_abort_map_done[current_level + 1] = root_elem->c->ptn[init_c] + 1;

        level_expecting_finished[0] = 0;
        level_sizes[0] = 1;
        level_reserved_sizes[0] = 1;

        level_expecting_finished[1] = sz;
        level_states[1] = new bfs_element<vertex_t>*[sz];
        level_reserved_sizes[1] = sz;
        level_sizes[1] = 0;

        finished_elems = new std::pair<bfs_element<vertex_t>*, int>[chunk_size * config.CONFIG_THREADS_REFINEMENT_WORKERS];
        finished_elems_sz = chunk_size * config.CONFIG_THREADS_REFINEMENT_WORKERS;
        PRINT("[BFS] BFS structure initialized, expecting " << sz << " on first level");
    }

    bool work_queues(int tolerance) {
        // no work left!
        if(current_level == target_level) {
            if(!done) {
                done = true;
                PRINT("[BFS] Finished BFS at " << current_level - 1 << " with " << level_sizes[current_level - 1] << " elements, maxweight " << level_maxweight[current_level - 1] << "");
            }
            return false;
        } else {
            done = false;
        }

        // dequeue and process on current level only
        size_t num = bfs_level_finished_elements[current_level].try_dequeue_bulk(finished_elems, finished_elems_sz);

        bool test = false;
        if(num == 0) {
            test = bfs_level_finished_elements[current_level].try_dequeue(finished_elems[0]);
            if(test)
                num = 1;
        }

        int todo, lvl, i;

        for(i = 0; i < num; ++i) {
            bfs_element<vertex_t> *elem = finished_elems[i].first;
            todo = finished_elems[i].second;
            lvl = current_level;

            if (elem == nullptr) {
                level_expecting_finished[lvl] -= todo;
            } else {
                level_states[lvl][level_sizes[lvl]] = elem;
                elem->id = level_sizes[lvl];
                if(elem->weight > level_maxweight[lvl])
                    level_maxweight[lvl] = elem->weight;
                if(elem->weight < level_minweight[lvl] && elem->weight >= 1)
                    level_minweight[lvl] = elem->weight;

                for(int j = 0; j < elem->base_sz; ++j)
                    assert(elem->base[j] >= 0 && elem->base[j] < domain_size);

                elem->level = lvl;
                assert(level_sizes[lvl] < level_reserved_sizes[lvl]);
                level_sizes[lvl] += 1;
                level_expecting_finished[lvl] -= 1;
                level_expecting_finished[lvl + 1] += todo;
            }
        }

        bool need_queue_fill = false;

        // advance level if possible
        if (level_expecting_finished[current_level] == 0) {
            int expected_size = level_expecting_finished[current_level + 1];

            PRINT("[BFS] Advancing to level " << current_level + 1 << " expecting " << level_sizes[current_level]
                                              << " -> " << expected_size << ", maxweight " << level_maxweight[current_level]
                                              << ", minweight " << level_minweight[current_level]);

            if(current_level == target_level - 1 && target_level <= base_size) {
                //if(expected_size < std::max(domain_size / 100, 1)) {
                if(level_sizes[current_level] == 1 && level_expecting_finished[current_level + 1] < chunk_size) {
                    //std::cout << "[B] Increasing target level (expected_size small), setting target level to " << current_level + 1 << std::endl;
                    //target_level += 1;
                }
            }

            if(expected_size < config.CONFIG_IR_SIZE_FACTOR * domain_size * tolerance || config.CONFIG_IR_FULL_BFS) {
                level_reserved_sizes[current_level + 1] = expected_size;
                level_states[current_level + 1] = new bfs_element<vertex_t> * [expected_size];
                level_sizes[current_level + 1] = 0;

                assert(expected_size > 0);

                need_queue_fill = true;
                current_level += 1;
            } else {
                PRINT("[BFS] Refusing to advance level (expected_size too large), setting target level to " << current_level + 1);

                level_reserved_sizes[current_level + 1] = expected_size;
                level_states[current_level + 1] = new bfs_element<vertex_t> * [expected_size];
                level_sizes[current_level + 1] = 0;

                target_level   = current_level + 1;
                reached_initial_target = false;
                current_level += 1;
            }
        }

        return need_queue_fill;
    }

    void reset_initial_target() {
        reached_initial_target = true;
    }
    void write_abort_map(int level, int pos, long val) {
        level_abort_map_mutex[level]->lock();
        level_abort_map[level].insert(std::pair<int, long>(pos, val));
        level_abort_map_done[level]--;
        level_abort_map_mutex[level]->unlock();
    }
    bool read_abort_map(int level, int pos, long val) {
        //std::cout << level_abort_map_done[level] << std::endl;
        if(level_abort_map_done[level] != 0) {
            //if(level_abort_map_done[level] < 0)
            //    std::cout << "bad" << level_abort_map_done[level] << std::endl;
            return true;
        }

        auto check = level_abort_map[level].find(std::pair<int, long>(pos, val));
        return !(check == level_abort_map[level].end());
    }
};

#endif //DEJAVU_BFS_H