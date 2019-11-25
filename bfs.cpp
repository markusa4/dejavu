#include <unordered_set>
#include "bfs.h"
#include "dejavu.h"

bfs_element::~bfs_element() {
    if(init_c)
        delete c;
    if(init_I)
        delete I;
    if(init_base)
        delete base;
}

bfs_workspace::bfs_workspace() {

}

void bfs_workspace::initialize(bfs_element* root_elem, int init_c, int domain_size, int base_size) {
    bfs_level_todo              = new moodycamel::ConcurrentQueue<std::tuple<bfs_element*, int, int>>[base_size + 2];
    bfs_level_finished_elements = new moodycamel::ConcurrentQueue<std::pair<bfs_element*, int>>[base_size + 2];

    this->domain_size = domain_size;
    this->base_size   = base_size;
    current_level = 1;
    target_level  = -1;
    level_states  = new bfs_element**[base_size + 2];
    level_sizes   = new int[base_size + 2];
    level_reserved_sizes = new int[base_size + 2];
    level_maxweight = new double[base_size + 2];
    level_minweight = new double[base_size + 2];
    level_abort_map_done  = new int[base_size + 2];
    level_abort_map_mutex = new std::mutex*[base_size + 2];
    level_abort_map = new std::unordered_set<std::pair<int, int>, pair_hash>[base_size + 2];

    abort_map_prune.store(0);

    level_expecting_finished = new int[base_size + 2];
    for(int i = 0; i < base_size + 2; ++i) {
        bfs_level_todo[i]              = moodycamel::ConcurrentQueue<std::tuple<bfs_element*, int, int>>(chunk_size, 0, config.CONFIG_THREADS_REFINEMENT_WORKERS + 1);
        bfs_level_finished_elements[i] = moodycamel::ConcurrentQueue<std::pair<bfs_element*, int>>(chunk_size, 0, config.CONFIG_THREADS_REFINEMENT_WORKERS + 1);
        level_expecting_finished[i] = 0;
        level_maxweight[i] = 1;
        level_minweight[i] = INT32_MAX;
        level_abort_map[i] = std::unordered_set<std::pair<int, int>, pair_hash>();
        level_abort_map_done[i] = -1;
        level_abort_map_mutex[i] = new std::mutex();
    }

    level_states[0]    = new bfs_element*[1];
    level_states[0][0] = root_elem;
    root_elem->weight = 1;
    root_elem->target_color = init_c;

    int sz = 0;
    for(int i = init_c; i < init_c + root_elem->c->ptn[init_c] + 1; ++i) {
        int next_v = root_elem->c->lab[i];
        sz += 1;
        bfs_level_todo[current_level].enqueue(std::tuple<bfs_element*, int, int>(root_elem, next_v, -1));
    }

    std::cout << "Abort map expecting: " <<  root_elem->c->ptn[init_c] + 1 << std::endl;
    level_abort_map_done[current_level + 1] = root_elem->c->ptn[init_c] + 1;

    level_expecting_finished[0] = 0;
    level_sizes[0] = 1;
    level_reserved_sizes[0] = 1;

    level_expecting_finished[1] = sz;
    level_states[1] = new bfs_element*[sz];
    level_reserved_sizes[1] = sz;
    level_sizes[1] = 0;

    finished_elems = new std::pair<bfs_element*, int>[chunk_size * config.CONFIG_THREADS_REFINEMENT_WORKERS];
    finished_elems_sz = chunk_size * config.CONFIG_THREADS_REFINEMENT_WORKERS;
    std::cout << "[B] BFS structure initialized, expecting " << sz << " on first level" << std::endl;
}

bool bfs_workspace::work_queues(int tolerance) {
    // no work left!
    if(current_level == target_level) {
        if(!done) {
            done = true;
            std::cout << "[B] Finished BFS at " << current_level - 1 << " with " << level_sizes[current_level - 1] << " elements, maxweight " << level_maxweight[current_level - 1] << "" << std::endl;
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
        bfs_element *elem = finished_elems[i].first;
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

        std::cout << "[B] BFS advancing to level " << current_level + 1 << " expecting " << level_sizes[current_level] << " -> " << expected_size << ", maxweight " << level_maxweight[current_level] << ", minweight " << level_minweight[current_level] << std::endl;

        if(current_level == target_level - 1 && target_level <= base_size) {
            //if(expected_size < std::max(domain_size / 100, 1)) {
            if(level_sizes[current_level] == 1 && level_expecting_finished[current_level + 1] < chunk_size) { // ToDo: this should be very efficient! make it efficient! (ToDos and back-and-forth between threads are probably the culprit)
                                                        // ToDo: once levels can be made cheaply high, prefer base points in canon such that no base points have to be fixed
                                                        // ToDo: or save skipperm for BW level
                 //std::cout << "[B] Increasing target level (expected_size small), setting target level to " << current_level + 1 << std::endl;
                 //target_level += 1;
            }
        }

        if(expected_size < config.CONFIG_IR_SIZE_FACTOR * domain_size * tolerance || config.CONFIG_IR_FULLBFS) {
            level_reserved_sizes[current_level + 1] = expected_size;
            level_states[current_level + 1] = new bfs_element * [expected_size];
            level_sizes[current_level + 1] = 0;

            assert(expected_size > 0);

            need_queue_fill = true;
            current_level += 1;
        } else {
            std::cout << "[B] Refusing to advance level (expected_size too large), setting target level to " << current_level + 1 << std::endl;

            level_reserved_sizes[current_level + 1] = expected_size;
            level_states[current_level + 1] = new bfs_element * [expected_size]; // maybe do this only if tolerance is increased?
            level_sizes[current_level + 1] = 0;

            target_level   = current_level + 1;
            reached_initial_target = false;
            current_level += 1;
        }
    }

    return need_queue_fill;
}

void bfs_workspace::reset_initial_target() {
    reached_initial_target = true;
}

void bfs_workspace::write_abort_map(int level, int pos, int val) {
    level_abort_map_mutex[level]->lock();
    level_abort_map[level].insert(std::pair<int, int>(pos, val));
    level_abort_map_done[level]--;
    level_abort_map_mutex[level]->unlock();
}

bool bfs_workspace::read_abort_map(int level, int pos, int val) {
    //std::cout << level_abort_map_done[level] << std::endl;
    if(level_abort_map_done[level] != 0) {
        //if(level_abort_map_done[level] < 0)
        //    std::cout << "bad" << level_abort_map_done[level] << std::endl;
        return true;
    }

    auto check = level_abort_map[level].find(std::pair<int, int>(pos, val));
    return !(check == level_abort_map[level].end());
}

bfs_workspace::~bfs_workspace() {
    for(int i = 0; i < base_size + 2; ++i)
        delete level_abort_map_mutex[i];

    for(int i = 0; i < current_level; ++i) {
        for(int j = 0; j < level_sizes[i]; ++j)
            delete level_states[i][j];
        delete[] level_states[i];
    }

    delete[] level_states;
    delete[] level_sizes;
    delete[] level_reserved_sizes;
    delete[] level_maxweight;
    delete[] level_minweight;
    delete[] level_abort_map_done;
    delete[] level_abort_map_mutex;
    delete[] level_abort_map;
    delete[] level_expecting_finished;
}
