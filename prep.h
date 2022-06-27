#ifndef DEJAVU_PREP_H
#define DEJAVU_PREP_H

#include "sgraph.h"
#include "refinement.h"
#include "schreier_shared.h"
#include "selector.h"
#include <vector>

struct recovery_map {
    int domain_size;
    mark_set del;
    mark_set del_e;

    std::vector<std::vector<int>> translation_layers;
    std::vector<std::vector<int>> backward_translation_layers;

    std::vector<int> backward_translation;
    std::vector<std::vector<int>> canonical_recovery_string;

    double base = 1;
    int    exp  = 0;
};

class preprocessor {
    work_list_t<std::vector<int>> add_edge_buff;
    work_list_t<int> worklist_deg0;
    work_list_t<int> worklist_deg1;
    mark_set add_edge_buff_act;

    std::vector<int> g_old_v;
    std::vector<int> g_old_e;
    std::vector<int> old_colmap;
    std::vector<int> g_old_d2;

    std::vector<int> translate_layer_fwd;
    std::vector<int> translate_layer_bwd;


private:
    void reset_automorphism(int* automorphism, int nsupp, int* supp) {
        for(int i = 0; i < nsupp; ++i) {
            automorphism[supp[i]] = supp[i];
        }
    };

    bool check_if_neighbour(sgraph* g, int v1, int v2) {
        const int v1_v = g->v[v1];
        const int v1_d = g->d[v1];
        for(int i = 0; i < v1_d; ++i) {
            const int neigh = g->e[v1_v + i];
            if(neigh == v2)
                return true;
        }
        return false;
    }

    void red_deg2_assume_cref(sgraph* g, int* colmap, recovery_map* rec, dejavu_consumer consume) {
        if(g->v_size <= 1)
            return;
        add_edge_buff_act.reset();
        rec->del.reset();

        worklist_deg1.reset();

        add_edge_buff.initialize(rec->domain_size);
        for(int i = 0; i < rec->domain_size; ++i)
            add_edge_buff.push_back(std::vector<int>()); // TODO: do this smarter... i know how many edges will end up here... should not allocate anything...
        add_edge_buff_act.initialize(rec->domain_size);

        std::vector<int> test_v_to_min_endpoint;
        test_v_to_min_endpoint.reserve(g->v_size);
        for(int i = 0; i < g->v_size; ++i)
            test_v_to_min_endpoint.push_back(-1);
        std::vector<int> test_v_to_max_endpoint;
        test_v_to_max_endpoint.reserve(g->v_size);
        for(int i = 0; i < g->v_size; ++i)
            test_v_to_max_endpoint.push_back(-1);

        work_list_t<int> smallest_endpoint_cnt;
        smallest_endpoint_cnt.initialize(g->v_size);
        for(int i = 0; i < g->v_size; ++i) {
            smallest_endpoint_cnt.push_back(0);
        }

        mark_set path_done;
        path_done.initialize(g->v_size);

        work_list_t<int> path;
        path.initialize(g->v_size);


        for(int i = 0; i < g->v_size; ++i) {
            switch(g->d[i]) {
                case 2:
                    worklist_deg1.push_back(i);
                    break;
                default:
                    break;
            }
        }

        // TODO: special, fast algorithm for length=1 paths? maybe some reduction that involves parallel-1 paths?
        // TODO: first pass to detect parallel paths and collpase them to single path?
        // detecting paths 1: check how many paths are connected to a node
        int num_paths1 = 0;
        double total_path_length = 0;
        while(!worklist_deg1.empty()) {
            const int v_child = worklist_deg1.pop_back();
            if (path_done.get(v_child))
                continue;
            assert(g->d[v_child] == 2);
            //if (g->d[v_child] != 2)
            //    continue;

            num_paths1 += 1;
            path.reset();
            path_done.set(v_child);
            path.push_back(v_child);

            int path_length = 1;

            int v_child_next = g->e[g->v[v_child]];
            int endpoint1, endpoint2 = -1;
            bool cycle = false;
            while(true) {
                if(g->d[v_child_next] != 2) {
                    // v_child_next is endpoint1
                    endpoint1 = v_child_next;
                    break;
                }
                assert(!path_done.get(v_child_next));
                path.push_back(v_child_next);
                path_done.set(v_child_next);
                ++path_length;
                assert(g->d[v_child_next] == 2);
                int v_child_next_next = g->e[g->v[v_child_next]];
                if(path_done.get(v_child_next_next)) {
                    // either went in reverse or cycle
                    v_child_next_next = g->e[g->v[v_child_next] + 1];
                    if(path_done.get(v_child_next_next)) {
                        cycle = true;
                        break;
                    }
                }
                assert(v_child_next_next != v_child);
                v_child_next = v_child_next_next;
            }

            if(!cycle) {
                v_child_next = g->e[g->v[v_child] + 1];
                while (true) {
                    if (g->d[v_child_next] != 2) {
                        // v_child_next is endpoint2
                        endpoint2 = v_child_next;
                        break;
                    }
                    assert(!path_done.get(v_child_next));
                    path.push_back(v_child_next);
                    path_done.set(v_child_next);
                    ++path_length;
                    int v_child_next_next = g->e[g->v[v_child_next]];
                    if (path_done.get(v_child_next_next)) {
                        // either went in reverse or cycle
                        v_child_next_next = g->e[g->v[v_child_next] + 1];
                        if (path_done.get(v_child_next_next)) {
                            cycle = true;
                            break;
                        }
                    }
                    assert(v_child_next_next != v_child);
                    v_child_next = v_child_next_next;
                }
                //std::cout << endpoint1 << "<-->" << endpoint2 << std::endl;
            } else {
                std::cout << "cycle" << std::endl;
            }

            if(!cycle) {
                assert(endpoint1 != endpoint2);
                assert(endpoint1 >= 0);
                assert(endpoint2 >= 0);
                smallest_endpoint_cnt.arr[endpoint1] += 1;
                smallest_endpoint_cnt.arr[endpoint2] += 1;
            }

            total_path_length += path_length;
        }

        // detecting paths 1: mark nodes for deletion and add edges to edge buffer
        for(int i = 0; i < g->v_size; ++i) {
            switch(g->d[i]) {
                case 2:
                    worklist_deg1.push_back(i);
                    break;
                default:
                    break;
            }
        }

        path_done.reset();
        int num_paths2 = 0;
        total_path_length = 0;
        int unique_smallest_endpoint_paths = 0;
        while(!worklist_deg1.empty()) {
            const int v_child = worklist_deg1.pop_back();
            if (path_done.get(v_child))
                continue;
            if (g->d[v_child] != 2)
                continue;

            num_paths2 += 1;
            path_done.set(v_child);

            int path_length = 1;
            path.reset();
            path.push_back(v_child); // probably need 2*n and set v_child into middle, then expand to left/right

            int v_child_next = g->e[g->v[v_child]];
            int endpoint1, endpoint2;
            bool cycle = false;
            while(true) {
                if(g->d[v_child_next] != 2) {
                    // v_child_next is endpoint1
                    endpoint1 = v_child_next;
                    break;
                }
                assert(!path_done.get(v_child_next));
                path_done.set(v_child_next);
                path.push_back(v_child_next);
                ++path_length;
                int v_child_next_next = g->e[g->v[v_child_next]];
                if(path_done.get(v_child_next_next)) {
                    // either went in reverse or cycle
                    v_child_next_next = g->e[g->v[v_child_next] + 1];
                    if(path_done.get(v_child_next_next)) {
                        cycle = true;
                        break;
                    }
                }
                assert(v_child_next_next != v_child);
                v_child_next = v_child_next_next;
            }

            if(!cycle) {
                v_child_next = g->e[g->v[v_child] + 1];
                while (true) {
                    if (g->d[v_child_next] != 2) {
                        // v_child_next is endpoint2
                        endpoint2 = v_child_next;
                        break;
                    }
                    assert(!path_done.get(v_child_next));
                    path_done.set(v_child_next);
                    path.push_back(v_child_next);
                    ++path_length;
                    int v_child_next_next = g->e[g->v[v_child_next]];
                    if (path_done.get(v_child_next_next)) {
                        // either went in reverse or cycle
                        v_child_next_next = g->e[g->v[v_child_next] + 1];
                        if (path_done.get(v_child_next_next)) {
                            cycle = true;
                            break;
                        }
                    }
                    assert(v_child_next_next != v_child);
                    v_child_next = v_child_next_next;
                }
            }

            ///////////////////////////////////////////////////
            // TESTCODE
            const int min_endpoint = std::min(endpoint1, endpoint2);
            const int max_endpoint = std::max(endpoint1, endpoint2);
            for(int i = 0; i < path.cur_pos; ++i) {
                test_v_to_min_endpoint[path.arr[i]] = min_endpoint;
                test_v_to_max_endpoint[path.arr[i]] = max_endpoint;
            }
            //////////////////////////////////////////////////

            total_path_length += path_length;

            if(cycle) {
                continue;
            }

            // unique endpoint reduction TODO make it more general
            if(colmap[endpoint1] != colmap[endpoint2]) {
                //if(smallest_endpoint_cnt.arr[endpoint1] == 1 && smallest_endpoint_cnt.arr[endpoint2] == 1) {
                //if((colmap[endpoint1] < colmap[endpoint2] && smallest_endpoint_cnt.arr[endpoint1] == 1)  ||
                //   (colmap[endpoint1] > colmap[endpoint2] && smallest_endpoint_cnt.arr[endpoint2] == 1)) { // TODO this might be dangerous, depends on ordering of colors? but should be fine no?
                if((smallest_endpoint_cnt.arr[endpoint1] == 1)  || smallest_endpoint_cnt.arr[endpoint2] == 1) { // TODO this might be dangerous
                    // TODO actually there is an underlying graph -- "DAC parts" of the graph can be used to attach paths to canonical endpoints?
                    // TODO: make sure endpoint1 is not neighbour of another vertex of color of endpoint2

                    if(check_if_neighbour(g, endpoint1, endpoint2)) {
                        //std::cout << "neighbours!" << std::endl;
                        continue;
                    }

                    const int col_endpoint2 = colmap[endpoint2];
                    bool col_cycle = false;
                    for(int f = 0; f < g->d[endpoint1]; ++f) {
                        const int col_other = colmap[g->e[g->v[endpoint1] + f]];
                        if(col_other == col_endpoint2) {
                            col_cycle = true;
                            break;
                        }
                    }
                    if(col_cycle) {
                        std::cout << "col cycle!" << std::endl;
                        continue;
                    }

                    add_edge_buff.arr[endpoint2].push_back(endpoint1);
                    add_edge_buff_act.set(endpoint2);
                    add_edge_buff.arr[endpoint1].push_back(endpoint2);
                    add_edge_buff_act.set(endpoint1);

                    assert(g->d[endpoint1] != 2);
                    assert(g->d[endpoint2] != 2);
                    unique_smallest_endpoint_paths += 1;
                    assert(path_length == path.cur_pos);
                    int deleted = 0;
                    for(int i = 0; i < path.cur_pos; ++i) {
                        ++deleted;
                        const int del_v = path.arr[i];
                        assert(g->d[del_v] == 2);
                        //assert(!rec->del2.get(del_v));
                        rec->del.set(del_v);
                    }

                    int unique_endpoint; // pick the unique endpoint to attach canonical information to
                    if((smallest_endpoint_cnt.arr[endpoint1] == 1) && (smallest_endpoint_cnt.arr[endpoint2] != 1)) {
                        unique_endpoint = endpoint1;
                    } else if((smallest_endpoint_cnt.arr[endpoint1] != 1) && (smallest_endpoint_cnt.arr[endpoint2] == 1)) {
                        unique_endpoint = endpoint2;
                    } else if(colmap[endpoint1] < colmap[endpoint2]) {
                        unique_endpoint = endpoint1;
                    } else {
                        assert(colmap[endpoint1] > colmap[endpoint2]);
                        unique_endpoint = endpoint2;
                    }

                    // translate_back
                    const int unique_endpoint_orig = translate_back(unique_endpoint);
                    // attach all represented vertices of path to unique_endpoint_orig in canonical fashion
                    path.sort_after_map(colmap);
                    for(int i = 0; i < path.cur_pos; ++i) {
                        if(path_length > 1) {
                            //std::cout << colmap[path.arr[i]] << " ";
                        }
                        const int path_v_orig = translate_back(path.arr[i]);
                        rec->canonical_recovery_string[unique_endpoint_orig].push_back(path_v_orig);
                        rec->canonical_recovery_string[unique_endpoint_orig].insert(rec->canonical_recovery_string[unique_endpoint_orig].end(),
                                                                                    rec->canonical_recovery_string[path_v_orig].begin(), rec->canonical_recovery_string[path_v_orig].end());
                    }
                    //if(path_length > 1)
                    //    std::cout<<std::endl;
                    path.reset();
                } else {
                    //std::cout << endpoint1 << "(" << colmap[endpoint1] << ")" << "<->" << endpoint2 << "(" << colmap[endpoint2] << ")" << std::endl;
                    // TODO: endpoint that is only connected to uniquely colored paths is also unique endpoint
                    // TODO: can reduce and canonize these by sorting according to color e.g.
                    std::vector<int> col_list;
                    for(int f = 0; f < g->d[endpoint1]; ++f) {
                        const int col_other = colmap[g->e[g->v[endpoint1] + f]];
                        if(g->d[g->e[g->v[endpoint1] + f]] == 2) {
                            col_list.push_back(col_other);
                        }
                    }
                    std::sort(col_list.begin(), col_list.end());
                    for(int i = 0; i < col_list.size(); ++i) {
                        //std::cout << col_list[i] << " ";
                    }
                    //std::cout << std::endl;
                }
            }
        }

        for(int i = 0; i < g->v_size; ++i) {
            if(g->d[i] == 2 && !rec->del.get(i)) {
                //std::cout << test_v_to_min_endpoint[i] << ":" << test_v_to_max_endpoint[i] << ":" << colmap[i] << std::endl;
            }
        }

        // TODO: multi-path removal, etc.
        assert(num_paths1 == num_paths2);
        if(num_paths1 != 0) {
            std::cout << "(prep-red) total paths: " << num_paths1 << " (avg length: " << total_path_length / num_paths1
                      << ")" << std::endl;
            std::cout << "(prep-red) unique endpoints: " << unique_smallest_endpoint_paths << std::endl;
        }
    }

    void write_canonical_recovery_string_to_automorphism(const int from, const int to) {
        assert(rec.canonical_recovery_string[from].size() == rec.canonical_recovery_string[to].size());
        for(int i = 0; i < rec.canonical_recovery_string[to].size(); ++i) {
            std::cout << i << std::endl;
            const int str_from = rec.canonical_recovery_string[from][i];
            const int str_to = rec.canonical_recovery_string[to][i];
            automorphism.arr[str_to]   = str_from;
            automorphism.arr[str_from] = str_to;
            automorphism_supp.push_back(str_from);
            automorphism_supp.push_back(str_to);
        }
        return;
    }

    void red_deg10_assume_cref(sgraph* g, const coloring<int>* c, recovery_map* rec, dejavu_consumer consume) {
        //std::cout << "(prep-red) deg1" << std::endl;
        rec->del.reset();

        worklist_deg0.reset();
        worklist_deg1.reset();

        mark_set is_parent;
        is_parent.initialize(g->v_size);

        g_old_d2.clear();
        g_old_d2.reserve(g->v_size);
        for(int i = 0; i < g->v_size; ++i) {
            g_old_d2.push_back(g->d[i]);
        }

        work_list_t<int> parentlist;
        parentlist.initialize(g->v_size);
        work_list_t<int> childlist;
        childlist.initialize(g->e_size);
        //for(int i = 0; i < g->e_size; ++i) {
        //    childlist.arr[i] = g->e[i];
        //}

        work_list_t<int> childcount;
        childcount.initialize(g->v_size);
        work_list_t<int> childcount_prev;
        childcount_prev.initialize(g->v_size);
        for(int i = 0; i < g->v_size; ++i) {
            childcount.push_back(0);
            childcount_prev.push_back(0);
        }

        work_list_t<std::pair<int, int>> stack1;
        work_list_t<int> map;
        stack1.initialize(g->v_size);
        map.initialize(g->v_size);

        for(int i = 0; i < g->v_size; ++i) {
            switch(g_old_d2[i]) {
                case 0:
                    worklist_deg0.push_back(i);
                    break;
                case 1:
                    worklist_deg1.push_back(i); // TODO: should not add if discrete, so that early out for empty worklists is used
                    // TODO: should also early-out with less overhead!
                    break;
                default:
                    break;
            }
        }

        while(!worklist_deg1.empty()) {
            const int v_child = worklist_deg1.pop_back();
            if(rec->del.get(v_child))
                continue;
            if(g_old_d2[v_child] != 1)
                continue;

            const int v_child_col = c->vertex_to_col[v_child];
            const int child_col_sz = c->ptn[v_child_col] + 1;
            //std::cout << "col " << v_child_col << "(" << child_col_sz << ")" << " from vertex " << v_child << std::endl;
            if(child_col_sz == 1) { // TODO this can never be mapped, this is bloat... should handle before?
                rec->del.set(v_child);
                continue;
            } else {
                parentlist.reset();
                is_parent.reset();
                bool is_pairs = false;
                for (int i = v_child_col; i < v_child_col + child_col_sz; ++i) {
                    int child = c->lab[i];
                    //assert(workspace_d[child] == 1);
                    // TODO: need special code for pairs in same color class

                    // search for parent
                    const int e_pos_child = g->v[child];
                    int parent = g->e[e_pos_child];

                    //if(is_pairs && rec->del1.get(parent))
                    //    continue;

                    int search_parent = 0;
                    while (rec->del.get(parent)) {
                        ++search_parent;
                        parent = g->e[e_pos_child + search_parent];
                    }

                    if(c->vertex_to_col[parent] == c->vertex_to_col[child]) { // dual-color pairs? not an issue...
                        is_pairs = true;
                        //std::cout << "pair" << std::endl;
                        //rec->del1.set(child);
                        //rec->del1.set(parent);
                        //std::cout << "yeah that probably doesn't work" << std::endl;
                        continue;
                    }

                    rec->del.set(child);

                    // save canonical info for parent
                    childlist.arr[g->v[parent] + childcount.arr[parent]] = child;
                    ++childcount.arr[parent];

                    if (!is_parent.get(parent)) {
                        is_parent.set(parent);
                        childcount_prev.arr[parent] = childcount.arr[parent] - 1;
                        parentlist.push_back(parent);
                    }

                    // adjust parent degree
                    g_old_d2[parent] -= 1;
                    if (g_old_d2[parent] == 1) {
                        worklist_deg1.push_back(parent);
                    } else if (g_old_d2[parent] == 0) {
                        worklist_deg0.push_back(parent);
                    }

                    assert(g_old_d2[parent] >= 0);
                }
            }

            while(!parentlist.empty()) {
                const int parent = parentlist.pop_back();
                const int childcount_from = childcount_prev.arr[parent];
                const int childcount_to   = childcount.arr[parent];
                // automorphism 1: long cycle (c1 ... cn)
                assert(childcount_to - childcount_from > 0);
                if(childcount_to - childcount_from == 1)
                    continue;
                int child_from = childlist.arr[g->v[parent] + childcount_from];

                // descending tree of child_from while writing map
                stack1.reset();
                map.reset(); // TODO: could be automorphism_supp, then reset automorphism_supp only half-way
                map.push_back(child_from);
                stack1.push_back(std::pair<int, int>(g->v[child_from], g->v[child_from] + childcount.arr[child_from]));
                while(!stack1.empty()) {
                    std::pair<int, int> from_to = stack1.pop_back();
                    int from       = from_to.first;
                    const int to   = from_to.second;
                    for(int f = from; f < to; ++f) {
                        const int next = childlist.arr[f];
                        const int from_next = g->v[next];
                        const int to_next = g->v[next] + childcount.arr[next];
                        map.push_back(next);
                        assert(next != parent);
                        if (from_next != to_next)
                            stack1.push_back(std::pair<int, int>(from_next, to_next));
                    }
                }
                int j = 2;
                for(int i = childcount_from + 1; i < childcount_to; ++i) {
                    multiply_to_group_size(j);
                    ++j;
                    const int child_to = childlist.arr[g->v[parent] + i];
                    assert(c->vertex_to_col[child_from] == c->vertex_to_col[child_to]);
                    assert(child_from != child_to);
                    int pos = 0;

                    automorphism_supp.reset();
                    // descending tree of child_to while writing automorphism
                    stack1.reset();
                    assert(map.arr[pos] != child_to);
                    const int to_1   = translate_back(child_to);
                    const int from_1 = translate_back(map.arr[pos]);
                    assert(automorphism.arr[to_1] == to_1);
                    assert(automorphism.arr[from_1] == from_1);
                    automorphism.arr[from_1] = to_1;
                    automorphism.arr[to_1]   = from_1;
                    automorphism_supp.push_back(from_1);
                    automorphism_supp.push_back(to_1);
                    write_canonical_recovery_string_to_automorphism(to_1, from_1);
                    ++pos;
                    // child_to and child_from could have canonical strings when translated back
                    assert(childcount.arr[child_to] == childcount.arr[child_from]);
                    stack1.push_back(std::pair<int, int>(g->v[child_to], g->v[child_to] + childcount.arr[child_to]));
                    while(!stack1.empty()) {
                        std::pair<int, int> from_to = stack1.pop_back();
                        int from = from_to.first;
                        const int to   = from_to.second;
                        for(int f = from; f < to; ++f) {
                            const int next      = childlist.arr[f];
                            const int from_next = g->v[next];
                            const int to_next   = g->v[next] + childcount.arr[next];
                            ++from;
                            assert(next >= 0);
                            assert(next < g->v_size);
                            assert(map.arr[pos] != next);

                            const int to_2   = translate_back(next);
                            const int from_2 = translate_back(map.arr[pos]);
                            assert(automorphism.arr[to_2] == to_2);
                            assert(automorphism.arr[from_2] == from_2);
                            automorphism.arr[from_2]   = to_2;
                            automorphism.arr[to_2]     = from_2;
                            automorphism_supp.push_back(from_2);
                            automorphism_supp.push_back(to_2);
                            write_canonical_recovery_string_to_automorphism(to_2, from_2);
                            ++pos;
                            if(from_next != to_next) // there was a semicolon here, should have been bug
                                stack1.push_back(std::pair<int, int>(from_next, to_next));
                        }
                    }

                    assert(pos == map.cur_pos);

                    assert(g_old_d2[child_to] == 1);
                    assert(g_old_d2[child_from] == 1);
                    assert(rec->del.get(child_to));
                    assert(rec->del.get(child_from));
                    consume(g->v_size, automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                    reset_automorphism(automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                    automorphism_supp.reset();
                }
            }
            parentlist.reset();
        }

        if(g->v_size > 1) {
            // search for parents that still remain in the graph, and rewrite childlist structure into canonical string
            for (int _i = 0; _i < g->v_size; ++_i) {
                if (!rec->del.get(_i) && childcount.arr[_i] > 0 && c->ptn[c->vertex_to_col[_i]] > 0) { // should it be childcount.arr[_i] > 0 or childcount.arr[_i] >= 0? was childcount.arr[_i] >= 0 but that doesnt make sense?
                    stack1.reset();
                    stack1.push_back(std::pair<int, int>(g->v[_i], g->v[_i] + childcount.arr[_i]));
                    const int orig_i = translate_back(_i);
                    rec->canonical_recovery_string[orig_i].reserve(childcount.arr[_i]); // TODO: segmentation fault! throw_length_error --- did i fix by switching orig_i to _i?
                    while (!stack1.empty()) {
                        std::pair<int, int> from_to = stack1.pop_back();
                        int from = from_to.first;
                        const int to = from_to.second;
                        if (from == to) {
                            continue;
                        } else {
                            const int next = childlist.arr[from];
                            const int from_next = g->v[next];
                            const int to_next = g->v[next] + childcount.arr[next];
                            ++from;
                            const int orig_next = translate_back(next);
                            rec->canonical_recovery_string[orig_i].push_back(orig_next);
                            // write canonical recovery string of orig_next into orig_i, since it is now represented by
                            // orig_i
                            assert(next != _i);
                            assert(orig_next != orig_i);
                            for(int j = 0; j < rec->canonical_recovery_string[orig_next].size(); ++j)
                                rec->canonical_recovery_string[orig_i].push_back(rec->canonical_recovery_string[orig_next][j]);
                            stack1.push_back(std::pair<int, int>(from, to));
                            stack1.push_back(std::pair<int, int>(from_next, to_next));
                        }
                    }
                }
            }
        }


        while(!worklist_deg0.empty()) {
            const int v_child = worklist_deg0.pop_back();
            if (rec->del.get(v_child))
                continue;
            assert(g_old_d2[v_child] == 0);

            is_parent.reset();
            const int v_child_col = c->vertex_to_col[v_child];
            const int child_col_sz = c->ptn[v_child_col] + 1;
            //std::cout << "col " << v_child_col << "(" << child_col_sz << ")" << " from vertex " << v_child << std::endl;
            const int parent_from = c->lab[v_child_col];
            rec->del.set(parent_from);

            if(child_col_sz == 1)
                continue;
            int j = 2;
            for (int i = v_child_col + 1; i < v_child_col + child_col_sz; ++i) {
                multiply_to_group_size(j);
                ++j;
                const int parent_to = c->lab[i];
                assert(g_old_d2[parent_to] == 0);
                rec->del.set(parent_to);

                const int orig_parent_from = translate_back(parent_from);
                const int orig_parent_to   = translate_back(parent_to);

                assert(rec->canonical_recovery_string[orig_parent_to].size() == rec->canonical_recovery_string[orig_parent_from].size());

                automorphism.arr[orig_parent_to]   = orig_parent_from;
                automorphism.arr[orig_parent_from] = orig_parent_to;
                automorphism_supp.push_back(orig_parent_from);
                automorphism_supp.push_back(orig_parent_to);
                for(int i = 0; i < rec->canonical_recovery_string[orig_parent_to].size(); ++i) {
                    const int str_from = rec->canonical_recovery_string[orig_parent_from][i];
                    const int str_to = rec->canonical_recovery_string[orig_parent_to][i];
                    automorphism.arr[str_to]   = str_from;
                    automorphism.arr[str_from] = str_to;
                    automorphism_supp.push_back(str_from);
                    automorphism_supp.push_back(str_to);
                }

                consume(g->v_size, automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                reset_automorphism(automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                automorphism_supp.reset();
            }
        }
    }

    void red_quotient_components(sgraph* g, int* colmap, recovery_map* rec, dejavu_consumer consume) {
        //std::cout << "(prep-red) deg1" << std::endl;
        rec->del.reset();

        worklist_deg0.reset();
        worklist_deg1.reset();

        mark_set found_match;
        found_match.initialize(g->v_size);

        mark_set connected_col;
        mark_set is_not_matched;
        connected_col.initialize(g->v_size);
        is_not_matched.initialize(g->v_size);

        int v_has_matching_color = 0;

        coloring<int> test_col;
        g->initialize_coloring(&test_col, colmap);

        int cnt = 0;
        int edge_flip_pot = 0;
        int edge_cnt = 0;

        std::vector<int> test_vec;
        for(int y = 0; y < g->v_size; ++y)
            test_vec.push_back(0);

        for(int i = 0; i < g->v_size; ) {
            //std::cout << i << ", deg: " << g->d[i] << ", col_sz: " << test_col.ptn[i] + 1 << std::endl;
            connected_col.reset();
            is_not_matched.reset();

            int connected_cols = 0;
            int connected_col_eq_sz = 0;
            int v = test_col.lab[i];
            for (int f = g->v[v]; f < g->v[v] + g->d[v]; ++f) {
                const int v_neigh = g->e[f];
                if (!connected_col.get(test_col.vertex_to_col[v_neigh])) {
                    assert(test_col.vertex_to_col[v_neigh] >= 0);
                    assert(test_col.vertex_to_col[v_neigh] < g->v_size);
                    assert(test_vec[test_col.vertex_to_col[v_neigh]] == 0);
                    connected_col.set(test_col.vertex_to_col[v_neigh]);
                    connected_cols += 1;
                    if (test_col.ptn[test_col.vertex_to_col[v_neigh]] == test_col.ptn[i])
                        connected_col_eq_sz += 1;
                } else {
                    assert(test_vec[test_col.vertex_to_col[v_neigh]] > 0);
                    is_not_matched.set(test_col.vertex_to_col[v_neigh]);
                }
                test_vec[test_col.vertex_to_col[v_neigh]] += 1;
            }

            for(int ii = 0; ii < test_col.ptn[i] + 1; ++ii) {
                const int vx = test_col.lab[i + ii];
                for (int f = g->v[vx]; f < g->v[vx] + g->d[vx]; ++f) {
                    const int v_neigh = g->e[f];

                    assert(test_vec[test_col.vertex_to_col[v_neigh]] >= 0);
                    assert(test_vec[test_col.vertex_to_col[v_neigh]] < g->v_size);
                    if (test_vec[test_col.vertex_to_col[v_neigh]] ==
                        test_col.ptn[test_col.vertex_to_col[v_neigh]] + 1) {
                        edge_cnt += 1;
                        rec->del_e.set(f); // mark edge for deletion (reverse edge is handled later automatically)
                    }

                    if(ii == 0) {
                        if(test_col.ptn[test_col.vertex_to_col[v_neigh]] == test_col.ptn[i] && test_vec[test_col.vertex_to_col[v_neigh]] ==
                            test_col.ptn[test_col.vertex_to_col[v_neigh]]) {
                            edge_flip_pot += (test_vec[test_col.vertex_to_col[v_neigh]] * (test_col.ptn[i] + 1)) - (test_col.ptn[i] + 1);
                        }
                    }
                }
            }

            for (int f = g->v[v]; f < g->v[v] + g->d[v]; ++f) {
                const int v_neigh = g->e[f];
                test_vec[test_col.vertex_to_col[v_neigh]] = 0;
            }

                i += test_col.ptn[i] + 1;
        }

        //std::cout << "flip_pot: " <<  edge_flip_pot << "/" << g->e_size << std::endl;
        //std::cout << "independent_con: " << edge_cnt << "/" << g->e_size << std::endl;
    }

    void red_quotient_matchings_(sgraph* g, int* colmap, recovery_map* rec, dejavu_consumer consume) {
        //std::cout << "(prep-red) deg1" << std::endl;
        rec->del.reset();

        worklist_deg0.reset();
        worklist_deg1.reset();

        mark_set found_match;
        found_match.initialize(g->v_size);

        mark_set connected_col;
        mark_set is_not_matched;
        connected_col.initialize(g->v_size);
        is_not_matched.initialize(g->v_size);

        int v_has_matching_color = 0;

        coloring<int> test_col;
        g->initialize_coloring(&test_col, colmap);

        int cnt = 0;
        int edge_cnt = 0;

        std::vector<int> test_vec;
        for(int y = 0; y < g->v_size; ++y)
            test_vec.push_back(0);

        for(int i = 0; i < g->v_size; ) {
            //std::cout << i << ", deg: " << g->d[i] << ", col_sz: " << test_col.ptn[i] + 1 << std::endl;
            if(!found_match.get(i) || true) {
                connected_col.reset();
                is_not_matched.reset();

                bool invalid = false;

                bool self_connected = false;

                int connected_cols = 0;
                int connected_col_eq_sz = 0;
                int v = test_col.lab[i];
                for (int f = g->v[v]; f < g->v[v] + g->d[v]; ++f) {
                    const int v_neigh = g->e[f];
                    if(test_col.vertex_to_col[v_neigh] == i) {
                        self_connected = true;
                    }
                    if (!connected_col.get(test_col.vertex_to_col[v_neigh])) {
                        assert(test_col.vertex_to_col[v_neigh] >= 0);
                        assert(test_col.vertex_to_col[v_neigh] < g->v_size);
                        assert(test_vec[test_col.vertex_to_col[v_neigh]] == 0);
                        connected_col.set(test_col.vertex_to_col[v_neigh]);
                        connected_cols += 1;
                        if(found_match.get(test_col.vertex_to_col[v_neigh]))
                            invalid = true;
                        if (test_col.ptn[test_col.vertex_to_col[v_neigh]] == test_col.ptn[i])
                            connected_col_eq_sz += 1;
                    } else {
                        assert(test_vec[test_col.vertex_to_col[v_neigh]] > 0);
                        is_not_matched.set(test_col.vertex_to_col[v_neigh]);
                    }
                    test_vec[test_col.vertex_to_col[v_neigh]] += 1;
                }

                //if(connected_cols == 2) { // TODO: can I do something here?
                //    std::cout << i << ", deg: " << g->d[i] << ", col_sz: " << test_col.ptn[i] + 1 << ", self_con: " << self_connected << std::endl;
                //}

                for (int f = g->v[v]; f < g->v[v] + g->d[v]; ++f) {
                    const int v_neigh = g->e[f];

                    assert(test_vec[test_col.vertex_to_col[v_neigh]] >= 0);
                    assert(test_vec[test_col.vertex_to_col[v_neigh]] < g->v_size);
                    if(test_vec[test_col.vertex_to_col[v_neigh]] == test_col.ptn[test_col.vertex_to_col[v_neigh]] + 1) {
                        // TODO: remove these edges!
                        std::cout << (test_col.ptn[test_col.vertex_to_col[v_neigh]] + 1) << " <-> " << (test_col.ptn[i] + 1) << std::endl;
                        edge_cnt += (test_col.ptn[test_col.vertex_to_col[v_neigh]] + 1) * (test_col.ptn[i] + 1);

                        //std::cout << "should remove edges" << std::endl;
                    }
                    //std::cout << "reset " << test_col.vertex_to_col[v_neigh] << std::endl;
                    test_vec[test_col.vertex_to_col[v_neigh]] = 0;

                    /*if (test_col.ptn[test_col.vertex_to_col[v_neigh]] == test_col.ptn[i]) {
                        if (!is_not_matched.get(test_col.vertex_to_col[v_neigh])) {
                            // TODO: test if a color connected to v_neigh is in connected_cols
                            bool indep_col = true;
                            for (int fx = g->v[v_neigh]; fx < g->v[v_neigh] + g->d[v_neigh]; ++fx) {
                                const int vx_neigh = g->e[fx];
                                if (connected_col.get(test_col.vertex_to_col[vx_neigh])) {
                                    indep_col = false;
                                    break;
                                }
                                if(found_match.get(test_col.vertex_to_col[vx_neigh]))
                                    invalid = true;
                            }
                            if (indep_col && !invalid) {
                                cnt += test_col.ptn[i] + 1;
                                found_match.set(i);
                                found_match.set(test_col.vertex_to_col[v_neigh]);
                            }
                            break;
                        }
                    }*/
                }

                // TODO: can I flatten color classes connected by matching? yes if they are only connected to independent color classes
            }
            i += test_col.ptn[i] + 1;
        }

        std::cout << "independent_con: " << edge_cnt << "/" << g->e_size << std::endl;
        std::cout << "has_matching_indep_color: " << cnt << "/" << g->v_size << std::endl;
        std::cout << "has_matching_color: " << v_has_matching_color << "/" << g->v_size << std::endl;
    }
    /*void store_graph_aspects(sgraph* g, recovery_map* rec) {
        rec->old_d.clear();
        rec->old_d.reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            rec->old_d.push_back(g->d[i]);
        }
    }*/

    void copy_coloring_to_colmap(const coloring<int>* c, int* colmap) {
        for(int i = 0; i < c->lab_sz; ++i) {
            colmap[i] = c->vertex_to_col[i];
        }
    }

    void perform_del(sgraph*g, int* colmap, recovery_map* rec) {
        // TODO: "sparse" deletion?

        // copy some stuff
        g_old_v.clear();
        g_old_e.clear();
        old_colmap.clear();
        g_old_d2.clear();
        translate_layer_fwd.clear();
        translate_layer_bwd.clear();

        for(int i = 0; i < rec->backward_translation_layers[rec->backward_translation_layers.size() - 1].size(); ++i)
            translate_layer_bwd.push_back(rec->backward_translation_layers[rec->backward_translation_layers.size() - 1][i]);

        // create translation array from old graph to new graph vertices
        int cnt = 0;
        int new_vsize = 0;
        //rec->translation_layers[fwd_ind].reserve(g->v_size);
        //rec->backward_translation_layers[back_ind].reserve(g->v_size);
        int pos_now = 0;
        for (int i = 0; i < g->v_size; ++i) {
            if (!rec->del.get(i)) {
                //rec->translation_layers[fwd_ind].push_back(cnt);
                translate_layer_fwd.push_back(cnt);
                //rec->backward_translation_layers[back_ind].push_back(rec->translation_layers[fwd_ind].size() - 1); // TODO at most need 2 arrays
                //translate_layer_bwd.push_back(translate_layer_fwd.size() - 1);
                //translate_layer_fwd.size() - 1
                const int translate_back = translate_layer_bwd[translate_layer_fwd.size() - 1];
                rec->backward_translation_layers[rec->backward_translation_layers.size() - 1][cnt] = translate_back;
                ++cnt;
                ++new_vsize;
            } else {
                //rec->translation_layers[fwd_ind].push_back(-1);
                translate_layer_fwd.push_back(-1);
            }
        }

        if(new_vsize == 0 || new_vsize == 1) {
            g->v_size = 0;
            g->e_size = 0;
            g->d_size = 0;
            return;
        }

        rec->backward_translation_layers[rec->backward_translation_layers.size() - 1].resize(cnt);

        g_old_v.reserve(g->v_size);
        g_old_e.reserve(g->e_size);
        g_old_d2.reserve(g->v_size);

        old_colmap.reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) { // TODO shrink graph in-place? at least should re-use some workspace arrays
            old_colmap.push_back(colmap[i]);
        }
        for (int i = 0; i < g->v_size; ++i) {
            g_old_v.push_back(g->v[i]);
        }
        for (int i = 0; i < g->v_size; ++i) {
            g_old_v.push_back(g->v[i]);
            g_old_d2.push_back(g->d[i]);
        }
        for (int i = 0; i < g->e_size; ++i) {
            g_old_e.push_back(g->e[i]);
        }

        // make graph smaller using the translation array
        int epos  = 0;
        for (int i = 0; i < g->v_size; ++i) {
            const int old_v = i;
            //const int new_v = rec->translation_layers[fwd_ind][i];
            const int new_v = translate_layer_fwd[i];

            if (new_v >= 0) {
                colmap[new_v] = old_colmap[old_v];
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + g_old_d2[old_v]; ++j) {
                    const int ve = g_old_e[j];
                    //const int new_ve = rec->translation_layers[fwd_ind][ve];
                    const int new_ve = translate_layer_fwd[ve];
                    if(new_ve >= 0) {
                        assert(new_ve < new_vsize);
                        assert(new_ve >= 0);
                        ++new_d;
                        g->e[epos] = new_ve;
                        ++epos;
                    }
                }
                g->d[new_v] = new_d;
            }
        }

        g->e_size = epos;

        g->v_size = new_vsize;
        g->d_size = new_vsize;
    }

    void perform_del_edge(sgraph*g, int* colmap, recovery_map* rec) {
        int pre_esize = g->e_size;
        // copy some stuff
        g_old_v.clear();
        g_old_e.clear();
        old_colmap.clear();
        g_old_d2.clear();

        g_old_v.reserve(g->v_size);
        g_old_e.reserve(g->e_size);
        g_old_d2.reserve(g->v_size);

        old_colmap.reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            old_colmap.push_back(colmap[i]);
        }
        for (int i = 0; i < g->v_size; ++i) {
            g_old_v.push_back(g->v[i]);
        }
        for (int i = 0; i < g->v_size; ++i) {
            g_old_v.push_back(g->v[i]);
            g_old_d2.push_back(g->d[i]);
        }
        for (int i = 0; i < g->e_size; ++i) {
            g_old_e.push_back(g->e[i]);
        }

        // create translation array from old graph to new graph vertices
        int cnt = 0;
        int new_vsize = g->v_size;

        // make graph smaller using the translation array
        int epos  = 0;
        for (int i = 0; i < g->v_size; ++i) {
            const int old_v = i;
            const int new_v = old_v;

            if (new_v >= 0) {
                colmap[new_v] = old_colmap[old_v];
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + g_old_d2[old_v]; ++j) {
                    const int ve = g_old_e[j];
                    const int new_ve = ve;
                    if(!rec->del_e.get(j)) {
                        assert(new_ve < new_vsize);
                        assert(new_ve >= 0);
                        ++new_d;
                        g->e[epos] = new_ve;
                        ++epos;
                    }
                }
                g->d[new_v] = new_d;
            }
        }

        g->e_size = epos;
        g->v_size = new_vsize;
        g->d_size = new_vsize;

        std::cout << "edges: " << pre_esize << " -> " << g->e_size << std::endl;
    }

    void perform_del2(sgraph*g, int* colmap, recovery_map* rec) {
        if(g->v_size <= 1)
            return;

        // copy some stuff
        g_old_v.clear();
        g_old_e.clear();
        old_colmap.clear();
        g_old_d2.clear();
        translate_layer_fwd.clear();
        translate_layer_bwd.clear();

        for(int i = 0; i < rec->backward_translation_layers[rec->backward_translation_layers.size() - 1].size(); ++i)
            translate_layer_bwd.push_back(rec->backward_translation_layers[rec->backward_translation_layers.size() - 1][i]);

        // create translation array from old graph to new graph vertices
        int cnt = 0;
        int new_vsize = 0;
        for (int i = 0; i < g->v_size; ++i) {
            if (!rec->del.get(i)) {
                translate_layer_fwd.push_back(cnt);
                const int translate_back = translate_layer_bwd[translate_layer_fwd.size() - 1];
                rec->backward_translation_layers[rec->backward_translation_layers.size() - 1][cnt] = translate_back;
                ++cnt;
                ++new_vsize;
            } else {
                //rec->translation_layers[fwd_ind].push_back(-1);
                translate_layer_fwd.push_back(-1);
            }
        }

        if(new_vsize == 0 || new_vsize == 1) {
            std::cout << "here" << std::endl;
            g->v_size = 0;
            g->e_size = 0;
            g->d_size = 0;
            return;
        }

        g_old_v.reserve(g->v_size);
        g_old_e.reserve(g->e_size);
        g_old_d2.reserve(g->v_size);
        if(colmap != nullptr) {
            old_colmap.reserve(g->v_size);
            for (int i = 0; i < g->v_size; ++i) {
                old_colmap.push_back(colmap[i]);
            }
        }
        for (int i = 0; i < g->v_size; ++i) {
            g_old_v.push_back(g->v[i]);
            g_old_d2.push_back(g->d[i]);
        }

        for (int i = 0; i < g->e_size; ++i) {
            g_old_e.push_back(g->e[i]);
        }

        rec->backward_translation_layers[rec->backward_translation_layers.size() - 1].resize(cnt);

        // make graph smaller using the translation array
        int epos  = 0;
        for (int i = 0; i < g->v_size; ++i) {
            const int old_v = i;
            //const int new_v = rec->translation_layers[fwd_ind][old_v];
            const int new_v = translate_layer_fwd[old_v];

            if (new_v >= 0) {
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + g_old_d2[old_v]; ++j) {
                    const int ve = g_old_e[j];
                    //const int new_ve = rec->translation_layers[fwd_ind][ve];
                    const int new_ve = translate_layer_fwd[ve];
                    if(new_ve >= 0) {
                        assert(new_ve < new_vsize);
                        assert(new_ve >= 0);
                        ++new_d;
                        g->e[epos] = new_ve;
                        ++epos;
                    }
                }
                if(add_edge_buff_act.get(old_v)) {
                    while(!add_edge_buff.arr[old_v].empty()) {
                        const int edge_to_old = add_edge_buff.arr[old_v].back();
                        add_edge_buff.arr[old_v].pop_back();
                        //const int edge_to_old = add_edge_buff.arr[old_v];
                        assert(add_edge_buff_act.get(edge_to_old));
                        //const int edge_to_new = rec->translation_layers[fwd_ind][edge_to_old];
                        const int edge_to_new = translate_layer_fwd[edge_to_old];
                        assert(edge_to_old >= 0);
                        assert(edge_to_old <= rec->domain_size);
                        assert(edge_to_new >= 0);
                        assert(edge_to_new <= new_vsize);
                        //std::cout << "adding edge " << old_v << "<->" << edge_to_old << std::endl;
                        ++new_d;
                        g->e[epos] = edge_to_new;
                        ++epos;
                    }
                }
                g->d[new_v] = new_d;
            }
        }

        //PRINT("(prep-red) shrinking graph esize "<< g->e_size << "->" << epos);
        g->e_size = epos;


        //PRINT("(prep-red) shrinking graph "<< g->v_size << "->" << new_vsize);

        // adapt colmap for remaining vertices
        if(colmap != nullptr) {
            for (int i = 0; i < g->v_size; ++i) {
                const int old_v = i;
                //const int new_v = rec->translation_layers[fwd_ind][i];
                const int new_v = translate_layer_fwd[i];
                assert(new_v < new_vsize);
                if (new_v >= 0) {
                    assert(old_colmap[old_v] >= 0);
                    assert(old_colmap[old_v] < rec->domain_size);
                    colmap[new_v] = old_colmap[old_v];
                }
            }
        }

        g->v_size = new_vsize;
        g->d_size = new_vsize;

        for(int i = 0; i < g->v_size; ++i) {
            assert(g->d[i] > 0?g->v[i] < g->e_size:true);
            assert(g->d[i] > 0?g->v[i] >= 0:true);
            assert(g->d[i] >= 0);
            assert(g->d[i] < g->v_size);
        }
        for(int i = 0; i < g->e_size; ++i) {
            assert(g->e[i] < g->v_size);
            assert(g->e[i] >= 0);
        }

        add_edge_buff_act.reset();
    }

    void mark_discrete_for_deletion(sgraph*g, int* colmap, recovery_map* rec) {
        int discrete_cnt = 0;
        worklist_deg0.reset();
        for(int i = 0; i < rec->domain_size; ++i) {
            worklist_deg0.push_back(0);
        }
        for(int i = 0; i < g->v_size; ++i) {
            assert(colmap[i] < rec->domain_size);
            worklist_deg0.arr[colmap[i]]++;
        }
        for(int i = 0; i < g->v_size; ++i) {
            if(worklist_deg0.arr[colmap[i]] == 1)
                rec->del.set(i);
        }
        worklist_deg0.reset();
    }

    void perform_del_discrete(sgraph*g, int* colmap, recovery_map* rec) {
        if(g->v_size <= 1)
            return;

        //work_list_t<int> work;
        //color_count.initialize(rec->domain_size);
        int discrete_cnt = 0;
        worklist_deg0.reset();
        for(int i = 0; i < rec->domain_size; ++i) {
            worklist_deg0.push_back(0);
        }
        for(int i = 0; i < g->v_size; ++i) {
            assert(colmap[i] < rec->domain_size);
            worklist_deg0.arr[colmap[i]]++;
        }
        for(int i = 0; i < g->v_size; ++i) {
            discrete_cnt += (worklist_deg0.arr[colmap[i]] == 1);
        }
        if(discrete_cnt == g->v_size) {
            g->v_size = 0;
            g->d_size = 0;
            g->e_size = 0;
            return;
        }
        if(discrete_cnt == 0) {
            return;
        }

        // add translation layers for this deletion
        /*rec->backward_translation_layers.emplace_back(std::vector<int>());
        const int back_ind = rec->backward_translation_layers.size()-1;
        rec->translation_layers.emplace_back(std::vector<int>());
        const int fwd_ind = rec->translation_layers.size()-1;*/

        // copy some stuff
        g_old_v.clear();
        g_old_e.clear();
        old_colmap.clear();
        g_old_d2.clear();
        translate_layer_fwd.clear();
        translate_layer_bwd.clear();

        for(int i = 0; i < rec->backward_translation_layers[rec->backward_translation_layers.size() - 1].size(); ++i)
            translate_layer_bwd.push_back(rec->backward_translation_layers[rec->backward_translation_layers.size() - 1][i]);

        g_old_v.reserve(g->v_size);
        g_old_e.reserve(g->e_size);
        g_old_d2.reserve(g->v_size);
        assert(colmap != nullptr);
        old_colmap.reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            old_colmap.push_back(colmap[i]);
        }

        for (int i = 0; i < g->v_size; ++i) {
            g_old_v.push_back(g->v[i]);
            g_old_d2.push_back(g->d[i]);
        }

        for (int i = 0; i < g->e_size; ++i) {
            g_old_e.push_back(g->e[i]);
        }

        // create translation array from old graph to new graph vertices
        int cnt = 0;
        int new_vsize = 0;
        /*rec->translation_layers[fwd_ind] = std::vector<int>();  // TODO: do I really need to save multiple of these?
        rec->backward_translation_layers[back_ind] = std::vector<int>();
        rec->translation_layers[fwd_ind].reserve(g->v_size);
        rec->backward_translation_layers[back_ind].reserve(g->v_size);*/
        for (int i = 0; i < g->v_size; ++i) {
            if (worklist_deg0.arr[colmap[i]] != 1) {
                translate_layer_fwd.push_back(cnt);
                //rec->backward_translation_layers[back_ind].push_back(rec->translation_layers[fwd_ind].size() - 1);
                //translate_layer_bwd.push_back(translate_layer_fwd.size() - 1);
                //translate_layer_fwd.size() - 1
                const int translate_back = translate_layer_bwd[translate_layer_fwd.size() - 1];
                rec->backward_translation_layers[rec->backward_translation_layers.size() - 1][cnt] = translate_back;
                ++cnt;
                ++new_vsize;
            } else {
                translate_layer_fwd.push_back(-1);
            }
        }

        rec->backward_translation_layers[rec->backward_translation_layers.size() - 1].resize(cnt);

        // make graph smaller using the translation array
        int epos  = 0;
        for (int i = 0; i < g->v_size; ++i) {
            const int old_v = i;
            assert(i < translate_layer_fwd.size());
            const int new_v = translate_layer_fwd[i];

            if (new_v >= 0) {
                // copy old colmap for remaining vertices
                colmap[new_v] = old_colmap[old_v];
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + g_old_d2[old_v]; ++j) {
                    const int ve = g_old_e[j];
                    assert(ve < translate_layer_fwd.size());
                    const int new_ve = translate_layer_fwd[ve];
                    if(new_ve >= 0) {
                        assert(new_ve < new_vsize);
                        assert(new_ve >= 0);
                        ++new_d;
                        g->e[epos] = new_ve;
                        ++epos;
                    }
                }
                if(new_d == 0) {
                    g->v[new_v] = 0;
                }
                g->d[new_v] = new_d;
            }
        }

        //PRINT("(prep-red) shrinking graph esize "<< g->e_size << "->" << epos);
        g->e_size = epos;


        //PRINT("(prep-red) shrinking graph "<< g->v_size << "->" << new_vsize);

        /*
        // copy old colmap for remaining vertices
        for (int i = 0; i < g->v_size; ++i) {
            const int old_v = i;
            const int new_v = rec->translate3[i];
            assert(new_v < new_vsize);
            if (new_v >= 0) {
                assert(old_colmap[old_v] >= 0);
                assert(old_colmap[old_v] < rec->domain_size);
                colmap[new_v] = old_colmap[old_v];
            }
        }*/

        g->v_size = new_vsize;
        g->d_size = new_vsize;

        for(int i = 0; i < g->v_size; ++i) {
            assert(g->d[i] > 0?g->v[i] < g->e_size:true);
            assert(g->d[i] > 0?g->v[i] >= 0:true);
            assert(g->d[i] >= 0);
            assert(g->d[i] < g->v_size);
        }
        for(int i = 0; i < g->e_size; ++i) {
            assert(g->e[i] < g->v_size);
            assert(g->e[i] >= 0);
        }
    }
public:
    recovery_map rec;
    coloring<int> c;
    work_list_t<int> automorphism;
    work_list_t<int> automorphism_supp;
    bool layers_melded = false;

    void meld_translation_layers() {
        if(layers_melded)
            return;
        const int layers = rec.translation_layers.size();
        rec.backward_translation.reserve(layers - 1);
        const int reduced_size = rec.backward_translation_layers[layers - 1].size();
        for(int i = 0; i < reduced_size; ++i)
            rec.backward_translation.push_back(i);
        for(int i = 0; i < reduced_size; ++i) {
            int next_v = i;
            for(int l = layers - 1; l >= 0; --l) {
                next_v = rec.backward_translation_layers[l][next_v];
            }
            rec.backward_translation[i] = next_v;
        }
        layers_melded = true;
    }

    int translate_back(int v) {
        const int layers = rec.translation_layers.size();
        for(int l = layers - 1; l >= 0; --l) {
            v = rec.backward_translation_layers[l][v];
        }
        return v;
    }

    void multiply_to_group_size(double n) {
        rec.base *= n;
        while(rec.base > 10) {
            rec.base = rec.base / 10;
            rec.exp += 1;
        }
    }

    void sparse_ir(sgraph*g, int* colmap, dejavu_consumer consume) {
        // TODO: ORBITS!
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() * ((1 * 5) * 5135235);
        int selector_seed = seed;

        coloring<int> c1;
        g->initialize_coloring(&c1, colmap); // could re-order to reduce overhead when not applicable
        coloring<int> c2;
        c2.copy(&c1);

        work_list_t<int> _automorphism;
        work_list_t<int> _automorphism_supp;

        _automorphism.initialize(g->v_size);
        _automorphism_supp.initialize(g->v_size);
        for(int i = 0; i < g->v_size; ++i)
            _automorphism.arr[i] = i;

        mark_set touched_color;
        work_list_t<int> touched_color_list;
        touched_color.initialize(g->v_size);
        touched_color_list.initialize(g->v_size);

        invariant I1, I2;
        I1.only_acc = true;
        I2.only_acc = true;

        refinement<int, int, int> R1;
        bool certify = true;
        int start_search_here = 0;
        while(certify) {
            // select a color class of size 2
            bool only_discrete_prev = true;
            int cell = -1;
            for (int i = start_search_here; i < c1.ptn_sz;) {
                if (c1.ptn[i] > 0 && only_discrete_prev) {
                    start_search_here = i;
                    only_discrete_prev = false;
                }
                if (c1.ptn[i] == 1) {
                    cell = i;
                    break;
                }
                i += c1.ptn[i] + 1;
            }
            if (cell == -1)
                return;

            touched_color.reset();
            touched_color_list.reset();

            const int ind_v1 = c1.lab[cell];
            const int ind_v2 = c1.lab[cell + 1];
            const int init_c1 = R1.individualize_vertex(&c1, ind_v1);

            touched_color.set(cell);
            touched_color.set(cell + 1);
            touched_color_list.push_back(cell);
            touched_color_list.push_back(cell + 1);

            R1.refine_coloring(g, &c1, &I1, init_c1, nullptr, -1, -1, nullptr, &touched_color, &touched_color_list);

            const int init_c2 = R1.individualize_vertex(&c2, ind_v2);
            R1.refine_coloring(g, &c2, &I2, init_c2, nullptr, -1, -1, nullptr, &touched_color, &touched_color_list);

            if(I1.acc != I2.acc) {
                for (int i = 0; i < g->v_size; ++i) {
                    colmap[i] = c1.vertex_to_col[i];
                }
                I2.acc = I1.acc;
                c2.copy(&c1);
                continue;
            }


            _automorphism_supp.reset();
            if(c1.cells != g->v_size) { // touched_colors doesn't work properly when early-out is used
                // read automorphism
                for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                    const int i = touched_color_list.arr[j];
                    if (c1.lab[i] != c2.lab[i]) {
                        //if(c1.ptn[i] == 1) {
                       //     std::cout << "(" << c1.lab[i] << ", " << c1.lab[i + 1] << ")" << "(" << c2.lab[i] << ", " << c2.lab[i + 1] << ")" << std::endl;
                       // }
                        //if (colmap[c1.lab[i]] != colmap[c2.lab[i]])
                        //    std::cout << "color mismatch" << std::endl;
                        _automorphism.arr[c1.lab[i]] = c2.lab[i];
                        _automorphism_supp.push_back(c1.lab[i]);
                    }
                }
            } else {
                for (int i = 0; i < g->v_size; ++i) {
                    if (c1.lab[i] != c2.lab[i]) {
                        //if (colmap[c1.lab[i]] != colmap[c2.lab[i]])
                        //    std::cout << "color mismatch" << std::endl;
                        _automorphism.arr[c1.lab[i]] = c2.lab[i];
                        _automorphism_supp.push_back(c1.lab[i]);
                    }
                }
            }

            certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
            assert(certify?R1.certify_automorphism(g, _automorphism.arr):true);
            //std::cout << "(pre-red) sparse-ir certify: " << certify << std::endl;
            if(certify) {
                pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr,
                                     consume);
                multiply_to_group_size(2);

                // reset c2 to c1
                if(c1.cells != g->v_size) {
                    for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                        const int i = touched_color_list.arr[j];
                        if (c1.lab[i] != c2.lab[i]) {
                            c2.lab[i] = c1.lab[i];
                            c2.vertex_to_col[c2.lab[i]] = c1.vertex_to_col[c1.lab[i]];
                            c2.vertex_to_lab[c2.lab[i]] = c1.vertex_to_lab[c1.lab[i]];
                        }
                    }
                } else {
                    for (int i = 0; i < g->v_size; ++i) {
                        if (c1.lab[i] != c2.lab[i]) {
                            c2.lab[i] = c1.lab[i];
                            c2.vertex_to_col[c2.lab[i]] = c1.vertex_to_col[c1.lab[i]];
                            c2.vertex_to_lab[c2.lab[i]] = c1.vertex_to_lab[c1.lab[i]];
                        }
                    }
                }
            } else {
                touched_color.reset();
                touched_color_list.reset();

                std::vector<int> save_colmap;
                save_colmap.reserve(g->v_size);
                for(int i = 0; i < g->v_size; ++i)
                    save_colmap.push_back(c1.vertex_to_col[i]);

                selector<int, int, int> S;
                strategy<int> strat;
                strat.cell_selector_type = SELECTOR_FIRST;
                S.empty_cache();
                int col = -1;
                while(true) {
                    col = S.select_color_dynamic(g, &c1, &strat);
                    if(col == -1) break;
                    const int rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c1.ptn[col] + 1));
                    const int v    = c1.lab[rpos];
                    const int init_color_class = R1.individualize_vertex(&c1, v);
                    R1.refine_coloring(g, &c1, &I1, init_color_class, nullptr, -1, -1,
                                       nullptr, nullptr, nullptr);
                }

                S.empty_cache();
                while(true) {
                    col = S.select_color_dynamic(g, &c2, &strat);
                    if(col == -1) break;
                    const int rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c2.ptn[col] + 1));
                    const int v    = c2.lab[rpos];
                    const int init_color_class = R1.individualize_vertex(&c2, v);
                    R1.refine_coloring(g, &c2, &I2, init_color_class, nullptr, -1, -1,
                                       nullptr, nullptr, nullptr);
                }

                std::cout << "(prep-red) sparse-ir random paths: "<< I1.acc << ", " << I2.acc << std::endl;
                reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                _automorphism_supp.reset();

                for (int i = 0; i < g->v_size; ++i) {
                    if (c1.lab[i] != c2.lab[i]) {
                        _automorphism.arr[c1.lab[i]] = c2.lab[i];
                        _automorphism_supp.push_back(c1.lab[i]);
                    }
                }

                certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                std::cout << "(prep-red) sparse-ir random path certify: " << certify << std::endl;
                if(certify) {
                    pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr,
                                         consume);
                    multiply_to_group_size(2);
                }
                reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                _automorphism_supp.reset();

                for(int i = 0; i < g->v_size; ++i) {
                    colmap[i] = save_colmap[i];
                }

                return;
                // TODO: could compute some type of invariant
            }
            reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
            automorphism_supp.reset();

            if(certify) {
                if(c1.cells == g->v_size) {
                    for (int i = 0; i < g->v_size; ++i) {
                        colmap[i] = c1.vertex_to_col[i];
                    }
                    break;
                } else {
                    for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                        const int c = touched_color_list.arr[j];
                        int f = 0;
                        while (f < c1.ptn[c] + 1) {
                            colmap[c1.lab[c + f]] = c1.vertex_to_col[c1.lab[c + f]];
                            ++f;
                        }
                    }
                }
            }
        }
    }

    void sparse_ir_but_better(sgraph*g, int* colmap, dejavu_consumer consume) {
        // TODO: orbits?
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() * ((1 * 5) * 5135235);
        int selector_seed = seed;

        std::vector<int> individualize_later; // collect vertices to individualize later

        coloring<int> c1;
        g->initialize_coloring(&c1, colmap); // could re-order to reduce overhead when not applicable
        coloring<int> c2;
        c2.copy(&c1);

        work_list_t<int> _automorphism;
        work_list_t<int> _automorphism_supp;

        _automorphism.initialize(g->v_size);
        _automorphism_supp.initialize(g->v_size);
        for(int i = 0; i < g->v_size; ++i)
            _automorphism.arr[i] = i;

        mark_set touched_color;
        work_list_t<int> touched_color_list;
        touched_color.initialize(g->v_size);
        touched_color_list.initialize(g->v_size);

        invariant I1, I2;
        I1.only_acc = true;
        I2.only_acc = true;

        refinement<int, int, int> R1;
        bool certify = true;
        int start_search_here = 0;
        while(certify) {
            // select a non-trivial color class
            bool only_discrete_prev = true;
            int cell = -1;
            for (int i = start_search_here; i < c1.ptn_sz;) {
                if (c1.ptn[i] > 0 && only_discrete_prev) {
                    start_search_here = i;
                    only_discrete_prev = false;
                }
                if (c1.ptn[i] > 0) {
                    cell = i;
                    break;
                }
                i += c1.ptn[i] + 1;
            }
            if (cell == -1)
                return;

            touched_color.reset();
            touched_color_list.reset();

            const int ind_v1 = c1.lab[cell];
            const int ind_v2 = c1.lab[cell + 1];
            const int init_c1 = R1.individualize_vertex(&c1, ind_v1);

            touched_color.set(cell);
            touched_color.set(cell + 1);
            touched_color_list.push_back(cell);
            touched_color_list.push_back(cell + 1);

            R1.refine_coloring(g, &c1, &I1, init_c1, nullptr, -1, -1, nullptr, &touched_color, &touched_color_list);

            const int init_c2 = R1.individualize_vertex(&c2, ind_v2);
            R1.refine_coloring(g, &c2, &I2, init_c2, nullptr, -1, -1, nullptr, &touched_color, &touched_color_list);

            if(I1.acc != I2.acc) {
                for (int i = 0; i < g->v_size; ++i) {
                    colmap[i] = c1.vertex_to_col[i];
                }
                I2.acc = I1.acc;
                c2.copy(&c1);
                continue;
            }


            _automorphism_supp.reset();
            if(c1.cells != g->v_size) { // touched_colors doesn't work properly when early-out is used
                // read automorphism
                for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                    const int i = touched_color_list.arr[j];
                    if (c1.lab[i] != c2.lab[i]) {
                        //if(c1.ptn[i] == 1) {
                        //     std::cout << "(" << c1.lab[i] << ", " << c1.lab[i + 1] << ")" << "(" << c2.lab[i] << ", " << c2.lab[i + 1] << ")" << std::endl;
                        // }
                        //if (colmap[c1.lab[i]] != colmap[c2.lab[i]])
                        //    std::cout << "color mismatch" << std::endl;
                        _automorphism.arr[c1.lab[i]] = c2.lab[i];
                        _automorphism_supp.push_back(c1.lab[i]);
                    }
                }
            } else {
                for (int i = 0; i < g->v_size; ++i) {
                    if (c1.lab[i] != c2.lab[i]) {
                        //if (colmap[c1.lab[i]] != colmap[c2.lab[i]])
                        //    std::cout << "color mismatch" << std::endl;
                        _automorphism.arr[c1.lab[i]] = c2.lab[i];
                        _automorphism_supp.push_back(c1.lab[i]);
                    }
                }
            }

            certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
            assert(certify?R1.certify_automorphism(g, _automorphism.arr):true);
            //std::cout << "(pre-red) sparse-ir certify: " << certify << std::endl;
            if(certify) {
                pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr,
                                     consume);
                multiply_to_group_size(2);

                // reset c2 to c1
                if(c1.cells != g->v_size) {
                    for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                        const int i = touched_color_list.arr[j];
                        if (c1.lab[i] != c2.lab[i]) {
                            c2.lab[i] = c1.lab[i];
                            c2.vertex_to_col[c2.lab[i]] = c1.vertex_to_col[c1.lab[i]];
                            c2.vertex_to_lab[c2.lab[i]] = c1.vertex_to_lab[c1.lab[i]];
                        }
                    }
                } else {
                    for (int i = 0; i < g->v_size; ++i) {
                        if (c1.lab[i] != c2.lab[i]) {
                            c2.lab[i] = c1.lab[i];
                            c2.vertex_to_col[c2.lab[i]] = c1.vertex_to_col[c1.lab[i]];
                            c2.vertex_to_lab[c2.lab[i]] = c1.vertex_to_lab[c1.lab[i]];
                        }
                    }
                }
            } else {
            }

            reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
            automorphism_supp.reset();

            if(certify) {
                if(c1.cells == g->v_size) {
                    for (int i = 0; i < g->v_size; ++i) {
                        colmap[i] = c1.vertex_to_col[i];
                    }
                    break;
                } else {
                    for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                        const int c = touched_color_list.arr[j];
                        int f = 0;
                        while (f < c1.ptn[c] + 1) {
                            colmap[c1.lab[c + f]] = c1.vertex_to_col[c1.lab[c + f]];
                            ++f;
                        }
                    }
                }
            }
        }
    }

    void pre_consumer_inplace(int _n, int* _automorphism, int _supp, int* _automorphism_supp, dejavu_consumer consume) {
        automorphism_supp.reset();

        for(int i = 0; i < _supp; ++i) {
            const int v_from      = _automorphism_supp[i];
            const int orig_v_from = translate_back(v_from);
            const int v_to        = _automorphism[v_from];
            assert(v_from != v_to);
            const int orig_v_to   = translate_back(v_to);
            assert(v_from >= 0);
            assert(v_to >= 0);
            assert(orig_v_from < rec.domain_size);
            assert(orig_v_from >= 0);
            assert(orig_v_to < rec.domain_size);
            assert(orig_v_to >= 0);
            assert(automorphism.arr[orig_v_from] == orig_v_from);
            automorphism.arr[orig_v_from] = orig_v_to;
            automorphism_supp.push_back(orig_v_from);

            assert(rec.canonical_recovery_string[orig_v_to].size() == rec.canonical_recovery_string[orig_v_from].size());

            for(int j = 0; j < rec.canonical_recovery_string[orig_v_to].size(); ++j) {
                const int v_from_t = rec.canonical_recovery_string[orig_v_from][j];
                const int v_to_t   = rec.canonical_recovery_string[orig_v_to][j];
                assert(automorphism.arr[v_from_t] == v_from_t);
                automorphism.arr[v_from_t] = v_to_t;
                automorphism_supp.push_back(v_from_t);
            }
        }

        consume(rec.domain_size, automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
        reset_automorphism(automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
        automorphism_supp.reset();
    }

    void pre_consumer(int _n, int* _automorphism, int _supp, int* _automorphism_supp, dejavu_consumer consume) {
        meld_translation_layers();
        automorphism_supp.reset();

        for(int i = 0; i < _supp; ++i) {
            const int v_from      = _automorphism_supp[i];
            const int orig_v_from = rec.backward_translation[v_from];
            const int v_to        = _automorphism[v_from];
            assert(v_from != v_to);
            const int orig_v_to   = rec.backward_translation[v_to];
            assert(v_from < rec.backward_translation.size());
            assert(v_from >= 0);
            assert(v_to < rec.backward_translation.size());
            assert(v_to >= 0);
            assert(orig_v_from < rec.domain_size);
            assert(orig_v_from >= 0);
            assert(orig_v_to < rec.domain_size);
            assert(orig_v_to >= 0);
            assert(automorphism.arr[orig_v_from] == orig_v_from);
            automorphism.arr[orig_v_from] = orig_v_to;
            automorphism_supp.push_back(orig_v_from);

            assert(rec.canonical_recovery_string[orig_v_to].size() == rec.canonical_recovery_string[orig_v_from].size());

            for(int j = 0; j < rec.canonical_recovery_string[orig_v_to].size(); ++j) {
                const int v_from_t = rec.canonical_recovery_string[orig_v_from][j];
                const int v_to_t   = rec.canonical_recovery_string[orig_v_to][j];
                assert(automorphism.arr[v_from_t] == v_from_t);
                automorphism.arr[v_from_t] = v_to_t;
                automorphism_supp.push_back(v_from_t);
            }
        }

        consume(rec.domain_size, automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
        reset_automorphism(automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
        automorphism_supp.reset();
    }

    void reduce(sgraph* g, int* colmap, dejavu_consumer consume) {
        {
            int deg0 = 0;
            int deg1 = 0;
            int deg2 = 0;
            for (int i = 0; i < g->v_size; ++i) {
                switch (g->d[i]) {
                    case 0:
                        ++deg0;
                        break;
                    case 1:
                        ++deg1;
                        break;
                    case 2:
                        ++deg2;
                        break;
                    default:
                        break;
                }
            }
            std::cout << "(pre-red) before reduction (0, 1, 2) " << deg0 << ", " << deg1 << ", " << deg2 << std::endl;
            std::cout << "(pre-red) before reduction (G, E) " << g->v_size << ", " << g->e_size << std::endl;
        }

        automorphism.initialize(g->v_size);
        for(int i = 0; i < g->v_size; ++i) {
            automorphism.push_back(i);
        }
        automorphism_supp.initialize(g->v_size);
        worklist_deg0.initialize(g->v_size);
        worklist_deg1.initialize(g->v_size);

        translate_layer_fwd.reserve(g->v_size);
        //translate_layer_bwd.reserve(g->v_size);
        rec.backward_translation_layers.emplace_back(std::vector<int>());
        const int back_ind = rec.backward_translation_layers.size()-1;
        rec.translation_layers.emplace_back(std::vector<int>());
        const int fwd_ind = rec.translation_layers.size()-1;
        for(int i = 0; i < g->v_size; ++i)
            rec.backward_translation_layers[back_ind].push_back(i);


        rec.domain_size = g->v_size;

        // assumes colmap is array of length g->v_size
        rec.del = mark_set();
        rec.del.initialize(g->v_size);
        rec.del_e = mark_set();
        rec.del_e.initialize(g->e_size);
        rec.canonical_recovery_string.reserve(g->v_size);
        for(int i = 0; i < rec.domain_size; ++i)
            rec.canonical_recovery_string.emplace_back(std::vector<int>());
        //store_graph_aspects(g, &rec);

        // refinement
        refinement<int, int, int> R;
        g->initialize_coloring(&c, colmap);
        R.refine_coloring_first(g, &c, -1);

        if(c.cells == g->v_size) {
            g->v_size = 0;
            g->e_size = 0;
            g->d_size = 0;
            return;
        }

        // eliminate degree 1 + 0
        red_deg10_assume_cref(g, &c, &rec, consume);
        copy_coloring_to_colmap(&c, colmap);
        mark_discrete_for_deletion(g, colmap, &rec);
        perform_del(g, colmap, &rec);

        // eliminate discrete colors
        //perform_del_discrete(g, colmap, &rec);

        // TODO: maybe degree 2, edge-colors

        /*red_quotient_matchings_(g, colmap, &rec, consume);*/

        // TODO: twins?


        // invariants: paths of length 2 for large regular components

        sparse_ir(g, colmap, consume);
        perform_del_discrete(g, colmap, &rec);

        sparse_ir(g, colmap, consume);
        perform_del_discrete(g, colmap, &rec);

        // eliminate degree 2, no edge-colors
        red_deg2_assume_cref(g, colmap, &rec, consume); // TODO: better more general algorithm
        perform_del2(g, colmap, &rec);
        rec.del.reset();

        red_deg2_assume_cref(g, colmap, &rec, consume);
        perform_del2(g, colmap, &rec);

        red_deg2_assume_cref(g, colmap, &rec, consume);
        perform_del2(g, colmap, &rec);

        //red_quotient_components(g, colmap, &rec, consume); // TODO: bug on highschool1-aigio.dimacs
        //perform_del_edge(g, colmap, &rec);

        if(g->v_size > 1) {
            int deg0 = 0;
            int deg1 = 0;
            int deg2 = 0;
            for (int i = 0; i < g->v_size; ++i) {
                switch (g->d[i]) {
                    case 0:
                        ++deg0;
                        break;
                    case 1:
                        ++deg1;
                        break;
                    case 2:
                        ++deg2;
                        break;
                    default:
                        break;
                }
            }
            std::cout << "(pre-red) after reduction (0, 1, 2) " << deg0 << ", " << deg1 << ", " << deg2 << std::endl;
            std::cout << "(pre-red) after reduction (G, E) " << g->v_size << ", " << g->e_size << std::endl;

            while(deg0 > 0 || deg1 > 0) {
                // refinement
                std::cout << "(prep-red) loop" << std::endl;
                coloring<int> c2;
                refinement<int, int, int> R2;
                g->initialize_coloring(&c2, colmap);
                R2.refine_coloring_first(g, &c2, -1);

                // eliminate degree 1 + 0
                red_deg10_assume_cref(g, &c2, &rec, consume);
                copy_coloring_to_colmap(&c2, colmap);
                mark_discrete_for_deletion(g, colmap, &rec);
                perform_del(g, colmap, &rec);
                //perform_del_discrete(g, colmap, &rec);
                //red_deg2_assume_cref(g, colmap, &rec, consume);
                //perform_del2(g, colmap, &rec);

                sparse_ir(g, colmap, consume);
                perform_del_discrete(g, colmap, &rec);

                {
                    deg0 = 0;
                    deg1 = 0;
                    deg2 = 0;
                    for (int i = 0; i < g->v_size; ++i) {
                        switch (g->d[i]) {
                            case 0:
                                ++deg0;
                                break;
                            case 1:
                                ++deg1;
                                break;
                            case 2:
                                ++deg2;
                                break;
                            default:
                                break;
                        }
                    }
                    std::cout << "(pre-red)[2] after reduction (0, 1, 2) " << deg0 << ", " << deg1 << ", " << deg2
                              << std::endl;
                    std::cout << "(pre-red)[2] after reduction (G, E) " << g->v_size << ", " << g->e_size << std::endl;
                }
                break;
            }
        }

        std::cout << "(prep-red) translation layers: " << rec.backward_translation_layers.size() << std::endl;

         //std::cout << "(pre-red) group size: " << rec.base << "*10^" << rec.exp << std::endl;

        g->sanity_check();


        // TODO "degree one random probing": make some random individualizations, search for deg1 after deleting discrete
        // TODO if color class of deg1 vertices = color class in original coloring, can remove those from graph?

        // TODO: multiple calls for now obviously independent components -- could use "buffer consumer" to translate domains
        // TODO: just use consumer for all the back-translation (just dont "restore"!): compactify translation layers here for this
        // TODO: the consumer just uses standardized "canonized" strings for the backward-reductions, no matter the reduction used

        // TODO: tell dejavu to skip first refinement

        // TODO: AL.bliss is lacking an automorphism!
    }

    void restore(sgraph* g, dejavu_stats* a, dejavu_consumer consume) {
        //sparse_group group;
        //group.domain_size = rec.domain_size;
        //group.initialize_from_automorphism_info(a, &rec);
    }
};

#endif //DEJAVU_PREP_H
