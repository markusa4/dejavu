#ifndef DEJAVU_PREP_H
#define DEJAVU_PREP_H

#include "sgraph.h"
#include "refinement.h"
#include "schreier_shared.h"
#include <vector>

struct recovery_map {
    int domain_size;
    mark_set del1;
    mark_set del2;

    std::vector<std::vector<int>> translation_layers;
    std::vector<std::vector<int>> backward_translation_layers;

    std::vector<int> backward_translation;

    std::vector<std::vector<int>> canonical_recovery_string;

    double base = 1;
    int    exp  = 0;
};

class preprocessor {
    work_list_t<int> add_edge_buff;
    mark_set add_edge_buff_act;
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
        rec->del2.reset();

        work_list_t<int> worklist_deg2;
        worklist_deg2.initialize(g->v_size);

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
                    worklist_deg2.push_back(i);
                    break;
                default:
                    break;
            }
        }

        // detecting paths 1: check how many paths are connected to a node
        int num_paths1 = 0;
        double total_path_length = 0;
        while(!worklist_deg2.empty()) {
            const int v_child = worklist_deg2.pop_back();
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

            v_child_next = g->e[g->v[v_child] + 1];
            while(true) {
                if(g->d[v_child_next] != 2) {
                    // v_child_next is endpoint2
                    endpoint2 = v_child_next;
                    break;
                }
                assert(!path_done.get(v_child_next));
                path.push_back(v_child_next);
                path_done.set(v_child_next);
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
                assert(endpoint1 != endpoint2);
                assert(endpoint1 >= 0);
                assert(endpoint2 >= 0);
            }

            total_path_length += path_length;

            smallest_endpoint_cnt.arr[endpoint1] += 1;
            smallest_endpoint_cnt.arr[endpoint2] += 1;
        }

        // detecting paths 1: mark nodes for deletion and add edges the edge buffer
        for(int i = 0; i < g->v_size; ++i) {
            switch(g->d[i]) {
                case 2:
                    worklist_deg2.push_back(i);
                    break;
                default:
                    break;
            }
        }

        path_done.reset();
        int num_paths2 = 0;
        total_path_length = 0;
        int unique_smallest_endpoint_paths = 0;
        while(!worklist_deg2.empty()) {
            const int v_child = worklist_deg2.pop_back();
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

            v_child_next = g->e[g->v[v_child] + 1];
            while(true) {
                if(g->d[v_child_next] != 2) {
                    // v_child_next is endpoint2
                    endpoint2 = v_child_next;
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

            total_path_length += path_length;

            if(cycle)
                continue;

            // unique endpoint reduction TODO make it more general
            if(colmap[endpoint1] != colmap[endpoint2]) {
                //if(smallest_endpoint_cnt.arr[endpoint1] == 1 && smallest_endpoint_cnt.arr[endpoint2] == 1) {
                //if((colmap[endpoint1] < colmap[endpoint2] && smallest_endpoint_cnt.arr[endpoint1] == 1)  ||
                //   (colmap[endpoint1] > colmap[endpoint2] && smallest_endpoint_cnt.arr[endpoint2] == 1)) { // TODO this might be dangerous, depends on ordering of colors? but should be fine no?
                if((smallest_endpoint_cnt.arr[endpoint1] == 1)  || smallest_endpoint_cnt.arr[endpoint2] == 1) { // TODO this might be dangerous
                    // TODO actually there is an underlying graph -- "DAC parts" of the graph can be used to attach paths to canonical endpoints?
                    // TODO: make sure endpoint1 is not neighbour of another vertex of color of endpoint2

                    if(check_if_neighbour(g, endpoint1, endpoint2)) {
                        std::cout << "neighbours!" << std::endl;
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
                    if(col_cycle)
                        continue;


                    assert(g->d[endpoint1] != 2);
                    assert(g->d[endpoint2] != 2);
                    unique_smallest_endpoint_paths += 1;
                    assert(path_length == path.cur_pos);
                    int deleted = 0;
                    for(int i = 0; i < path.cur_pos; ++i) {
                        ++deleted;
                        const int del_v = path.arr[i];
                        assert(g->d[del_v] == 2);
                        assert(!rec->del2.get(del_v));
                        rec->del2.set(del_v);
                    }

                    int unique_endpoint; // pick the unique endpoint to attach canonical information to
                    if((smallest_endpoint_cnt.arr[endpoint1] == 1) && (smallest_endpoint_cnt.arr[endpoint2] != 1)) {
                        unique_endpoint = endpoint1;
                        add_edge_buff.arr[endpoint1] = endpoint2;
                        add_edge_buff_act.set(endpoint1);
                    } else if((smallest_endpoint_cnt.arr[endpoint1] != 1) && (smallest_endpoint_cnt.arr[endpoint2] == 1)) {
                        unique_endpoint = endpoint2;
                        add_edge_buff.arr[endpoint2] = endpoint1;
                        add_edge_buff_act.set(endpoint2);
                    } else if(colmap[endpoint1] < colmap[endpoint2]) {
                        unique_endpoint = endpoint1;
                        add_edge_buff.arr[endpoint1] = endpoint2;
                        add_edge_buff_act.set(endpoint1);
                    } else {
                        assert(colmap[endpoint1] > colmap[endpoint2]);
                        unique_endpoint = endpoint2;
                        add_edge_buff.arr[endpoint2] = endpoint1;
                        add_edge_buff_act.set(endpoint2);
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
                }
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

    void red_deg10_assume_cref(sgraph* g, const coloring<int>* c, recovery_map* rec, dejavu_consumer consume) {
        std::cout << "(prep-red) deg1" << std::endl;
        rec->del1.reset();

        work_list_t<int> worklist_deg0;
        work_list_t<int> worklist_deg1;

        mark_set is_parent;
        is_parent.initialize(g->v_size);

        std::vector<int> workspace_d;
        workspace_d.reserve(g->v_size);
        for(int i = 0; i < g->v_size; ++i) {
            workspace_d.push_back(g->d[i]);
        }

        work_list_t<int> parentlist;
        parentlist.initialize(g->v_size);
        work_list_t<int> childlist;
        childlist.initialize(g->e_size);
        for(int i = 0; i < g->e_size; ++i) {
            childlist.arr[i] = g->e[i];
        }

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

        worklist_deg0.initialize(g->v_size);
        worklist_deg1.initialize(g->v_size);

        for(int i = 0; i < g->v_size; ++i) {
            switch(workspace_d[i]) {
                case 0:
                    worklist_deg0.push_back(i);
                    break;
                case 1:
                    worklist_deg1.push_back(i);
                    break;
                default:
                    break;
            }
        }

        while(!worklist_deg1.empty()) {
            const int v_child = worklist_deg1.pop_back();
            if(rec->del1.get(v_child))
                continue;
            if(workspace_d[v_child] != 1)
                continue;

            const int v_child_col = c->vertex_to_col[v_child];
            const int child_col_sz = c->ptn[v_child_col] + 1;
            //std::cout << "col " << v_child_col << "(" << child_col_sz << ")" << " from vertex " << v_child << std::endl;
            if(child_col_sz == 1) { // TODO this can never be mapped, this is bloat... should handle before?
                /*rec->del1.set(v_child);
                const int e_pos_child = g->v[v_child];
                int parent = g->e[e_pos_child];
                int search_parent = 0;
                while (rec->del1.get(parent)) {
                    ++search_parent;
                    parent = g->e[e_pos_child + search_parent];
                }

                // save canonical info for parent
                childlist.arr[g->v[parent] + childcount.arr[parent]] = v_child;
                ++childcount.arr[parent];

                // adjust parent degree
                workspace_d[parent] -= 1;
                if (workspace_d[parent] == 1) {
                    worklist_deg1.push_back(parent);
                } else if (workspace_d[parent] == 0) {
                    worklist_deg0.push_back(parent);
                }*/
                continue;
            } else {
                parentlist.reset();
                is_parent.reset();
                for (int i = v_child_col; i < v_child_col + child_col_sz; ++i) {
                    int child = c->lab[i];
                    assert(workspace_d[child] == 1);
                    // TODO: surely need special code for pairs in same color class
                    rec->del1.set(child);

                    // search for parent
                    const int e_pos_child = g->v[child];
                    int parent = g->e[e_pos_child];
                    int search_parent = 0;
                    while (rec->del1.get(parent)) {
                        ++search_parent;
                        parent = g->e[e_pos_child + search_parent];
                    }

                    // save canonical info for parent
                    childlist.arr[g->v[parent] + childcount.arr[parent]] = child;
                    ++childcount.arr[parent];

                    if (!is_parent.get(parent)) {
                        is_parent.set(parent);
                        childcount_prev.arr[parent] = childcount.arr[parent] - 1;
                        parentlist.push_back(parent);
                    }

                    // adjust parent degree
                    workspace_d[parent] -= 1;
                    if (workspace_d[parent] == 1) {
                        worklist_deg1.push_back(parent);
                    } else if (workspace_d[parent] == 0) {
                        worklist_deg0.push_back(parent);
                    }

                    assert(workspace_d[parent] >= 0);
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
                map.reset();
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
                    automorphism.arr[map.arr[pos]] = child_to;
                    automorphism.arr[child_to]     = map.arr[pos];
                    automorphism_supp.push_back(map.arr[pos]);
                    automorphism_supp.push_back(child_to);
                    ++pos;
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
                            assert(automorphism.arr[map.arr[pos]] == map.arr[pos]);
                            assert(automorphism.arr[next] == next);
                            automorphism.arr[map.arr[pos]] = next;
                            automorphism.arr[next]     = map.arr[pos];
                            automorphism_supp.push_back(map.arr[pos]);
                            automorphism_supp.push_back(next);
                            ++pos;
                            if(from_next != to_next);
                                stack1.push_back(std::pair<int, int>(from_next, to_next));
                        }
                    }

                    assert(pos == map.cur_pos);

                    assert(workspace_d[child_to] == 1);
                    assert(workspace_d[child_from] == 1);
                    assert(rec->del1.get(child_to));
                    assert(rec->del1.get(child_from));
                    consume(g->v_size, automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                    reset_automorphism(automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                    automorphism_supp.reset();
                }
            }
            parentlist.reset();
        }

        if(g->v_size > 1) {
            // search for parents that still remain in the graph, and rewrite childlist structure into canonical string
            for (int i = 0; i < g->v_size; ++i) {
                if (!rec->del1.get(i) && childcount.arr[i] >= 0) {
                    stack1.reset();
                    stack1.push_back(std::pair<int, int>(g->v[i], g->v[i] + childcount.arr[i]));
                    rec->canonical_recovery_string[i].reserve(childcount.arr[i]);
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
                            rec->canonical_recovery_string[i].push_back(next);
                            stack1.push_back(std::pair<int, int>(from, to));
                            stack1.push_back(std::pair<int, int>(from_next, to_next));
                        }
                    }
                }
            }
        }


        while(!worklist_deg0.empty()) {
            const int v_child = worklist_deg0.pop_back();
            if (rec->del1.get(v_child))
                continue;
            assert(workspace_d[v_child] == 0);

            is_parent.reset();
            const int v_child_col = c->vertex_to_col[v_child];
            const int child_col_sz = c->ptn[v_child_col] + 1;
            //std::cout << "col " << v_child_col << "(" << child_col_sz << ")" << " from vertex " << v_child << std::endl;
            const int parent_from = c->lab[v_child_col];
            rec->del1.set(parent_from);

            if(child_col_sz == 1)
                continue;
            int j = 2;
            for (int i = v_child_col + 1; i < v_child_col + child_col_sz; ++i) {
                multiply_to_group_size(j);
                ++j;
                const int parent_to = c->lab[i];
                assert(workspace_d[parent_to] == 0);
                rec->del1.set(parent_to);

                assert(rec->canonical_recovery_string[parent_to].size() == rec->canonical_recovery_string[parent_from].size());
                automorphism.arr[parent_to]   = parent_from;
                automorphism.arr[parent_from] = parent_to;
                automorphism_supp.push_back(parent_from);
                automorphism_supp.push_back(parent_to);
                for(int i = 0; i < rec->canonical_recovery_string[parent_to].size(); ++i) {
                    const int str_from = rec->canonical_recovery_string[parent_from][i];
                    const int str_to = rec->canonical_recovery_string[parent_to][i];
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
        // add translation layers for this deletion
        rec->backward_translation_layers.emplace_back(std::vector<int>());
        const int back_ind = rec->backward_translation_layers.size()-1;
        rec->translation_layers.emplace_back(std::vector<int>());
        const int fwd_ind = rec->translation_layers.size()-1;

        // copy some stuff
        std::vector<int> g_old_v;
        std::vector<int> g_old_e;
        std::vector<int> old_colmap;
        std::vector<int> g_old_d2;
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
        int new_vsize = 0;
        rec->translation_layers[fwd_ind].reserve(g->v_size);
        rec->backward_translation_layers[back_ind].reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            if (!rec->del1.get(i)) {
                rec->translation_layers[fwd_ind].push_back(cnt);
                rec->backward_translation_layers[back_ind].push_back(rec->translation_layers[fwd_ind].size() - 1);
                ++cnt;
                ++new_vsize;
            } else {
                rec->translation_layers[fwd_ind].push_back(-1);
            }
        }

        // make graph smaller using the translation array
        int epos  = 0;
        for (int i = 0; i < g->v_size; ++i) {
            const int old_v = i;
            const int new_v = rec->translation_layers[fwd_ind][i];

            if (new_v >= 0) {
                colmap[new_v] = old_colmap[old_v];
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + g_old_d2[old_v]; ++j) {
                    const int ve = g_old_e[j];
                    const int new_ve = rec->translation_layers[fwd_ind][ve];
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


        PRINT("(prep-red) shrinking graph "<< g->v_size << "->" << new_vsize);

        g->v_size = new_vsize;
        g->d_size = new_vsize;

        // TODO: remove discrete parts of graph
    }

    void perform_del2(sgraph*g, int* colmap, recovery_map* rec) {
        if(g->v_size <= 1)
            return;

        // add translation layers for this deletion
        rec->backward_translation_layers.emplace_back(std::vector<int>());
        const int back_ind = rec->backward_translation_layers.size()-1;
        rec->translation_layers.emplace_back(std::vector<int>());
        const int fwd_ind = rec->translation_layers.size()-1;

        // copy some stuff
        std::vector<int> g_old_v;
        std::vector<int> g_old_e;
        std::vector<int> g_old_d2;
        std::vector<int> old_colmap;
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

        // create translation array from old graph to new graph vertices
        int cnt = 0;
        int new_vsize = 0;
        rec->translation_layers[fwd_ind].clear();
        rec->backward_translation_layers[back_ind].clear();
        rec->translation_layers[fwd_ind].reserve(g->v_size);
        rec->backward_translation_layers[back_ind].reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            if (!rec->del2.get(i)) {
                rec->translation_layers[fwd_ind].push_back(cnt);
                rec->backward_translation_layers[back_ind].push_back(rec->translation_layers[fwd_ind].size() - 1);
                ++cnt;
                ++new_vsize;
            } else {
                rec->translation_layers[fwd_ind].push_back(-1);
            }
        }

        // make graph smaller using the translation array
        int epos  = 0;
        for (int i = 0; i < g->v_size; ++i) {
            const int old_v = i;
            const int new_v = rec->translation_layers[fwd_ind][i];

            if (new_v >= 0) {
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + g_old_d2[old_v]; ++j) {
                    const int ve = g_old_e[j];
                    const int new_ve = rec->translation_layers[fwd_ind][ve];
                    if(new_ve >= 0) {
                        assert(new_ve < new_vsize);
                        assert(new_ve >= 0);
                        ++new_d;
                        g->e[epos] = new_ve;
                        ++epos;
                    }
                }
                if(add_edge_buff_act.get(old_v)) {
                    const int edge_to_old = add_edge_buff.arr[old_v];
                    assert(add_edge_buff_act.get(edge_to_old));
                    assert(add_edge_buff.arr[edge_to_old] >= 0);
                    assert(add_edge_buff.arr[edge_to_old] == old_v);
                    // TODO can i insert this edge for to_old as well here?
                    const int edge_to_new = rec->translation_layers[fwd_ind][edge_to_old];
                    assert(edge_to_old >= 0);
                    assert(edge_to_old <= rec->domain_size);
                    assert(edge_to_new >= 0);
                    assert(edge_to_new <= new_vsize);
                    //std::cout << "adding edge " << old_v << "<->" << edge_to_old << std::endl;
                    ++new_d;
                    g->e[epos] = edge_to_new;
                    ++epos;
                }
                g->d[new_v] = new_d;
            }
        }

        PRINT("(prep-red) shrinking graph esize "<< g->e_size << "->" << epos);
        g->e_size = epos;


        PRINT("(prep-red) shrinking graph "<< g->v_size << "->" << new_vsize);

        // adapt colmap for remaining vertices
        if(colmap != nullptr) {
            for (int i = 0; i < g->v_size; ++i) {
                const int old_v = i;
                const int new_v = rec->translation_layers[fwd_ind][i];
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
            assert(g->v[i] < g->e_size);
            assert(g->v[i] >= 0);
            assert(g->d[i] >= 0);
            assert(g->d[i] < g->v_size);
        }
        for(int i = 0; i < g->e_size; ++i) {
            assert(g->e[i] < g->v_size);
            assert(g->e[i] >= 0);
        }

        add_edge_buff_act.reset();
        // TODO: remove discrete parts of graph
    }

    void perform_del_discrete(sgraph*g, int* colmap, recovery_map* rec) {
        if(g->v_size <= 1)
            return;

        // add translation layers for this deletion
        rec->backward_translation_layers.emplace_back(std::vector<int>());
        const int back_ind = rec->backward_translation_layers.size()-1;
        rec->translation_layers.emplace_back(std::vector<int>());
        const int fwd_ind = rec->translation_layers.size()-1;

        work_list_t<int> color_count;
        color_count.initialize(rec->domain_size);
        int discrete_cnt = 0;
        for(int i = 0; i < rec->domain_size; ++i) {
            color_count.push_back(0);
        }
        for(int i = 0; i < g->v_size; ++i) {
            assert(colmap[i] < rec->domain_size);
            color_count.arr[colmap[i]]++;
        }
        /*for(int i = 0; i < g->v_size; ++i) {
            if(color_count.arr[colmap[i]] == 1) {
                ++discrete_cnt;
            }
        }
        std::cout << "(prep-red) discrete vertices: " <<discrete_cnt << std::endl;
        */

        // copy some stuff
        std::vector<int> g_old_v;
        std::vector<int> g_old_e;
        std::vector<int> g_old_d2;
        std::vector<int> old_colmap;
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
        rec->translation_layers[fwd_ind] = std::vector<int>();
        rec->backward_translation_layers[back_ind] = std::vector<int>();
        rec->translation_layers[fwd_ind].reserve(g->v_size);
        rec->backward_translation_layers[back_ind].reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            if (color_count.arr[colmap[i]] != 1) {
                rec->translation_layers[fwd_ind].push_back(cnt);
                rec->backward_translation_layers[back_ind].push_back(rec->translation_layers[fwd_ind].size() - 1);
                ++cnt;
                ++new_vsize;
            } else {
                rec->translation_layers[fwd_ind].push_back(-1);
            }
        }

        // make graph smaller using the translation array
        int epos  = 0;
        for (int i = 0; i < g->v_size; ++i) {
            const int old_v = i;
            assert(i < rec->translation_layers[fwd_ind].size());
            const int new_v = rec->translation_layers[fwd_ind][i];

            if (new_v >= 0) {
                // copy old colmap for remaining vertices
                colmap[new_v] = old_colmap[old_v];
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + g_old_d2[old_v]; ++j) {
                    const int ve = g_old_e[j];
                    assert(ve < rec->translation_layers[fwd_ind].size());
                    const int new_ve = rec->translation_layers[fwd_ind][ve];
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

        PRINT("(prep-red) shrinking graph esize "<< g->e_size << "->" << epos);
        g->e_size = epos;


        PRINT("(prep-red) shrinking graph "<< g->v_size << "->" << new_vsize);

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
            assert(g->v[i] < g->e_size);
            assert(g->v[i] >= 0);
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
        automorphism.initialize(g->v_size);
        for(int i = 0; i < g->v_size; ++i) {
            automorphism.push_back(i);
        }
        automorphism_supp.initialize(g->v_size);

        int deg0 = 0;
        int deg1 = 0;
        int deg2 = 0;
        for(int i = 0; i < g->v_size; ++i) {
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
        std::cout << "(pre-red) before reduction " << deg0 << ", "  << deg1 << ", "  << deg2 << std::endl;

        rec.domain_size = g->v_size;
        add_edge_buff.initialize(rec.domain_size);
        for(int i = 0; i < rec.domain_size; ++i)
            add_edge_buff.push_back(-1);
        add_edge_buff_act.initialize(rec.domain_size);
        // assumes colmap is array of length g->v_size
        rec.del1 = mark_set();
        rec.del1.initialize(g->v_size);
        rec.del2 = mark_set();
        rec.del2.initialize(g->v_size);
        rec.canonical_recovery_string.reserve(g->v_size);
        for(int i = 0; i < rec.domain_size; ++i)
            rec.canonical_recovery_string.emplace_back(std::vector<int>());
        //store_graph_aspects(g, &rec);
        // singleton-only refinement, then cut graph
        refinement<int, int, int> R;
        g->initialize_coloring(&c, colmap);
        R.refine_coloring_first(g, &c, -1);

        // elimination 1: degree 1 + 0
        //red_deg10_assume_cref(g, &c, &rec, consume);
        copy_coloring_to_colmap(&c, colmap);
        //perform_del(g, colmap, &rec);

        // elimination 2: degree 0

        // invariants 1: paths of length 2 for large regular components

        // elimination 4: discrete colors
        //perform_del_discrete(g, colmap, &rec);

        // elimination 3: degree 2
        red_deg2_assume_cref(g, colmap, &rec, consume);
        perform_del2(g, colmap, &rec);
        //rec.del2.reset();

        //red_deg2_assume_cref(g, colmap, &rec, consume); // deletes an automorphism for some reason
        //perform_del2(g, colmap, &rec);

        deg0 = 0;
        deg1 = 0;
        deg2 = 0;
        for(int i = 0; i < g->v_size; ++i) {
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
        std::cout << "(pre-red) after reduction " << deg0 << ", "  << deg1 << ", "  << deg2 << std::endl;
        std::cout << "(pre-red) group size: " << rec.base << "*10^" << rec.exp << std::endl;

        g->sanity_check();

        // TODO: multiple calls for now obviously independent components -- could use "buffer consumer" to translate domains
        // TODO: just use consumer for all the back-translation (just dont "restore"!): compactify translation layers here for this
        // TODO: the consumer just uses standardized "canonized" strings for the backward-reductions, no matter the reduction used
    }

    void restore(sgraph* g, automorphism_info* a, dejavu_consumer consume) {
        //sparse_group group;
        //group.domain_size = rec.domain_size;
        //group.initialize_from_automorphism_info(a, &rec);
    }
};

#endif //DEJAVU_PREP_H
