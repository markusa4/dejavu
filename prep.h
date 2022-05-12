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

    std::vector<int> translate1;
    std::vector<int> backwards_translate1;

    std::vector<int> translate2;
    std::vector<int> backwards_translate2;

    std::vector<int> translate3;
    std::vector<int> backwards_translate3;

    std::vector<std::vector<int>> canonical_recovery_string;
};

class sparse_automorphism {
private:
    int cnt_supp = 0;
    int support  = 0;
    std::vector<int> perm_from;
    std::vector<int> perm_to;
public:
    sparse_automorphism(int supp) {
        support = supp;
        perm_from.reserve(supp);
        perm_to.reserve(supp);
    }
    void add_map(int from, int to) {
        assert(cnt_supp < support);
        perm_from.push_back(from);
        perm_to.push_back(to);
        ++cnt_supp;
    }
};

class sparse_group {
public:
    int domain_size = 0;
    std::vector<int> orbit;
    std::vector<sparse_automorphism> gens_sparse;

    void add_gen_dense(int* p) {

    }

    void add_gen_sparse(sparse_automorphism* gen) {

    }

    void initialize_from_automorphism_info(automorphism_info* a, recovery_map* rec) {
        shared_permnode* perms_first = *a->gens;
        shared_permnode* perms_next = perms_first;
        double average_supp = 0;
        int gens = 0;
        do {
            ++gens;
            int supp = 0;
            for(int i = 0; i < perms_next->nalloc; ++i) {
                if(i != perms_next->p[i])
                    ++supp;
            }
            average_supp += supp;
            perms_next = perms_next->next;
        } while(perms_next != perms_first);
        std::cout << "(prep-res) gens: " << gens << std::endl;
        std::cout << "(prep-res) average support: " << average_supp / gens << std::endl;
    }
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

    void red_deg2_assume_cref(sgraph* g, int* colmap, recovery_map* rec, dejavu_consumer consume) {
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
                if(smallest_endpoint_cnt.arr[endpoint1] == 1 && smallest_endpoint_cnt.arr[endpoint2] == 1) {
                    assert(g->d[endpoint1] != 2);
                    assert(g->d[endpoint2] != 2);
                    unique_smallest_endpoint_paths += 1;
                    add_edge_buff.arr[endpoint1] = endpoint2;
                    add_edge_buff_act.set(endpoint1);
                    add_edge_buff.arr[endpoint2] = endpoint1;
                    add_edge_buff_act.set(endpoint2);
                    assert(path_length == path.cur_pos);
                    int deleted = 0;
                    while(!path.empty()) {
                        ++deleted;
                        const int del_v = path.pop_back();
                        assert(g->d[del_v] == 2);
                        assert(!rec->del2.get(del_v));
                        rec->del2.set(del_v);
                    }
                }
            }
        }

        // TODO: maybe add multi-path removal, etc.
        assert(num_paths1 == num_paths2);
        if(num_paths1 != 0) {
            std::cout << "(prep-red) total paths: " << num_paths1 << " (avg length: " << total_path_length / num_paths1
                      << ")" << std::endl;
            std::cout << "(prep-red) unique smallest endpoints: " << unique_smallest_endpoint_paths << std::endl;
        }
    }

    void red_deg10_assume_cref(sgraph* g, const coloring<int>* c, recovery_map* rec, dejavu_consumer consume) {
        std::cout << "(prep-red) deg1" << std::endl;
        rec->del1.reset();
        rec->translate1.clear();
        rec->backwards_translate1.clear();

        work_list_t<int> worklist_deg0;
        work_list_t<int> worklist_deg1;

        mark_set is_parent;
        is_parent.initialize(g->v_size);

        work_list_t<int> automorphism;
        automorphism.initialize(g->v_size);
        for(int i = 0; i < g->v_size; ++i) {
            automorphism.push_back(i);
        }

        std::vector<int> workspace_d;
        workspace_d.reserve(g->v_size);
        for(int i = 0; i < g->v_size; ++i) {
            workspace_d.push_back(g->d[i]);
        }

        work_list_t<int> automorphism_supp;
        automorphism_supp.initialize(g->v_size);

        work_list_t<int> parentlist;
        parentlist.initialize(g->v_size);
        work_list_t<int> childlist;
        childlist.initialize(g->e_size);
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

            is_parent.reset();
            const int v_child_col = c->vertex_to_col[v_child];
            const int child_col_sz = c->ptn[v_child_col] + 1;
            //std::cout << "col " << v_child_col << "(" << child_col_sz << ")" << " from vertex " << v_child << std::endl;
            for(int i = v_child_col; i < v_child_col + child_col_sz; ++i) {
                const int child  = c->lab[i];
                assert(workspace_d[child] == 1);
                // TODO: surely need special code for pairs in same color class

                // search for parent
                int parent = g->e[g->v[child]];
                int search_parent = 0;
                while(rec->del1.get(parent)) {
                    ++search_parent;
                    parent = g->e[g->v[child] + search_parent];
                }

                // remove children and save canonical info for parent
                childlist.arr[g->v[parent] + childcount.arr[parent]] = child;
                ++childcount.arr[parent];

                if(!is_parent.get(parent)) {
                    is_parent.set(parent);
                    childcount_prev.arr[parent] = childcount.arr[parent] - 1;
                    parentlist.push_back(parent);
                }

                rec->del1.set(child);

                // adjust parent degree
                workspace_d[parent] -= 1;
                if(workspace_d[parent] == 1) {
                    worklist_deg1.push_back(parent);
                } else if(workspace_d[parent] == 0) {
                    worklist_deg0.push_back(parent);
                }

                assert(workspace_d[parent] >= 0);
            }

            while(!parentlist.empty()) {
                const int parent = parentlist.pop_back();
                const int childcount_from = childcount_prev.arr[parent];
                const int childcount_to   = childcount.arr[parent];
                // automorphism 1: long cycle (c1 ... cn)
                if(childcount_to - childcount_from == 1)
                    continue;
                int child_from = childlist.arr[g->v[parent]];

                // descending tree of child_from while writing map
                stack1.reset();
                map.push_back(child_from);
                stack1.push_back(std::pair<int, int>(g->v[child_from], g->v[child_from] + childcount.arr[child_from]));
                while(!stack1.empty()) {
                    std::pair<int, int> from_to = stack1.pop_back();
                    int from = from_to.first;
                    const int to   = from_to.second;
                    if(from == to) {
                        continue;
                    } else {
                        const int next      = childlist.arr[from];
                        const int from_next = g->v[next];
                        const int to_next   = g->v[next] + childcount.arr[next];
                        ++from;
                        map.push_back(next);
                        stack1.push_back(std::pair<int, int>(from, to));
                        stack1.push_back(std::pair<int, int>(from_next, to_next));
                    }
                }

                for(int i = childcount_from + 1; i < childcount_to; ++i) {
                    const int child_to = childlist.arr[g->v[parent] + i];
                    int pos = 0;

                    // descending tree of child_to while writing automorphism
                    stack1.reset();
                    assert(map.arr[pos] != child_to);
                    automorphism.arr[map.arr[pos]] = child_to;
                    automorphism.arr[child_to]     = map.arr[pos];
                    automorphism_supp.push_back(map.arr[pos]);
                    automorphism_supp.push_back(child_to);
                    ++pos;
                    stack1.push_back(std::pair<int, int>(g->v[child_from], g->v[child_from] + childcount.arr[child_from]));
                    while(!stack1.empty()) {
                        std::pair<int, int> from_to = stack1.pop_back();
                        int from = from_to.first;
                        const int to   = from_to.second;
                        if(from == to) {
                            continue;
                        } else {
                            const int next      = childlist.arr[from];
                            const int from_next = g->v[next];
                            const int to_next   = g->v[next] + childcount.arr[next];
                            ++from;
                            assert(map.arr[pos] != next);
                            assert(automorphism.arr[map.arr[pos]] == map.arr[pos]);
                            assert(automorphism.arr[next] == next);
                            automorphism.arr[map.arr[pos]] = next;
                            automorphism.arr[next]     = map.arr[pos];
                            automorphism_supp.push_back(map.arr[pos]);
                            automorphism_supp.push_back(next);
                            ++pos;
                            stack1.push_back(std::pair<int, int>(from, to));
                            stack1.push_back(std::pair<int, int>(from_next, to_next));
                        }
                    }

                    assert(workspace_d[child_to] == 1);
                    assert(workspace_d[child_from] == 1);
                    assert(rec->del1.get(child_to));
                    assert(rec->del1.get(child_from));
                    consume(g->v_size, automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                    reset_automorphism(automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                }
            }
        }

        // search for parents that still remain in the graph, and rewrite childlist structure into canonical string
        for(int i = 0; i < g->v_size; ++i) {
            if(!rec->del1.get(i) && childcount.arr[i] >= 0) {
                stack1.reset();
                stack1.push_back(std::pair<int, int>(g->v[i], g->v[i] + childcount.arr[i]));
                rec->canonical_recovery_string[i].reserve(childcount.arr[i]);
                while(!stack1.empty()) {
                    std::pair<int, int> from_to = stack1.pop_back();
                    int from = from_to.first;
                    const int to   = from_to.second;
                    if(from == to) {
                        continue;
                    } else {
                        const int next      = childlist.arr[from];
                        const int from_next = g->v[next];
                        const int to_next   = g->v[next] + childcount.arr[next];
                        ++from;
                        rec->canonical_recovery_string[i].push_back(next);
                        stack1.push_back(std::pair<int, int>(from, to));
                        stack1.push_back(std::pair<int, int>(from_next, to_next));
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

            for (int i = v_child_col + 1; i < v_child_col + child_col_sz; ++i) {
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
        // copy some stuff
        std::vector<int> g_old_v;
        std::vector<int> g_old_e;
        std::vector<int> old_colmap;
        std::vector<int> g_old_d2;
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
        rec->translate1.reserve(g->v_size);
        rec->backwards_translate1.reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            if (!rec->del1.get(i)) {
                rec->translate1.push_back(cnt);
                rec->backwards_translate1.push_back(rec->translate1.size() - 1);
                ++cnt;
                ++new_vsize;
            } else {
                rec->translate1.push_back(-1);
            }
        }

        // make graph smaller using the translation array
        int epos  = 0;
        for (int i = 0; i < g->v_size; ++i) {
            const int old_v = i;
            const int new_v = rec->translate1[i];

            if (new_v >= 0) {
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + g_old_d2[old_v]; ++j) {
                    const int ve = g_old_e[j];
                    const int new_ve = rec->translate1[ve];
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

        // adapt colmap for remaining vertices
        if(colmap != nullptr) {
            for (int i = 0; i < g->v_size; ++i) {
                const int old_v = i;
                const int new_v = rec->translate1[i];
                assert(new_v < new_vsize);
                if (new_v >= 0) {
                    colmap[new_v] = old_colmap[old_v];
                }
            }
        }

        g->v_size = new_vsize;
        g->d_size = new_vsize;

        // TODO: remove discrete parts of graph
    }

    void perform_del2(sgraph*g, int* colmap, recovery_map* rec) {
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
        rec->translate2.reserve(g->v_size);
        rec->backwards_translate2.reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            if (!rec->del2.get(i)) {
                rec->translate2.push_back(cnt);
                rec->backwards_translate2.push_back(rec->translate2.size() - 1);
                ++cnt;
                ++new_vsize;
            } else {
                rec->translate2.push_back(-1);
            }
        }

        // make graph smaller using the translation array
        int epos  = 0;
        for (int i = 0; i < g->v_size; ++i) {
            const int old_v = i;
            const int new_v = rec->translate2[i];

            if (new_v >= 0) {
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + g_old_d2[old_v]; ++j) {
                    const int ve = g_old_e[j];
                    const int new_ve = rec->translate2[ve];
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
                    const int edge_to_new = rec->translate2[edge_to_old];
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
                const int new_v = rec->translate2[i];
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
        work_list_t<int> color_count;
        color_count.initialize(rec->domain_size);
        int discrete_cnt = 0;
        for(int i = 0; i < rec->domain_size; ++i) {
            color_count.push_back(0);
        }
        for(int i = 0; i < g->v_size; ++i) {
            color_count.arr[colmap[i]]++;
        }
        for(int i = 0; i < rec->domain_size; ++i) {
            if(color_count.arr[colmap[i]] == 1) {
                ++discrete_cnt;
            }
        }
        std::cout << "(prep-red) discrete vertices: " <<discrete_cnt << std::endl;

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
        rec->translate3 = std::vector<int>();
        rec->backwards_translate3 = std::vector<int>();
        rec->translate3.reserve(g->v_size);
        rec->backwards_translate3.reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            if (color_count.arr[colmap[i]] != 1) {
                rec->translate3.push_back(cnt);
                rec->backwards_translate3.push_back(rec->translate3.size() - 1);
                ++cnt;
                ++new_vsize;
            } else {
                rec->translate3.push_back(-1);
            }
        }

        // make graph smaller using the translation array
        int epos  = 0;
        for (int i = 0; i < g->v_size; ++i) {
            const int old_v = i;
            assert(i < rec->translate3.size());
            const int new_v = rec->translate3[i];

            if (new_v >= 0) {
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + g_old_d2[old_v]; ++j) {
                    const int ve = g_old_e[j];
                    assert(ve < rec->translate3.size());
                    const int new_ve = rec->translate3[ve];
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

        PRINT("(prep-red) shrinking graph esize "<< g->e_size << "->" << epos);
        g->e_size = epos;


        PRINT("(prep-red) shrinking graph "<< g->v_size << "->" << new_vsize);

        // copy old colmap for remaining vertices
        if(colmap != nullptr) {
            for (int i = 0; i < g->v_size; ++i) {
                const int old_v = i;
                const int new_v = rec->translate3[i];
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
    }
public:
    recovery_map rec;
    coloring<int> c;
    void reduce(sgraph* g, int* colmap, dejavu_consumer consume) {
        rec.domain_size = g->v_size;
        add_edge_buff.initialize(rec.domain_size);
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
        red_deg10_assume_cref(g, &c, &rec, consume);
        copy_coloring_to_colmap(&c, colmap);
        perform_del(g, colmap, &rec);

        // elimination 2: degree 0

        // elimination 3: degree 2
        red_deg2_assume_cref(g, colmap, &rec, consume);
        perform_del2(g, colmap, &rec);
        // invariants 1: paths of length 2 for large regular components

        // elimination 4: discrete colors
        perform_del_discrete(g, colmap, &rec);

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
        std::cout << "(pre-red) after reduction " << deg0 << ", "  << deg1 << ", "  << deg2 << std::endl;

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
