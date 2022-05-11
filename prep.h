#ifndef DEJAVU_PREP_H
#define DEJAVU_PREP_H

#include "sgraph.h"
#include "refinement.h"
#include "schreier_shared.h"
#include <vector>

struct recovery_map {
    mark_set del;
    std::vector<int> translate1;
    std::vector<int> backwards_translate1;
    std::vector<int> translate2;
    std::vector<int> backwards_translate2;
    std::vector<int> old_d;

    std::vector<int> deg1_to_parent_ptn;
    std::vector<int> elim_order;
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
private:
    int initial_vsize = 0;
    void red_singleton_only_refinement() {

    }

    void res_singleton_only_refinement() {

    }

    void reset_automorphism(int* automorphism, int nsupp, int* supp) {
        for(int i = 0; i < nsupp; ++i) {
            automorphism[supp[i]] = supp[i];
        }
    };

    void red_deg2_assume_cref(sgraph* g, int* colmap, recovery_map* rec, dejavu_consumer consume) {
        work_list_t<int> worklist_deg2;
        worklist_deg2.initialize(g->v_size);
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

        int num_paths = 0;
        double total_path_length = 0;
        while(!worklist_deg2.empty()) {
            const int v_child = worklist_deg2.pop_back();
            if (path_done.get(v_child))
                continue;
            if (g->d[v_child] != 2)
                continue;

            num_paths += 1;
            path_done.set(v_child);

            int path_length = 1;
            path.initialize(v_child); // probably need 2*n and set v_child into middle, then expand to left/right

            int v_child_next = g->e[g->v[v_child]];
            int endpoint1, endpoint2;
            bool cycle = false;
            while(true) {
                if(g->d[v_child_next] != 2) {
                    // v_child_next is endpoint1
                    endpoint1 = v_child_next;
                    break;
                }
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
                v_child_next = v_child_next_next;
            }

            v_child_next = g->e[g->v[v_child] + 1];
            while(true) {
                if(g->d[v_child_next] != 2) {
                    // v_child_next is endpoint2
                    endpoint2 = v_child_next;
                    break;
                }
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
                v_child_next = v_child_next_next;
            }

            total_path_length += path_length;

            // TODO: some practical assumptions that could be used: colors of endpoints distinct, endpoints only connected to 1 path, path has length 1

            // UNIQUE ENDPOINT REDUCTION
            // TODO: "unique path-node representation": if an endpoint is only connected to 1 path, then if the endpoint is moved to another, the paths are interchanged as well
            //  -- hence, the path can just be represented by one of its endpoints, e.g., the one with the smaller color
            if(path_length == 1) {
                if(colmap[endpoint1] != colmap[endpoint2]) {
                    std::cout << "distinct color path from " << endpoint1 << " -> " << endpoint2 << " length " << path_length << std::endl;
                }
            }

            if(!cycle) {
                //std::cout << "path from " << endpoint1 << " -> " << endpoint2 << " length " << path_length << std::endl;
            } else {

                //std::cout << "cycle length " << path_length << std::endl;
            }

            // TODO: shrink multi-paths to one? then just leave out the paths?
        }
        if(num_paths != 0)
            std::cout << "(prep-red) total paths: " << num_paths << " (avg length: " << total_path_length / num_paths << ")" << std::endl;

    }

    void red_deg1_assume_cref(sgraph* g, coloring<int>* c, recovery_map* rec, dejavu_consumer consume) {
        work_list_t<int> worklist_deg1;

        mark_set is_parent;
        is_parent.initialize(g->v_size);

        work_list_t<int> automorphism;
        automorphism.initialize(g->v_size);
        for(int i = 0; i < g->v_size; ++i) {
            automorphism.push_back(i);
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

        worklist_deg1.initialize(g->v_size);

        rec->deg1_to_parent_ptn.reserve(g->v_size);
        rec->elim_order.reserve(g->v_size);

        for(int i = 0; i < g->v_size; ++i) {
            rec->deg1_to_parent_ptn.push_back(i);
        }

        for(int i = 0; i < g->v_size; ++i) {
            switch(g->d[i]) {
                case 1:
                    worklist_deg1.push_back(i);
                    break;
                default:
                    break;
            }
        }

        while(!worklist_deg1.empty()) {
            const int v_child = worklist_deg1.pop_back();
            if(rec->del.get(v_child))
                continue;
            if(g->d[v_child] != 1)
                continue;

            is_parent.reset();
            const int v_child_col = c->vertex_to_col[v_child];
            const int child_col_sz = c->ptn[v_child_col] + 1;
            //std::cout << "col " << v_child_col << "(" << child_col_sz << ")" << " from vertex " << v_child << std::endl;
            for(int i = v_child_col; i < v_child_col + child_col_sz; ++i) {
                const int child  = c->lab[i];
                assert(g->d[child] == 1);
                // TODO: surely need special code for pairs in same color class

                // search for parent
                int parent = g->e[g->v[child]];
                int search_parent = 0;
                while(rec->del.get(parent)) {
                    ++search_parent;
                    parent = g->e[g->v[child] + search_parent];
                }
                //

                // remove children and save info for parent
                // TODO: always call consumer for "obvious" automorphisms eagerly?
                // TODO: later on only the automorphisms that map unresolved parents must be expanded, i.e., the tree structures represented by them expanded
                // TODO: but I also need this "expansion" here... expand parent -> a canonical ordering of represented vertices that is the automorphism expansion? maybe save this canonical ordering in upward fashion / a readable fashion? one n-array suffices?
                rec->deg1_to_parent_ptn[child] = parent; // TODO: use g->v to determine space in the array to write children in canonical order to -- it is automatically canonical because we are proceeding color-wise?
                                                         // TODO: can I somehow use colors / individualization instead?
                childlist.arr[g->v[parent] + childcount.arr[parent]] = child;
                ++childcount.arr[parent];

                if(!is_parent.get(parent)) {
                    is_parent.set(parent);
                    childcount_prev.arr[parent] = childcount.arr[parent] - 1;
                    parentlist.push_back(parent);
                }

                rec->del.set(child);

                // adjust parent degree
                g->d[parent] -= 1;
                if(g->d[parent] == 1) {
                    worklist_deg1.push_back(parent);
                }

                assert(g->d[parent] >= 0);
            }
            while(!parentlist.empty()) {
                const int parent = parentlist.pop_back();
                const int childcount_from = childcount_prev.arr[parent];
                const int childcount_to   = childcount.arr[parent];
                // automorphism 1: long cycle (c1 ... cn)
                if(childcount_to - childcount_from == 1)
                    continue;
                int child_from = childlist.arr[g->v[parent] + childcount_to - 1];
                for(int i = childcount_from; i < childcount_to; ++i) {
                    const int child_to = childlist.arr[g->v[parent] + i];
                    assert(g->d[child_to] == 1);
                    assert(g->d[child_from] == 1);
                    assert(rec->del.get(child_to));
                    assert(rec->del.get(child_from));
                    // TODO: map child_from -> child_to in automorphism
                    consume(g->v_size, automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                    reset_automorphism(automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                }
            }
        }
    }

    void store_graph_aspects(sgraph* g, recovery_map* rec) {
        rec->old_d.clear();
        rec->old_d.reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            rec->old_d.push_back(g->d[i]);
        }
    }

    void perform_del(sgraph*g, int* colmap, recovery_map* rec) {
        initial_vsize = g->v_size;

        // copy some stuff
        std::vector<int> g_old_v;
        std::vector<int> g_old_e;
        std::vector<int> old_colmap;
        g_old_v.reserve(g->v_size);
        g_old_e.reserve(g->e_size);
        if(colmap != nullptr) {
            old_colmap.reserve(g->v_size);
            for (int i = 0; i < g->v_size; ++i) {
                old_colmap.push_back(colmap[i]);
            }
        }
        for (int i = 0; i < g->v_size; ++i) {
            g_old_v.push_back(g->v[i]);
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
            if (!rec->del.get(i)) {
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
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + rec->old_d[old_v]; ++j) {
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
                    colmap[new_v] = c.vertex_to_col[old_v];
                }
            }
        }

        g->v_size = new_vsize;
        g->d_size = new_vsize;

        // TODO: remove discrete parts of graph
    }

    void perform_del_discrete(sgraph*g, int* colmap, recovery_map* rec) {
        work_list_t<int> color_count;
        color_count.initialize(initial_vsize);
        int discrete_cnt = 0;
        for(int i = 0; i < initial_vsize; ++i) {
            color_count.push_back(0);
        }
        for(int i = 0; i < g->v_size; ++i) {
            color_count.arr[colmap[i]]++;
        }
        for(int i = 0; i < initial_vsize; ++i) {
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
        rec->translate2.reserve(g->v_size);
        rec->backwards_translate2.reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            if (color_count.arr[colmap[i]] != 1) {
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
            assert(i < rec->translate2.size());
            const int new_v = rec->translate2[i];

            if (new_v >= 0) {
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + g_old_d2[old_v]; ++j) {
                    const int ve = g_old_e[j];
                    assert(ve < rec->translate2.size());
                    const int new_ve = rec->translate2[ve];
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
                const int new_v = rec->translate2[i];
                assert(new_v < new_vsize);
                if (new_v >= 0) {
                    assert(old_colmap[old_v] >= 0);
                    assert(old_colmap[old_v] < initial_vsize);
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

    void res_low_deg() {
        // deg1:
        // restoration: each generator that permutes parents is extended to also permute X_1 -> X_2
        // restoration: add generators for symmetrical group for vertices of each parent
    }
public:
    recovery_map rec;
    coloring<int> c;
    void reduce(sgraph* g, int* colmap, dejavu_consumer consume) {
        // assumes colmap is array of length g->v_size
        rec.del = mark_set();
        rec.del.initialize(g->v_size);
        store_graph_aspects(g, &rec);
        // singleton-only refinement, then cut graph
        refinement<int, int, int> R;
        g->initialize_coloring(&c, colmap);
        R.refine_coloring_first(g, &c, -1);

        // elimination 1: degree 1
        red_deg1_assume_cref(g, &c, &rec, consume);
        perform_del(g, colmap, &rec);
        // elimination 2: degree 0

        // elimination 3: degree 2
        red_deg2_assume_cref(g, colmap, &rec, consume);

        // invariants 1: paths of length 2 for large regular components

        // elimination 4: discrete colors
        perform_del_discrete(g, colmap, &rec);

        // TODO: multiple calls for now obviously independent components -- could use "buffer consumer" to translate domains
        // TODO: just use consumer for all the back-translation (just dont "restore"!): meld translation layers here for this
        // TODO: for the consumer could use also use "canonization" strings for all kinds of reductions that represent parts of the graph as one node or edge
    }

    void restore(sgraph* g, automorphism_info* a, dejavu_consumer consume) {
        sparse_group group;
        group.domain_size = rec.old_d.size();
        group.initialize_from_automorphism_info(a, &rec);
    }
};

#endif //DEJAVU_PREP_H
