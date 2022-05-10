#ifndef DEJAVU_PREP_H
#define DEJAVU_PREP_H

#include "sgraph.h"
#include "refinement.h"
#include "schreier_shared.h"
#include <vector>

struct deg1_record {
    int parent;
    std::vector<int> children;
};

struct recovery_map {
    mark_set del;
    std::vector<int> translate;
    std::vector<int> backwards_translate;
    std::vector<deg1_record> deg1_children;
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
    }
};


class preprocessor {
private:
    void red_singleton_only_refinement() {

    }

    void res_singleton_only_refinement() {

    }

    void red_low_deg_assume_cref(sgraph* g, coloring<int>* c, recovery_map* rec) {
        work_list_t<int> worklist_deg0;
        work_list_t<int> worklist_deg1;
        work_list_t<int> worklist_deg2;

        worklist_deg0.initialize(g->v_size);
        worklist_deg1.initialize(g->v_size);
        worklist_deg2.initialize(g->v_size);

        rec->deg1_to_parent_ptn.reserve(g->v_size);
        rec->elim_order.reserve(g->v_size);

        for(int i = 0; i < g->v_size; ++i) {
            rec->deg1_to_parent_ptn.push_back(i);
        }

        for(int i = 0; i < g->v_size; ++i) {
            switch(g->d[i]) {
                case 0:
                    worklist_deg0.push_back(i);
                    break;
                case 1:
                    worklist_deg1.push_back(i);
                    break;
                case 2:
                    worklist_deg2.push_back(i);
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
                rec->deg1_to_parent_ptn[child] = parent;
                rec->elim_order.push_back(child);
                /*workspace1.push_back(child);
                if(!workset1.get(parent)) {
                    workset1.set(parent);
                    rec->deg1_children.emplace_back(deg1_record());
                    deg1_record *last_rec = &rec->deg1_children[rec->deg1_children.size()-1];
                    workspace2.arr[parent] = last_rec;

                    last_rec->children.reserve(child_col_sz);  // TODO: this should be one large array with pointers into it (like graph)
                    last_rec->children.push_back(child);
                } else {
                    deg1_record *last_rec = workspace2.arr[parent];
                    assert(last_rec != nullptr);
                    last_rec->children.push_back(child);
                }*/
                rec->del.set(child);

                // adjust parent degree
                g->d[parent] -= 1;
                if(g->d[parent] == 2) {
                    worklist_deg2.push_back(parent);
                } else if(g->d[parent] == 1) {
                    worklist_deg1.push_back(parent);
                } else if(g->d[parent] == 0) {
                    worklist_deg0.push_back(parent);
                }

                assert(g->d[parent] >= 0);
            }
            //std::cout << std::endl;
        }
        /*for(int i = 0; i < g->v_size; ++i) {
            std::cout << rec->del.get(i);
        }
        std::cout << std::endl;*/
    }

    void store_graph_aspects(sgraph* g, recovery_map* rec) {
        rec->old_d.clear();
        rec->old_d.reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            rec->old_d.push_back(g->d[i]);
        }
    }

    void perform_del(sgraph*g, int* colmap, recovery_map* rec) {
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
        rec->translate.reserve(g->v_size);
        rec->backwards_translate.reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            if (!rec->del.get(i)) {
                rec->translate.push_back(cnt);
                rec->backwards_translate.push_back(rec->translate.size() - 1);
                ++cnt;
                ++new_vsize;
            } else {
                rec->translate.push_back(-1);
            }
        }

        // make graph smaller using the translation array
        int epos  = 0;
        for (int i = 0; i < g->v_size; ++i) {
            const int old_v = i;
            const int new_v = rec->translate[i];

            if (new_v >= 0) {
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + rec->old_d[old_v]; ++j) {
                    const int ve = g_old_e[j];
                    const int new_ve = rec->translate[ve];
                    if(ve >= 0) {
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

        // copy old colmap for remaining vertices
        if(colmap != nullptr) {
            for (int i = 0; i < g->v_size; ++i) {
                const int old_v = i;
                const int new_v = rec->translate[i];
                if (new_v >= 0) {
                    colmap[new_v] = old_colmap[old_v];
                }
            }
        }

    }

    void res_low_deg() {
        // deg1:
        // restoration: each generator that permutes parents is extended to also permute X_1 -> X_2
        // restoration: add generators for symmetrical group for vertices of each parent
    }
public:
    recovery_map rec;
    void reduce(sgraph* g, int* colmap) {
        rec.del = mark_set();
        rec.del.initialize(g->v_size);
        store_graph_aspects(g, &rec);
        coloring<int> c;
        // singleton-only refinement, then cut graph
        refinement<int, int, int> R;
        g->initialize_coloring(&c, colmap);
        R.refine_coloring_first(g, &c, -1);
        red_low_deg_assume_cref(g, &c, &rec);

        perform_del(g, colmap, &rec);
    }

    void restore(sgraph* g, automorphism_info* a, dejavu_consumer) {
        sparse_group group;
        group.domain_size = rec.old_d.size();
        group.initialize_from_automorphism_info(a, &rec);
    }
};

#endif //DEJAVU_PREP_H
