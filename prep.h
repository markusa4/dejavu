#ifndef DEJAVU_PREP_H
#define DEJAVU_PREP_H

#include "sgraph.h"
#include "refinement.h"
#include "schreier_shared.h"
#include "selector.h"
#include <vector>

class tiny_orbit {
    mark_set          touched;
    work_list_t<int>  reset_arr;
    work_list_t<int>  map_arr;
public:
    int find_and_cut_orbit(const int v1) {
        int orbit1 = map_arr.arr[v1];
        while(orbit1 != map_arr.arr[orbit1])
            orbit1 = map_arr.arr[orbit1];
        map_arr.arr[v1] = orbit1;
        return orbit1;
    }

    void combine_orbits(const int v1, const int v2) {
        if(v1 != v2) {
            if(!touched.get(v1))
                reset_arr.push_back(v1);
            if(!touched.get(v2))
                reset_arr.push_back(v2);
            touched.set(v1);
            touched.set(v2);
            int orbit1 = find_and_cut_orbit(v1);
            int orbit2 = find_and_cut_orbit(v2);
            if(orbit1 < orbit2) {
                map_arr.arr[orbit2] = orbit1;
            } else {
                map_arr.arr[orbit1] = orbit2;
            }
        }
    }

    bool are_in_same_orbit(const int v1, const int v2) {
        if(v1 == v2)
            return true;
        int orbit1 = find_and_cut_orbit(v1);
        int orbit2 = find_and_cut_orbit(v2);
        return (orbit1 == orbit2);
    }

    void reset() {
        while(!reset_arr.empty()) {
            const int v = reset_arr.pop_back();
            map_arr.arr[v] = v;
        }
        touched.reset();
    }

    void initialize(int domain_size) {
        touched.initialize(domain_size);
        reset_arr.initialize(domain_size);
        map_arr.initialize(domain_size);
        for(int i = 0; i < domain_size; ++i) {
            map_arr.push_back(i);
        }
    }
};

class preprocessor {
    int domain_size;
    mark_set del;
    mark_set del_e;

    tiny_orbit orbit;

    std::vector<std::vector<int>> translation_layers;
    std::vector<std::vector<int>> backward_translation_layers;

    std::vector<int> backward_translation;
    std::vector<std::vector<int>> canonical_recovery_string;
public:
    double base = 1;
    int    exp  = 0;
    int    avg_support_sparse_ir        = 0;
    int    avg_reached_end_of_component = 0;
    double avg_end_of_comp = 0.0;

    work_list_t<std::vector<int>> add_edge_buff;
    work_list_t<int> worklist_deg0;
    work_list_t<int> worklist_deg1;
    mark_set add_edge_buff_act;
    refinement<int, int, int> R1;

    std::vector<int> g_old_v;
    std::vector<int> old_colmap;
    std::vector<int> g_old_d2;

    work_list_t<int> edge_scratch;

    std::vector<int> translate_layer_fwd;
    std::vector<int> translate_layer_bwd;

    std::vector<int> quotient_component_worklist_col;
    std::vector<int> quotient_component_worklist_col_sz;
    std::vector<int> quotient_component_worklist_v;
    std::vector<std::pair<int,int>> quotient_component_worklist_boundary;
    std::vector<std::pair<int,int>> quotient_component_worklist_boundary_swap;
    std::vector<int> quotient_component_workspace;

    std::vector<int> quotient_component_touched;
    std::vector<int> quotient_component_touched_swap;


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

    void red_deg2_path_size_1(sgraph* g, int* colmap, dejavu_consumer consume) {
        if(g->v_size <= 1)
            return;

        del.reset();

        int found_match = 0;

        coloring<int> test_col;
        g->initialize_coloring(&test_col, colmap);

        mark_set touched_endpoints;
        touched_endpoints.initialize(g->v_size);

        work_list_t<int> already_matched_to;
        work_list_t<int> already_matched_to_sz;
        already_matched_to.initialize(g->e_size);
        already_matched_to_sz.initialize(g->v_size);

        std::vector<int> underlying_match;
        underlying_match.reserve(g->v_size);
        for(int i = 0; i < g->v_size; ++i) {
            underlying_match.push_back(-1);
            already_matched_to_sz.push_back(0);
        }

        add_edge_buff_act.reset();
        add_edge_buff.reset();
        for(int i = 0; i < g->v_size; ++i)
            add_edge_buff.push_back(std::vector<int>()); // TODO: do this smarter... i know how many edges will end up here... should not allocate anything...

        for(int i = 0; i < g->v_size;) {
            const int test_v = test_col.lab[i];
            const int path_col_sz = test_col.ptn[i];
            if(g->d[test_v] == 2) {
                const int n1 = g->e[g->v[test_v] + 0];
                const int n2 = g->e[g->v[test_v] + 1];
                if(g->d[n1] != 2 && g->d[n2] != 2) {
                    // relevant path of length 1
                    const int col_n1 = test_col.vertex_to_col[n1];
                    const int col_n2 = test_col.vertex_to_col[n2];
                    const int col_sz_n1 = test_col.ptn[test_col.vertex_to_col[n1]];
                    const int col_sz_n2 = test_col.ptn[test_col.vertex_to_col[n2]];

                    if(col_sz_n1 != col_sz_n2 || col_sz_n1 != path_col_sz || col_n1 == col_n2) {
                        i += test_col.ptn[i] + 1;
                        continue;
                    }
                    //std::cout << "path col_sz " << test_col.ptn[i] + 1 << " neighbour sizes : " << col_sz_n1 + 1 << ", " << col_sz_n2 + 1 << std::endl;

                    bool already_matched_n1_n2 = false;

                    const int already_match_pt1 = g->v[test_col.lab[col_n1]];
                    const int already_match_pt2 = g->v[test_col.lab[col_n2]];

                    if(touched_endpoints.get(col_n1) && touched_endpoints.get(col_n2)) {
                        for(int j = 0; j < already_matched_to_sz.arr[col_n1]; ++j) {
                            if(already_matched_to.arr[already_match_pt1 + j] == col_n2)
                                already_matched_n1_n2 = true;
                        }
                    }

                    if(touched_endpoints.get(col_n1) && touched_endpoints.get(col_n2) && already_matched_n1_n2 && already_matched_to_sz.arr[col_n1] == 1 && already_matched_to_sz.arr[col_n2] == 1) {
                        //std::cout << "already touched endpoint 1 and 2 check inf matches perfectly" << std::endl;

                        const bool matching_same_cols1 = test_col.vertex_to_col[underlying_match[n1]] == test_col.vertex_to_col[n2];
                        const bool matching_same_cols2 = test_col.vertex_to_col[underlying_match[n2]] == test_col.vertex_to_col[n1];
                        const bool match_same_cols = matching_same_cols1 && matching_same_cols2;

                        bool check_if_match = true;
                        for(int f = 0; f < test_col.ptn[i] + 1; ++f) {
                            const int _test_v = test_col.lab[i + f];
                            const int _n1 = g->e[g->v[_test_v] + 0];
                            const int _n2 = g->e[g->v[_test_v] + 1];
                            check_if_match = check_if_match && (underlying_match[_n1] == _n2);
                            if(!check_if_match)
                                break;
                        }
                        //std::cout << check_if_match << std::endl;

                        if(check_if_match) {
                            found_match += test_col.ptn[i] + 1;
                            for(int f = 0; f < test_col.ptn[i] + 1; ++f) {
                                const int _test_v = test_col.lab[i + f];
                                del.set(_test_v);
                            }
                            for(int f = 0; f < test_col.ptn[i] + 1; ++f) {
                                const int _test_v = test_col.lab[i + f];
                                const int _n1 = g->e[g->v[_test_v] + 0];
                                const int _n2 = g->e[g->v[_test_v] + 1];
                                int can_n;
                                if(test_col.vertex_to_col[_n1] < test_col.vertex_to_col[_n2])
                                    can_n = _n1;
                                else
                                    can_n = _n2;
                                const int orig_test_v = translate_back(_test_v);
                                const int orig_n1 = translate_back(can_n);
                                for (int s = 0; s < canonical_recovery_string[orig_test_v].size(); ++s)
                                    canonical_recovery_string[orig_n1].push_back(
                                            canonical_recovery_string[orig_test_v][s]);
                            }
                            //std::cout << "matched previous matching, removing paths" << std::endl;
                        }


                        i += test_col.ptn[i] + 1;
                        continue;
                    }

                    if(touched_endpoints.get(col_n1) && touched_endpoints.get(col_n2) && already_matched_n1_n2) {
                        //std::cout << "already touched endpoint 1 and 2, already matched multiple" << std::endl;
                        i += test_col.ptn[i] + 1;
                        continue;
                    }
                    /*if(touched_endpoints.get(col_sz_n2)) {
                        std::cout << "already touched endpoint 2" << std::endl;
                        i += test_col.ptn[i] + 1;
                        continue;
                    }*/

                    const int col_endpoint2 = colmap[n2]; // colmap???
                    bool col_cycle = false;
                    for(int f = 0; f < g->d[n1]; ++f) {
                        const int col_other = colmap[g->e[g->v[n1] + f]];
                        if(col_other == col_endpoint2) {
                            col_cycle = true;
                            break;
                        }
                    }

                    if(col_cycle) {
                        PRINT("neighbours :-(");
                        i += test_col.ptn[i] + 1;
                        continue;
                    }

                    //std::cout << "matching " << col_n1 << " <-> " << col_n2 << std::endl;

                    already_matched_to.arr[already_match_pt1 + already_matched_to_sz.arr[col_n1]] = col_n2; // overwrites itself, need canonical vertex for color
                    already_matched_to.arr[already_match_pt2 + already_matched_to_sz.arr[col_n2]] = col_n1;
                    ++already_matched_to_sz.arr[col_n1];
                    ++already_matched_to_sz.arr[col_n2];

                    touched_endpoints.set(col_n1);
                    touched_endpoints.set(col_n2);

                    found_match += test_col.ptn[i] + 1;

                    for(int f = 0; f < test_col.ptn[i] + 1; ++f) {
                        const int _test_v = test_col.lab[i + f];
                        assert(g->d[_test_v] == 2);
                        const int _n1 = g->e[g->v[_test_v] + 0];
                        const int _n2 = g->e[g->v[_test_v] + 1];
                        underlying_match[_n1] = _n2;
                        underlying_match[_n2] = _n1;
                        for(int t = 0; t < add_edge_buff.arr[_n2].size(); ++t)
                            assert(add_edge_buff.arr[_n2][t] != _n1);
                        add_edge_buff.arr[_n2].push_back(_n1);
                        add_edge_buff_act.set(_n2);
                        add_edge_buff.arr[_n1].push_back(_n2);
                        add_edge_buff_act.set(_n1);
                        del.set(_test_v);

                        int can_n;
                        if(test_col.vertex_to_col[_n1] < test_col.vertex_to_col[_n2])
                            can_n = _n1;
                        else
                            can_n = _n2;
                        const int orig_test_v = translate_back(_test_v);
                        const int orig_n1 = translate_back(can_n);
                        for(int s = 0; s < canonical_recovery_string[orig_test_v].size(); ++s)
                            canonical_recovery_string[orig_n1].push_back(canonical_recovery_string[orig_test_v][s]);

                    }
                }
            }
            i += test_col.ptn[i] + 1;
        }
        PRINT("found matching vertices: " << found_match);
    }


    void red_deg2_assume_cref(sgraph* g, int* colmap, dejavu_consumer consume) {
        if(g->v_size <= 1)
            return;
        add_edge_buff_act.reset();
        del.reset();

        worklist_deg1.reset();

        add_edge_buff.initialize(domain_size);
        for(int i = 0; i < domain_size; ++i)
            add_edge_buff.push_back(std::vector<int>()); // TODO: do this smarter... i know how many edges will end up here... should not allocate anything...
        add_edge_buff_act.initialize(domain_size);

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
                        if(endpoint1 == endpoint2)
                            cycle = true;
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

            total_path_length += path_length;

            if(cycle) {
                continue;
            }

            // unique endpoint reduction TODO make it more general
            if(colmap[endpoint1] != colmap[endpoint2]) {
                //if(smallest_endpoint_cnt.arr[endpoint1] == 1 && smallest_endpoint_cnt.arr[endpoint2] == 1) {
                if((smallest_endpoint_cnt.arr[endpoint1] == 1)  || smallest_endpoint_cnt.arr[endpoint2] == 1) {
                    // TODO actually there is an underlying graph -- "DAC parts" of the graph can be used to attach paths to canonical endpoints?

                    /*if(check_if_neighbour(g, endpoint1, endpoint2)) {
                        //std::cout << "neighbours!" << std::endl;
                        continue;
                    }*/

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
                        //PRINT("col neighbours!");
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
                        //assert(!del2.get(del_v));
                        del.set(del_v);
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
                        canonical_recovery_string[unique_endpoint_orig].push_back(path_v_orig);
                        canonical_recovery_string[unique_endpoint_orig].insert(canonical_recovery_string[unique_endpoint_orig].end(),
                                                                                    canonical_recovery_string[path_v_orig].begin(), canonical_recovery_string[path_v_orig].end());
                    }
                    //if(path_length > 1)
                    //    std::cout<<std::endl;
                    path.reset();
                } else {
                    // TODO: endpoint that is only connected to uniquely colored paths is also unique endpoint
                    // TODO: can reduce and canonize these by sorting according to color e.g.
                }
            }
        }

        assert(num_paths1 == num_paths2);
        if(num_paths1 != 0) {
            PRINT("(prep-red) total paths: " << num_paths1 << " (avg length: " << total_path_length / num_paths1
                      << ")" << std::endl << "(prep-red) unique endpoints: " << unique_smallest_endpoint_paths);
        }
    }

    void write_canonical_recovery_string_to_automorphism(const int from, const int to) {
        assert(canonical_recovery_string[from].size() == canonical_recovery_string[to].size());
        for(int i = 0; i < canonical_recovery_string[to].size(); ++i) {
            //std::cout << i << std::endl;
            const int str_from = canonical_recovery_string[from][i];
            const int str_to = canonical_recovery_string[to][i];
            automorphism.arr[str_to]   = str_from;
            automorphism.arr[str_from] = str_to;
            automorphism_supp.push_back(str_from);
            automorphism_supp.push_back(str_to);
        }
        return;
    }

    void write_deg1_parent_to_map(sgraph* g, work_list_t<std::pair<int, int>>* stack1, work_list_t<int>* map, work_list_t<int>* childlist, work_list_t<int>* childcount, const int parent, const int child_from) {
        stack1->reset();
        map->reset();
        map->push_back(child_from);
        stack1->push_back(std::pair<int, int>(g->v[child_from], g->v[child_from] + childcount->arr[child_from]));
        while(!stack1->empty()) {
            std::pair<int, int> from_to = stack1->pop_back();
            int from       = from_to.first;
            const int to   = from_to.second;
            for(int f = from; f < to; ++f) {
                const int next = childlist->arr[f];
                const int from_next = g->v[next];
                const int to_next = g->v[next] + childcount->arr[next];
                map->push_back(next);
                assert(next != parent);
                if (from_next != to_next)
                    stack1->push_back(std::pair<int, int>(from_next, to_next));
            }
        }
    }

    void red_deg10_assume_cref(sgraph* g, const coloring<int>* c, dejavu_consumer consume) {
        //del.reset();

        worklist_deg0.reset();
        worklist_deg1.reset();

        mark_set is_parent;
        is_parent.initialize(g->v_size);

        g_old_d2.clear();
        g_old_d2.reserve(g->v_size);
        for(int i = 0; i < g->v_size; ++i) {
            g_old_d2.push_back(g->d[i]);
        }

        work_list_t<int> pair_match;
        pair_match.initialize(g->v_size);

        work_list_t<int> parentlist;
        parentlist.initialize(g->v_size);

        work_list_t<int> childcount;
        childcount.initialize(g->v_size);
        for(int i = 0; i < g->v_size; ++i)
            childcount.push_back(0);

        work_list_t<int> childcount_prev;
        childcount_prev.initialize(g->v_size);
        for(int i = 0; i < g->v_size; ++i)
            childcount_prev.push_back(0);

        work_list_t<std::pair<int, int>> stack1;
        work_list_t<int> map;
        stack1.initialize(g->v_size);
        map.initialize(g->v_size);


        for(int i = 0; i < c->ptn_sz;) {
            const int v = c->lab[i];
            switch(g_old_d2[v]) {
                case 0:
                    worklist_deg0.push_back(v);
                    break;
                case 1:
                    if(c->ptn[c->vertex_to_col[v]] > 0)
                        worklist_deg1.push_back(v);
                    // TODO: should early-out with less overhead!
                    break;
                default:
                    break;
            }
            i += c->ptn[i] + 1;
        }

        while(!worklist_deg1.empty()) {
            const int v_child = worklist_deg1.pop_back();
            if(del.get(v_child))
                continue;
            if(g_old_d2[v_child] != 1)
                continue;

            const int v_child_col = c->vertex_to_col[v_child];
            const int child_col_sz = c->ptn[v_child_col] + 1;
            bool is_pairs = false;
            bool permute_parents_instead = false;
            //std::cout << "col " << v_child_col << "(" << child_col_sz << ")" << " from vertex " << v_child << std::endl;
            if(child_col_sz == 1) { // TODO this can never be mapped, this is bloat... should handle before?
                del.set(v_child);
                continue;
            } else {
                parentlist.reset();
                is_parent.reset();
                for (int i = v_child_col; i < v_child_col + child_col_sz; ++i) {
                    int child = c->lab[i];

                    // search for parent
                    const int e_pos_child = g->v[child];
                    int parent = g->e[e_pos_child];

                    if(is_pairs && del.get(child))
                        continue;

                    int search_parent = 0;
                    while (del.get(parent)) {
                        ++search_parent;
                        parent = g->e[e_pos_child + search_parent];
                    }

                    assert(is_pairs?c->vertex_to_col[parent] == c->vertex_to_col[child]:true);
                    if(c->vertex_to_col[parent] == c->vertex_to_col[child]) {
                        is_pairs = true;
                        del.set(child);
                        del.set(parent);
                        pair_match.arr[child] = parent;
                        pair_match.arr[parent] = child;
                        if(parent < child) {
                            parentlist.push_back(parent);
                        } else {
                            parentlist.push_back(child);
                        }
                        continue;
                    }

                    del.set(child);

                    // save canonical info for parent
                    edge_scratch.arr[g->v[parent] + childcount.arr[parent]] = child;
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

            if(is_pairs) {
                for(int j = 0; j < parentlist.cur_pos; ++j) {
                    const int first_pair_parent = parentlist.arr[j];
                    const int pair_from = first_pair_parent; // need to use both childlist and canonical recovery again
                    const int pair_to = pair_match.arr[pair_from];

                    stack1.reset();
                    map.reset();
                    map.push_back(pair_from);
                    stack1.push_back(std::pair<int, int>(g->v[pair_from], g->v[pair_from] + childcount.arr[pair_from]));
                    while (!stack1.empty()) {
                        std::pair<int, int> from_to = stack1.pop_back();
                        int from = from_to.first;
                        const int to = from_to.second;
                        for (int f = from; f < to; ++f) {
                            const int next = edge_scratch.arr[f];
                            const int from_next = g->v[next];
                            const int to_next = g->v[next] + childcount.arr[next];
                            map.push_back(next);
                            assert(next != pair_to);
                            if (from_next != to_next)
                                stack1.push_back(std::pair<int, int>(from_next, to_next));
                        }
                    }

                    multiply_to_group_size(2);

                    assert(c->vertex_to_col[pair_from] == c->vertex_to_col[pair_to]);
                    assert(pair_from != pair_to);
                    int pos = 0;

                    automorphism_supp.reset();
                    // descending tree of child_to while writing automorphism
                    stack1.reset();
                    assert(map.arr[pos] != pair_to);
                    const int to_1 = translate_back(pair_to);
                    const int from_1 = translate_back(map.arr[pos]);
                    assert(automorphism.arr[to_1] == to_1);
                    assert(automorphism.arr[from_1] == from_1);
                    automorphism.arr[from_1] = to_1;
                    automorphism.arr[to_1] = from_1;
                    automorphism_supp.push_back(from_1);
                    automorphism_supp.push_back(to_1);
                    write_canonical_recovery_string_to_automorphism(to_1, from_1);
                    ++pos;
                    // child_to and child_from could have canonical strings when translated back
                    assert(childcount.arr[pair_to] == childcount.arr[pair_from]);
                    stack1.push_back(std::pair<int, int>(g->v[pair_to], g->v[pair_to] + childcount.arr[pair_to]));
                    while (!stack1.empty()) {
                        std::pair<int, int> from_to = stack1.pop_back();
                        int from = from_to.first;
                        const int to = from_to.second;
                        for (int f = from; f < to; ++f) {
                            const int next = edge_scratch.arr[f];
                            const int from_next = g->v[next];
                            const int to_next = g->v[next] + childcount.arr[next];
                            ++from;
                            assert(next >= 0);
                            assert(next < g->v_size);
                            assert(map.arr[pos] != next);

                            const int to_2 = translate_back(next);
                            const int from_2 = translate_back(map.arr[pos]);
                            assert(automorphism.arr[to_2] == to_2);
                            assert(automorphism.arr[from_2] == from_2);
                            automorphism.arr[from_2] = to_2;
                            automorphism.arr[to_2] = from_2;
                            automorphism_supp.push_back(from_2);
                            automorphism_supp.push_back(to_2);
                            write_canonical_recovery_string_to_automorphism(to_2, from_2);
                            ++pos;
                            if (from_next != to_next) // there was a semicolon here, should have been bug
                                stack1.push_back(std::pair<int, int>(from_next, to_next));
                        }
                    }

                    assert(pos == map.cur_pos);

                    assert(g_old_d2[pair_to] == 1);
                    assert(g_old_d2[pair_from] == 1);
                    assert(del.get(pair_to));
                    assert(del.get(pair_from));
                    consume(g->v_size, automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                    reset_automorphism(automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                    automorphism_supp.reset();
                }

                for(int j = 0; j < parentlist.cur_pos; ++j) {
                    const int first_pair_parent = parentlist.arr[j];
                    const int pair_from = first_pair_parent; // need to use both childlist and canonical recovery again
                    const int pair_to = pair_match.arr[pair_from];

                    const int original_parent = translate_back(first_pair_parent);
                    canonical_recovery_string[original_parent].push_back(translate_back(pair_to));
                }

                permute_parents_instead = true;
            }

            while(!parentlist.empty()) {
                int parent, childcount_from, childcount_to, child_from;
                if(!permute_parents_instead) {
                    parent = parentlist.pop_back();
                    childcount_from = childcount_prev.arr[parent];
                    childcount_to = childcount.arr[parent];
                } else {
                    parent = -1;
                    childcount_from = 0;
                    childcount_to = parentlist.cur_pos;
                }
                // automorphism 1: long cycle (c1 ... cn)
                assert(childcount_to - childcount_from > 0);
                if(childcount_to - childcount_from == 1) {
                    if(permute_parents_instead)
                        break;
                    continue;
                }
                if(!permute_parents_instead) {
                    child_from = edge_scratch.arr[g->v[parent] + childcount_from];
                } else {
                    child_from = parentlist.arr[0];
                }

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
                        const int next = edge_scratch.arr[f];
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
                    int child_to;
                    if(!permute_parents_instead) {
                        child_to = edge_scratch.arr[g->v[parent] + i];
                    } else {
                        child_to = parentlist.arr[i];
                    }
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
                            const int next      = edge_scratch.arr[f];
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
                    assert(del.get(child_to));
                    assert(del.get(child_from));
                    consume(g->v_size, automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                    reset_automorphism(automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
                    automorphism_supp.reset();
                }
                if(permute_parents_instead) {
                    break;      // :-)
                }
            }
            parentlist.reset();
        }

        if(g->v_size > 1) {
            // search for parents that still remain in the graph, and rewrite childlist structure into canonical string
            for (int _i = 0; _i < g->v_size; ++_i) {
                if (!del.get(_i) && childcount.arr[_i] > 0 && c->ptn[c->vertex_to_col[_i]] > 0) { // should it be childcount.arr[_i] > 0 or childcount.arr[_i] >= 0? was childcount.arr[_i] >= 0 but that doesnt make sense?
                    stack1.reset();
                    stack1.push_back(std::pair<int, int>(g->v[_i], g->v[_i] + childcount.arr[_i]));
                    const int orig_i = translate_back(_i);
                    canonical_recovery_string[orig_i].reserve(childcount.arr[_i]);
                    while (!stack1.empty()) {
                        std::pair<int, int> from_to = stack1.pop_back();
                        int from = from_to.first;
                        const int to = from_to.second;
                        if (from == to) {
                            continue;
                        } else {
                            const int next = edge_scratch.arr[from];
                            const int from_next = g->v[next];
                            const int to_next = g->v[next] + childcount.arr[next];
                            ++from;
                            const int orig_next = translate_back(next);
                            canonical_recovery_string[orig_i].push_back(orig_next);
                            // write canonical recovery string of orig_next into orig_i, since it is now represented by
                            // orig_i
                            assert(next != _i);
                            assert(orig_next != orig_i);
                            for(int j = 0; j < canonical_recovery_string[orig_next].size(); ++j)
                                canonical_recovery_string[orig_i].push_back(canonical_recovery_string[orig_next][j]);
                            stack1.push_back(std::pair<int, int>(from, to));
                            stack1.push_back(std::pair<int, int>(from_next, to_next));
                        }
                    }
                }
            }
        }


        while(!worklist_deg0.empty()) {
            const int v_child = worklist_deg0.pop_back();
            if (del.get(v_child))
                continue;
            assert(g_old_d2[v_child] == 0);

            is_parent.reset();
            const int v_child_col = c->vertex_to_col[v_child];
            const int child_col_sz = c->ptn[v_child_col] + 1;
            //std::cout << "col " << v_child_col << "(" << child_col_sz << ")" << " from vertex " << v_child << std::endl;
            const int parent_from = c->lab[v_child_col];
            del.set(parent_from);

            if(child_col_sz == 1)
                continue;
            int j = 2;
            for (int i = v_child_col + 1; i < v_child_col + child_col_sz; ++i) {
                multiply_to_group_size(j);
                ++j;
                const int parent_to = c->lab[i];
                assert(g_old_d2[parent_to] == 0);
                del.set(parent_to);

                const int orig_parent_from = translate_back(parent_from);
                const int orig_parent_to   = translate_back(parent_to);

                assert(canonical_recovery_string[orig_parent_to].size() == canonical_recovery_string[orig_parent_from].size());

                automorphism.arr[orig_parent_to]   = orig_parent_from;
                automorphism.arr[orig_parent_from] = orig_parent_to;
                automorphism_supp.push_back(orig_parent_from);
                automorphism_supp.push_back(orig_parent_to);
                for(int i = 0; i < canonical_recovery_string[orig_parent_to].size(); ++i) {
                    const int str_from = canonical_recovery_string[orig_parent_from][i];
                    const int str_to = canonical_recovery_string[orig_parent_to][i];
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

    void red_quotient_components(sgraph* g, int* colmap, dejavu_consumer consume) {
        if(g->v_size <= 1)
            return;

        //std::cout << "(prep-red) deg1" << std::endl;
        del_e.reset();

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
                        del_e.set(f); // mark edge for deletion (reverse edge is handled later automatically)
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

    void red_quotient_matchings_(sgraph* g, int* colmap, dejavu_consumer consume) {
        //std::cout << "(prep-red) deg1" << std::endl;
        del.reset();

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
                        //std::cout << (test_col.ptn[test_col.vertex_to_col[v_neigh]] + 1) << " <-> " << (test_col.ptn[i] + 1) << std::endl;
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

        //std::cout << "independent_con: " << edge_cnt << "/" << g->e_size << std::endl;
        //std::cout << "has_matching_indep_color: " << cnt << "/" << g->v_size << std::endl;
        //std::cout << "has_matching_color: " << v_has_matching_color << "/" << g->v_size << std::endl;
    }
    /*void store_graph_aspects(sgraph* g, recovery_map* rec) {
        old_d.clear();
        old_d.reserve(g->v_size);
        for (int i = 0; i < g->v_size; ++i) {
            old_d.push_back(g->d[i]);
        }
    }*/

    void copy_coloring_to_colmap(const coloring<int>* c, int* colmap) {
        for(int i = 0; i < c->lab_sz; ++i) {
            colmap[i] = c->vertex_to_col[i];
        }
    }

    void perform_del(sgraph*g, int* colmap) {
        // TODO: "sparse" deletion?

        // copy some stuff
        g_old_v.clear();
        //g_old_e.clear();
        edge_scratch.reset();
        old_colmap.clear();
        g_old_d2.clear();
        translate_layer_fwd.clear();
        translate_layer_bwd.clear();

        for(int i = 0; i < backward_translation_layers[backward_translation_layers.size() - 1].size(); ++i)
            translate_layer_bwd.push_back(backward_translation_layers[backward_translation_layers.size() - 1][i]);

        // create translation array from old graph to new graph vertices
        int cnt = 0;
        int new_vsize = 0;
        //translation_layers[fwd_ind].reserve(g->v_size);
        //backward_translation_layers[back_ind].reserve(g->v_size);
        int pos_now = 0;
        for (int i = 0; i < g->v_size; ++i) {
            if (!del.get(i)) {
                //translation_layers[fwd_ind].push_back(cnt);
                translate_layer_fwd.push_back(cnt);
                //backward_translation_layers[back_ind].push_back(translation_layers[fwd_ind].size() - 1); // TODO at most need 2 arrays
                //translate_layer_bwd.push_back(translate_layer_fwd.size() - 1);
                //translate_layer_fwd.size() - 1
                const int translate_back = translate_layer_bwd[translate_layer_fwd.size() - 1];
                backward_translation_layers[backward_translation_layers.size() - 1][cnt] = translate_back;
                ++cnt;
                ++new_vsize;
            } else {
                //translation_layers[fwd_ind].push_back(-1);
                translate_layer_fwd.push_back(-1);
            }
        }

        if(new_vsize == 0 || new_vsize == 1) {
            g->v_size = 0;
            g->e_size = 0;
            g->d_size = 0;
            return;
        }

        backward_translation_layers[backward_translation_layers.size() - 1].resize(cnt);

        g_old_v.reserve(g->v_size);
        //g_old_e.reserve(g->e_size);
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
            edge_scratch.push_back(g->e[i]);
        }

        // make graph smaller using the translation array
        int epos  = 0;
        for (int i = 0; i < g->v_size; ++i) { // TODO: maybe iterate over new_vsize and use backward array?
            const int old_v = i;
            //const int new_v = translation_layers[fwd_ind][i];
            const int new_v = translate_layer_fwd[i];

            if (new_v >= 0) {
                colmap[new_v] = old_colmap[old_v];
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + g_old_d2[old_v]; ++j) {
                    const int ve = edge_scratch.arr[j];
                    //const int new_ve = translation_layers[fwd_ind][ve];
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
        del.reset();
    }

    void perform_del_edge(sgraph*g, int* colmap) {
        if(g->v_size <= 1)
            return;

        int pre_esize = g->e_size;
        // copy some stuff
        g_old_v.clear();
        //g_old_e.clear();
        edge_scratch.reset();
        old_colmap.clear();
        g_old_d2.clear();

        g_old_v.reserve(g->v_size);
        //g_old_e.reserve(g->e_size);
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
            edge_scratch.push_back(g->e[i]);
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
                    const int ve = edge_scratch.arr[j];
                    const int new_ve = ve;
                    if(!del_e.get(j)) {
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

        //std::cout << "edges: " << pre_esize << " -> " << g->e_size << std::endl;
    }

    void perform_del2(sgraph*g, int* colmap) {
        if(g->v_size <= 1)
            return;

        // copy some stuff
        g_old_v.clear();
        old_colmap.clear();
        g_old_d2.clear();
        translate_layer_fwd.clear();
        translate_layer_bwd.clear();

        edge_scratch.reset();

        for(int i = 0; i < backward_translation_layers[backward_translation_layers.size() - 1].size(); ++i)
            translate_layer_bwd.push_back(backward_translation_layers[backward_translation_layers.size() - 1][i]);

        // create translation array from old graph to new graph vertices
        int cnt = 0;
        int new_vsize = 0;
        for (int i = 0; i < g->v_size; ++i) {
            if (!del.get(i)) {
                translate_layer_fwd.push_back(cnt);
                const int translate_back = translate_layer_bwd[translate_layer_fwd.size() - 1];
                backward_translation_layers[backward_translation_layers.size() - 1][cnt] = translate_back;
                ++cnt;
                ++new_vsize;
            } else {
                //translation_layers[fwd_ind].push_back(-1);
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
        //g_old_e.reserve(g->e_size);
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
            edge_scratch.push_back(g->e[i]);
        }

        backward_translation_layers[backward_translation_layers.size() - 1].resize(cnt);

        // make graph smaller using the translation array
        int epos  = 0;
        for (int i = 0; i < g->v_size; ++i) {
            const int old_v = i;
            //const int new_v = translation_layers[fwd_ind][old_v];
            const int new_v = translate_layer_fwd[old_v];

            if (new_v >= 0) {
                int new_d = 0;
                g->v[new_v] = epos;
                for(int j = g_old_v[old_v]; j < g_old_v[old_v] + g_old_d2[old_v]; ++j) {
                    const int ve = edge_scratch.arr[j];
                    //const int new_ve = translation_layers[fwd_ind][ve];
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
                        //const int edge_to_new = translation_layers[fwd_ind][edge_to_old];
                        const int edge_to_new = translate_layer_fwd[edge_to_old];
                        assert(edge_to_old >= 0);
                        assert(edge_to_old <= domain_size);
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
                //const int new_v = translation_layers[fwd_ind][i];
                const int new_v = translate_layer_fwd[i];
                assert(new_v < new_vsize);
                if (new_v >= 0) {
                    assert(old_colmap[old_v] >= 0);
                    assert(old_colmap[old_v] < domain_size);
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
        del.reset();
    }

    void mark_discrete_for_deletion(sgraph*g, int* colmap) {
        int discrete_cnt = 0;
        worklist_deg0.reset();
        for(int i = 0; i < domain_size; ++i) {
            worklist_deg0.push_back(0);
        }
        for(int i = 0; i < g->v_size; ++i) {
            assert(colmap[i] < domain_size);
            worklist_deg0.arr[colmap[i]]++;
        }
        for(int i = 0; i < g->v_size; ++i) {
            if(worklist_deg0.arr[colmap[i]] == 1)
                del.set(i);
        }
        worklist_deg0.reset();
    }

    void perform_del_discrete(sgraph*g, int* colmap) {
        if(g->v_size <= 1)
            return;

        int discrete_cnt = 0;
        worklist_deg0.reset();
        for(int i = 0; i < domain_size; ++i) {
            worklist_deg0.push_back(0);
        }
        for(int i = 0; i < g->v_size; ++i) {
            assert(colmap[i] < domain_size);
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

        // copy some stuff
        g_old_v.clear();
        old_colmap.clear();
        g_old_d2.clear();
        translate_layer_fwd.clear();
        translate_layer_bwd.clear();

        for(int i = 0; i < backward_translation_layers[backward_translation_layers.size() - 1].size(); ++i)
            translate_layer_bwd.push_back(backward_translation_layers[backward_translation_layers.size() - 1][i]);

        g_old_v.reserve(g->v_size);
        g_old_d2.reserve(g->v_size);

        edge_scratch.reset();

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
            edge_scratch.push_back(g->e[i]);
        }

        // create translation array from old graph to new graph vertices
        int cnt = 0;
        int new_vsize = 0;
        for (int i = 0; i < g->v_size; ++i) {
            if (worklist_deg0.arr[colmap[i]] != 1) {
                translate_layer_fwd.push_back(cnt);
                const int translate_back = translate_layer_bwd[translate_layer_fwd.size() - 1];
                backward_translation_layers[backward_translation_layers.size() - 1][cnt] = translate_back;
                ++cnt;
                ++new_vsize;
            } else {
                translate_layer_fwd.push_back(-1);
            }
        }

        backward_translation_layers[backward_translation_layers.size() - 1].resize(cnt);

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
                    const int ve = edge_scratch.arr[j];
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

        g->e_size = epos;
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
        del.reset();
    }
public:
    coloring<int> c;
    work_list_t<int> automorphism;
    work_list_t<int> automorphism_supp;
    bool layers_melded = false;

    void meld_translation_layers() {
        if(layers_melded)
            return;
        const int layers = translation_layers.size();
        backward_translation.reserve(layers - 1);
        const int reduced_size = backward_translation_layers[layers - 1].size();
        for(int i = 0; i < reduced_size; ++i)
            backward_translation.push_back(i);
        for(int i = 0; i < reduced_size; ++i) {
            int next_v = i;
            for(int l = layers - 1; l >= 0; --l) {
                next_v = backward_translation_layers[l][next_v];
            }
            backward_translation[i] = next_v;
        }
        layers_melded = true;
    }

    int translate_back(int v) {
        const int layers = translation_layers.size();
        for(int l = layers - 1; l >= 0; --l) {
            v = backward_translation_layers[l][v];
        }
        return v;
    }

    void multiply_to_group_size(double n) {
        base *= n;
        while(base > 10) {
            base = base / 10;
            exp += 1;
        }
    }

    int select_color_component(sgraph*g, coloring<int>* c1, int component_start_pos) {
        //const int my_component = v_to_component.arr[c1->lab[quotient_component_worklist_col[component_start_pos]]];
        int cell = -1;
        bool only_discrete_prev = true;
        for (int _i = component_start_pos; _i < quotient_component_worklist_col.size(); ++_i) {
            const int col    = quotient_component_worklist_col[_i];
            const int col_sz = quotient_component_worklist_col_sz[_i];

            if (col == -1) {// reached end of component
                break;
            }

            for(int i = col; i < col + col_sz + 1; ++i) {
                if(c1->vertex_to_col[c1->lab[i]] != i)
                    continue; // not a color
                if (c1->ptn[i] > 0 && only_discrete_prev) {
                    //start_search_here = i;
                    only_discrete_prev = false;
                }
                if (c1->ptn[i] >= 1) {
                    cell = i;
                    break;
                }
            }
            if(cell != -1)
                break;
        }
        return cell;
    }

    int select_color_component(sgraph*g, coloring<int>* c1, int component_start_pos, int component_end_pos) {
        int cell = -1;
        bool only_discrete_prev = true;
        for (int _i = component_start_pos; _i < component_end_pos; ++_i) {
            const int col    = quotient_component_worklist_col[_i];
            const int col_sz = quotient_component_worklist_col_sz[_i];

            if (col == -1) {// reached end of component
                break;
            }

            for(int i = col; i < col + col_sz + 1; ++i) {
                if(c1->vertex_to_col[c1->lab[i]] != i)
                    continue; // not a color
                if (c1->ptn[i] > 0 && only_discrete_prev) {
                    //start_search_here = i;
                    only_discrete_prev = false;
                }
                if (c1->ptn[i] >= 1) {
                    cell = i;
                    break;
                }
            }
            if(cell != -1)
                break;
        }
        return cell;
    }

    work_list_t<int> _automorphism;
    work_list_t<int> _automorphism_supp;
    std::vector<int> save_colmap;
    mark_set touched_color;
    work_list_t<int> touched_color_list;
    bool ir_quotient_component_init = false;

    int sparse_ir_col_sz_2_quotient_components(sgraph*g, int* colmap, dejavu_consumer consume) {
        //::cout << "component algo" << std::endl;

        int automorphisms_found = 0;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() * ((1 * 5) * 5135235);
        int selector_seed = seed;

        //work_list_t<int> _automorphism;
        //work_list_t<int> _automorphism_supp;

        save_colmap.reserve(g->v_size);
        save_colmap.clear();

        if(!ir_quotient_component_init) {
            _automorphism.initialize(g->v_size);
            _automorphism_supp.initialize(g->v_size);
            for (int i = 0; i < g->v_size; ++i)
                _automorphism.arr[i] = i;
            touched_color.initialize(g->v_size);
            touched_color_list.initialize(g->v_size);
            orbit.initialize(g->v_size);
            ir_quotient_component_init = true;
        }

        //mark_set touched_color;
        //work_list_t<int> touched_color_list;
        touched_color.reset();
        touched_color_list.reset();

        invariant I1, I2;
        I1.only_acc = true;
        I2.only_acc = true;

        coloring<int> c1;
        g->initialize_coloring(&c1, colmap);
        for(int i = 0; i < g->v_size; ++i)
            colmap[i] = c1.vertex_to_col[i];

        coloring<int> c2;
        c2.copy(&c1);
        c2.copy_ptn(&c1);

        //refinement<int, int, int> R1;
        bool certify = true;
        int quotient_component_start_pos   = 0;
        int quotient_component_start_pos_v = 0;

        while(quotient_component_start_pos < quotient_component_worklist_col.size()) {
            assert(quotient_component_start_pos_v < quotient_component_worklist_v.size());
            // additional loop here for components?
            //std::cout << "component " << quotient_component_start_pos << ", " << quotient_component_start_pos_v << std::endl;
            int start_search_here = quotient_component_start_pos;
            assert(quotient_component_worklist_col[quotient_component_start_pos] >= 0);
            assert(quotient_component_worklist_v[quotient_component_start_pos_v] >= 0);
            certify = true;

            while(certify) {
                // select a color class of size 2
                bool only_discrete_prev = true;
                int cell = -1;
                for (int _i = start_search_here; _i < quotient_component_worklist_col.size(); ++_i) { // TODO
                    const int v = quotient_component_worklist_col[_i];
                    if (v == -1) {// reached end of component
                        break;
                    }
                    const int i = v;
                    if (c1.ptn[i] > 0 && only_discrete_prev) {
                        start_search_here = i;
                        only_discrete_prev = false;
                    }
                    if (c1.ptn[i] == 1) {
                        cell = i;
                        break;
                    }
                }

                if (cell == -1)
                    break;

                touched_color.reset();
                touched_color_list.reset();

                const int ind_v1 = c1.lab[cell];
                const int ind_v2 = c1.lab[cell + 1];
                const int init_c1 = R1.individualize_vertex(&c1, ind_v1);

                touched_color.set(cell);
                touched_color.set(init_c1);
                touched_color_list.push_back(cell);
                touched_color_list.push_back(init_c1);

                R1.refine_coloring(g, &c1, &I1, init_c1, nullptr, -1, -1, nullptr, &touched_color, &touched_color_list);

                const int init_c2 = R1.individualize_vertex(&c2, ind_v2);
                R1.refine_coloring(g, &c2, &I2, init_c2, nullptr, -1, -1, nullptr, &touched_color, &touched_color_list);

                if (I1.acc != I2.acc) {
                    for (int i = 0; i < g->v_size; ++i) {
                        colmap[i] = c1.vertex_to_col[i];
                    }
                    I2.acc = I1.acc;
                    c2.copy(&c1);
                    continue;
                }


                _automorphism_supp.reset();
                if (c1.cells != g->v_size) { // touched_colors doesn't work properly when early-out is used
                    // read automorphism
                    for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                        const int _c = touched_color_list.arr[j];
                        int f = 0;
                        while (f < c1.ptn[_c] + 1) {
                            const int i = _c + f;
                            ++f;
                            if (c1.lab[i] != c2.lab[i]) {
                                _automorphism.arr[c1.lab[i]] = c2.lab[i];
                                _automorphism_supp.push_back(c1.lab[i]);
                            }
                        }
                    }
                } else {
                    for (int i = 0; i < g->v_size; ++i) {
                        if (c1.lab[i] != c2.lab[i]) {
                            _automorphism.arr[c1.lab[i]] = c2.lab[i];
                            _automorphism_supp.push_back(c1.lab[i]);
                        }
                    }
                }

                certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr, _automorphism_supp.cur_pos,
                                                         _automorphism_supp.arr);
                assert(certify ? R1.certify_automorphism(g, _automorphism.arr) : true);
                //std::cout << "(pre-red) sparse-ir certify: " << certify << std::endl;
                if (certify) {
                    ++automorphisms_found;
                    pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr,
                                         consume);
                    multiply_to_group_size(2);

                    // reset c2 to c1
                    if (c1.cells != g->v_size) {
                        for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                            const int _c = touched_color_list.arr[j];
                            int f = 0;
                            c2.cells = c1.cells;
                            while (f < c1.ptn[_c] + 1) {
                                const int i = _c + f;
                                ++f;
                                c2.ptn[i] = c1.ptn[i];
                                c2.lab[i] = c1.lab[i];
                                c2.vertex_to_col[c2.lab[i]] = c1.vertex_to_col[c1.lab[i]];
                                c2.vertex_to_lab[c2.lab[i]] = c1.vertex_to_lab[c1.lab[i]];
                            }
                        }
                    } else {
                        // discrete
                        //std::cout << "this here" << std::endl;
                    }
                } else {
                    touched_color.reset();
                    touched_color_list.reset();

                    save_colmap.clear();
                    for (int _i = quotient_component_start_pos_v; _i < quotient_component_worklist_v.size() &&
                                                                  quotient_component_worklist_v[_i] != -1; ++_i) {
                        const int i = quotient_component_worklist_v[_i];
                        save_colmap.push_back(c1.vertex_to_col[i]);
                    }

                    //selector<int, int, int> S;
                    //strategy<int> strat;
                    //strat.cell_selector_type = SELECTOR_FIRST;
                    //S.empty_cache();
                    int col = -1;
                    //std::cout << "first path" << std::endl;

                    while (true) {
                        // stay within component!
                        col = select_color_component(g, &c1, quotient_component_start_pos);
                        if (col == -1) break;
                        const int rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c1.ptn[col] + 1));
                        const int v = c1.lab[rpos];
                        const int init_color_class = R1.individualize_vertex(&c1, v);
                        R1.refine_coloring(g, &c1, &I1, init_color_class, nullptr, -1, -1,
                                           nullptr, nullptr, nullptr);
                    }

                    while (true) {
                        col = select_color_component(g, &c2, quotient_component_start_pos);
                        if (col == -1) break;
                        const int rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c2.ptn[col] + 1));
                        const int v = c2.lab[rpos];
                        const int init_color_class = R1.individualize_vertex(&c2, v);
                        R1.refine_coloring(g, &c2, &I2, init_color_class, nullptr, -1, -1,
                                           nullptr, nullptr, nullptr);
                    }

                    //std::cout << "(prep-red) sparse-ir random paths: " << I1.acc << ", " << I2.acc << std::endl;
                    reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                    _automorphism_supp.reset();


                    for(int _i = quotient_component_start_pos_v; _i < quotient_component_worklist_v.size() &&
                                                                   quotient_component_worklist_v[_i] != -1; ++_i) { // TODO: only consider component! can use touched as long as not final component
                        const int i = quotient_component_worklist_v[_i];
                        const int lab_p = c1.vertex_to_lab[i];
                        assert(i >= 0);
                        assert(i < g->v_size);
                        if (c1.lab[lab_p] != c2.lab[lab_p]) {
                            _automorphism.arr[c1.lab[lab_p]] = c2.lab[lab_p];
                            _automorphism_supp.push_back(c1.lab[lab_p]);
                        }
                    }

                    certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr, _automorphism_supp.cur_pos,
                                                             _automorphism_supp.arr);
                    //std::cout << "(prep-red) sparse-ir random path certify: " << certify << std::endl;
                    if (certify) {
                        ++automorphisms_found;
                        pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos,
                                             _automorphism_supp.arr,
                                             consume);
                        multiply_to_group_size(2);
                    }
                    reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                    _automorphism_supp.reset();

                    if (certify) {
                        /*for (int i = 0; i < g->v_size; ++i) {
                            colmap[i] = save_colmap[i];
                        }*/
                        int f = 0;
                        for (int _i = quotient_component_start_pos_v; _i < quotient_component_worklist_v.size() &&
                                                                      quotient_component_worklist_v[_i] != -1; ++_i) {
                            const int i = quotient_component_worklist_v[_i];
                            colmap[i] = save_colmap[f];
                            ++f;
                        }
                    }

                    break;
                }
                reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                automorphism_supp.reset();

                if (certify) {
                    if (c1.cells == g->v_size) {
                        for (int i = 0; i < g->v_size; ++i) {
                            colmap[i] = c1.vertex_to_col[i];
                        }
                    } else {
                        for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                            const int _c = touched_color_list.arr[j];
                            int f = 0;
                            while (f < c1.ptn[_c] + 1) {
                                assert(c1.vertex_to_col[c1.lab[_c + f]] == _c);
                                colmap[c1.lab[_c + f]] = c1.vertex_to_col[c1.lab[_c + f]];
                                ++f;
                            }
                        }
                    }
                }
            }
            while(quotient_component_start_pos < quotient_component_worklist_col.size() &&
                  quotient_component_worklist_col[quotient_component_start_pos] != -1)
                ++quotient_component_start_pos;
            ++quotient_component_start_pos;

            while(quotient_component_start_pos_v < quotient_component_worklist_v.size() &&
                  quotient_component_worklist_v[quotient_component_start_pos_v] != -1)
                ++quotient_component_start_pos_v;
            ++quotient_component_start_pos_v;
        }
        return automorphisms_found;
    }


    int sparse_ir_col_sz_2(sgraph*g, int* colmap, dejavu_consumer consume) {
        // TODO: ORBITS!
        int found_automorphisms = 0;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() * ((1 * 5) * 5135235);
        int selector_seed = seed;

        coloring<int> c1;
        g->initialize_coloring(&c1, colmap); // could re-order to reduce overhead when not applicable
        for(int i = 0; i < g->v_size; ++i)
            colmap[i] = c1.vertex_to_col[i];

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

        //refinement<int, int, int> R1;
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
                return found_automorphisms;

            touched_color.reset();
            touched_color_list.reset();

            const int ind_v1 = c1.lab[cell];
            const int ind_v2 = c1.lab[cell + 1];
            const int init_c1 = R1.individualize_vertex(&c1, ind_v1);

            touched_color.set(cell);
            touched_color.set(init_c1);
            touched_color_list.push_back(cell);
            touched_color_list.push_back(init_c1);

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
                    const int _c = touched_color_list.arr[j];
                    int f = 0;
                    while (f < c1.ptn[_c] + 1) {
                        const int i = _c + f;
                        ++f;
                        if (c1.lab[i] != c2.lab[i]) {
                            // TODO: what if color is larger than 1?
                            //if(c1.ptn[i] == 1) {
                            //     std::cout << "(" << c1.lab[i] << ", " << c1.lab[i + 1] << ")" << "(" << c2.lab[i] << ", " << c2.lab[i + 1] << ")" << std::endl;
                            // }
                            //if (colmap[c1.lab[i]] != colmap[c2.lab[i]])
                            //    std::cout << "color mismatch" << std::endl;
                            _automorphism.arr[c1.lab[i]] = c2.lab[i];
                            _automorphism_supp.push_back(c1.lab[i]);
                        }
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
                ++found_automorphisms;
                pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr,
                                     consume);
                multiply_to_group_size(2);

                // reset c2 to c1
                if(c1.cells != g->v_size) {
                    for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                        const int _c = touched_color_list.arr[j];
                        int f = 0;
                        while (f < c1.ptn[_c] + 1) {
                            const int i = _c + f;
                            ++f;
                            if (c1.lab[i] != c2.lab[i]) {
                                c2.lab[i] = c1.lab[i];
                                c2.vertex_to_col[c2.lab[i]] = c1.vertex_to_col[c1.lab[i]];
                                c2.vertex_to_lab[c2.lab[i]] = c1.vertex_to_lab[c1.lab[i]];
                            }
                        }
                    }
                } else {
                    //std::cout << "this here" << std::endl;
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

                //std::cout << "(prep-red) sparse-ir random paths: "<< I1.acc << ", " << I2.acc << std::endl;
                reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                _automorphism_supp.reset();

                for (int i = 0; i < g->v_size; ++i) {
                    if (c1.lab[i] != c2.lab[i]) {
                        _automorphism.arr[c1.lab[i]] = c2.lab[i];
                        _automorphism_supp.push_back(c1.lab[i]);
                    }
                }

                certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                //std::cout << "(prep-red) sparse-ir random path certify: " << certify << std::endl;
                if(certify) {
                    ++found_automorphisms;
                    pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr,
                                         consume);
                    multiply_to_group_size(2);
                }
                reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                _automorphism_supp.reset();

                if(certify) {
                    for (int i = 0; i < g->v_size; ++i) {
                        colmap[i] = save_colmap[i];
                    }
                }

                return found_automorphisms;
                // TODO: could compute some type of invariant
            }
            reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
            automorphism_supp.reset();

            if(certify) {
                if(c1.cells == g->v_size) {
                    for (int i = 0; i < g->v_size; ++i) {
                        colmap[i] = c1.vertex_to_col[i];
                    }
                } else {
                    for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                        const int _c = touched_color_list.arr[j];
                        int f = 0;
                        while (f < c1.ptn[_c] + 1) {
                            assert(c1.vertex_to_col[c1.lab[_c + f]] == _c);
                            colmap[c1.lab[_c + f]] = c1.vertex_to_col[c1.lab[_c + f]];
                            ++f;
                        }
                    }
                }
            }
        }
        return found_automorphisms;
    }

    void sparse_ir(sgraph*g, int* colmap, dejavu_consumer consume, selector_type sel_type) {
        if(g->v_size <= 1)
            return;

        // TODO: ORBITS!
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() * ((1 * 5) * 5135235);
        int selector_seed = seed;

        std::vector<int> individualize_later;

        coloring<int> c1;
        g->initialize_coloring(&c1, colmap); // could re-order to reduce overhead when not applicable
        for(int i = 0; i < g->v_size; ++i)
            colmap[i] = c1.vertex_to_col[i];

        coloring<int> c2;
        c2.copy(&c1);
        c2.copy_ptn(&c1);

        for(int i = 0; i < g->v_size; ++i) {
            assert(c2.vertex_to_col[i] == c1.vertex_to_col[i]);
            assert(c2.vertex_to_lab[i] == c1.vertex_to_lab[i]);
            assert(c2.lab[i] == c1.lab[i]);
            assert(c2.ptn[i] == c1.ptn[i]);
        }

        coloring<int> original_c;
        original_c.copy(&c1);
        original_c.copy_ptn(&c1);

        coloring<int> color_cache;
        color_cache.copy(&c1);

        selector<int, int, int> S;
        strategy<int> m;
        m.cell_selector_type = sel_type;

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
        I1.acc = 0;
        I2.acc = 0;

        //refinement<int, int, int> R1;
        bool certify = true;
        int start_search_here = 0;

        S.empty_cache();

        while(true) {
            // select a color class
            const int cell = S.select_color_dynamic(g, &c1, &m);
            /*bool only_discrete_prev = true;
            int cell = -1;
            for (int i = start_search_here; i < c1.ptn_sz;) {
                if (c1.ptn[i] > 0 && only_discrete_prev) {
                    start_search_here = i;
                    only_discrete_prev = false;
                }
                if (c1.ptn[i] >= 1) {
                    cell = i;
                    break;
                }
                i += c1.ptn[i] + 1;
            }*/
            if (cell == -1)
                break;

            const int cell_sz = c1.ptn[cell];

            bool is_orig_color    = original_c.vertex_to_col[original_c.lab[cell]] == cell;
            bool is_orig_color_sz = original_c.ptn[cell] == cell_sz;

            //std::cout << cell << " sz " << c1.ptn[cell] + 1 << " is_orig_col " << is_orig_color << " orig_sz " << original_c.ptn[cell] + 1 << " is_orig_sz " << is_orig_color_sz << std::endl;

            touched_color.reset();
            touched_color_list.reset();

            for(int i = 0; i < g->v_size; ++i) {
                assert(c2.vertex_to_col[i] == c1.vertex_to_col[i]);
                assert(c2.vertex_to_lab[i] == c1.vertex_to_lab[i]);
                assert(c2.lab[i] == c1.lab[i]);
                assert(c2.ptn[i] == c1.ptn[i]);
            }
            assert(c1.cells == c2.cells);

            const int ind_v1 = c1.lab[cell];
            assert(c1.vertex_to_col[ind_v1] == cell);
            assert(c1.ptn[cell] > 0);
            const long acc_prev = I1.acc;
            I2.acc = acc_prev;
            I1.acc = acc_prev;
            const int init_c1 = R1.individualize_vertex(&c1, ind_v1);

            bool all_certified = true;

            //std::cout << cell << " ind " << init_c1 << std::endl;
            touched_color.set(cell);
            touched_color.set(init_c1);
            touched_color_list.push_back(cell);
            touched_color_list.push_back(init_c1);

            R1.refine_coloring(g, &c1, &I1, init_c1, nullptr, -1, -1, nullptr, &touched_color, &touched_color_list);

            // color cache to reset c2 back "before" individualization
            for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                const int _c = touched_color_list.arr[j];
                int f = 0;
                while (f < c1.ptn[_c] + 1) {
                    const int i = _c + f;
                    ++f;
                    color_cache.ptn[i] = c2.ptn[i];
                    color_cache.lab[i] = c2.lab[i];
                    color_cache.vertex_to_col[color_cache.lab[i]] = c2.vertex_to_col[c2.lab[i]];
                    color_cache.vertex_to_lab[color_cache.lab[i]] = c2.vertex_to_lab[c2.lab[i]];
                }
            }
            color_cache.cells = c2.cells;

            if(!is_orig_color || !is_orig_color_sz) {// no point in looking for automorphisms
                all_certified = false;
            }

            bool hard_reset = false;

            // TODO: REPEAT FOR ALL v2
            for(int it_c = 1; (it_c < cell_sz + 1) && all_certified; ++it_c) { // TODO should be it_c = 1
                const int ind_v2 = c2.lab[cell + it_c];
                //assert(ind_v1 == ind_v2);
                assert(c2.lab[cell + it_c] == color_cache.lab[cell + it_c]);
                //std::cout << "ind " << ind_v2 << std::endl;
                assert(c2.vertex_to_col[ind_v2] == cell);
                assert(c2.ptn[cell] > 0);
                assert(c2.ptn[cell] == cell_sz);
                I2.acc = acc_prev;
                const int init_c2 = R1.individualize_vertex(&c2, ind_v2);
                assert(init_c1 == init_c2);
                R1.refine_coloring(g, &c2, &I2, init_c2, nullptr, -1, -1, nullptr, &touched_color, &touched_color_list);

                //std::cout << I1.acc << ", "  << I2.acc << std::endl;
                if (I1.acc != I2.acc) {
                    //std::cout << I1.acc << ", "  << I2.acc << std::endl;
                    //std::cout << "acc mismatch " << cell_sz << ", " << c1.cells << ", " << g->v_size << std::endl;
                    //std::cout << "acc mismatch " << cell_sz << ", " << c2.cells << ", " << g->v_size << std::endl;
                    if(cell_sz == 1)
                        individualize_later.push_back(ind_v1);
                    all_certified = false;
                    hard_reset = true;
                    break;
                }

                _automorphism_supp.reset();
                if (c1.cells != g->v_size) { // touched_colors doesn't work properly when early-out is used
                    // read automorphism
                    for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                        const int _c = touched_color_list.arr[j];
                        int f = 0;
                        while (f < c1.ptn[_c] + 1) {
                            const int i = _c + f;
                            ++f;
                            if (c1.lab[i] != c2.lab[i]) {
                                _automorphism.arr[c1.lab[i]] = c2.lab[i];
                                _automorphism_supp.push_back(c1.lab[i]);
                            }
                        }
                    }
                } else {
                    for (int i = 0; i < g->v_size; ++i) {
                        if (c1.lab[i] != c2.lab[i]) {
                            _automorphism.arr[c1.lab[i]] = c2.lab[i];
                            _automorphism_supp.push_back(c1.lab[i]);
                        }
                    }
                }

                certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr, _automorphism_supp.cur_pos,
                                                         _automorphism_supp.arr);
                assert(certify ? R1.certify_automorphism(g, _automorphism.arr) : true);
                //std::cout << "(pre-red) " << it_c << " sparse-ir certify: " << certify << " all_certified: " << all_certified << std::endl;
                all_certified = certify && all_certified;
                if (certify) {
                    pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos,
                                         _automorphism_supp.arr,
                                         consume);
                    // reset c2 to color_cache
                    c2.cells = color_cache.cells;
                    if (c2.cells != g->v_size) {
                        for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                            const int _c = touched_color_list.arr[j];
                            int f = 0;
                            while (f < c1.ptn[_c] + 1) {
                                const int i = _c + f;
                                ++f;
                                c2.ptn[i] = color_cache.ptn[i];
                                c2.lab[i] = color_cache.lab[i];
                                c2.vertex_to_col[c2.lab[i]] = color_cache.vertex_to_col[color_cache.lab[i]];
                                c2.vertex_to_lab[c2.lab[i]] = color_cache.vertex_to_lab[color_cache.lab[i]];
                            }
                        }
                    } else {
                        //std::cout << "this here" << std::endl;
                    }
                }

                reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                automorphism_supp.reset();
            }

            // reset c2 to c1
            c2.cells = c1.cells;
            if (c1.cells != g->v_size && !hard_reset) {
                for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                    const int _c = touched_color_list.arr[j];
                    int f = 0;
                    while (f < c1.ptn[_c] + 1) {
                        const int i = _c + f;
                        ++f;
                        c2.ptn[i] = c1.ptn[i];
                        c2.lab[i] = c1.lab[i];
                        c2.vertex_to_col[c2.lab[i]] = c1.vertex_to_col[c1.lab[i]];
                        c2.vertex_to_lab[c2.lab[i]] = c1.vertex_to_lab[c1.lab[i]];
                    }
                }

                for(int i = 0; i < g->v_size; ++i) {
                    assert(c2.vertex_to_col[i] == c1.vertex_to_col[i]);
                    assert(c2.vertex_to_lab[i] == c1.vertex_to_lab[i]);
                    assert(c2.lab[i] == c1.lab[i]);
                    assert(c2.ptn[i] == c1.ptn[i]);
                }
            } else {
                c2.copy_force(&c1);
                c2.copy_ptn(&c1);

                for(int i = 0; i < g->v_size; ++i) {
                    assert(c2.vertex_to_col[i] == c1.vertex_to_col[i]);
                    assert(c2.vertex_to_lab[i] == c1.vertex_to_lab[i]);
                    assert(c2.lab[i] == c1.lab[i]);
                    assert(c2.ptn[i] == c1.ptn[i]);
                }
            }

            if(all_certified) {
                //std::cout << "all certified, found orbit" << std::endl;
                assert(cell_sz == original_c.ptn[cell]);
                multiply_to_group_size(cell_sz + 1);
                individualize_later.push_back(ind_v1);
            }
        }

        if(!individualize_later.empty()) {
            PRINT("(pre-red) sparse-ir completed orbits: " << individualize_later.size());
            for (int i = 0; i < individualize_later.size(); ++i) {
                const int ind_vert = individualize_later[i];
                const int init_c = R1.individualize_vertex(&original_c, ind_vert);
                R1.refine_coloring(g, &original_c, &I1, init_c, nullptr, -1, -1, nullptr, nullptr, nullptr);
            }
            for (int i = 0; i < g->v_size; ++i) {
                colmap[i] = original_c.vertex_to_col[i];
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
            assert(orig_v_from < domain_size);
            assert(orig_v_from >= 0);
            assert(orig_v_to < domain_size);
            assert(orig_v_to >= 0);
            assert(automorphism.arr[orig_v_from] == orig_v_from);
            automorphism.arr[orig_v_from] = orig_v_to;
            automorphism_supp.push_back(orig_v_from);

            assert(canonical_recovery_string[orig_v_to].size() == canonical_recovery_string[orig_v_from].size());

            for(int j = 0; j < canonical_recovery_string[orig_v_to].size(); ++j) {
                const int v_from_t = canonical_recovery_string[orig_v_from][j];
                const int v_to_t   = canonical_recovery_string[orig_v_to][j];
                assert(automorphism.arr[v_from_t] == v_from_t);
                automorphism.arr[v_from_t] = v_to_t;
                automorphism_supp.push_back(v_from_t);
            }
        }

        consume(domain_size, automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
        reset_automorphism(automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
        automorphism_supp.reset();
    }

    void pre_consumer(int _n, int* _automorphism, int _supp, int* _automorphism_supp, dejavu_consumer consume) {
        meld_translation_layers();
        automorphism_supp.reset();

        for(int i = 0; i < _supp; ++i) {
            const int v_from      = _automorphism_supp[i];
            const int orig_v_from = backward_translation[v_from];
            const int v_to        = _automorphism[v_from];
            assert(v_from != v_to);
            const int orig_v_to   = backward_translation[v_to];
            assert(v_from < backward_translation.size());
            assert(v_from >= 0);
            assert(v_to < backward_translation.size());
            assert(v_to >= 0);
            assert(orig_v_from < domain_size);
            assert(orig_v_from >= 0);
            assert(orig_v_to < domain_size);
            assert(orig_v_to >= 0);
            assert(automorphism.arr[orig_v_from] == orig_v_from);
            automorphism.arr[orig_v_from] = orig_v_to;
            automorphism_supp.push_back(orig_v_from);

            assert(canonical_recovery_string[orig_v_to].size() == canonical_recovery_string[orig_v_from].size());

            for(int j = 0; j < canonical_recovery_string[orig_v_to].size(); ++j) {
                const int v_from_t = canonical_recovery_string[orig_v_from][j];
                const int v_to_t   = canonical_recovery_string[orig_v_to][j];
                assert(automorphism.arr[v_from_t] == v_from_t);
                automorphism.arr[v_from_t] = v_to_t;
                automorphism_supp.push_back(v_from_t);
            }
        }

        consume(domain_size, automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
        reset_automorphism(automorphism.arr, automorphism_supp.cur_pos, automorphism_supp.arr);
        automorphism_supp.reset();
    }

    mark_set seen_vertex;
    mark_set seen_color;
    bool init_quotient_arrays = false;
    std::vector<int> worklist;

    void compute_quotient_graph_components_update(sgraph* g, coloring<int>* c1, dejavu_consumer consume) {
        if(!init_quotient_arrays) {
            seen_vertex.initialize(g->v_size);
            seen_color.initialize(g->v_size);

            worklist.reserve(g->v_size);

            quotient_component_worklist_v.reserve(g->v_size + 32);
            quotient_component_worklist_col.reserve(g->v_size + 32);
            quotient_component_worklist_col_sz.reserve(g->v_size + 32);
            init_quotient_arrays = true;
        }

        seen_vertex.reset();
        seen_color.reset();

        int touched_vertices = 0;

        quotient_component_workspace.clear();
        quotient_component_worklist_col.clear();
        quotient_component_worklist_col_sz.clear();
        quotient_component_worklist_boundary_swap.swap(quotient_component_worklist_boundary);
        quotient_component_worklist_boundary.clear();
        quotient_component_touched_swap.clear();
        quotient_component_touched_swap.swap(quotient_component_touched);

        if(quotient_component_touched_swap.empty()) {
            quotient_component_worklist_v.clear();
            for (int vs = 0; vs < g->v_size; ++vs) {
                if (c1->ptn[c1->vertex_to_col[vs]] == 0) {// ignore discrete vertices
                    seen_vertex.set(vs);
                }
            }

            int component = 0;
            int current_component_sz = 0;
            for (int vs = 0; vs < g->v_size; ++vs) {
                if (seen_vertex.get(vs))
                    continue;
                worklist.push_back(vs);

                while (!worklist.empty()) {
                    const int next_v = worklist.back();
                    worklist.pop_back();
                    if (seen_vertex.get(next_v))
                        continue;
                    /*if (c1->ptn[c1->vertex_to_col[next_v]] == 0) {// ignore discrete vertices
                        seen_vertex.set(next_v);
                        continue;
                    }*/

                    ++touched_vertices;
                    current_component_sz += 1;
                    seen_vertex.set(next_v);
                    quotient_component_worklist_v.push_back(next_v);
                    for (int i = 0; i < g->d[next_v]; ++i) {
                        if (!seen_vertex.get(g->e[g->v[next_v] + i]))
                            worklist.push_back(g->e[g->v[next_v] + i]); // neighbours
                    }
                    const int col = c1->vertex_to_col[next_v];
                    assert(next_v < g->v_size);
                    assert(col < g->v_size);
                    if (!seen_color.get(col)) {
                        quotient_component_worklist_col.push_back(col);
                        quotient_component_worklist_col_sz.push_back(c1->ptn[col]);
                        seen_color.set(col);
                        for (int i = 0; i < c1->ptn[col] + 1; ++i) {
                            assert(col + i < g->v_size);
                            assert(c1->vertex_to_col[c1->lab[col + i]] == c1->vertex_to_col[next_v]);
                            if (c1->lab[col + i] != next_v) {
                                if (!seen_vertex.get(c1->lab[col + i]))
                                    worklist.push_back(c1->lab[col + i]);
                            }
                        }
                    }
                }
                quotient_component_worklist_boundary.emplace_back(
                        std::pair<int, int>(quotient_component_worklist_col.size(),
                                            quotient_component_worklist_v.size()));
                //std::cout << current_component_sz << " ";
                current_component_sz = 0;
                ++component;
            }
            //std::cout << "[" << component  << "]" << std::endl;
        } else {
            // go component by component and only check old touched components
            int old_component = 0;
            int new_component = 0;
            int touched_comp_i = 0;
            quotient_component_worklist_col.clear();
            quotient_component_worklist_col_sz.clear();

            while(old_component < quotient_component_worklist_boundary_swap.size()) {
                int next_touched_old_component = -1;
                if(touched_comp_i < quotient_component_touched_swap.size()) {
                    next_touched_old_component = quotient_component_touched_swap[touched_comp_i];
                }
                //std::cout << "next touched " << next_touched_old_component << std::endl;

                if (old_component == next_touched_old_component) {
                    //std::cout << "computing old_component " << old_component << std::endl;
                    ++touched_comp_i;

                    int component_vstart = 0;
                    int component_cstart = 0;
                    if(old_component > 0) {
                        component_vstart = quotient_component_worklist_boundary_swap[old_component-1].second;
                        component_cstart = quotient_component_worklist_boundary_swap[old_component-1].first;
                    }
                    const int component_vend   = quotient_component_worklist_boundary_swap[old_component].second;
                    const int component_cend   = quotient_component_worklist_boundary_swap[old_component].first;

                    int v_write_pos = component_vstart;
                    quotient_component_workspace.clear();
                    int discrete_vertices = 0;
                    for (int vi = component_vstart; vi < component_vend; ++vi) {
                        const int vs = quotient_component_worklist_v[vi]; // need to copy this away
                        if (c1->ptn[c1->vertex_to_col[vs]] == 0) { // set up fake discrete vertex component, not gonna use it
                            assert(c1->vertex_to_col[vs] == c1->vertex_to_lab[vs]);
                            seen_vertex.set(vs);
                            ++v_write_pos;
                            ++discrete_vertices;
                        } else {
                            quotient_component_workspace.push_back(vs);
                        }
                    }

                    if (discrete_vertices > 0) {
                        quotient_component_worklist_boundary.emplace_back( // discrete vertex component
                                std::pair<int, int>(quotient_component_worklist_col.size(),
                                                    v_write_pos));
                        ++new_component;
                    }

                    int current_component_sz = 0;
                    for (int vi = 0; vi < quotient_component_workspace.size(); ++vi) {
                        const int vs = quotient_component_workspace[vi];
                        if (seen_vertex.get(vs))
                            continue;
                        worklist.push_back(vs);

                        while (!worklist.empty()) {
                            const int next_v = worklist.back();
                            worklist.pop_back();
                            if (seen_vertex.get(next_v))
                                continue;
                            if (c1->ptn[c1->vertex_to_col[next_v]] == 0) {// ignore discrete vertices
                                assert(c1->vertex_to_col[next_v] == c1->vertex_to_lab[next_v]);
                                seen_vertex.set(next_v);
                                continue;
                            }

                            ++touched_vertices;
                            current_component_sz += 1;
                            seen_vertex.set(next_v);
                            quotient_component_worklist_v[v_write_pos] = next_v;
                            ++v_write_pos;
                            for (int i = 0; i < g->d[next_v]; ++i) {
                                if (!seen_vertex.get(g->e[g->v[next_v] + i]))
                                    worklist.push_back(g->e[g->v[next_v] + i]); // neighbours
                            }
                            const int col = c1->vertex_to_col[next_v];
                            assert(next_v < g->v_size);
                            assert(col < g->v_size);
                            if (!seen_color.get(col)) {
                                quotient_component_worklist_col.push_back(col);
                                quotient_component_worklist_col_sz.push_back(c1->ptn[col]);
                                seen_color.set(col);
                                for (int i = 0; i < c1->ptn[col] + 1; ++i) {
                                    assert(col + i < g->v_size);
                                    assert(c1->vertex_to_col[c1->lab[col + i]] == c1->vertex_to_col[next_v]);
                                    if (c1->lab[col + i] != next_v) {
                                        if (!seen_vertex.get(c1->lab[col + i]))
                                            worklist.push_back(c1->lab[col + i]);
                                    }
                                }
                            }
                        }
                        quotient_component_touched.push_back(new_component);
                        quotient_component_worklist_boundary.emplace_back(
                                std::pair<int, int>(quotient_component_worklist_col.size(),
                                                    v_write_pos));
                        //std::cout << current_component_sz << " ";
                        current_component_sz = 0;
                        ++new_component;
                    }

                    assert(v_write_pos == quotient_component_worklist_boundary_swap[old_component].second);
                    ++old_component;
                } else {
                    //std::cout << "skipping old_component " << old_component << std::endl;
                    ++new_component;
                    quotient_component_worklist_boundary.emplace_back(
                            std::pair<int, int>(quotient_component_worklist_col.size(),quotient_component_worklist_boundary_swap[old_component].second));
                    ++old_component; // ignore component, just push old boundaries
                }
            }
            //std::cout << "[" << new_component  << "]" << std::endl;
            //std::cout << "[" << quotient_component_touched.size() << ", " << touched_vertices << "]" << std::endl;
        }
    }

    int sparse_ir_col_sz_2_quotient_components(sgraph* g, int* colmap, dejavu_consumer consume, int num_paths) {
        if(g->v_size <= 1)
            return 0;

        quotient_component_touched.clear();
        quotient_component_touched_swap.clear();

        coloring<int> c1;
        g->initialize_coloring(&c1, colmap);

        for (int i = 0; i < g->v_size; ++i)
            colmap[i] = c1.vertex_to_col[i];

        int global_automorphisms_found = 0;

        std::vector<int> save_colmap_v_to_col;
        std::vector<int> save_colmap_v_to_lab;
        std::vector<int> save_colmap_ptn;
        int save_cell_number = -1;

        if (!ir_quotient_component_init) {
            _automorphism.initialize(g->v_size);
            _automorphism_supp.initialize(g->v_size);
            for (int i = 0; i < g->v_size; ++i)
                _automorphism.arr[i] = i;
            touched_color.initialize(g->v_size);
            touched_color_list.initialize(g->v_size);
            orbit.initialize(g->v_size);
            ir_quotient_component_init = true;
        }

        coloring<int> c2;
        c2.copy(&c1);

        for(int x = 0; x < num_paths; ++x) {
            if(x > 0 && quotient_component_touched.empty()) {
                break;
            }

            compute_quotient_graph_components_update(g, &c1, consume);

            int automorphisms_found = 0;

            save_colmap.reserve(g->v_size);
            save_colmap.clear();

            quotient_component_touched_swap.clear();
            quotient_component_touched_swap.swap(quotient_component_touched);

            touched_color.reset();
            touched_color_list.reset();

            invariant I1, I2;
            I1.only_acc = true;
            I2.only_acc = true;

            bool certify = true;
            int quotient_component_start_pos = 0;
            int quotient_component_start_pos_v = 0;

            int component = 0;
            int next_touched_component_i = 0;
            int next_touched_component = -1;

            int quotient_component_end_pos = quotient_component_worklist_boundary[component].first;
            int quotient_component_end_pos_v = quotient_component_worklist_boundary[component].second;

            while (quotient_component_start_pos < quotient_component_worklist_col.size()) { // TODO: get rid of the color array?
                assert(quotient_component_start_pos_v < quotient_component_worklist_v.size());
                // additional loop here for components?
                int start_search_here = quotient_component_start_pos;

                //if(quotient_component_start_pos == quotient_component_end_pos)
                //    std::cout << "skipped" << std::endl;

                assert(quotient_component_worklist_col[quotient_component_start_pos] >= 0);
                assert(quotient_component_worklist_v[quotient_component_start_pos_v] >= 0);
                certify = true;
                bool touched_current_component = false;

                if (next_touched_component_i < quotient_component_touched_swap.size())
                    next_touched_component = quotient_component_touched_swap[next_touched_component_i];

                while (certify) { // (component == next_touched_component || quotient_component_touched_swap.empty())
                    ++next_touched_component_i;
                    // select a color class of size 2
                    bool only_discrete_prev = true;
                    int cell = -1;
                    for (int _i = start_search_here; _i < quotient_component_end_pos; ++_i) {
                        const int col    = quotient_component_worklist_col[_i];
                        const int col_sz = quotient_component_worklist_col_sz[_i];

                        if (only_discrete_prev)
                            start_search_here = _i;

                        if (col == -1) {// reached end of component
                            break;
                        }
                        for (int i = col; i < col + col_sz + 1; ++i) {
                            if (c1.vertex_to_col[c1.lab[i]] != i)
                                continue; // not a color
                            if (c1.ptn[i] > 0 && only_discrete_prev) {
                                //start_search_here = i;
                                only_discrete_prev = false;
                            }
                            //if((c1.ptn[i] > 0) && (g->d[c1.lab[i]] == 0 || g->d[c1.lab[i]] == 1)) {
                            //    cell = i;
                            //    break;
                            //}
                            if (c1.ptn[i] == 1) {
                                cell = i;
                                break;
                            }
                        }

                        if(cell >= 0 && (g->d[c1.lab[cell]] == 0 || g->d[c1.lab[cell]] == 1)) {
                            break;
                        }

                        if(cell>= 0)
                            break;
                    }

                    //if(cell >= 0 && (g->d[c1.lab[cell]] == 0 || g->d[c1.lab[cell]] == 1)) {
                    //    cell = -1;
                    //}

                    if (cell == -1) {
                        break;
                    } else {
                    }

                    touched_color.reset();
                    touched_color_list.reset();

                    I1.acc = 0;
                    I2.acc = 0;

                    const int ind_v1 = c1.lab[cell];
                    const int ind_v2 = c1.lab[cell + 1];
                    const int init_c1 = R1.individualize_vertex(&c1, ind_v1);

                    touched_color.set(cell);
                    touched_color.set(init_c1);
                    touched_color_list.push_back(cell);
                    touched_color_list.push_back(init_c1);

                    R1.refine_coloring(g, &c1, &I1, init_c1, nullptr, -1, -1, nullptr, &touched_color,
                                       &touched_color_list);

                    const int init_c2 = R1.individualize_vertex(&c2, ind_v2);
                    R1.refine_coloring(g, &c2, &I2, init_c2, nullptr, -1, -1, nullptr, &touched_color,
                                       &touched_color_list);

                    if (I1.acc != I2.acc) {
                        std::cout << "probably broken code 0" << std::endl;
                        for (int i = 0; i < g->v_size; ++i) { // TODO only use component
                            colmap[i] = c1.vertex_to_col[i];
                        }
                        I2.acc = I1.acc;
                        c2.copy(&c1);
                        continue;
                    }

                    _automorphism_supp.reset();
                    if (c1.cells != g->v_size) { // touched_colors doesn't work properly when early-out is used
                        // read automorphism
                        for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                            const int _c = touched_color_list.arr[j];
                            int f = 0;
                            while (f < c1.ptn[_c] + 1) {
                                const int i = _c + f;
                                ++f;
                                if (c1.lab[i] != c2.lab[i]) {
                                    _automorphism.arr[c1.lab[i]] = c2.lab[i];
                                    _automorphism_supp.push_back(c1.lab[i]);
                                }
                            }
                        }
                    } else {
                        //for (int i = 0; i < g->v_size; ++i) {
                        for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                            const int v = quotient_component_worklist_v[_i];
                            const int i = c2.vertex_to_lab[v];
                            if (c1.lab[i] != c2.lab[i]) {
                                assert(i >= 0);
                                assert(i < g->v_size);
                                _automorphism.arr[c1.lab[i]] = c2.lab[i];
                                _automorphism_supp.push_back(c1.lab[i]);
                            }
                        }
                    }

                    touched_current_component = true;
                    certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr, _automorphism_supp.cur_pos,
                                                             _automorphism_supp.arr);
                    assert(certify ? R1.certify_automorphism(g, _automorphism.arr) : true);
                    //std::cout << "(pre-red) sparse-ir certify: " << certify << std::endl;
                    if (certify) {
                        ++automorphisms_found;
                        pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos,
                                             _automorphism_supp.arr,
                                             consume);
                        multiply_to_group_size(2);

                        // reset c2 to c1
                        if (c1.cells != g->v_size) {
                            for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                                const int _c = touched_color_list.arr[j];
                                int f = 0;
                                c2.cells = c1.cells;
                                while (f < c1.ptn[_c] + 1) {
                                    const int i = _c + f;
                                    ++f;
                                    c2.ptn[i] = c1.ptn[i];
                                    c2.lab[i] = c1.lab[i];
                                    c2.vertex_to_col[c2.lab[i]] = c1.vertex_to_col[c1.lab[i]];
                                    c2.vertex_to_lab[c2.lab[i]] = c1.vertex_to_lab[c1.lab[i]];
                                }
                            }
                        } else {
                            // discrete
                            for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                                const int v = quotient_component_worklist_v[_i];
                                colmap[v] = c1.vertex_to_col[v];
                            }
                            touched_current_component = false;
                            break;
                        }
                    } else { // save entire coloring, including lab, ptn, v_to_lab
                        //touched_color.reset();
                        //touched_color_list.reset();
                        save_colmap.clear();
                        save_colmap_v_to_col.clear();
                        save_colmap_v_to_lab.clear();
                        save_colmap_ptn.clear();
                        for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                            const int v = quotient_component_worklist_v[_i];
                            save_colmap_v_to_col.push_back(c1.vertex_to_col[v]);
                            const int lab_pos = c1.vertex_to_lab[v];
                            save_colmap_v_to_lab.push_back(lab_pos);
                            save_colmap_ptn.push_back(c1.ptn[lab_pos]);
                            save_cell_number = c1.cells;
                            assert(c1.lab[lab_pos] == v);
                        }

                        int col = -1;

                        certify = false;

                        reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                        _automorphism_supp.reset();

                        while (true) {
                            // stay within component!
                            col = select_color_component(g, &c1, quotient_component_start_pos,
                                                         quotient_component_end_pos);
                            if (col == -1) break;
                            const int rpos = col + (0 % (c1.ptn[col] + 1));
                            const int v1 = c1.lab[rpos];
                            const int init_color_class1 = R1.individualize_vertex(&c1, v1);

                            if(!touched_color.get(init_color_class1)) {
                               touched_color.set(init_color_class1);
                               touched_color_list.push_back(init_color_class1);
                            }
                            if(!touched_color.get(col)) {
                                touched_color.set(col);
                                touched_color_list.push_back(col);
                            }

                            R1.refine_coloring(g, &c1, &I1, init_color_class1, nullptr, -1, -1,
                                               nullptr, &touched_color, &touched_color_list);

                            //const int rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c2.ptn[col] + 1));
                            const int v2 = c2.lab[rpos];
                            const int init_color_class2 = R1.individualize_vertex(&c2, v2);
                            R1.refine_coloring(g, &c2, &I2, init_color_class2, nullptr, -1, -1,
                                               nullptr, &touched_color, &touched_color_list);

                            /*col = select_color_component(g, &c1, quotient_component_start_pos,
                                                         quotient_component_end_pos);

                            if (col == -1) break;*/

                            for (int j = 0; j < touched_color_list.cur_pos; ++j) { // check incremental singleton automorphisms
                                const int _c = touched_color_list.arr[j];
                                int f = 0;
                                if(c1.ptn[_c] == 0) {
                                    while (f < c1.ptn[_c] + 1) {
                                        const int i = _c + f;
                                        ++f;
                                        if (c1.lab[i] != c2.lab[i]) {
                                            _automorphism.arr[c1.lab[i]] = c2.lab[i];
                                            _automorphism_supp.push_back(c1.lab[i]);
                                        }
                                    }
                                }
                            }

                            certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr,
                                                                     _automorphism_supp.cur_pos,
                                                                     _automorphism_supp.arr);

                            //reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                            //_automorphism_supp.reset();
                            touched_color.reset();
                            touched_color_list.reset();

                            if(certify)    break;
                            //if (col == -1) break;
                        }

                        if(!certify) {
                            reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                            _automorphism_supp.reset();

                            for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                                const int i = quotient_component_worklist_v[_i];
                                const int lab_p = c1.vertex_to_lab[i];
                                assert(i >= 0);
                                assert(i < g->v_size);
                                if (c1.lab[lab_p] != c2.lab[lab_p]) {
                                    _automorphism.arr[c1.lab[lab_p]] = c2.lab[lab_p];
                                    _automorphism_supp.push_back(c1.lab[lab_p]);
                                }
                            }
                        }

                        certify = certify || R1.certify_automorphism_sparse(g, colmap, _automorphism.arr,
                                                                 _automorphism_supp.cur_pos,
                                                                 _automorphism_supp.arr);
                        //std::cout << "(prep-red) sparse-ir random path certify: " << certify << std::endl;
                        if (certify) {
                            assert(R1.certify_automorphism(g, _automorphism.arr));
                            touched_current_component = true;
                            ++automorphisms_found;
                            pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos,
                                                 _automorphism_supp.arr,
                                                 consume);
                            multiply_to_group_size(2);
                            /*const int test_col = select_color_component(g, &c1, quotient_component_start_pos, quotient_component_end_pos);
                            if(test_col == -1) {
                                touched_current_component = false;
                            }*/
                        }
                        reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                        _automorphism_supp.reset();

                        if (certify) {
                            //std::cout << "certified normal, reading coloring" << std::endl;
                            /*for (int i = 0; i < g->v_size; ++i) {
                                colmap[i] = save_colmap[i];
                            }*/
                            int f = 0;
                            for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) { // TODO could use touched...
                                const int v = quotient_component_worklist_v[_i];
                                colmap[v] = save_colmap_v_to_col[f]; // TODO reset c1 / c2 instead such that paths can just be restarted

                                c1.vertex_to_col[v] = save_colmap_v_to_col[f];
                                const int lab_pos   = save_colmap_v_to_lab[f];
                                const int ptn_at_lab_pos = save_colmap_ptn[f];
                                c1.vertex_to_lab[v] = lab_pos;
                                c1.lab[lab_pos]     = v;
                                c1.ptn[lab_pos]     = ptn_at_lab_pos;
                                c1.cells = save_cell_number;

                                c2.vertex_to_col[v] = save_colmap_v_to_col[f];
                                c2.vertex_to_lab[v] = lab_pos;
                                c2.lab[lab_pos]     = v;
                                c2.ptn[lab_pos]     = ptn_at_lab_pos;
                                c2.cells = save_cell_number;

                                ++f;
                            }
                            for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                                const int v = quotient_component_worklist_v[_i];
                                const int list_index = _i - quotient_component_start_pos_v;
                                assert(c1.ptn[c1.vertex_to_col[v]] >= 0);
                                assert(c1.lab[c1.vertex_to_lab[v]] == v);
                                assert(c1.vertex_to_col[c1.lab[c1.vertex_to_col[v]]] == c1.vertex_to_col[v]);

                                assert(c2.ptn[c2.vertex_to_col[v]] >= 0);
                                assert(c2.lab[c2.vertex_to_lab[v]] == v);
                                assert(c2.vertex_to_col[c2.lab[c2.vertex_to_col[v]]] == c2.vertex_to_col[v]);
                            }
                        } else {
                            std::cout << "did not certifyX" << (c1.cells == g->v_size) << std::endl;
                            // TODO hmmmm
                            // TODO: mark this component such that it is never touched again?
                            // c1/c2 are now 'broken' for this component, never touch it again
                            // can't reliably find automorphisms here
                            touched_current_component = false;
                        }

                        break;
                    }

                    reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                    automorphism_supp.reset();

                    if (certify) {
                        if (c1.cells == g->v_size) {
                            //std::cout << "probably broken code 2" << std::endl;
                            for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                                const int i = quotient_component_worklist_v[_i];
                                assert(i >= 0);
                                assert(i < g->v_size);
                                colmap[i] = c1.vertex_to_col[i];
                            }
                        } else {
                            for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                                const int _c = touched_color_list.arr[j];
                                int f = 0;
                                while (f < c1.ptn[_c] + 1) {
                                    assert(c1.vertex_to_col[c1.lab[_c + f]] == _c);
                                    colmap[c1.lab[_c + f]] = c1.vertex_to_col[c1.lab[_c + f]]; // should be c1 or c2?
                                    ++f;
                                }
                            }
                        }
                    }
                }

                if(c1.cells == g->v_size)
                    break;

                if (touched_current_component) {
                    quotient_component_touched.push_back(component);
                    del_discrete_edges_inplace_component(g, &c1, quotient_component_start_pos);
                }

                quotient_component_start_pos = quotient_component_worklist_boundary[component].first;
                quotient_component_start_pos_v = quotient_component_worklist_boundary[component].second;

                ++component;

                quotient_component_end_pos = quotient_component_worklist_boundary[component].first;
                quotient_component_end_pos_v = quotient_component_worklist_boundary[component].second;
            }
            global_automorphisms_found += automorphisms_found;

            if(c1.cells == g->v_size)
                break;
            if(automorphisms_found == 0)
                break;
        }
        //for(int i = 0; i < quotient_component_touched.size(); ++i)
        //    std::cout << quotient_component_touched[i] << " ";
        //std::cout << std::endl;
        return global_automorphisms_found;
    }

    void add_automorphism_to_orbit(tiny_orbit* orbit, int* automorphism, int nsupp, int* supp) {
        for(int i = 0; i < nsupp; ++i) {
            orbit->combine_orbits(automorphism[supp[i]], supp[i]);
        }
    }

    int select_color_component_min_cost(sgraph*g, coloring<int>* c1, int component_start_pos, int component_end_pos, int max_cell_size) {
        bool only_discrete_prev = true;
        int cell = -1;
        int cell_d = 1;
        for (int _i = component_start_pos; _i < component_end_pos; ++_i) {
            const int col = quotient_component_worklist_col[_i];
            const int col_sz = quotient_component_worklist_col_sz[_i];

            //if (only_discrete_prev)
            //    start_search_here = _i;

            if (col == -1) {// reached end of component
                break;
            }
            for (int i = col; i < col + col_sz + 1; ++i) {
                if (c1->vertex_to_col[c1->lab[i]] != i)
                    continue; // not a color
                if (c1->ptn[i] > 0 && only_discrete_prev) {
                    //start_search_here = i;
                    only_discrete_prev = false;
                }
                const int i_d = g->d[c1->lab[i]];
                //if(c1->ptn[i] >= 1 && (i_d == 0 || i_d == 1)) {
                //    cell = i;
                //    break;
                //}
                if(c1->ptn[i] >= 1 && (i_d == 0)) {
                    cell = i;
                    break;
                }
                if (c1->ptn[i] >= 1 && c1->ptn[i] + 1 <= max_cell_size && (i_d != 1)) { // && (i_d != 1)
                    if(cell == -1 || c1->ptn[i] < c1->ptn[cell]) {
                        cell = i;
                        cell_d = i_d;
                    }
                }
            }
            //if(cell != -1 && (g->d[c1->lab[cell]] == 0 || g->d[c1->lab[cell]] == 1)) {
            //    cell = -1;
            //    break;
            //}
            if(cell != -1 && (g->d[c1->lab[cell]] == 0)) {
                cell = -1;
                break;
            }

            if(cell != -1 && c1->ptn[cell] < 3)
                break;
        }

        return cell;
    }

    int sparse_ir_quotient_components(sgraph* g, int* colmap, dejavu_consumer consume, int max_col_size, int num_paths) {
        if(g->v_size <= 1 || num_paths <= 0)
            return 0;

        // TODO this needs to use exploit orbits at least a little bit

        avg_support_sparse_ir     = 0;
        avg_reached_end_of_component = 0;
        int avg_support_sparse_ir_num = 0;
        int avg_support_sparse_ir_num2 = 0;

        //std::cout << "sparse_ir_quotient_components" << std::endl;

        quotient_component_touched.clear();
        quotient_component_touched_swap.clear();

        mark_set touched_color_cache;
        work_list_t<int> touched_color_list_cache;
        touched_color_cache.initialize(g->v_size);
        touched_color_list_cache.initialize(g->v_size);

        coloring<int> c1;
        g->initialize_coloring(&c1, colmap);

        for (int i = 0; i < g->v_size; ++i)
            colmap[i] = c1.vertex_to_col[i];

        int global_automorphisms_found = 0;

        //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() * ((1 * 5) * 5135235);
        //int selector_seed = seed;

        if (!ir_quotient_component_init) {
            _automorphism.initialize(g->v_size);
            _automorphism_supp.initialize(g->v_size);
            for (int i = 0; i < g->v_size; ++i)
                _automorphism.arr[i] = i;
            touched_color.initialize(g->v_size);
            touched_color_list.initialize(g->v_size);
            ir_quotient_component_init = true;
            orbit.initialize(g->v_size);
        }

        reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
        _automorphism_supp.reset();

        coloring<int> c2;
        c2.copy(&c1);
        c2.copy_ptn(&c1);

        coloring<int> c3;
        c3.copy(&c1);
        c3.copy_ptn(&c1);

        for(int x = 0; x < num_paths; ++x) {
            if(x > 0 && quotient_component_touched.empty()) {
                break;
            }

            //if(x == 4) {
            //    del_discrete_edges_inplace(g, &c1); // TODO: do this per-component
            //}
            //std::cout << "path: " << x << std::endl;
            compute_quotient_graph_components_update(g, &c1, consume);

            int automorphisms_found = 0;

            quotient_component_touched_swap.clear();
            quotient_component_touched_swap.swap(quotient_component_touched);

            touched_color.reset();
            touched_color_list.reset();

            invariant I1, I2;
            I1.only_acc = true;
            I2.only_acc = true;

            bool certify = true;
            int quotient_component_start_pos = 0;
            int quotient_component_start_pos_v = 0;

            int component = 0;
            int next_touched_component_i = 0;
            int next_touched_component = -1;

            int quotient_component_end_pos = quotient_component_worklist_boundary[component].first;
            int quotient_component_end_pos_v = quotient_component_worklist_boundary[component].second;

            while (quotient_component_start_pos < quotient_component_worklist_col.size()) { // TODO: get rid of the color array?
                assert(quotient_component_start_pos_v < quotient_component_worklist_v.size());
                int start_search_here = quotient_component_start_pos;

                assert(quotient_component_worklist_col[quotient_component_start_pos] >= 0);
                assert(quotient_component_worklist_v[quotient_component_start_pos_v] >= 0);
                certify = true;
                bool touched_current_component = false;

                if (next_touched_component_i < quotient_component_touched_swap.size())
                    next_touched_component = quotient_component_touched_swap[next_touched_component_i];

                while (certify) { // (component == next_touched_component || quotient_component_touched_swap.empty())
                    ++next_touched_component_i;
                    // select a color class of size 2
                    int cell = select_color_component_min_cost(g, &c1, quotient_component_start_pos, quotient_component_end_pos, max_col_size);

                    orbit.reset();

                    //if(cell >= 0 && (g->d[c1.lab[cell]] == 0 || g->d[c1.lab[cell]] == 1)) {
                    //    cell = -1;
                    //}

                    if (cell == -1) {
                        //std::cout << "skipping component " << quotient_component_start_pos << std::endl;
                        break;
                    }

                    const int cell_sz = c1.ptn[cell];

                    //std::cout << "cell " << cell << " cell_sz " << cell_sz << std::endl;

                    touched_color.reset();
                    touched_color_list.reset();


                    I1.acc = 0;
                    I2.acc = 0;

                    const int ind_v1 = c1.lab[cell];
                    const int ind_v2 = c1.lab[cell + 1];
                    const int init_c1 = R1.individualize_vertex(&c1, ind_v1);

                    touched_color.set(cell);
                    touched_color.set(init_c1);
                    touched_color_list.push_back(cell);
                    touched_color_list.push_back(init_c1);
                    assert(cell != init_c1);

                    R1.refine_coloring(g, &c1, &I1, init_c1, nullptr, -1, -1, nullptr, &touched_color,
                                       &touched_color_list); // TODO: save difference of c1 to before individualization now using touched and c2? or have c3 that is only updated once c1 is confirmed?

                    const int init_c2 = R1.individualize_vertex(&c2, ind_v2);
                    R1.refine_coloring(g, &c2, &I2, init_c2, nullptr, -1, -1, nullptr, &touched_color,
                                       &touched_color_list);

                    if (I1.acc != I2.acc) {
                        std::cout << "probably broken code 0" << std::endl; // TODO: I am not recovering both colorings?
                        /*for (int i = 0; i < g->v_size; ++i) { // TODO only use component
                            colmap[i] = c1.vertex_to_col[i];
                        }*/
                        I2.acc = I1.acc;
                        c1.copy(&c3); // TODO not this
                        touched_current_component = false;
                        break;
                    }

                    reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                    _automorphism_supp.reset();
                    if (c1.cells != g->v_size) { // touched_colors doesn't work properly when early-out is used
                        // read automorphism
                        for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                            const int _c = touched_color_list.arr[j];
                            int f = 0;
                            while (f < c1.ptn[_c] + 1) {
                                const int i = _c + f;
                                ++f;
                                if (c1.lab[i] != c2.lab[i]) {
                                    _automorphism.arr[c1.lab[i]] = c2.lab[i];
                                    _automorphism_supp.push_back(c1.lab[i]);
                                }
                            }
                        }
                    } else {
                        //for (int i = 0; i < g->v_size; ++i) {
                        for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                            assert(_i >= 0);
                            assert(_i < g->v_size);
                            const int v = quotient_component_worklist_v[_i];
                            assert(v >= 0);
                            assert(v < g->v_size);
                            const int i = c2.vertex_to_lab[v];
                            assert(i >= 0);
                            assert(i < g->v_size);
                            if (c1.lab[i] != c2.lab[i]) {
                                _automorphism.arr[c1.lab[i]] = c2.lab[i];
                                _automorphism_supp.push_back(c1.lab[i]);
                            }
                        }
                    }

                    certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr, _automorphism_supp.cur_pos,
                                                             _automorphism_supp.arr);
                    assert(certify ? R1.certify_automorphism(g, _automorphism.arr) : true);
                    //std::cout << "(pre-red) sparse-ir certify: " << certify << std::endl;
                    bool all_certified = certify;
                    if (certify) {
                        //std::cout << "certified sparse" << std::endl;
                        avg_support_sparse_ir += _automorphism_supp.cur_pos;
                        avg_support_sparse_ir_num += 1;
                        pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos,
                                             _automorphism_supp.arr,
                                             consume);
                        assert(_automorphism.arr[ind_v1] == ind_v2);
                        add_automorphism_to_orbit(&orbit, _automorphism.arr, _automorphism_supp.cur_pos,
                                                  _automorphism_supp.arr);

                        for(int k = 2; (k < cell_sz + 1) && all_certified; ++k) {
                            // reset c2 to before individualization
                            assert(touched_color_list.cur_pos <= quotient_component_end_pos_v - quotient_component_start_pos_v);
                            const int ind_vk = c3.lab[cell + k];
                            if(orbit.are_in_same_orbit(ind_v1, ind_vk)) {
                                //std::cout << "skip" << std::endl;
                                continue;
                            }

                            if(g->v_size == c1.cells || touched_color_list.cur_pos >= (0.5*(quotient_component_end_pos_v - quotient_component_start_pos_v))) {
                                for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                                    const int v = quotient_component_worklist_v[_i];
                                    const int lab_pos = c3.vertex_to_lab[v];
                                    c2.vertex_to_col[v] = c3.vertex_to_col[v];
                                    c2.vertex_to_lab[v] = c3.vertex_to_lab[v];
                                    c2.lab[lab_pos] = v;
                                    c2.ptn[lab_pos] = c3.ptn[lab_pos];
                                    c2.cells = c3.cells;
                                }
                            } else {
                                for (int _i = 0; _i < touched_color_list.cur_pos; ++_i) {
                                    const int reset_col = touched_color_list.arr[_i];
                                    if(c3.vertex_to_col[c3.lab[reset_col]] == reset_col ) {
                                        for(int l = 0; l < c3.ptn[reset_col] + 1; ++l) {
                                            const int v = c3.lab[reset_col + l];
                                            const int lab_pos = c3.vertex_to_lab[v];
                                            assert(lab_pos == reset_col + l);
                                            c2.vertex_to_col[v] = c3.vertex_to_col[v];
                                            c2.vertex_to_lab[v] = c3.vertex_to_lab[v];
                                            c2.lab[lab_pos] = v;
                                            c2.ptn[lab_pos] = c3.ptn[lab_pos];
                                        }
                                        c2.cells = c3.cells;
                                    }
                                }
                            }

                            // individualize and test cell + k
                            const int init_ck = R1.individualize_vertex(&c2, ind_vk);
                            I2.acc = 0;
                            R1.refine_coloring(g, &c2, &I2, init_ck, nullptr, -1, -1, nullptr, &touched_color,
                                               &touched_color_list);
                            reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                            _automorphism_supp.reset();

                            if(I1.acc == I2.acc && c1.cells == c2.cells) {
                                if (c1.cells !=
                                    g->v_size) { // touched_colors doesn't work properly when early-out is used
                                    // read automorphism
                                    for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                                        const int _c = touched_color_list.arr[j];
                                        int f = 0;
                                        while (f < c1.ptn[_c] + 1) {
                                            const int i = _c + f;
                                            ++f;
                                            if (c1.lab[i] != c2.lab[i]) {
                                                _automorphism.arr[c1.lab[i]] = c2.lab[i];
                                                _automorphism_supp.push_back(c1.lab[i]);
                                            }
                                        }
                                    }
                                } else {
                                    //for (int i = 0; i < g->v_size; ++i) {
                                    for (int _i = quotient_component_start_pos_v;
                                         _i < quotient_component_end_pos_v; ++_i) {
                                        const int v = quotient_component_worklist_v[_i];
                                        const int i = c2.vertex_to_lab[v];
                                        if (c1.lab[i] != c2.lab[i]) {
                                            assert(i >= 0);
                                            assert(i < g->v_size);
                                            _automorphism.arr[c1.lab[i]] = c2.lab[i];
                                            _automorphism_supp.push_back(c1.lab[i]);
                                        }
                                    }
                                }
                                certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr,
                                                                         _automorphism_supp.cur_pos,
                                                                         _automorphism_supp.arr);
                                if(certify) {
                                    avg_support_sparse_ir += _automorphism_supp.cur_pos;
                                    avg_support_sparse_ir_num += 1;
                                    pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos,
                                                         _automorphism_supp.arr,
                                                         consume);
                                    assert(_automorphism.arr[ind_v1] == ind_vk);
                                    add_automorphism_to_orbit(&orbit, _automorphism.arr, _automorphism_supp.cur_pos,
                                                              _automorphism_supp.arr);
                                }
                            } else {
                                certify = false;
                            }
                            all_certified = certify && all_certified;
                        }

                        if(all_certified) {
                            automorphisms_found += cell_sz ;
                            touched_current_component = true;
                            avg_support_sparse_ir_num2   += 1;
                            //std::cout << "all_certified early" << std::endl;
                            multiply_to_group_size(cell_sz + 1);
                            // reset c2 and c3 to c1
                            if (c1.cells != g->v_size) {
                                for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                                    const int _c = touched_color_list.arr[j];
                                    int f = 0;
                                    c2.cells = c1.cells;
                                    c3.cells = c1.cells;
                                    while (f < c1.ptn[_c] + 1) {
                                        const int i = _c + f;
                                        ++f;
                                        c2.ptn[i] = c1.ptn[i];
                                        c2.lab[i] = c1.lab[i];
                                        c2.vertex_to_col[c2.lab[i]] = c1.vertex_to_col[c1.lab[i]];
                                        c2.vertex_to_lab[c2.lab[i]] = c1.vertex_to_lab[c1.lab[i]];

                                        c3.ptn[i] = c1.ptn[i];
                                        c3.lab[i] = c1.lab[i];
                                        c3.vertex_to_col[c3.lab[i]] = c1.vertex_to_col[c1.lab[i]];
                                        c3.vertex_to_lab[c3.lab[i]] = c1.vertex_to_lab[c1.lab[i]];
                                    }
                                }
                            }
                        }
                    }

                    // TODO: DEBUG
                    //if(!all_certified)
                    //    touched_current_component = false;

                    if(!all_certified) {
                        int col = -1;
                        certify = false;
                        reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                        _automorphism_supp.reset();

                        // reset c2 to before individualization
                        assert(touched_color_list.cur_pos <= quotient_component_end_pos_v - quotient_component_start_pos_v);
                        if(g->v_size == c1.cells || touched_color_list.cur_pos >= (0.5*(quotient_component_end_pos_v - quotient_component_start_pos_v))) {
                            for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                                const int v = quotient_component_worklist_v[_i];
                                const int lab_pos = c3.vertex_to_lab[v];
                                c2.vertex_to_col[v] = c3.vertex_to_col[v];
                                c2.vertex_to_lab[v] = c3.vertex_to_lab[v];
                                c2.lab[lab_pos] = v;
                                c2.ptn[lab_pos] = c3.ptn[lab_pos];
                                c2.cells = c3.cells;
                            }
                        } else {
                            for (int _i = 0; _i < touched_color_list.cur_pos; ++_i) {
                                const int reset_col = touched_color_list.arr[_i];
                                if(c3.vertex_to_col[c3.lab[reset_col]] == reset_col ) {
                                    for(int l = 0; l < c3.ptn[reset_col] + 1; ++l) {
                                        const int v = c3.lab[reset_col + l];
                                        const int lab_pos = c3.vertex_to_lab[v];
                                        assert(lab_pos == reset_col + l);
                                        c2.vertex_to_col[v] = c3.vertex_to_col[v];
                                        c2.vertex_to_lab[v] = c3.vertex_to_lab[v];
                                        c2.lab[lab_pos] = v;
                                        c2.ptn[lab_pos] = c3.ptn[lab_pos];
                                    }
                                    c2.cells = c3.cells;
                                }
                            }
                        }

                        I2.acc = 0;
                        const int init_c2 = R1.individualize_vertex(&c2, ind_v2);
                        R1.refine_coloring(g, &c2, &I2, init_c2, nullptr, -1, -1, nullptr, &touched_color,
                                           &touched_color_list);


                        int num_inds = 0;
                        int touched_start_pos = 0;

                        touched_color_cache.reset();
                        touched_color_list_cache.reset();

                        for(int k = 0; k < touched_color_list.cur_pos; ++k) {
                            touched_color_cache.set(touched_color_list.arr[k]);
                            touched_color_list_cache.arr[k] = touched_color_list.arr[k];
                        }
                        touched_color_list_cache.cur_pos = touched_color_list.cur_pos;


                        while (true) {
                            // stay within component!
                            col = select_color_component(g, &c1, quotient_component_start_pos,
                                                         quotient_component_end_pos);
                            if (col == -1) break;
                            const int rpos = col + (0 % (c1.ptn[col] + 1));
                            const int v1 = c1.lab[rpos];
                            const int init_color_class1 = R1.individualize_vertex(&c1, v1);

                            if(!touched_color_cache.get(init_color_class1)) {
                                touched_color_cache.set(init_color_class1);
                                touched_color_list_cache.push_back(init_color_class1);
                            }
                            if(!touched_color_cache.get(col)) {
                                touched_color_cache.set(col);
                                touched_color_list_cache.push_back(col);
                            }

                            R1.refine_coloring(g, &c1, &I1, init_color_class1, nullptr, -1, -1,
                                               nullptr, &touched_color_cache, &touched_color_list_cache);

                            //const int rpos = col + (intRand(0, INT32_MAX, selector_seed) % (c2.ptn[col] + 1));
                            const int v2 = c2.lab[rpos];
                            const int init_color_class2 = R1.individualize_vertex(&c2, v2);
                            R1.refine_coloring(g, &c2, &I2, init_color_class2, nullptr, -1, -1,
                                               nullptr, &touched_color_cache, &touched_color_list_cache);

                            if(c1.cells == g->v_size) {
                                //std::cout << "filling colorsA" << std::endl;
                                for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                                    const int v = quotient_component_worklist_v[_i];
                                    const int v_col1 = c1.vertex_to_col[v];
                                    if (!touched_color_cache.get(v_col1) && !touched_color.get(v_col1)) {
                                        touched_color_cache.set(v_col1);
                                        touched_color_list_cache.push_back(v_col1);
                                    }
                                }
                            }

                            //assert(I1.acc == I2.acc);

                            for (int j = 0; j < touched_color_list_cache.cur_pos; ++j) { // check incremental singleton automorphisms
                                const int _c = touched_color_list_cache.arr[j];
                                if(!touched_color.get(_c)) {
                                    touched_color.set(_c);
                                    touched_color_list.push_back(_c);
                                }
                                int f = 0;
                                assert(c1.cells != g->v_size || (c1.ptn[_c] == 0));
                                if(c1.ptn[_c] == 0 && c2.ptn[_c] == 0) {
                                    while (f < c1.ptn[_c] + 1) {
                                        const int i = _c + f;
                                        ++f;
                                        if (c1.lab[i] != c2.lab[i]) {
                                            _automorphism.arr[c1.lab[i]] = c2.lab[i];
                                            _automorphism_supp.push_back(c1.lab[i]);
                                        }
                                    }
                                }
                            }

                            ++num_inds;

                            certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr,
                                                                     _automorphism_supp.cur_pos,
                                                                     _automorphism_supp.arr);

                            touched_color_cache.reset();
                            touched_color_list_cache.reset();

                            if(certify)    break;
                        }

                        if(!certify && (g->v_size == c1.cells)) {
                            reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                            _automorphism_supp.reset();
                            for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                                const int _c = touched_color_list.arr[j];
                                int f = 0;
                                assert(c1.cells != g->v_size || (c1.ptn[_c] == 0));
                                if(c1.ptn[_c] == 0) {
                                    while (f < c1.ptn[_c] + 1) {
                                        const int i = _c + f;
                                        ++f;
                                        if (c1.lab[i] != c2.lab[i]) {
                                            _automorphism.arr[c1.lab[i]] = c2.lab[i];
                                            _automorphism_supp.push_back(c1.lab[i]);
                                        }
                                    }
                                }
                            }
                            certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr,
                                                                     _automorphism_supp.cur_pos,
                                                                     _automorphism_supp.arr);
                            //if(certify)
                            //    std::cout << "recheck success " << (g->v_size == c1.cells) << std::endl;
                        }

                        bool certify_before = certify;

                        if (certify) {
                            assert(R1.certify_automorphism(g, _automorphism.arr));
                            avg_support_sparse_ir += _automorphism_supp.cur_pos;
                            avg_support_sparse_ir_num += 1;
                            pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos,
                                                 _automorphism_supp.arr,
                                                 consume);
                            add_automorphism_to_orbit(&orbit, _automorphism.arr, _automorphism_supp.cur_pos,
                                                      _automorphism_supp.arr);
                           // multiply_to_group_size(2);
                           all_certified = true;

                            reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                            _automorphism_supp.reset();

                            for(int k = 2; (k < cell_sz + 1) && all_certified; ++k) {
                                const int ind_vk = c3.lab[cell + k];
                                if(orbit.are_in_same_orbit(ind_v1, ind_vk)) {
                                    //std::cout << "skip" << std::endl;
                                    continue;
                                }

                                // TODO: reset c2 using only touched?
                                // reset c2 to before individualization
                                assert(touched_color_list.cur_pos <= quotient_component_end_pos_v - quotient_component_start_pos_v);
                                if(g->v_size == c1.cells || touched_color_list.cur_pos >= (0.5*(quotient_component_end_pos_v - quotient_component_start_pos_v))) {
                                    for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                                        const int v = quotient_component_worklist_v[_i];
                                        const int lab_pos = c3.vertex_to_lab[v];
                                        c2.vertex_to_col[v] = c3.vertex_to_col[v];
                                        c2.vertex_to_lab[v] = c3.vertex_to_lab[v];
                                        c2.lab[lab_pos] = v;
                                        c2.ptn[lab_pos] = c3.ptn[lab_pos];
                                        c2.cells = c3.cells;
                                    }
                                } else {
                                    for (int _i = 0; _i < touched_color_list.cur_pos; ++_i) {
                                        const int reset_col = touched_color_list.arr[_i];
                                        if(c3.vertex_to_col[c3.lab[reset_col]] == reset_col ) {
                                            for(int l = 0; l < c3.ptn[reset_col] + 1; ++l) {
                                                const int v = c3.lab[reset_col + l];
                                                const int lab_pos = c3.vertex_to_lab[v];
                                                assert(lab_pos == reset_col + l);
                                                c2.vertex_to_col[v] = c3.vertex_to_col[v];
                                                c2.vertex_to_lab[v] = c3.vertex_to_lab[v];
                                                c2.lab[lab_pos] = v;
                                                c2.ptn[lab_pos] = c3.ptn[lab_pos];
                                            }
                                            c2.cells = c3.cells;
                                        }
                                    }
                                }

                                // individualize and test cell + k
                                const int init_ck = R1.individualize_vertex(&c2, ind_vk);
                                if(!touched_color.get(cell)) {
                                    touched_color.set(cell);
                                    touched_color_list.push_back(cell);
                                }
                                if(!touched_color.get(init_ck)) {
                                    touched_color.set(init_ck);
                                    touched_color_list.push_back(init_ck);
                                }
                                I2.acc = 0;
                                R1.refine_coloring(g, &c2, &I2, init_ck, nullptr, -1, -1, nullptr, &touched_color,
                                                   &touched_color_list);
                                reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                                _automorphism_supp.reset();

                                for(int i = 0; i < g->v_size; ++i) {
                                    assert(_automorphism.arr[i] == i);
                                }

                                int repeat_num_inds = 0;

                                bool should_add_ind = false;
                                do {
                                    while (repeat_num_inds < num_inds) {
                                        // stay within component!
                                        col = select_color_component(g, &c2, quotient_component_start_pos,
                                                                     quotient_component_end_pos);
                                        if (col == -1) {
                                            break;
                                        }
                                        const int rpos =
                                                col + (0 % (c2.ptn[col] + 1));
                                        const int v2 = c2.lab[rpos];
                                        const int init_color_class2 = R1.individualize_vertex(&c2, v2);
                                        if (!touched_color.get(init_color_class2)) {
                                            touched_color.set(init_color_class2);
                                            touched_color_list.push_back(init_color_class2);
                                        }
                                        if(!touched_color.get(col)) {
                                            touched_color.set(col);
                                            touched_color_list.push_back(col);
                                        }
                                        R1.refine_coloring(g, &c2, &I2, init_color_class2, nullptr, -1, -1,
                                                           nullptr, &touched_color, &touched_color_list);
                                        ++repeat_num_inds;
                                    }

                                    should_add_ind = false;

                                   // std::cout << I1.acc << ", " << I2.acc << std::endl;

                                    if (I1.acc == I2.acc && c1.cells == c2.cells) {
                                        for (int j = 0; j < touched_color_list.cur_pos; ++j) { // check the singleton automorphism
                                            const int _c = touched_color_list.arr[j];
                                            int f = 0;
                                            if(c1.ptn[_c] == 0) {
                                                while (f < c1.ptn[_c] + 1) {
                                                    const int i = _c + f;
                                                    ++f;
                                                    if (c1.lab[i] != c2.lab[i]) {
                                                        //assert(_automorphism.arr[c1.lab[i]] == c1.lab[i]);
                                                        _automorphism.arr[c1.lab[i]] = c2.lab[i];
                                                        _automorphism_supp.push_back(c1.lab[i]);
                                                    }
                                                }
                                            }
                                        }

                                        certify = R1.certify_automorphism_sparse(g, colmap, _automorphism.arr,
                                                                                 _automorphism_supp.cur_pos,
                                                                                 _automorphism_supp.arr);
                                        should_add_ind = !certify;
                                    } else {
                                        certify = false;
                                    }
                                    if(should_add_ind) {
                                        col = select_color_component(g, &c1, quotient_component_start_pos,
                                                                     quotient_component_end_pos);
                                        if (col == -1) { should_add_ind = false; break; }
                                        const int rpos = col + (0 % (c1.ptn[col] + 1));
                                        const int v1 = c1.lab[rpos];
                                        const int init_color_class1 = R1.individualize_vertex(&c1, v1);
                                        R1.refine_coloring(g, &c1, &I1, init_color_class1, nullptr, -1, -1,
                                                           nullptr, &touched_color, &touched_color_list);
                                        num_inds += 1;

                                        reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                                        _automorphism_supp.reset();
                                    }
                                } while(should_add_ind);

                                all_certified = certify && all_certified;
                                if(certify) {
                                    avg_support_sparse_ir += _automorphism_supp.cur_pos;
                                    avg_support_sparse_ir_num += 1;
                                    pre_consumer_inplace(g->v_size, _automorphism.arr, _automorphism_supp.cur_pos,
                                                         _automorphism_supp.arr,
                                                         consume);
                                    assert(_automorphism.arr[ind_v1] == ind_vk);
                                    add_automorphism_to_orbit(&orbit, _automorphism.arr, _automorphism_supp.cur_pos,
                                                              _automorphism_supp.arr);
                                }
                                reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                                _automorphism_supp.reset();
                            }

                           if(all_certified) {
                               automorphisms_found += cell_sz + 1;
                               multiply_to_group_size(cell_sz + 1);
                               touched_current_component = true;

                               avg_support_sparse_ir_num2   += 1;
                               const int test_col = select_color_component(g, &c1, quotient_component_start_pos, quotient_component_end_pos);
                               if(test_col == -1) {
                                   avg_reached_end_of_component += 1;
                                   touched_current_component = false;
                               }
                           }
                        }

                        certify = all_certified;

                        if (certify) {
                            int f = 0;
                            // TODO i think these operations are the bottleneck
                            I1.acc = 0;
                            assert(touched_color.get(c3.vertex_to_col[ind_v1]));
                            const int init_c3 = R1.individualize_vertex(&c3, ind_v1);
                            assert(touched_color.get(init_c3));
                            const int touched_col_prev = touched_color_list.cur_pos;
                            R1.refine_coloring(g, &c3, &I1, init_c3, nullptr, -1, -1,
                                               nullptr, &touched_color, &touched_color_list);
                            assert(touched_col_prev == touched_color_list.cur_pos);

                            assert(touched_color_list.cur_pos <= quotient_component_end_pos_v - quotient_component_start_pos_v);
                            if(g->v_size == c1.cells || touched_color_list.cur_pos >= (0.5*(quotient_component_end_pos_v - quotient_component_start_pos_v))) {
                                for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                                    const int v = quotient_component_worklist_v[_i];
                                    colmap[v] = c3.vertex_to_col[v];

                                    c1.vertex_to_col[v] = c3.vertex_to_col[v];
                                    c1.vertex_to_lab[v] = c3.vertex_to_lab[v];
                                    const int lab_pos = c3.vertex_to_lab[v];
                                    c1.lab[lab_pos] = v;
                                    c1.ptn[lab_pos] = c3.ptn[lab_pos];
                                    c1.cells = c3.cells;

                                    c2.vertex_to_col[v] = c3.vertex_to_col[v];
                                    c2.vertex_to_lab[v] = c3.vertex_to_lab[v];
                                    c2.lab[lab_pos] = v;
                                    c2.ptn[lab_pos] = c3.ptn[lab_pos];
                                    c2.cells = c3.cells;
                                }
                            } else {
                                for (int _i = 0; _i < touched_color_list.cur_pos; ++_i) {
                                    const int reset_col = touched_color_list.arr[_i];
                                    if(c3.vertex_to_col[c3.lab[reset_col]] == reset_col ) {
                                        for(int k = 0; k < c3.ptn[reset_col] + 1; ++k) {
                                            const int v = c3.lab[reset_col + k];
                                            colmap[v] = c3.vertex_to_col[v];

                                            c1.vertex_to_col[v] = c3.vertex_to_col[v];
                                            c1.vertex_to_lab[v] = c3.vertex_to_lab[v];
                                            const int lab_pos = c3.vertex_to_lab[v];
                                            assert(lab_pos == reset_col + k);
                                            c1.lab[lab_pos] = v;
                                            c1.ptn[lab_pos] = c3.ptn[lab_pos];

                                            c2.vertex_to_col[v] = c3.vertex_to_col[v];
                                            c2.vertex_to_lab[v] = c3.vertex_to_lab[v];
                                            c2.lab[lab_pos] = v;
                                            c2.ptn[lab_pos] = c3.ptn[lab_pos];
                                        }
                                        c1.cells = c3.cells;
                                        c2.cells = c3.cells;
                                    }
                                }
                            }


                            for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                                const int v = quotient_component_worklist_v[_i];
                                const int list_index = _i - quotient_component_start_pos_v;
                                assert(c1.ptn[c1.vertex_to_col[v]] >= 0);
                                assert(c1.lab[c1.vertex_to_lab[v]] == v);
                                assert(c1.vertex_to_col[c1.lab[c1.vertex_to_col[v]]] == c1.vertex_to_col[v]);

                                assert(c2.ptn[c2.vertex_to_col[v]] >= 0);
                                assert(c2.lab[c2.vertex_to_lab[v]] == v);
                                assert(c2.vertex_to_col[c2.lab[c2.vertex_to_col[v]]] == c2.vertex_to_col[v]);

                                assert(c1.vertex_to_col[v] == c3.vertex_to_col[v]);
                                assert(c2.vertex_to_col[v] == c3.vertex_to_col[v]);
                                assert(c1.ptn[c1.vertex_to_col[v]] == c3.ptn[c1.vertex_to_col[v]]);
                                assert(c2.ptn[c2.vertex_to_col[v]] == c3.ptn[c1.vertex_to_col[v]]);
                            }
                        } else {
                            std::cout << "did not certifyY" << (c1.cells == g->v_size) << ", " << cell_sz << std::endl;
                            int f = 0;
                            for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) { // TODO: why is this necessary?
                                const int v = quotient_component_worklist_v[_i];

                                // reset c1 / c2 such that paths can be restarted in next outer loop iteration
                                c1.vertex_to_col[v] = c3.vertex_to_col[v];
                                c1.vertex_to_lab[v] = c3.vertex_to_lab[v];
                                const int lab_pos = c3.vertex_to_lab[v];
                                c1.lab[lab_pos] = v;
                                c1.ptn[lab_pos] = c3.ptn[lab_pos];
                                c1.cells = c3.cells;

                                c2.vertex_to_col[v] = c3.vertex_to_col[v];
                                c2.vertex_to_lab[v] = c3.vertex_to_lab[v];
                                c2.lab[lab_pos] = v;
                                c2.ptn[lab_pos] = c3.ptn[lab_pos];
                                c2.cells = c3.cells;
                            }
                            // we did not find an automorphism
                            // c1/c2 are now 'broken' for this component, i.e., we can not find automorphisms here
                            // so we just never touch it again
                            //std::cout << "touched_current_component: " << touched_current_component << std::endl;
                            touched_current_component = false;
                        }

                        break;
                    }
                    reset_automorphism(_automorphism.arr, _automorphism_supp.cur_pos, _automorphism_supp.arr);
                    automorphism_supp.reset();

                    if (certify) {
                        if (c1.cells == g->v_size) {
                            //std::cout << "probably broken code 2" << std::endl;
                            for (int _i = quotient_component_start_pos_v; _i < quotient_component_end_pos_v; ++_i) {
                                const int i = quotient_component_worklist_v[_i];
                                assert(i >= 0);
                                assert(i < g->v_size);
                                colmap[i] = c1.vertex_to_col[i];
                            }
                        } else {
                            for (int j = 0; j < touched_color_list.cur_pos; ++j) {
                                const int _c = touched_color_list.arr[j];
                                int f = 0;
                                while (f < c1.ptn[_c] + 1) {
                                    assert(c1.vertex_to_col[c1.lab[_c + f]] == _c);
                                    colmap[c1.lab[_c + f]] = c1.vertex_to_col[c1.lab[_c + f]]; // should be c1 or c2?
                                    ++f;
                                }
                            }
                        }
                    }
                }

                if (touched_current_component) {
                    quotient_component_touched.push_back(component);
                    del_discrete_edges_inplace_component(g, &c1, quotient_component_start_pos);
                }
                quotient_component_start_pos = quotient_component_worklist_boundary[component].first;
                quotient_component_start_pos_v = quotient_component_worklist_boundary[component].second;

                ++component;

                quotient_component_end_pos = quotient_component_worklist_boundary[component].first;
                quotient_component_end_pos_v = quotient_component_worklist_boundary[component].second;
            }
            global_automorphisms_found += automorphisms_found;
            if(automorphisms_found == 0)
                break;
        }
        //for(int i = 0; i < quotient_component_touched.size(); ++i)
        //    std::cout << quotient_component_touched[i] << " ";
        //std::cout << std::endl;
        if(avg_support_sparse_ir > 0 && avg_support_sparse_ir_num > 0) {
            avg_support_sparse_ir = avg_support_sparse_ir / avg_support_sparse_ir_num;
        }
        std::cout << avg_support_sparse_ir << "/" << g->v_size << std::endl;
        std::cout << (1.0*avg_reached_end_of_component) / (avg_support_sparse_ir_num2*1.0) << std::endl;
        if(avg_support_sparse_ir_num2 > 0) {
            avg_end_of_comp = (1.0 * avg_reached_end_of_component) / (avg_support_sparse_ir_num2 * 1.0);
        } else {
            avg_end_of_comp = 0;
        }
        return global_automorphisms_found;
    }


    void compute_quotient_graph_components(sgraph* g, int* colmap, dejavu_consumer consume) {
        coloring<int> c1;
        g->initialize_coloring(&c1, colmap);

        if(!init_quotient_arrays) {
            seen_vertex.initialize(g->v_size);
            seen_color.initialize(g->v_size);

            worklist.reserve(g->v_size);

            quotient_component_worklist_v.reserve(g->v_size + 32);
            quotient_component_worklist_col.reserve(g->v_size + 32);
            quotient_component_worklist_col_sz.reserve(g->v_size + 32);
            init_quotient_arrays = true;
        }

        seen_vertex.reset();
        seen_color.reset();

        std::vector<int> worklist;
        worklist.reserve(g->v_size);

        quotient_component_worklist_v.reserve(g->v_size + 32);
        quotient_component_worklist_col.reserve(g->v_size + 32);
        quotient_component_worklist_col_sz.reserve(g->v_size + 32);

        quotient_component_worklist_v.clear();
        quotient_component_worklist_col.clear();
        quotient_component_worklist_col_sz.clear();

        int component = 0;
        int current_component_sz = 0;
        for(int vs = 0; vs < g->v_size; ++vs) {
            if(seen_vertex.get(vs))
                continue;
            worklist.push_back(vs);

            while(!worklist.empty()) {
                const int next_v = worklist.back();
                worklist.pop_back();
                if(seen_vertex.get(next_v))
                    continue;
                current_component_sz += 1;
                seen_vertex.set(next_v);
                quotient_component_worklist_v.push_back(next_v);
                for(int i = 0; i < g->d[next_v]; ++i) {
                    if(!seen_vertex.get(g->e[g->v[next_v] + i]))
                        worklist.push_back(g->e[g->v[next_v] + i]); // neighbours
                }
                const int col = c1.vertex_to_col[next_v];
                assert(next_v < g->v_size);
                assert(col < g->v_size);
                if(!seen_color.get(col)) {
                    quotient_component_worklist_col.push_back(col);
                    quotient_component_worklist_col_sz.push_back(c1.ptn[col]);
                    seen_color.set(col);
                    for(int i = 0; i < c1.ptn[col] + 1; ++i) {
                        assert(col + i < g->v_size);
                        assert(c1.vertex_to_col[c1.lab[col + i]] == c1.vertex_to_col[next_v]);
                        if(c1.lab[col + i] != next_v) {
                            if(!seen_vertex.get(c1.lab[col + i]))
                                worklist.push_back(c1.lab[col + i]);
                        }
                    }
                }
            }

            //std::cout << "current_component_sz: " << current_component_sz << std::endl;

            quotient_component_worklist_col_sz.push_back(-1);
            quotient_component_worklist_col.push_back(-1); // border of component
            quotient_component_worklist_v.push_back(-1); // border of component
            std::cout << current_component_sz << " ";
            current_component_sz = 0;
            ++component;
        }
        std::cout << "[" << component  << "]" << std::endl;
    }

    void del_discrete_edges_inplace(sgraph* g, coloring<int>*c) {
        int rem_edges = 0;
        int discrete_vert = 0;
        del.reset();
        for(int i = 0; i < c->lab_sz;) {
            const int col_sz = c->ptn[i];
            if(col_sz == 0) {
                ++discrete_vert;
                del.set(c->lab[i]);
            }
            i += col_sz + 1;
        }

        for(int v = 0; v < g->v_size; ++v) {
            if(del.get(v)) {
                rem_edges += g->d[v];
                g->d[v] = 0;
                continue;
            }
            int write_pt_front = g->v[v];
            int write_pt_back  = g->v[v] + g->d[v] - 1;
            for(int n = g->v[v]; n < g->v[v] + g->d[v];) {
                const int neigh = g->e[n];
                if(del.get(neigh)) {
                    const int swap_neigh = g->e[g->v[v] + g->d[v] - 1];
                    g->e[g->v[v] + g->d[v] - 1] = neigh;
                    g->e[n] = swap_neigh;
                    --g->d[v];
                    ++rem_edges;
                } else {
                    ++n;
                }
            }
        }
        //std::cout << discrete_vert << ", " << rem_edges << std::endl;
    }

    void del_discrete_edges_inplace_component(sgraph* g, coloring<int>*c, int component_start_pos) {
        int rem_edges = 0;
        int discrete_vert = 0;
        for (int _i = component_start_pos; _i < quotient_component_worklist_col.size(); ++_i) {
            const int col = quotient_component_worklist_col[_i];
            const int col_sz = quotient_component_worklist_col_sz[_i];

            if (col == -1) {// reached end of component
                break;
            }

            for (int i = col; i < col + col_sz + 1; ++i) {
                if (c->vertex_to_col[c->lab[i]] != i)
                    continue; // not a color
                if (c->ptn[i] == 0) {
                    ++discrete_vert;
                    del.set(c->lab[i]);
                }
            }
        }

        for (int _i = component_start_pos; _i < quotient_component_worklist_col.size(); ++_i) {
            const int col = quotient_component_worklist_col[_i];
            const int col_sz = quotient_component_worklist_col_sz[_i];

            if (col == -1) {// reached end of component
                break;
            }

            for (int i = col; i < col + col_sz + 1; ++i) {
                const int v = c->lab[i];
                if (del.get(v)) {
                    rem_edges += g->d[v];
                    g->d[v] = 0;
                    continue;
                }
                int write_pt_front = g->v[v];
                int write_pt_back = g->v[v] + g->d[v] - 1;
                for (int n = g->v[v]; n < g->v[v] + g->d[v];) {
                    const int neigh = g->e[n];
                    if (del.get(neigh)) {
                        const int swap_neigh = g->e[g->v[v] + g->d[v] - 1];
                        g->e[g->v[v] + g->d[v] - 1] = neigh;
                        g->e[n] = swap_neigh;
                        --g->d[v];
                        ++rem_edges;
                    } else {
                        ++n;
                    }
                }
            }
        }
    }

    void reduce(sgraph* g, int* colmap, dejavu_consumer consume) {
        {
            /*int deg0 = 0;
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
            std::cout << "(pre-red) before reduction (0, 1, 2) " << deg0 << ", " << deg1 << ", " << deg2 << std::endl;*/
            std::cout << "(pre-red) before reduction (G, E) " << g->v_size << ", " << g->e_size << std::endl;
        }
        g->initialize_coloring(&c, colmap);
        R1.refine_coloring_first(g, &c, -1);

        if(c.cells == g->v_size) {
            g->v_size = 0;
            g->e_size = 0;
            g->d_size = 0;
            return;
        }

        automorphism.initialize(g->v_size);
        for(int i = 0; i < g->v_size; ++i) {
            automorphism.push_back(i);
        }
        automorphism_supp.initialize(g->v_size);
        worklist_deg0.initialize(g->v_size);
        worklist_deg1.initialize(g->v_size);
        edge_scratch.initialize(g->e_size);

        translate_layer_fwd.reserve(g->v_size);

        backward_translation_layers.emplace_back(std::vector<int>());
        const int back_ind = backward_translation_layers.size()-1;
        translation_layers.emplace_back(std::vector<int>());
        const int fwd_ind = translation_layers.size()-1;
        backward_translation_layers[back_ind].reserve(g->v_size);
        for(int i = 0; i < g->v_size; ++i)
            backward_translation_layers[back_ind].push_back(i);

        domain_size = g->v_size;

        // assumes colmap is array of length g->v_size
        del = mark_set();
        del.initialize(g->v_size);

        canonical_recovery_string.reserve(g->v_size);
        for(int i = 0; i < domain_size; ++i)
            canonical_recovery_string.emplace_back(std::vector<int>());

        // eliminate degree 1 + 0
        // TODO delete edges to discrete vertices in-place here?
        del_discrete_edges_inplace(g, &c);
        red_deg10_assume_cref(g, &c, consume);
        copy_coloring_to_colmap(&c, colmap);
        //mark_discrete_for_deletion(g, colmap, &rec);
        perform_del(g, colmap);


        if(g->v_size == 0)
            return;

        del_e = mark_set();
        del_e.initialize(g->e_size);
        add_edge_buff.initialize(g->v_size);
        add_edge_buff_act.initialize(g->v_size);

        red_quotient_components(g, colmap, consume);
        perform_del_edge(g, colmap);

        red_deg2_path_size_1(g, colmap, consume);
        perform_del2(g, colmap);

        const int auto_found = sparse_ir_col_sz_2_quotient_components(g, colmap, consume, 16);
        if(auto_found > 0) {
            perform_del_discrete(g, colmap);
        }
        std::cout << "sparse_ir_qco_2 " << auto_found << std::endl;

        sparse_ir(g, colmap, consume, SELECTOR_LARGEST);
        perform_del_discrete(g, colmap);

        red_deg2_path_size_1(g, colmap, consume);
        perform_del2(g, colmap);

        // TODO: look at average support of automorphisms to decide whether its worth?
        // TODO: active support vs. "implied" support - amalgam of number of individualizations and support is crucial factor
        // TODO: path length vs. base?
        // TODO: If I'm just badly emulating the IR solver without proper Schreier-Sims, should immediately stop
        const int test = sparse_ir_quotient_components(g, colmap, consume, 16, 4);
        std::cout << "sparse_ir_qco_X " << test << std::endl;
        if(test > 0) {
            perform_del_discrete(g, colmap);
        }

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

            if(deg0 == 0 && deg1 == 0 && deg2 != 0) {
                red_deg2_assume_cref(g, colmap, consume);
                perform_del2(g, colmap);
            }

            std::cout << "(pre-red) after reduction (0, 1, 2) " << deg0 << ", " << deg1 << ", " << deg2 << std::endl;
            std::cout << "(pre-red) after reduction (G, E) " << g->v_size << ", " << g->e_size << std::endl;

            int back_off = 1;
            int add_on   = 0;

            while(deg0 > 0 || deg1 > 0) {
                if(avg_end_of_comp > 0.5 && (avg_support_sparse_ir * 1.0 / g->v_size * 1.0) > 0.1) {
                    back_off = back_off * 16; // TODO this is not "fine" enough
                }
                if(avg_end_of_comp < 0.1 || (avg_support_sparse_ir * 1.0 / g->v_size * 1.0) < 0.01) {
                    add_on += 4;
                }
                // refinement
                std::cout << "(prep-red) loop [collected grp_sz " << base << "*10^" << exp << "] [bo, ao ] "<< back_off << ", " << add_on << std::endl;
                coloring<int> c2;
                //refinement<int, int, int> R2;
                g->initialize_coloring(&c2, colmap);
                //R2.refine_coloring_first(g, &c2, -1);

                // eliminate degree 1 + 0
                red_deg10_assume_cref(g, &c2, consume);
                copy_coloring_to_colmap(&c2, colmap);
                mark_discrete_for_deletion(g, colmap);
                perform_del(g, colmap);
                //perform_del_discrete(g, colmap, &rec);
                red_deg2_assume_cref(g, colmap, consume);
                perform_del2(g, colmap);

                //red_quotient_components(g, colmap, consume);
                //perform_del_edge(g, colmap);

                red_deg2_path_size_1(g, colmap, consume);
                perform_del2(g, colmap);

                const int test2 = sparse_ir_col_sz_2_quotient_components(g, colmap, consume, 16);
                std::cout << "sparse_ir_qco_2 " << test2 << std::endl;
                perform_del_discrete(g, colmap);

                sparse_ir(g, colmap, consume, SELECTOR_LARGEST);
                perform_del_discrete(g, colmap);

                red_deg2_path_size_1(g, colmap, consume);
                perform_del2(g, colmap);

                const int test3 = sparse_ir_quotient_components(g, colmap, consume, 16, 5 - back_off + add_on);
                std::cout << "sparse_ir_qco_X " << test3 << std::endl;
                if(test3 > 0) {
                    perform_del_discrete(g, colmap);
                }
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
                //break;
            }
        }

        std::cout << "(prep-red) translation layers: " << backward_translation_layers.size() << std::endl;

         //std::cout << "(pre-red) group size: " << base << "*10^" << exp << std::endl;

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
        //group.domain_size = domain_size;
        //group.initialize_from_automorphism_info(a, &rec);
    }
};

#endif //DEJAVU_PREP_H
