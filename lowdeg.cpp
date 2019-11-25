// currently not used

#include <iostream>
#include <cmath>
#include "lowdeg.h"
#include "refinement.h"
#include "group_shared.h"


void split_with_respect_to_array(coloring* c, int col, int col_sz, refinement* R, work_set_int* A) {
    int i, v, acc, val, _col, h, v_new_color;
    R->degrees_worklist.reset();
    R->neighbour_sizes.reset();
    // collect degrees...
    int first_index = A->get(c->lab[col]) + 1;
    R->degrees_worklist.push_back(first_index);
    int total = 0;
    for(i = col; i < col + col_sz; ++i) {
        v = c->lab[i];
        int index = A->get(v) + 1;
        if(index == first_index)
            continue;
        total += 1;
        if(R->neighbour_sizes.inc(index) == 0)
            R->degrees_worklist.push_back(index);
    }
    R->neighbour_sizes.inc(first_index);
    R->neighbour_sizes.set(first_index, col_sz - total - 1);


    // if we need to split...
    if(R->degrees_worklist.cur_pos != 1) {
        // do accumulative counting...
        R->degrees_worklist.sort();
        // enrich neighbour_sizes to accumulative counting array
        acc = 0;
        while(!R->degrees_worklist.empty()) {
            i   = R->degrees_worklist.pop_back();
            val = R->neighbour_sizes.get(i) + 1;
            if(val >= 1) {
                R->neighbour_sizes.set(i, val + acc);
                acc += val;
                _col = col + col_sz - (R->neighbour_sizes.get(i)); // use new color here to remove color_workset?
                c->ptn[_col] = -1; // this is val - 1, actually...
            }
        }

        memcpy(R->vertex_worklist.arr, c->lab + col, col_sz * sizeof(int));
        R->vertex_worklist.cur_pos = col_sz;

        // determine colors and rearrange
        // split color classes according to count in counting array
        while(!R->vertex_worklist.empty()) {
            v = R->vertex_worklist.pop_back();
            v_new_color = col + col_sz - (R->neighbour_sizes.get(A->get(v) + 1));

            assert(v_new_color >= col);
            assert(v_new_color < col + col_sz);

            c->lab[v_new_color + c->ptn[v_new_color] + 1] = v;
            c->vertex_to_col[v] = v_new_color;
            c->vertex_to_lab[v] = v_new_color + c->ptn[v_new_color] + 1;
            c->ptn[v_new_color] += 1;
        }

        for(h = col; h < col + col_sz;) {
            assert(h >= 0 && h < c->ptn_sz);
            assert(c->ptn[h] + 1 > 0);
            if(h != 0)
                c->ptn[h - 1] = 0;
            h += c->ptn[h] + 1;
        }
    }
}

std::pair<sgraph*, coloring*> lowdeg::preprocess(coloring* c, sgraph* g, refinement* R) {
    this->g = g;
    this->c = new coloring;
    this->c->copy_force(c);

    domain_size = c->lab_sz;
    R->assure_initialized(g);

    sgraph* reduced_graph = new sgraph;
    reduced_graph->v = new int[g->v_size];
    reduced_graph->e = new int[g->e_size];
    reduced_graph->d = new int[g->d_size];
    reduced_graph->v_size = g->v_size;
    reduced_graph->e_size = g->e_size;
    reduced_graph->d_size = g->d_size;
    reduced_graph->max_degree = g->max_degree;

    std::cout << "dejavu_lowdeg_pre-------------------------------------------------" << std::endl;

    int deg1_vert = 0;
    int deg2_vert = 0;
    int pairs = 0;
    int paths = 0;

    //coloring* c = new coloring;
    coloring* reduced_c = new coloring;
    //g->initialize_coloring(c);

    int v, vpos, epos, reduce_vpos;


    // ToDo: need to do deg1_counter "with respect to" every degree 1 class, then its correct...
    // ToDo: keep separate degree1 counter for every degree1 class?

    // ToDo: firs find degree1 classes? then make counters,
    // ToDo: map colors -> counter_id++ while discovering -> counters

    // discover degree1 classes
    bool seen_deg_1 = false;
    for(int i = 0; i < c->lab_sz;) {
        if(g->d[c->lab[i]] == 1) {
            seen_deg_1 = true;
            break;
        }
        i += c->ptn[i] + 1;
    }
    if(!seen_deg_1) {
        abort = true;
        return std::pair<sgraph*, coloring*>(g, c);
    }


    deg1_counter.initialize(g->v_size);
    g_to_reduced_g.initialize(g->v_size);
    reduced_g_to_g.initialize(g->v_size);
    // work_set deg2_marker_visited;
    work_set_int deg2_path_length;
    //deg2_marker_visited.initialize(g->v_size);
    //deg2_path_length.initialize(g->v_size);

    work_set color_workset;
    color_workset.initialize(g->v_size);


    reduce_vpos = 0;
    for(int i = 0; i < c->lab_sz; ++i) {
        v = c->lab[i];
        if(g->d[v] != 1) {
            g_to_reduced_g.set(v, reduce_vpos);
            reduced_g_to_g.set(reduce_vpos++, v);
            assert(reduced_g_to_g.get(g_to_reduced_g.get(v)) == v);
            continue;
        }
        if(g->d[v] == 1) {
            deg1_counter.inc_nr(g->e[g->v[v]]); // ToDo: collect more information here, worklists?
            //counters[deg1_to_counter.get(c->vertex_to_col[v])].inc_nr(g->e[g->v[v]]);
            // ToDo: also for enriching generators in the end
            color_workset.set(c->vertex_to_col[g->e[g->v[v]]]);
            deg1_vert += 1;
            if(deg1_counter.get(v) == 0) { // detect pairs // ToDo: use specific counter for current color class...
                pairs += 1;
            }
        }
    }

    reduced_domain_size = reduce_vpos;
    std::cout << "reduce_vpos: " << reduce_vpos << std::endl;
    std::cout << "deg1_vert: " << deg1_vert << std::endl;
    std::cout << "deg2_vert: " << deg2_vert << std::endl;
    std::cout << "paths: " << paths << std::endl;
    std::cout << "pairs: "     << pairs << std::endl; // properly handle pairs
    std::cout << "cycle: " << std::endl; // properly handle cycles

    if(pairs > 0) {
        abort = true;
        return std::pair<sgraph*, coloring*>(g, c);
    }

    std::cout << "Writing graph..." << std::endl;
    vpos = 0;
    epos = 0;

    int con_v, map_v, col, col_sz, h, j, acc, val, _col, v_new_color;

    reduced_c->initialize(reduce_vpos);
    int labpos = 0;

    for(int i = 0; i < c->lab_sz;) {
        v = c->lab[i];
        if(g->d[v] == 1) {
            //i += c->ptn[i] + 1;
            i++;
            continue;
        }

        col = i;
        col_sz = c->ptn[i] + 1;

        if(!color_workset.get(col)) {
            // std::cout << "untouched" << std::endl;
            // copy color to reduced graph
            int first_labpos = labpos;
            for (h = col; h < col + col_sz; ++h) {
                v = c->lab[h];
                map_v = g_to_reduced_g.get(v);
                assert(map_v >= 0);
                assert(map_v < reduce_vpos);

                reduced_c->lab[labpos]          = map_v;
                reduced_c->ptn[labpos]          = c->ptn[h];
                reduced_c->vertex_to_lab[map_v] = labpos++;
                assert(first_labpos >= 0 && first_labpos < reduce_vpos);
                reduced_c->vertex_to_col[map_v] = first_labpos; // last seen color

                reduced_graph->v[map_v] = epos;
                int d = 0;
                for (j = g->v[v]; j < g->v[v] + g->d[v]; ++j) {
                    con_v = g->e[j];
                    if (g->d[con_v] == 1)
                        continue;
                    reduced_graph->e[epos++] = g_to_reduced_g.get(con_v);
                    assert(g_to_reduced_g.get(con_v) >= 0);
                    ++d;
                }

                reduced_graph->d[map_v] = d;
            }
        } else {
            // ToDo: do this with "all counters"
            split_with_respect_to_array(c, col, col_sz, R, &deg1_counter);

            bool is_col = true;
            // carbon copy cell and coloring now...
            int first_labpos = labpos;
            for (h = col; h < col + col_sz; ++h) {
                if(is_col) {
                    if(labpos > 0)
                        reduced_c->ptn[labpos - 1] = 0;
                    _col = labpos;
                }
                is_col = (c->ptn[h] == 0);

                v = c->lab[h];
                map_v = g_to_reduced_g.get(v);

                // copy lab and ptn
                reduced_c->lab[labpos]          = map_v;
                reduced_c->ptn[labpos]          = c->ptn[h];
                reduced_c->vertex_to_lab[map_v] = labpos++;
                assert(reduced_c->lab[reduced_c->vertex_to_lab[map_v]] == map_v);
                assert(_col >= 0 && _col < reduce_vpos);
                reduced_c->vertex_to_col[map_v] = _col; // last seen color
                assert(_col >= first_labpos && _col <= first_labpos + col_sz);

                assert(map_v >= 0);
                assert(map_v < reduce_vpos);

                reduced_graph->v[map_v] = epos;
                int d = 0;
                for (j = g->v[v]; j < g->v[v] + g->d[v]; ++j) {
                    con_v = g->e[j];
                    if (g->d[con_v] == 1)
                        continue;
                    reduced_graph->e[epos++] = g_to_reduced_g.get(con_v);
                    assert(g_to_reduced_g.get(con_v) >= 0);
                    ++d;
                }

                assert(reduced_graph->d[map_v] <= g->d[v]);
                reduced_graph->d[map_v] = d;
            }
        }
        i = h;
    }

    color_workset.reset();

    assert(labpos == reduce_vpos);

    for(int i = 0; i < labpos;++i) {
        assert(reduced_c->lab[i] >= 0 && reduced_c->lab[i] < reduce_vpos);
        assert(reduced_c->lab[reduced_c->vertex_to_lab[i]] == i);
    }

    // assert proper ptn
    for(int i = 0; i < labpos;) {
        if(i != 0)
            assert(reduced_c->ptn[i - 1] == 0);
        i += reduced_c->ptn[i] + 1;
    }

    reduced_graph->v_size = reduce_vpos;
    reduced_graph->e_size = epos;
    reduced_graph->d_size = reduce_vpos;

    // assert graph stuff
    for(int i = 0; i < reduced_graph->v_size; ++i) {
        assert(reduced_graph->v[i] <= reduced_graph->e_size);
        assert(reduced_graph->d[i] >= 0 && reduced_graph->d[i] <= reduced_graph->max_degree);
        for(int j = reduced_graph->v[i]; j < reduced_graph->v[i] + reduced_graph->d[i]; ++j)
            assert(reduced_graph->e[j] >= 0 && reduced_graph->e[j] < reduced_graph->v_size);
    }

    std::cout << "Reduced (" << g->v_size << ", " << g->e_size << ") -> (" << reduced_graph->v_size << ", " << reduced_graph->e_size << ")" << std::endl;


    std::cout << "------------------------------------------------------------------" << std::endl;


    if(g->v_size == reduced_graph->v_size) {
        return std::pair<sgraph *, coloring *>(reduced_graph, reduced_c);
    } else {
        // recursive...
        //R->assure_initialized(g);
        R->refine_coloring_first(reduced_graph, reduced_c, -1); // check if discrete and abort?
        int i;
        for(i = 0; i < reduced_domain_size; ++i)
            if(reduced_c->ptn[i] != 0) break;
        if(i == reduced_domain_size) {
            std::cout << "Reduced discrete..." << std::endl;
            return std::pair<sgraph *, coloring *>(reduced_graph, reduced_c);
        }
        nest = new lowdeg;
        return nest->preprocess(reduced_c, reduced_graph, R);
    }
}

std::pair<sgraph*, coloring*> lowdeg::preprocess2(coloring* c, sgraph* g, refinement* R) {
    domain_size = c->lab_sz;

    sgraph* reduced_graph = new sgraph;
    reduced_graph->v = new int[g->v_size];
    reduced_graph->e = new int[g->e_size];
    reduced_graph->d = new int[g->d_size];
    reduced_graph->v_size = g->v_size;
    reduced_graph->e_size = g->e_size;
    reduced_graph->d_size = g->d_size;
    reduced_graph->max_degree = g->max_degree;

    std::cout << "dejavu_lowdeg_pre-------------------------------------------------" << std::endl;
    deg1_counter.initialize(g->v_size);
    g_to_reduced_g.initialize(g->v_size);
    reduced_g_to_g.initialize(g->v_size);
    work_set deg2_marker_visited;
    work_set_int deg2_path_length;
    deg2_marker_visited.initialize(g->v_size);
    deg2_path_length.initialize(g->v_size);

    work_set color_workset;
    color_workset.initialize(g->v_size);
    color_workset.reset();

    int deg1_vert = 0;
    int deg2_vert = 0;
    int pairs = 0;
    int paths = 0;

    //coloring* c = new coloring;
    coloring* reduced_c = new coloring;
    //g->initialize_coloring(c);

    int v, vpos, epos, reduce_vpos;


    // ToDo: need to do deg1_counter "with respect to" every degree 1 class, then its correct...
    // ToDo: keep separate degree1 counter for every degree1 class?

    reduce_vpos = 0;
    for(int i = 0; i < c->lab_sz; ++i) {
        v = c->lab[i];
        if(g->d[v] != 1 && g->d[v] != 2) {
            g_to_reduced_g.set(v, reduce_vpos);
            reduced_g_to_g.set(reduce_vpos++, v);
            assert(reduced_g_to_g.get(g_to_reduced_g.get(v)) == v);
            continue;
        }
        if(g->d[v] == 1) {
            deg1_counter.inc_nr(g->e[g->v[v]]); // ToDo: collect more information here, worklists?
                                                // ToDo: also for enriching generators in the end
            color_workset.set(c->vertex_to_col[g->e[g->v[v]]]);
            deg1_vert += 1;
            if(deg1_counter.get(v) == 0) { // detect pairs
                pairs += 1;
            }
        }
        if(g->d[v] == 2) {
            deg2_vert += 1;
            if(deg2_marker_visited.get(v))
                continue;
            deg2_marker_visited.set(v);
            g_to_reduced_g.set(v, reduce_vpos);
            reduced_g_to_g.set(reduce_vpos++, v);
            assert(reduced_g_to_g.get(g_to_reduced_g.get(v)) == v);
            int end1_v = -2, end2_v = -2, length = 0;
            int last1_v, last2_v, temp;
            // explore path in both directions
            // v then connects both ends and carries as weight the length of the path

            // a path ends if...
            // 1) d > 2
            // 2) deg1 neighbours != 0 (= 1 if 1) does not hold)
            // 3) deg2_marker already set (cycle!)

            // search for end1_v
            if(end1_v != v) {
                end1_v = g->e[g->v[v]];
                last1_v = v;
                while(true) {
                    if(g->d[end1_v] != 2) {
                        break;
                    }
                    deg2_marker_visited.set(end1_v);
                    // advance v
                    temp = (g->e[g->v[end1_v]] == last1_v)?(g->e[g->v[end1_v] + 1]):g->e[g->v[end1_v]];
                    last1_v = end1_v;
                    end1_v = temp;
                    ++length;
                }
            }

            // search for end2_v
            if(end2_v != v) {
                end2_v = g->e[g->v[v] + 1];
                last2_v = v;
                while(true) {
                    if(g->d[end2_v] != 2) {
                        break;
                    }
                    deg2_marker_visited.set(end2_v);
                    // advance v
                    temp = (g->e[g->v[end2_v]] == last2_v)?(g->e[g->v[end2_v] + 1]):g->e[g->v[end2_v]];
                    last2_v = end2_v;
                    end2_v = temp;
                    ++length;
                }
            }

            deg2_path_length.set(v, length); // to split deg2 class later...
            ++paths;
            //std::cout << "path from " << end1_v << "(" << g->d[end1_v] << ")" << " to " << end2_v << "(" << g->d[end2_v] << ")" << " l: " << length << std::endl;

            // connect end1_v to end2_v via v
            if(length > 0) {
                g->e[g->v[v]]     = end1_v;
                g->e[g->v[v] + 1] = end2_v;

                for(int j = g->v[end1_v]; j < g->v[end1_v] + g->d[end1_v]; ++j) {
                    if(g->e[j] == last1_v) {
                        g->e[j] = v;
                        break;
                    }
                }

                for(int j = g->v[end2_v]; j < g->v[end2_v] + g->d[end2_v]; ++j) {
                    if(g->e[j] == last2_v) {
                        g->e[j] = v;
                        break;
                    }
                }
                // save whether v has 0, 1, 2 deg1 neighbours and use later to distinguish more paths...
                deg1_counter.set(v, (g->d[end1_v] == 1) + (g->d[end2_v] == 1) - 1);
            }
        }
    }

    reduced_domain_size = reduce_vpos;
    std::cout << "reduce_vpos: " << reduce_vpos << std::endl;
    std::cout << "deg1_vert: " << deg1_vert << std::endl;
    std::cout << "deg2_vert: " << deg2_vert << std::endl;
    std::cout << "paths: " << paths << std::endl;
    std::cout << "pairs: "     << pairs << std::endl; // properly handle pairs
    std::cout << "cycle: " << std::endl; // properly handle cycles

    std::cout << "Writing graph..." << std::endl;
    vpos = 0;
    epos = 0;

    int con_v, map_v, col, col_sz, h, j, acc, val, _col, v_new_color;

    reduced_c->initialize(reduce_vpos);
    int labpos = 0;

    for(int i = 0; i < c->lab_sz;) {
        v = c->lab[i];
        if(g->d[v] == 1) {
            //i += c->ptn[i] + 1;
            i++;
            continue;
        }

        col = i;
        col_sz = c->ptn[i] + 1;

        if(!R->color_workset.get(col) && g->d[c->lab[col]] != 2) {
           // std::cout << "untouched" << std::endl;
            // copy color to reduced graph
            int first_labpos = labpos;
            for (h = col; h < col + col_sz; ++h) {
                v = c->lab[h];
                map_v = g_to_reduced_g.get(v);
                assert(map_v >= 0);
                assert(map_v < reduce_vpos);

                reduced_c->lab[labpos]          = map_v;
                reduced_c->ptn[labpos]          = c->ptn[h];
                reduced_c->vertex_to_lab[map_v] = labpos++;
                assert(first_labpos >= 0 && first_labpos < reduce_vpos);
                reduced_c->vertex_to_col[map_v] = first_labpos; // last seen color

                reduced_graph->v[map_v] = epos;
                int d = 0;
                for (j = g->v[v]; j < g->v[v] + g->d[v]; ++j) {
                    con_v = g->e[j];
                    if (g->d[con_v] == 1)
                        continue;
                    reduced_graph->e[epos++] = g_to_reduced_g.get(con_v);
                    assert(g_to_reduced_g.get(con_v) >= 0);
                    ++d;
                }

                reduced_graph->d[map_v] = d;
            }
        } else if (g->d[c->lab[col]] != 2) {
            split_with_respect_to_array(c, col, col_sz, R, &deg1_counter);

            bool is_col = true;
            // carbon copy cell and coloring now...
            int first_labpos = labpos;
            for (h = col; h < col + col_sz; ++h) {
                if(is_col) {
                    if(labpos > 0)
                        reduced_c->ptn[labpos - 1] = 0;
                    _col = labpos;
                }
                is_col = (c->ptn[h] == 0);

                v = c->lab[h];
                map_v = g_to_reduced_g.get(v);

                // copy lab and ptn
                reduced_c->lab[labpos]          = map_v;
                reduced_c->ptn[labpos]          = c->ptn[h];
                reduced_c->vertex_to_lab[map_v] = labpos++;
                assert(reduced_c->lab[reduced_c->vertex_to_lab[map_v]] == map_v);
                assert(_col >= 0 && _col < reduce_vpos);
                reduced_c->vertex_to_col[map_v] = _col; // last seen color
                assert(_col >= first_labpos && _col <= first_labpos + col_sz);

                assert(map_v >= 0);
                assert(map_v < reduce_vpos);

                reduced_graph->v[map_v] = epos;
                int d = 0;
                for (j = g->v[v]; j < g->v[v] + g->d[v]; ++j) {
                    con_v = g->e[j];
                    if (g->d[con_v] == 1)
                        continue;
                    reduced_graph->e[epos++] = g_to_reduced_g.get(con_v);
                    assert(g_to_reduced_g.get(con_v) >= 0);
                    ++d;
                }

                assert(reduced_graph->d[map_v] <= g->d[v]);
                reduced_graph->d[map_v] = d;
            }
        } else {
            // split col according to path length and degree 1 counter
            // only write necessary vertices to graph
            int first_labpos = labpos;
            int actual_size = 0;

            split_with_respect_to_array(c, col, col_sz, R, &deg1_counter);
            for(h = col; h < col + col_sz;) {
                int h_size = c->ptn[h] + 1;
                split_with_respect_to_array(c, h, h_size, R, &deg2_path_length);
                h += h_size;
            }
            if(labpos > 0) {
                c->ptn[labpos - 1] = 0;
            }

            //std::cout << "deg2" << std::endl;

            bool is_col = true;
            // carbon copy cell and coloring now...
            for (h = col; h < col + col_sz; ++h) {
                v = c->lab[h];
                if(deg2_path_length.get(v) >= 0) {
                    if (is_col) {
                        if (labpos > 0)
                            reduced_c->ptn[labpos - 1] = 0;
                        _col = labpos;
                    }
                    is_col = (c->ptn[h] == 0);
                    map_v = g_to_reduced_g.get(v);

                    // copy lab and ptn
                    reduced_c->lab[labpos] = map_v;
                    reduced_c->ptn[labpos] = c->ptn[h];
                    reduced_c->vertex_to_lab[map_v] = labpos++;
                    assert(reduced_c->lab[reduced_c->vertex_to_lab[map_v]] == map_v);
                    assert(_col >= 0 && _col < reduce_vpos);
                    reduced_c->vertex_to_col[map_v] = _col; // last seen color
                    assert(_col >= first_labpos && _col <= first_labpos + col_sz);

                    assert(map_v >= 0);
                    assert(map_v < reduce_vpos);

                    reduced_graph->v[map_v] = epos;
                    int d = 0;
                    for (j = g->v[v]; j < g->v[v] + g->d[v]; ++j) {
                        con_v = g->e[j];
                        if (g->d[con_v] == 1)
                            continue;
                        reduced_graph->e[epos++] = g_to_reduced_g.get(con_v);
                        assert(g_to_reduced_g.get(con_v) >= 0);
                        ++d;
                    }
                    reduced_graph->d[map_v] = d;
                    assert(reduced_graph->d[map_v] <= g->d[v]);
                }
            }

            //std::cout << "col_sz: " << col_sz << " actual_sz: " << actual_size << std::endl;
            //reduced_c->ptn[first_labpos] = actual_size - 1;
            //reduced_c->ptn[h - 1] = 0;
        }
        i = h;
    }

    R->degrees_worklist.reset();
    color_workset.reset();
    R->vertex_worklist.reset();
    R->neighbour_sizes.reset();

    /*for(int i = 0; i < labpos; ++i)
        std::cout << reduced_c->lab[i] << " ";
    std::cout << std::endl;
    for(int i = 0; i < labpos; ++i)
        std::cout << reduced_c->ptn[i] << " ";
    std::cout << std::endl;*/

    assert(labpos == reduce_vpos);

    for(int i = 0; i < labpos;++i) {
        assert(reduced_c->lab[i] >= 0 && reduced_c->lab[i] < reduce_vpos);
        assert(reduced_c->lab[reduced_c->vertex_to_lab[i]] == i);
    }

    // assert proper ptn
    for(int i = 0; i < labpos;) {
        if(i != 0)
            assert(reduced_c->ptn[i - 1] == 0);
        i += reduced_c->ptn[i] + 1;
    }

    reduced_graph->v_size = reduce_vpos;
    reduced_graph->e_size = epos;
    reduced_graph->d_size = reduce_vpos;


    // assert graph stuff
    for(int i = 0; i < reduced_graph->v_size; ++i) {
        assert(reduced_graph->v[i] <= reduced_graph->e_size);
        assert(reduced_graph->d[i] >= 0 && reduced_graph->d[i] <= reduced_graph->max_degree);
        for(int j = reduced_graph->v[i]; j < reduced_graph->v[i] + reduced_graph->d[i]; ++j)
            assert(reduced_graph->e[j] >= 0 && reduced_graph->e[j] < reduced_graph->v_size);
    }

    std::cout << "Reduced (" << g->v_size << ", " << g->e_size << ") -> (" << reduced_graph->v_size << ", " << reduced_graph->e_size << ")" << std::endl;


    std::cout << "------------------------------------------------------------------" << std::endl;


    if(g->v_size == reduced_graph->v_size) {
        R->refine_coloring_first(reduced_graph, reduced_c, -1); // check if discrete and abort?
        return std::pair<sgraph *, coloring *>(reduced_graph, reduced_c);
    } else {
        // recursive...
        nest = new lowdeg;
        return nest->preprocess(reduced_c, reduced_graph, R);
    }
}

double factorial(int n) {
    double res = n;
    while(--n > 0) res *= n;
    return (res == 0)?1:res;
}

long double lowdeg::postprocess(group_shared* G) {
    //std::cout << "dejavu_lowdeg_post------------------------------------------------" << std::endl;
    double grpsize1 = 1;
    long double addsize = 1;
    int grpsize2 = 0;

    if(abort) {
        if(G == nullptr)
            return 1;
        mgrouporder(G->b, G->base_size, G->gp, &G->gens, &grpsize1, &grpsize2, G->domain_size);
        return grpsize1 * pow(10, grpsize2);
    }

    int i, v, deg1_neighbours;
    // determine add size
    work_set_int calculator;
    work_list    colors_seen;
    colors_seen.initialize(domain_size);
    calculator.initialize(domain_size);

    for(i = 0; i < reduced_domain_size; ++i) {
        v = reduced_g_to_g.get(i);

        if(deg1_counter.get(v) <= 0)
            continue;

        for(int j = g->v[v]; j < g->v[v] + g->d[v]; ++j) {
            int neighb_v = g->e[j];
            if(g->d[neighb_v] == 1) {
                int col = c->vertex_to_col[neighb_v];
                if(calculator.get(col) == -1) {
                    colors_seen.push_back(col);
                }
                calculator.inc(col);
            }
        }

        while(!colors_seen.empty()) {
            int col        = colors_seen.pop_back();
            int col_neighb = calculator.get(col) + 1;
            addsize *= factorial(col_neighb);
        }

        calculator.reset();
        colors_seen.reset();

        //for(int j = 0; j < counters_sz; ++j) {
        //    deg1_neighbours = counters[j].get(v) + 1; // ToDo: wrong: can only flip them if same color!
        //
        //}
    }

    if(nest == nullptr) {
        if(G != nullptr) {
            mgrouporder(G->b, G->base_size, G->gp, &G->gens, &grpsize1, &grpsize2, G->domain_size);
            return addsize * grpsize1 * pow(10, grpsize2);
        } else {
            return addsize;
        }
    } else {
        grpsize1 = nest->postprocess(G);
    }
    std::cout << "Group size: ";
    std::cout << grpsize1 * addsize << std::endl;
    return grpsize1 * addsize;

    std::cout << "------------------------------------------------------------------" << std::endl;
}
