#include "dejavu_api.h"

configstruct config;
volatile int dejavu_kill_request = 0;
thread_local int numnodes;
thread_local int colorcost;

struct graph_constr {
    std::vector<int> labels;
    std::vector<std::pair<int, int>> edges;
    std::vector<int> edge_labels;
    bool has_edge_labels;

    graph_constr() {
        labels = std::vector<int>();
        edges = std::vector<std::pair<int, int>>();
        edge_labels = std::vector<int>();
        has_edge_labels = false;
    }
};

class path_node {
public:
    std::vector<int> vertex_to_col;
    long invariant;
    std::vector<int> base_points;
    long get_invariant() {
        return invariant;
    }
    void set_vertex_to_col(int* c, int len) {
        vertex_to_col = std::vector<int>(c, c + len);
    }
    void set_base_points(std::vector<int>& path_vec) {
        base_points.swap(path_vec);
    }
    void set_invariant(long inv) {
        invariant = inv;
    }
};

struct path_constr {
    std::vector<path_node> paths;
};

struct coloring_constr {
    std::vector<int> coloring;
};

static std::vector<graph_constr*> graphs;
static std::vector<path_constr*>  paths;
static std::vector<coloring_constr*>  colorings;

void initialize() {
    graphs = std::vector<graph_constr*>();
    paths = std::vector<path_constr*>();
    colorings = std::vector<coloring_constr*>();
}

void clean() {
    for(int i = 0; i < graphs.size(); ++i) {
        if(graphs[i] != nullptr) {
            delete graphs[i];
            graphs[i] = nullptr;
        }
    }
    graphs.clear();
    for(int i = 0; i < paths.size(); ++i) {
        if(paths[i] != nullptr) {
            delete paths[i];
            paths[i] = nullptr;
        }
    }
    paths.clear();
    for(int i = 0; i < colorings.size(); ++i) {
        if(colorings[i] != nullptr) {
            delete colorings[i];
            colorings[i] = nullptr;
        }
    }
    colorings.clear();
}

void make_coloring_from_graph_constr(int* vertex_to_col, graph_constr* graph) {
    for(int i = 0; i < graph->labels.size(); ++i) {
        vertex_to_col[i] = graph->labels[i];
    }
    return;
}

void make_sgraph_from_graph_constr(sgraph* g, int** vertex_to_col, graph_constr* graph, bool directed_dimacs) {
    std::vector<std::vector<int>> incidence_list;
    const int nv = graph->labels.size();
    int ne;
    if(directed_dimacs) {
        //PRINT("[api] Parsing with directed DIMACS option...")
        assert(edges.size() % 2 == 0);
        ne = graph->edges.size() / 2;
        g->initialize(nv, ne * 2);
    } else {
        ne = graph-> edges.size();
        g->initialize(nv, ne * 2);
    }
    const int avg_deg = ceil(ne / nv) + 1;
    incidence_list.reserve(nv);
    for(int i = 0; i < nv; ++i) {
        incidence_list.emplace_back(std::vector<int>());
        incidence_list[incidence_list.size() - 1].reserve(avg_deg);
    }

    for(int i = 0; i < graph->edges.size(); ++i) {
        const int nv1 = graph->edges[i].first;
        const int nv2 = graph->edges[i].second;
        incidence_list[nv1].push_back(nv2);
        if(!directed_dimacs) {
            incidence_list[nv2].push_back(nv1);
        }
    }

    int epos = 0;
    int vpos = 0;
    int maxd = 0;

    for(size_t i = 0; i < incidence_list.size(); ++i) {
        g->v[vpos] = epos;
        g->d[vpos] = incidence_list[i].size();
        if(g->d[vpos] > maxd)
            maxd = g->d[vpos];
        vpos += 1;
        for(size_t j = 0; j < incidence_list[i].size(); ++j) {
            g->e[epos] = incidence_list[i][j];
            epos += 1;
        }
    }

    g->v_size = nv;
    g->d_size = nv;
    g->e_size = 2 * ne;

    g->max_degree = maxd;

    if(vertex_to_col != nullptr) {
        *vertex_to_col = new int[g->v_size];
        make_coloring_from_graph_constr(*vertex_to_col, graph);
    }
    return;
}

int graph_create(int size) {
    auto new_graph = new graph_constr();
    new_graph->edges.reserve(size);
    new_graph->labels.resize(size);
    graphs.push_back(new_graph);
    return graphs.size() - 1;
}

void graph_delete(int graph_handle) {
    graph_constr* graph = graphs[graph_handle];
    if(graph != nullptr) {
        delete graph;
        graphs[graph_handle] = nullptr;
    }
}

void _graph_add_edge(graph_constr* graph, int v1, int v2) {
    graph->edges.emplace_back(std::pair<int, int>(v1, v2));
}


void graph_add_edge(int graph_handle, int v1, int v2) {
    graph_constr* graph = graphs[graph_handle];
    _graph_add_edge(graph, v1, v2);
}

void graph_add_edge_labelled(int graph_handle, int v1, int v2, int l) {
    graph_constr* graph = graphs[graph_handle];
    graph->has_edge_labels = true;
    graph->edges.emplace_back(std::pair<int, int>(v1, v2));
    graph->edge_labels.emplace_back(l);
    assert(graph->edges.size() == graph->edge_labels.size());
}

void _graph_label(graph_constr* graph, int v, int l) {
    graph->labels[v] = l;
}
void graph_label(int graph_handle, int v, int l) {
    graph_constr* graph = graphs[graph_handle];
    _graph_label(graph, v, l);
}

void make_sgraph_from_handle(sgraph* g, int** vertex_to_col, int graph_handle, bool directed_dimacs) {
    graph_constr* graph = graphs[graph_handle];
    const int nv = graph->labels.size();

    if(graph->has_edge_labels) {
        // make new graph and call again
        graph_constr subdivision_graph;

        subdivision_graph = *graph;
        subdivision_graph.has_edge_labels = false;
        subdivision_graph.edges.clear();

        int highest_label = *std::max_element(subdivision_graph.labels.begin(), subdivision_graph.labels.end());
        subdivision_graph.labels.resize(graph->labels.size() + graph->edges.size());
        subdivision_graph.edges.reserve(graph->edges.size() * 2);

        int edge_vertex_id = graph->labels.size();
        for(int i = 0; i < graph->edges.size(); ++i) {
            if(directed_dimacs) {
                if(graph->edges[i].first > graph->edges[i].second) {
                    continue;
                }
            }
            _graph_add_edge(&subdivision_graph, graph->edges[i].first, edge_vertex_id);
            _graph_add_edge(&subdivision_graph, edge_vertex_id, graph->edges[i].second);
            _graph_label(&subdivision_graph, edge_vertex_id, 1 + highest_label + graph->edge_labels[i]);
            ++edge_vertex_id;
        }
        /*for(int i = 0; i < subdivision_graph.labels.size(); ++i) {
            std::cout << "v" << i << " l: " << subdivision_graph.labels[i] << std::endl;
        }
        for(int i = 0; i < subdivision_graph.edges.size(); ++i) {
            std::cout << "v" << subdivision_graph.edges[i].first << "--" << "v" << subdivision_graph.edges[i].second << std::endl;
        }*/
        make_sgraph_from_graph_constr(g, vertex_to_col, &subdivision_graph, directed_dimacs);
    } else {
        make_sgraph_from_graph_constr(g, vertex_to_col, graph, directed_dimacs);
    }
}

int _graph_size_from_handle(int graph_handle) {
    return graphs[graph_handle]->labels.size();
}

int random_paths(int graph_handle, int max_length, int num, bool fill_paths) {
    sgraph g;
    int** vertex_to_col = new int*;
    make_sgraph_from_handle(&g, vertex_to_col, graph_handle, false);
    if(graphs[graph_handle]->has_edge_labels) {
        config.CONFIG_IR_SELECTOR_FORBIDDEN_TAIL = graphs[graph_handle]->labels.size();
    } else {
        config.CONFIG_IR_SELECTOR_FORBIDDEN_TAIL = INT32_MAX - 1;
    }
    dejavu_api v;

    std::set<std::tuple<int*, int, int*, long>> paths_result;

    v.random_paths(&g, *vertex_to_col, max_length, num, &paths_result);
    delete[] *vertex_to_col;
    delete vertex_to_col;

    path_constr* new_paths = new path_constr;
    std::set<std::tuple<int*, int, int*, long>>::iterator it;
    //std::vector<node> _nodes;
    mark_set used_color;
    used_color.initialize(g.v_size);
    std::vector<int> path_vec;
    for (it = paths_result.begin(); it != paths_result.end(); ++it) {
        std::tuple<int*, int, int*, long> f = *it;
        used_color.reset();
        path_vec.clear();
        const int path_length = std::get<1>(f);;
        const int* path = std::get<0>(f);
        const int* path_vertex_to_col = std::get<2>(f);
        for(int i = 0; i < path_length; ++i) {
            used_color.set(path_vertex_to_col[path[i]]);
            path_vec.push_back(path[i]);
        }
        if(fill_paths) {
            for(int i = 0; i < g.v_size && i < config.CONFIG_IR_SELECTOR_FORBIDDEN_TAIL; ++i) {
                if(path_vec.size() >= max_length)
                    break;
                if(!used_color.get(path_vertex_to_col[i])) {
                    used_color.set(path_vertex_to_col[i]);
                    path_vec.push_back(i);
                }
            }
        }

        new_paths->paths.emplace_back(path_node());
        new_paths->paths[new_paths->paths.size()-1].set_invariant(std::get<3>(f));
        new_paths->paths[new_paths->paths.size()-1].set_vertex_to_col(std::get<2>(f), g.v_size);
        new_paths->paths[new_paths->paths.size()-1].set_base_points(path_vec);
        delete[] std::get<0>(f);
        delete[] std::get<2>(f);
    }
    paths.push_back(new_paths);
    config.CONFIG_IR_SELECTOR_FORBIDDEN_TAIL = INT32_MAX - 1;
    return paths.size() - 1;
}


int path_get_num(int path_handle) {
    return paths[path_handle]->paths.size();
}

int path_get_size(int path_handle, int path_id) {
    return paths[path_handle]->paths[path_id].base_points.size();
}

int path_get_inv(int path_handle, int path_id) {
    return paths[path_handle]->paths[path_id].invariant;
}

int path_get_point(int path_handle, int path_id, int path_pos) {
    return paths[path_handle]->paths[path_id].base_points[path_pos];
}

int path_get_vertex_color(int path_handle, int path_id, int v) {
    return paths[path_handle]->paths[path_id].vertex_to_col[v];
}

volatile bool are_isomorphic(int graph_handle1, int graph_handle2, int err) {
    sgraph g1;
    sgraph g2;

    make_sgraph_from_handle(&g1, nullptr, graph_handle1, false);
    make_sgraph_from_handle(&g2, nullptr, graph_handle2, false);

    config.CONFIG_RAND_ABORT = err;
    dejavu_iso solver;
    bool is_iso = solver.iso(&g1, &g2);
    return is_iso;
}

extern std::vector<bijection<int>> automorphisms(int graph_handle) {
    //dejavu_auto solver;
    //solver.automorphisms(g, vertex_to_col, nullptr);

    // ToDo: return results
}