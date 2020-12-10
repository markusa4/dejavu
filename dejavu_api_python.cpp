#include <boost/python.hpp>
#include "dejavu_api.h"

configstruct config;
volatile int dejavu_kill_request = 0;
thread_local int numnodes;
thread_local int colorcost;

template<class T>
boost::python::list std_vector_to_py_list(const std::vector<T>& v)
{
    boost::python::list l;
    typename std::vector<T>::const_iterator it;
    for (it = v.begin(); it != v.end(); ++it)
        l.append(*it);
    return l;
}

class pygraph {
    int  vertices = 0;
public:
    bool directed_dimacs = false;
    std::vector<std::pair<int, int>> edges;
    void set_size(int vertices) {
        this->vertices = vertices;
    }

    int get_size() {
        return this->vertices;
    }

    void set_directed_dimacs(bool directed_dimacs) {
        this->directed_dimacs = directed_dimacs;
    }

    void add_edge(int from, int to) {
        edges.emplace_back(std::pair<int, int>(from, to));
    }

    void add_edges(boost::python::list edges_from, boost::python::list edges_to) {
        boost::python::ssize_t len = boost::python::len(edges_from);
        for(int i=0; i<len; i++) {
            const int from = boost::python::extract<int>(edges_from[i]);
            const int to = boost::python::extract<int>(edges_to[i]);
            edges.emplace_back(std::pair<int, int>(from, to));
        }
    }

    void to_sgraph(sgraph* g) {
        std::vector<std::vector<int>> incidence_list;
        const int nv = vertices;
        int ne;
        if(directed_dimacs) {
            //PRINT("[api] Parsing with directed DIMACS option...")
            assert(edges.size() % 2 == 0);
            ne = edges.size() / 2;
            g->initialize(nv, ne * 2);
        } else {
            ne = edges.size();
            g->initialize(nv, ne * 2);
        }
        const int avg_deg = ceil(ne / nv) + 1;
        incidence_list.reserve(nv);
        for(int i = 0; i < nv; ++i) {
            incidence_list.emplace_back(std::vector<int>());
            incidence_list[incidence_list.size() - 1].reserve(avg_deg);
        }

        for(int i = 0; i < edges.size(); ++i) {
            const int nv1 = edges[i].first;
            const int nv2 = edges[i].second;
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
        return;
    }
};

class pycoloring {
public:
    std::vector<int> colors;
    void add_color(int color) {
        colors.push_back(color);
    }
    void add_colors(boost::python::list vertex_to_col) {
        boost::python::ssize_t len = boost::python::len(vertex_to_col);
        for(int i=0; i<len; i++) {
            const int color = boost::python::extract<int>(vertex_to_col[i]);
            colors.push_back(color);
        }
    }
    int get_color(int node) {
        return colors[node];
    }
};

class node {
    std::vector<int> vertex_to_col;
    long invariant;
    std::vector<int> base_points;
public:
    boost::python::list get_vertex_to_col() {
        return std_vector_to_py_list(vertex_to_col);
    }
    boost::python::list get_base_points() {
        return std_vector_to_py_list(base_points);
    }
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

class nodes {
    int size = 0;
    std::vector<node> nodes;
public:
    int get_size() {
        return size;
    }
    node get_node(int i) {
        return nodes[i];
    }
    node* push_back() {
        nodes.emplace_back(node());
        size += 1;
        return &nodes[nodes.size() - 1];
    }
};


boost::python::list _random_paths(pygraph g, pycoloring c, int max_length, int num, bool fill_paths) {
    sgraph sg;
    g.to_sgraph(&sg);
    sg.print();
    int* vertex_to_col = new int[sg.v_size];
    for(int i = 0; i < sg.v_size; ++i)
        vertex_to_col[i] = c.get_color(i);
    std::set<std::tuple<int*, int, int*, long>> paths;
    random_paths(&sg, vertex_to_col, max_length, num, &paths);
    std::set<std::tuple<int*, int, int*, long>>::iterator it;
    std::vector<node> _nodes;
    mark_set used_color;
    used_color.initialize(sg.v_size);
    std::vector<int> path_vec;
    for (it = paths.begin(); it != paths.end(); ++it) {
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
            for(int i = 0; i < sg.v_size && i < config.CONFIG_IR_SELECTOR_FORBIDDEN_TAIL; ++i) {
                if(path_vec.size() >= max_length)
                    break;
                if(!used_color.get(path_vertex_to_col[i])) {
                    used_color.set(path_vertex_to_col[i]);
                    path_vec.push_back(i);
                }
            }
        }

        _nodes.emplace_back(node());
        _nodes[_nodes.size()-1].set_invariant(std::get<3>(f));
        _nodes[_nodes.size()-1].set_vertex_to_col(std::get<2>(f), sg.v_size);
        _nodes[_nodes.size()-1].set_base_points(path_vec);
        delete[] std::get<0>(f);
        delete[] std::get<2>(f);
    }

    delete[] vertex_to_col;
    return std_vector_to_py_list(_nodes);
}

boost::python::list _random_paths_edge_colored(pygraph g, pycoloring c, pycoloring ce, int max_length, int num, bool fill_paths) {
    pygraph ge;
    pycoloring c_ce;

    for(int i = 0; i < g.get_size(); ++i) {
        c_ce.add_color(c.get_color(i));
    }

    int edge_vertex_id = g.get_size();
    for(int i = 0; i < g.edges.size(); ++i) {
        if(g.directed_dimacs) {
            if(g.edges[i].first > g.edges[i].second) {
                continue;
            }
        }
        ge.add_edge(g.edges[i].first, edge_vertex_id);
        ge.add_edge(edge_vertex_id, g.edges[i].second);
        c_ce.add_color(ce.get_color(i));
        ++edge_vertex_id;
    }

    ge.set_directed_dimacs(false);
    ge.set_size(edge_vertex_id);
    config.CONFIG_IR_SELECTOR_FORBIDDEN_TAIL = g.get_size();
    boost::python::list l = _random_paths(g, c_ce, max_length, num, fill_paths);
    config.CONFIG_IR_SELECTOR_FORBIDDEN_TAIL = INT32_MAX - 1;
    return l;
}

BOOST_PYTHON_MODULE(libdejavu_python) {
        using namespace boost::python;
        def("_random_paths", _random_paths);
        def("_random_paths_edge_colored", _random_paths_edge_colored);
        class_<nodes>("nodes")
            .def("get_size", &nodes::get_size)
            .def("get_node", &nodes::get_node)
            ;
        class_<node>("node")
            .def("get_vertex_to_col", &node::get_vertex_to_col)
            .def("get_base_points", &node::get_base_points)
            .def("get_invariant", &node::get_invariant)
            ;
        class_<pygraph>("pygraph")
                .def("add_edge", &pygraph::add_edge)
                .def("add_edges", &pygraph::add_edges)
                .def("set_size", &pygraph::set_size)
                .def("set_directed_dimacs", &pygraph::set_directed_dimacs)
                ;
        class_<pycoloring>("pycoloring")
                .def("add_color", &pycoloring::add_color)
                .def("add_colors", &pycoloring::add_colors)
                .def("get_color", &pycoloring::get_color)
                ;
}
