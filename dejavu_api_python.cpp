#include <boost/python.hpp>
#include "dejavu_api.h"

configstruct config;
volatile int dejavu_kill_request = 0;
thread_local int numnodes;
thread_local int colorcost;

class graph {
    int vertices = 0;
    std::vector<std::pair<int, int>> edges;
public:
    void set_size(int vertices) {
        this->vertices = vertices;
    }

    void add_edge(int from, int to) {
        edges.emplace_back(std::pair<int, int>(from, to));
    }

    sgraph to_sgraph() {
        sgraph_t<int, int, int> g;
        std::vector<std::vector<int>> incidence_list;
        const int nv = vertices;
        const int ne = edges.size();
        g.v = new int[nv];
        g.d = new int[nv];
        g.e = new int[ne * 2];
        for(int i = 0; i < nv; ++i) {
            incidence_list.emplace_back(std::vector<int>());
        }

        for(int i = 0; i < edges.size(); ++i) {
            const int nv1 = edges[i].first;
            const int nv2 = edges[i].second;
            incidence_list[nv1 - 1].push_back(nv2 - 1);
            incidence_list[nv2 - 1].push_back(nv1 - 1);
        }

        int epos = 0;
        int vpos = 0;
        int maxd = 0;

        for(size_t i = 0; i < incidence_list.size(); ++i) {
            g.v[vpos] = epos;
            g.d[vpos] = incidence_list[i].size();
            if(g.d[vpos] > maxd)
                maxd = g.d[vpos];
            vpos += 1;
            for(size_t j = 0; j < incidence_list[i].size(); ++j) {
                g.e[epos] = incidence_list[i][j];
                epos += 1;
            }
        }

        g.v_size = nv;
        g.d_size = nv;
        g.e_size = 2 * ne;

        g.max_degree = maxd;
        return g;
    }
};

class node {
    std::vector<int> vertex_to_col;
    long invariant;
public:
    std::vector<int> get_vertex_to_col() {
        return vertex_to_col;
    }
    long get_invariant() {
        return invariant;
    }
    void set_vertex_to_col(int* c, int len) {
        for (int i = 0; i < len; i++) {
            vertex_to_col.push_back(c[i]);
        }
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
        return &nodes[nodes.size() - 1];
    }
};


nodes _random_paths(graph g, int max_length, int num) {
    sgraph sg = g.to_sgraph();
    std::set<std::pair<int*, long>> paths;
    random_paths(&sg, max_length, num, &paths);
    std::set<std::pair<int*, long>>::iterator it;
    nodes n;
    for (it = paths.begin(); it != paths.end(); ++it) {
        std::pair<int*, long> f = *it;
        node* new_node = n.push_back();
        new_node->set_invariant(f.second);
        new_node->set_vertex_to_col(f.first, sg.v_size);
        delete[] f.first;
    }
    // ToDo delete sg
    return n;
}

BOOST_PYTHON_MODULE(libdejavu_python) {
        using namespace boost::python;
        def("random_paths", _random_paths);
        class_<nodes>("nodes")
            .def("get_size", &nodes::get_size)
            .def("get_nodes", &nodes::get_node)
            ;
        class_<node>("node")
            .def("get_vertex_to_col", &node::get_vertex_to_col)
            .def("get_invariant", &node::get_invariant)
            ;
    class_<graph>("graph")
            .def("add_edge", &graph::add_edge)
            .def("set_size", &graph::set_size)
            ;
}
