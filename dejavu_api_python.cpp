#include <boost/python.hpp>
#include "dejavu_api.h"

configstruct config;
volatile int dejavu_kill_request = 0;
thread_local int numnodes;
thread_local int colorcost;

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


nodes _random_paths(sgraph* g, int max_length, int num) {
    std::set<std::pair<int*, long>> paths;
    random_paths(g, max_length, num, &paths);
    std::set<std::pair<int*, long>>::iterator it;
    nodes n;
    for (it = paths.begin(); it != paths.end(); ++it) {
        std::pair<int*, long> f = *it;
        node* new_node = n.push_back();
        new_node->set_invariant(f.second);
        new_node->set_vertex_to_col(f.first, g->v_size);
        delete[] f.first;
    }
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
}
