#include "dejavu_api.h"

configstruct config;
volatile int dejavu_kill_request = 0;
thread_local int numnodes;
thread_local int colorcost;

int test_xyz(int a) {
    return 5;
}

void random_paths(sgraph *g, int *vertex_to_col, int max_length, int num,
                         std::set<std::tuple<int *, int, int *, long>> *paths) {
    dejavu_api v;
    v.random_paths(g, vertex_to_col, max_length, num, paths);
}

bijection<int> are_isomorphic(sgraph_t<int, int, int> *g1, int *vertex_to_col1,
        sgraph_t<int, int, int> *g2, int *vertex_to_col2) {
//dejavu_iso solver;
//solver.iso(g1, g2);

// ToDo: return results
}

extern std::vector<bijection<int>> automorphisms(sgraph_t<int, int, int> *g, int *vertex_to_col) {
    dejavu_auto solver;
    solver.automorphisms(g, vertex_to_col, nullptr);

    // ToDo: return results
}