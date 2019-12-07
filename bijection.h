#ifndef DEJAVU_BIJECTION_H
#define DEJAVU_BIJECTION_H

#include "coloring.h"

class bijection {
public:
    bool init = false;
    bool mark = false;
    int* map;
    int map_sz;
    bool non_uniform  = false;
    bool foreign_base = false;
    bool certified    = false;
    int map_vertex(int v);
    bijection();
    ~bijection();
    void read_from_coloring(coloring *c);
    void compose(bijection* p);
    void inverse();
    static void random_bijection(bijection* p, int n, unsigned seed);

    void print();
    void not_deletable();
    void deletable();
};


#endif //DEJAVU_BIJECTION_H
