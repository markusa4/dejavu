//
// Created by markus on 19.09.19.
//

#ifndef BRUTUS_BIJECTION_H
#define BRUTUS_BIJECTION_H


#include <vector>
#include "coloring.h"

class bijection {
public:
    bool init = false;
    bool mark = false;
    int* map;
    int map_sz;
    bool non_uniform  = false;
    bool foreign_base = false;
    int map_vertex(int v);
    bijection();
    ~bijection();
    void read_from_coloring(coloring *c);
    void compose(bijection* p);
    void inverse();
    static void random_bijection(bijection* p, int n);

    void print();
    void not_deletable();
    void deletable();
};


#endif //BRUTUS_BIJECTION_H
