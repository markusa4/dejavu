//
// Created by markus on 19.09.19.
//

#ifndef BRUTUS_BIJECTION_H
#define BRUTUS_BIJECTION_H


#include <vector>
#include "coloring.h"

class bijection {
public:
    std::vector<int> map;
    int map_vertex(int v);
    bijection();
    ~bijection();
    void read_from_coloring(coloring *c);
    void compose(bijection p);
    void inverse();
    static bijection random_bijection(int n);
};


#endif //BRUTUS_BIJECTION_H
