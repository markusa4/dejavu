//
// Created by markus on 20.09.19.
//

#ifndef BRUTUS_GROUP_H
#define BRUTUS_GROUP_H


#include "bijection.h"
// #include "nauty/traces.h"
//#include "nauty/schreier.h"

#include "my_schreier.h"


class sequential_group {
public:
    int domain_size;
    int base_size;
    int* b;
    int added;
    schreier *gp;
    permnode *gens;
    sequential_group(int domain_size, bijection* base_points);
    ~sequential_group();
    bool add_permutation(bijection* p);
    void print_group_size();
};


#endif //BRUTUS_GROUP_H
