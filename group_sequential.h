//
// Created by markus on 20.09.19.
//

#ifndef BRUTUS_GROUP_H
#define BRUTUS_GROUP_H


#include "bijection.h"
// #include "nauty/traces.h"
//#include "nauty/schreier.h"

#include "schreier_shared.h"
#include "concurrentqueue.h"

class group_sequential {
public:
    int domain_size;
    int base_size;
    int* b;
    int added;
    mschreier *gp;
    mpermnode *gens;
    group_sequential(int domain_size, bijection* base_points);
    ~group_sequential();
    bool add_permutation(bijection* p);
    void print_group_size();
};


#endif //BRUTUS_GROUP_H
