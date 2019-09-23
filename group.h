//
// Created by markus on 20.09.19.
//

#ifndef BRUTUS_GROUP_H
#define BRUTUS_GROUP_H


#include "bijection.h"
#include "nauty/traces.h"
#include "nauty/naugroup.h"

class group {
public:
    int domain_size;
    int* b;
    schreier *gp;
    permnode *gens;
    group(int domain_size);
    ~group();
    bool add_permutation(bijection* p);
    void print_group_size();
};


#endif //BRUTUS_GROUP_H
