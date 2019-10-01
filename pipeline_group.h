//
// Created by markus on 01.10.19.
//

#ifndef DEJAVU_PIPELINE_GROUP_H
#define DEJAVU_PIPELINE_GROUP_H


#include "my_schreier.h"
#include "bijection.h"

class pipeline_group {
public:
    int stages;
    int domain_size;
    int base_size;
    int* b;
    int added;
    schreier *gp;
    permnode *gens;
    pipeline_group(int domain_size, bijection* base_points);
    ~pipeline_group();
    bool add_permutation(bijection* p);
    void print_group_size();
    void control_pipeline();
    void pipeline_stage(int n);
};


#endif //DEJAVU_PIPELINE_GROUP_H
