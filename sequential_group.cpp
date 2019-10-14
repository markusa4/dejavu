//
// Created by markus on 20.09.19.
//

#include <assert.h>
#include <iostream>
#include "sequential_group.h"

bool sequential_group::add_permutation(bijection *p) {
    // copy to proper array
    assert(p->map_sz == domain_size);
    int* _p = p->map;//new int[domain_size];
    //std::cout << "perm: ";
    //for(int k = 0; k < p->map.size(); ++k) {
   //     _p[k] = p->map[k];
   //     //std::cout << p->map[k] << " ";
   // }
    //std::cout << std::endl;
    mschreier_fails(1);
    mexpandschreier(gp, &gens, domain_size);
    bool was_added = maddgenerator(&gp, &gens, _p, domain_size);
    added += 1;
    delete[] _p;
    if(was_added) {
        // search for base points that are not fixed
    }
    return was_added;
}

sequential_group::sequential_group(int domain_size, bijection* base_points) {
    mschreier_fails(0);
    added = 0;
    this->domain_size = domain_size;
    this->base_size = 0;
    std::cout << "Creating new group... " << std::endl;
    mnewgroup(&gp, &gens, domain_size);
    b = new int[domain_size];
    for(int i = 0; i < base_points->map_sz; ++i) {
        b[i] = base_points->map[i];
        base_size += 1;
    }
    mgetorbits(b, base_size, gp, &gens, domain_size);
}

sequential_group::~sequential_group() {
    mfreeschreier(&gp, &gens);
}

void sequential_group::print_group_size() {
    double grpsize1;
    int grpsize2;
    mgrouporder(b, base_size, gp,&gens, &grpsize1, &grpsize2, domain_size);
    std::cout << grpsize1 << " * 10^" << (grpsize2) << std::endl; // 10 ^ grpsize2 actually
    //deleteunmarked(&gens);
    std::cout << "Generators: " << mschreier_gens(gens) << std::endl;
}

