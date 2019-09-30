//
// Created by markus on 20.09.19.
//

#include <assert.h>
#include <iostream>
#include "group.h"

bool group::add_permutation(bijection *p) {
    // copy to proper array
    assert(p->map.size() == domain_size);
    int* _p = new int[domain_size];
    //std::cout << "perm: ";
    for(int k = 0; k < p->map.size(); ++k) {
        _p[k] = p->map[k];
        //std::cout << p->map[k] << " ";
    }
    //std::cout << std::endl;
    schreier_fails(1);
    expandschreier(gp, &gens, domain_size);
    bool was_added = addgenerator(&gp, &gens, _p, domain_size);
    added += 1;
    delete[] _p;
    if(was_added) {
        // search for base points that are not fixed
       /* bool fixed = false;
        int last_it = 0;
        while(!fixed) {
            fixed = true;
            int *orbits = getorbits(b, base_size, gp, &gens, domain_size);
            for (int i = last_it; i < domain_size; i++) {
                if (orbits[i] != i) {
                    //std::cout << "extended base: " << base_size << std::endl;
                    assert(base_size < domain_size);
                    b[base_size] = orbits[i];
                    base_size += 1;
                    fixed = false;
                    last_it = i;
                    break;
                }
            }
        }*/
    }
    return was_added;
}

group::group(int domain_size, bijection* base_points) {
    schreier_fails(0);
    added = 0;
    this->domain_size = domain_size;
    this->base_size = 0;
    std::cout << "Creating new group... " << std::endl;
    newgroup(&gp, &gens, domain_size);
    b = new int[domain_size];
    for(int i = 0; i < base_points->map.size(); ++i) {
        b[i] = base_points->map[i];
        base_size += 1;
    }
    getorbits(b, base_size, gp, &gens, domain_size);
}

group::~group() {
    freeschreier(&gp, &gens);
}

void group::print_group_size() {
    double grpsize1;
    int grpsize2;
    grouporder(b, base_size, gp,&gens, &grpsize1, &grpsize2, domain_size);
    std::cout << grpsize1 << " * 10^" << (grpsize2) << std::endl; // 10 ^ grpsize2 actually
    //deleteunmarked(&gens);
    std::cout << "Generators: " << schreier_gens(gens) << std::endl;
}

