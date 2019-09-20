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
    std::cout << "perm: ";
    for(int k = 0; k < p->map.size(); ++k) {
        _p[k] = p->map[k];
        std::cout << p->map[k] << " ";
    }
    std::cout << std::endl;
    bool was_added = condaddgenerator(&gp, &gens, _p, domain_size);
    delete[] _p;
    expandschreier(gp, &gens, domain_size);
    return was_added;
}

group::group(int domain_size) {
    schreier_fails(10);
    this->domain_size = domain_size;
    std::cout << "Creating new group... " << std::endl;
    newgroup(&gp, &gens, domain_size);
    b = new int[domain_size];
    for(int i = 0; i < domain_size; ++i) {
        b[i] = domain_size - i - 1;
    }
    getorbits(b, domain_size - 1, gp, &gens, domain_size);
}

void group::print_group_size() {
    double grpsize1;
    int grpsize2;
    grouporder(b, domain_size - 1, gp, &gens, &grpsize1, &grpsize2, domain_size);
    std::cout << grpsize1 << " * 10^" << (grpsize2) << std::endl; // 10 ^ grpsize2 actually
    //deleteunmarked(&gens);
    std::cout << "gens: " << schreier_gens(gens) << std::endl;
}
