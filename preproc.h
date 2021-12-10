#include <vector>
#include "sgraph.h"

#ifndef DEJAVU_PREPROC_H
#define DEJAVU_PREPROC_H


class automorphism_sparse {
    std::vector<int> cycles;
};

class automorphism_dense {
    std::vector<int> v_to_v;
};

class generating_set_sparse {
    std::vector<automorphism_sparse> generators;
    void add_automorphism(automorphism_sparse& _auto) {
        generators.push_back(_auto);
    }
};

class generating_set_dense {
    std::vector<automorphism_dense> generators;
    void add_automorphism(automorphism_dense& _auto) {
        generators.push_back(_auto);
    }
};

class restore_info {
    restore_info* next = nullptr;
};

std::tuple<sgraph, restore_info> reduce(sgraph graph) {
    // degree 0
}

generating_set_sparse restore(generating_set_dense g, restore_info r) {
    // degree 0
}

#endif //DEJAVU_PREPROC_H
