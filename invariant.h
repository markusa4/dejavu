#ifndef BRUTUS_INVARIANT_H
#define BRUTUS_INVARIANT_H


#include <stack>
#include <vector>
#include "assert.h"

class alignas(64) invariant {
public:
    std::vector<int>* vec_invariant = nullptr;
    invariant*        compareI;
    std::vector<int>* compare_vec;
    bool has_compare = false;
    bool no_write = false;
    int cur_pos = -1;
    inline bool write_top_and_compare(int i) {
        if(no_write) {
            return (*compare_vec)[++cur_pos] == i;
        } else {
            vec_invariant->push_back(i);
            cur_pos += 1;
           // assert(cur_pos == vec_invariant.size() - 1);
            if (has_compare) {
                if ((compareI->vec_invariant)->size() < vec_invariant->size())
                    return false;
                return (*vec_invariant)[cur_pos] == (*compareI->vec_invariant)[cur_pos];
            } else {
                return true;
            }
        }
    }
    void set_compare_invariant(invariant *I);

    void create_vector() {
        vec_invariant = new std::vector<int>();
    }
};


#endif //BRUTUS_INVARIANT_H
