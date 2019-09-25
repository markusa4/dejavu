//
// Created by markus on 24.09.19.
//

#ifndef DEJAVU_BADTREE_H
#define DEJAVU_BADTREE_H


#include <set>
#include <vector>
#include <unordered_map>

class badtree {
    std::unordered_map<int, badtree*> choices;
    bool traverse_down(int choice);
    void mark_bad(std::vector<int>* path);
};


#endif //DEJAVU_BADTREE_H
