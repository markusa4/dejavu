//
// Created by markus on 19.09.19.
//

#ifndef BRUTUS_COLORING_H
#define BRUTUS_COLORING_H


#include <vector>

class coloring { // broken copy
public:
    void rewrite_ptn(coloring* c);
    void copy(coloring* c);
    int* lab;
    int* ptn;
    int lab_sz;
    int ptn_sz;
    bool init = false;
    std::vector<int>  vertex_to_col;
    std::vector<int>  vertex_to_lab;
    ~coloring();

    void copy_force(coloring *c);
};

#endif //BRUTUS_COLORING_H
