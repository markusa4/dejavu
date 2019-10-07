//
// Created by markus on 19.09.19.
//

#ifndef BRUTUS_COLORING_H
#define BRUTUS_COLORING_H


#include <vector>

class coloring {
public:
    void rewrite_ptn(coloring* c);
    std::vector<int>  lab;
    std::vector<int>  ptn;
    std::vector<int>  vertex_to_col;
    std::vector<int>  vertex_to_lab;
};

#endif //BRUTUS_COLORING_H
