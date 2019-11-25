#ifndef DEJAVU_COLORING_H
#define DEJAVU_COLORING_H


#include <vector>

class coloring {
public:
    void copy(coloring* c);
    int* lab;
    int* ptn;
    int lab_sz;
    int ptn_sz;
    bool init = false;
    int* vertex_to_col;
    int* vertex_to_lab;

    std::vector<std::pair<int, int>> color_choices;

    ~coloring();

    void copy_force(coloring *c);

    void initialize(int domain_size);

    bool check();
};

#endif //DEJAVU_COLORING_H
