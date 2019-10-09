//
// Created by markus on 09.10.19.
//

#ifndef DEJAVU_COLORING_BUCKET_H
#define DEJAVU_COLORING_BUCKET_H


class coloring_bucket {
public:
    int* lab;
    int* ptn;
    int numcells = 1;
    int lab_sz;
    int ptn_sz;

    bool init = false;
    ~coloring_bucket();
    coloring_bucket();
    void copy(coloring_bucket* c);
    void print();
};


#endif //DEJAVU_COLORING_BUCKET_H
