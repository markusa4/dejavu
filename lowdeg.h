//
// Created by markus on 11/11/2019.
//

#ifndef DEJAVU_LOWDEG_H
#define DEJAVU_LOWDEG_H


#include "sgraph.h"
#include "pipeline_schreier.h"

class lowdeg {
public:
    std::pair<sgraph*, coloring*> preprocess(sgraph* g);
    void postprocess(mschreier* gp, mpermnode* gens);
};


#endif //DEJAVU_LOWDEG_H
