//
// Created by markus on 19.09.19.
//

#ifndef BRUTUS_SELECTOR_H
#define BRUTUS_SELECTOR_H


#include "coloring.h"
#include "graph.h"

class selector {
public:
    int select_color(graph *g, coloring *c);
};


#endif //BRUTUS_SELECTOR_H
