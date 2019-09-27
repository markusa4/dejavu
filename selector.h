//
// Created by markus on 19.09.19.
//

#ifndef BRUTUS_SELECTOR_H
#define BRUTUS_SELECTOR_H


#include "coloring.h"
#include "sgraph.h"

class selector {
public:
    int select_color(sgraph *g, coloring *c);
    int select_color2(sgraph *g, coloring *c);

    int select_color3(sgraph *g, coloring *c);
};


#endif //BRUTUS_SELECTOR_H
