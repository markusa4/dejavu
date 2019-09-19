//
// Created by markus on 19.09.19.
//

#include "selector.h"

// "first non-trivial color", -1 if coloring is discrete
int selector::select_color(graph *g, coloring *c) {
    for(int i = 0; i < c->ptn.size();){
        if(c->ptn[i] != 0) {
            return i;
        }
        i += c->ptn[i] + 1;
    }
    return - 1;
}
