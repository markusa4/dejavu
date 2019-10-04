//
// Created by markus on 24.09.19.
//

#ifndef DEJAVU_CONFIGURATION_H
#define DEJAVU_CONFIGURATION_H

struct configstruct {
    bool CONFIG_IR_BACKTRACK = true;
    bool CONFIG_IR_BACKTRACK_RANDOM  = true;
    int CONFIG_IR_CELL_SELECTOR = 0;
    int CONFIG_IR_INVARIANT = 0;

    int CONFIG_THREADS_NO_PIPELINE = 2;
    int CONFIG_THREADS_REFINEMENT_WORKERS = 2;
    int CONFIG_THREADS_PIPELINE_DEPTH = 2;
};

extern configstruct config;

#endif //DEJAVU_CONFIGURATION_H
