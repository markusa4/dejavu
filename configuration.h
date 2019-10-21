//
// Created by markus on 24.09.19.
//

#ifndef DEJAVU_CONFIGURATION_H
#define DEJAVU_CONFIGURATION_H

struct configstruct {
    bool CONFIG_IR_BACKTRACK = false;
    bool CONFIG_IR_BACKTRACK_RANDOM  = false;
    int  CONFIG_IR_CELL_SELECTOR = 1;
    int  CONFIG_IR_INVARIANT = 0;
    int  CONFIG_IR_REFINEMENT = 0;
    bool CONFIG_IR_FAST_AUTOPRE = true; // ToDo: option to stop this dynamically from group

    int CONFIG_RAND_ABORT = 5;
    int CONFIG_RAND_ABORT_RAND = 5;

    int CONFIG_THREADS_NO_PIPELINE = 1;
    int CONFIG_THREADS_REFINEMENT_WORKERS = 3;
    int CONFIG_THREADS_PIPELINE_DEPTH = 1;
    int CONFIG_THREADS_PIPELINE_STAGE_MIN = 10;
    bool CONFIG_THREADS_COPYG = false;
    bool CONFIG_THREADS_COLLABORATE = false;
};

extern configstruct config;

#endif //DEJAVU_CONFIGURATION_H
