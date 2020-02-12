#ifndef DEJAVU_CONFIGURATION_H
#define DEJAVU_CONFIGURATION_H

struct configstruct {
    int  CONFIG_IR_CELL_SELECTOR  = 3;     // not used
    bool CONFIG_IR_DENSE          = true;  // automatically set by the solver
    int  CONFIG_IR_SIZE_FACTOR    = 10;    // trade off between restarts and allowed breadth-first width
    bool CONFIG_IR_FULL_INVARIANT = false; // uses a complete invariant and no certification if enabled
    bool CONFIG_IR_FULLBFS = false;        // enforces full traversal of the search tree (maybe good for asymmetric)
    bool CONFIG_IR_FORCE_SELECTOR = false;

    bool CONFIG_WRITE_AUTOMORPHISMS = false;

    bool CONFIG_PREPROCESS = false;

    int CONFIG_BULK_ALLOCATION = 10;

    int CONFIG_RAND_ABORT      = 6;        // determines error probability (higher value means lower error probability)
                                           // error is at most (1 / 2)^(i-1), where i is the given value
                                           // 6 equals an error of roughly 3.2%
    int CONFIG_RAND_ABORT_RAND = -1;

    int CONFIG_THREADS_REFINEMENT_WORKERS = 1; // number of threads to use
};

extern configstruct config;
extern volatile int dejavu_kill_request;

#endif //DEJAVU_CONFIGURATION_H
