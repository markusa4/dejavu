#ifndef DEJAVU_CONFIGURATION_H
#define DEJAVU_CONFIGURATION_H

// hook for automorphisms found in dejavu
// void dejavu_hook(int n, const int *perm, int nsupp, const int *support)
// int n:        domain size of the graph and automorphism group
// int* perm:    permutation in one-line notation, i.e., perm[i] = j if and only if the permutation maps i to j
// int nsupp:    number of vertices moved by the permutation, i.e., support of the permutation and length of int* support array
// int* support: vertices moved by the permutation
// IMPORTANT NOTE: Try to avoid sequential reads of perm, rather use the support array to only access those parts of the
// permutation that are non-trivial.
typedef void dejavu_hook(int, const int *, int, const int *);

struct configstruct {
    int  CONFIG_IR_CELL_SELECTOR  = 3;     // selector type, if CONFIG_IR_FORCE_SELECTOR is set
    int  CONFIG_IR_SIZE_FACTOR    = 10;    // tradeoff between restarts and allowed breadth-first width (default 10)
    bool CONFIG_IR_FULL_INVARIANT = false; // uses a complete invariant and no certification if enabled
    bool CONFIG_IR_FULL_BFS       = false; // enforces full traversal of the search tree (maybe good for asymmetric)
    bool CONFIG_IR_FORCE_SELECTOR = false;
    bool CONFIG_IR_IDLE_SKIP      = true;
    bool CONFIG_IR_FAST_TOLERANCE_INC = true;
    bool CONFIG_IR_INDIVIDUALIZE_EARLY = false;
    int  CONFIG_IR_EXPAND_DEVIATION = 5;  // default 5
    bool CONFIG_IR_FORCE_EXPAND_DEVIATION = false;
    int  CONFIG_IR_LEAF_STORE_LIMIT = 64; // default 64
    int  CONFIG_IR_SELECTOR_FORBIDDEN_TAIL = INT32_MAX - 1;
    bool CONFIG_IR_SELECTOR_FILL_GRAPH = false;
    bool CONFIG_IR_ALWAYS_STORE = false;
    bool config_IR_SKIP_FIRST_REFINEMENT = false; // use if initial coloring has already been refined using >=1-WL

    bool CONFIG_BULK_ALLOCATOR = true;

    bool CONFIG_IR_WRITE_GROUPORDER      = false;
    bool CONFIG_WRITE_AUTOMORPHISMS      = false;
    bool CONFIG_PREPROCESS_COMPRESS      = false;
    bool CONFIG_PREPROCESS_EDGELIST_SORT = false;
    bool CONFIG_PREPROCESS               = false;
    bool CONFIG_SOLVE_ISO                = false;

    bool ONLY_COLOR_REF_INVARIANT       = false;

    int CONFIG_RAND_ABORT = 8;             // determines error probability (higher value means lower error probability)
                                           // error is at most (1 / 2)^(i-1), where i is the given value
                                           // 6 equals an error of roughly 3.2%
                                           // 8 lower than 0.78%

    int CONFIG_THREADS_REFINEMENT_WORKERS = 0; // number of (additional) threads to use
};

extern configstruct config;
extern volatile int dejavu_kill_request;

extern thread_local int numnodes;
extern thread_local int colorcost;

#endif //DEJAVU_CONFIGURATION_H
