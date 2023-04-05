#ifndef DEJAVU_CONFIGURATION_H
#define DEJAVU_CONFIGURATION_H

#include <functional>

// hook for automorphisms found in dejavu
// void dejavu_hook(int n, const int *perm, int nsupp, const int *support)
// int n:        domain size of the graph and automorphism group
// int* perm:    permutation in one-line notation, i.e., perm[i] = j if and only if the permutation maps i to j
// int nsupp:    number of vertices moved by the permutation, i.e., support of the permutation and length of int* support array
// int* support: vertices moved by the permutation
// IMPORTANT NOTE: Try to avoid sequential reads of perm, rather use the support array to only access those parts of the
// permutation that are non-trivial.
typedef const std::function<void(int, const int *, int, const int *)> dejavu_hook;

struct configstruct {
    bool CONFIG_PREPROCESS_COMPRESS = true;
    bool CONFIG_IR_FORCE_SELECTOR = false;
    int  CONFIG_IR_CELL_SELECTOR  = 3;     // selector type, if CONFIG_IR_FORCE_SELECTOR is set
    int  CONFIG_IR_SIZE_FACTOR    = 10;    // tradeoff between restarts and allowed breadth-first width
    bool CONFIG_IR_FULL_INVARIANT = false; // uses a complete invariant and no certification if enabled
    bool CONFIG_IR_FULL_BFS       = false; // enforces full traversal of the search tree
    bool CONFIG_IR_IDLE_SKIP      = true;  // blueprints
    bool CONFIG_IR_FAST_TOLERANCE_INC  = true;  // tolerance of solver is raised more quickly if search tree appears difficult
    bool CONFIG_IR_INDIVIDUALIZE_EARLY = false; // experimental feature, based on an idea by Adolfo Piperno
    int  CONFIG_IR_EXPAND_DEVIATION = 5;   // additional cells processed after deviation is detected
    bool CONFIG_IR_FORCE_EXPAND_DEVIATION = false;
    int  CONFIG_IR_LEAF_STORE_LIMIT = 64;  // amount of leaf colorings stored in full at any given point
    int  CONFIG_IR_SELECTOR_FORBIDDEN_TAIL = INT32_MAX - 1; // selector may not select colors above this limit
    bool CONFIG_IR_SELECTOR_FILL_GRAPH = false;
    bool CONFIG_IR_ALWAYS_STORE = false;
    bool CONFIG_IR_SKIP_FIRST_REFINEMENT = false; // use if initial coloring has already been refined using >=1-WL
    bool CONFIG_BULK_ALLOCATOR      = true; // bulk allocation of multiple colorings, use to relief heap allocator
    bool CONFIG_PREPROCESS          = true;  // enable sassy
    bool CONFIG_PREP_DEACT_PROBE    = false; // sassy: no probing
    bool CONFIG_PREP_DEACT_DEG01    = false; // sassy: no degree 0,1 processing
    bool CONFIG_PREP_DEACT_DEG2     = false; // sassy: no degree 2   processing
    bool CONFIG_PREP_ALT_SCHEDULE   = false;
    bool CONFIG_IR_REFINE_EARLYOUT_LATE  = false;
    bool CONFIG_IR_WRITE_GROUPORDER      = false; // print grouporder
    bool CONFIG_WRITE_AUTOMORPHISMS      = false; // print automorphisms
    bool CONFIG_WRITE_AUTOMORPHISMS_GAP  = false; // print automorphisms in GAP format
    bool CONFIG_PREPROCESS_EDGELIST_SORT = false; // sort edgelists
    bool CONFIG_SOLVE_ISO                = false;
    bool CONFIG_ONLY_COLOR_REF_INVARIANT = false;
    int  CONFIG_RAND_ABORT = 8;            // error probability (higher value means lower error probability)
                                           // error is <(1/2)^(i-1), where i is the given value
                                           // 6 equals an error of roughly 3.2%
                                           // 8 lower than 0.78%
                                           // error is one-sided, probability 0.78% means a generator is missed with
                                           // probability at most 0.78%
                                           // when multiple threads are used, error probability only valid under
                                           // assumption that there is no bias in finding already found automorphisms
                                           // more quickly than other automorphisms
    int CONFIG_THREADS_REFINEMENT_WORKERS = 0; // max number of (additional) threads
};

extern configstruct config;
extern volatile int dejavu_kill_request;

extern thread_local int numnodes;
extern thread_local int colorcost;

#endif //DEJAVU_CONFIGURATION_H
