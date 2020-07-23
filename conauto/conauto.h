#ifndef DEJAVU_CONAUTO_H
#define DEJAVU_CONAUTO_H

#define NOT_ADJ 0
#define ARC_IN 1
#define ARC_OUT 2
#define ARC_IO 3

struct adjgraph {
    uint64_t *deg;
    int8_t **adj;
    uint16_t num_vert;
    uint32_t num_arc;
};

int are_isomorphic ( const struct adjgraph *g, const struct adjgraph *h );

#endif //DEJAVU_CONAUTO_H
