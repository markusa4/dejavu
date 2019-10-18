//
// Created by markus on 18/10/2019.
//

#include "bfs.h"

bfs_element::~bfs_element() {
    if(init_c)
        delete c;
    if(init_I)
        delete I;
    if(init_base)
        delete base;
}

bfs::bfs() {

}

void bfs::initialize(bfs_element *root_node, int domain_size) {

}
