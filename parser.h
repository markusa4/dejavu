//
// Created by markus on 19.09.19.
//

#ifndef BRUTUS_PARSER_H
#define BRUTUS_PARSER_H
#include "graph.h"
#include <string>

class parser {
public:
    void parse_dimacs_file(std::string filename, graph* g);
};


#endif //BRUTUS_PARSER_H
