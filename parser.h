//
// Created by markus on 19.09.19.
//

#ifndef BRUTUS_PARSER_H
#define BRUTUS_PARSER_H
#include "sgraph.h"
#include <string>

class parser {
public:
    void parse_dimacs_file(std::string filename, sgraph* g);

    void parse_dimacs_file_digraph(std::string filename, sgraph *g);
};


#endif //BRUTUS_PARSER_H
