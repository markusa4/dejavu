#ifndef DEJAVU_PARSER_H
#define DEJAVU_PARSER_H
#include "sgraph.h"
#include <string>

class parser {
public:
    void parse_dimacs_file_dynamic(std::string filename, dynamic_sgraph* g);
    void parse_dimacs_file(std::string filename, sgraph_t<int, int, int>* g, int** colmap);
    void parse_dimacs_file_g(std::string filename, sgraph_t<int, int, int>* g);
    void parse_dimacs_file_digraph(std::string filename, sgraph_t<int, int, int>* g);
};


#endif //DEJAVU_PARSER_H
