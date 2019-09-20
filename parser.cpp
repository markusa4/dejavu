//
// Created by markus on 19.09.19.
//

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "parser.h"

using std::string;
using std::vector;

void parser::parse_dimacs_file(std::string filename, sgraph* g) {
    std::ifstream infile(filename);
    vector<vector<int>> incidence_list;
    string line;
    int nv, ne;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        char m;
        if (!(iss >> m)) break;
        switch (m) {
            case 'p':
                iss.ignore(6);
                iss >> nv >> ne;
                g->v.reserve(nv);
                g->d.reserve(nv);
                g->e.reserve(ne * 2);
                for(int i = 0; i < nv; ++i) {
                    incidence_list.emplace_back(vector<int>());
                }
                break;
            case 'e':
                int nv1, nv2;
                iss >> nv1 >> nv2;
                incidence_list[nv1 - 1].push_back(nv2 - 1);
                incidence_list[nv2 - 1].push_back(nv1 - 1);
                break;
            default:
                break;
        }
    }

    for(int i = 0; i < incidence_list.size(); ++i) {
        g->v.push_back(g->e.size());
        g->d.push_back(incidence_list[i].size());
        for(int j = 0; j < incidence_list[i].size(); ++j) {
            g->e.push_back(incidence_list[i][j]);
        }
    }

    std::cout << "Vertices: " << g->v.size() << std::endl;
    std::cout << "Edges: " << g->e.size() << std::endl;

    assert(nv == g->v.size());
    assert(nv == g->d.size());
    assert(2 * ne == g->e.size());
}
