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
    std::cout << "Graph: \t\t" << filename << std::endl;
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
                g->v = new int[nv];
                g->d = new int[nv];
                g->e = new int[ne * 2];
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

    int epos = 0;
    int vpos = 0;
    for(int i = 0; i < incidence_list.size(); ++i) {
        g->v[vpos] = epos;
        g->d[vpos] = incidence_list[i].size();
        vpos += 1;
        for(int j = 0; j < incidence_list[i].size(); ++j) {
            g->e[epos] = incidence_list[i][j];
            epos += 1;
        }
    }

    g->v_size = nv;
    g->d_size = nv;
    g->e_size = 2 * ne;

    std::cout << "Vertices: \t" << g->v_size << std::endl;
    std::cout << "Edges: \t\t" << g->e_size << std::endl;

    assert(nv == g->v_size);
    assert(nv == g->d_size);
    assert(2 * ne == g->e_size);
}
