//
// Created by markus on 19.09.19.
//

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <set>
#include "parser.h"

using std::string;
using std::vector;

void parser::parse_dimacs_file_g(std::string filename, sgraph* g) {
    std::cout << "Graph: \t\t" << filename << std::endl;
    std::ifstream infile(filename);
    vector<vector<int>> incidence_list;
    vector<std::set<int>>    incidence_set;
    std::set<int> degrees;
    string line;
    int nv, ne;
    bool first = true;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if(line.empty())
            break;
        switch (first) {
            case true:
                iss >> nv >> ne;
                g->v = new int[nv];
                g->d = new int[nv];
                g->e = new int[ne * 2];
                for(int i = 0; i < nv; ++i) {
                    incidence_list.emplace_back(vector<int>());
                }
                first = false;
                break;
            case false:
                int nv1, nv2;
                iss >> nv1 >> nv2;
                incidence_list[nv1].push_back(nv2);
                incidence_list[nv2].push_back(nv1);
                break;
            default:
                break;
        }
    }

    int epos = 0;
    int vpos = 0;

    int maxd = 0;

    for(int i = 0; i < incidence_list.size(); ++i) {
        g->v[vpos] = epos;
        g->d[vpos] = incidence_list[i].size();
        degrees.insert(g->d[vpos]);
        if(g->d[vpos] > maxd)
            maxd = g->d[vpos];
        vpos += 1;
        for(int j = 0; j < incidence_list[i].size(); ++j) {
            g->e[epos] = incidence_list[i][j];
            epos += 1;
        }
    }

    g->v_size = nv;
    g->d_size = nv;
    g->e_size = 2 * ne;

    g->max_degree = maxd;

    std::cout << "Vertices: \t" << g->v_size << std::endl;
    std::cout << "Edges: \t\t" << g->e_size << std::endl;
    std::cout << "Degrees: \t";
    for(auto it = degrees.begin(); it != degrees.end(); ++it)
        std::cout << *it << ", ";
    std::cout << std::endl;

    assert(nv == g->v_size);
    assert(nv == g->d_size);
    assert(2 * ne == g->e_size);
}
/*
void parser::parse_dimacs_file_cnf(std::string filename, sgraph* g) {
    std::cout << "Graph: \t\t" << filename << std::endl;
    std::ifstream infile(filename);
    vector<vector<int>> incidence_list;
    vector<std::set<int>>    incidence_set;
    std::set<int> degrees;
    string line;
    int nv, ne;
    bool first = true;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if(line.empty())
            break;
        switch (first) {
            case true:
                iss >> nv >> ne;

                for(int i = 0; i < nv; ++i) {
                    incidence_list.emplace_back(vector<int>());
                }
                first = false;
                break;
            case false:
                int nv1, nv2;
                iss >> nv1 >> nv2;
                incidence_list[nv1].push_back(nv2);
                incidence_list[nv2].push_back(nv1);
                break;
            default:
                break;
        }
    }

    int epos = 0;
    int vpos = 0;

    int maxd = 0;

    for(int i = 0; i < incidence_list.size(); ++i) {
        g->v[vpos] = epos;
        g->d[vpos] = incidence_list[i].size();
        degrees.insert(g->d[vpos]);
        if(g->d[vpos] > maxd)
            maxd = g->d[vpos];
        vpos += 1;
        for(int j = 0; j < incidence_list[i].size(); ++j) {
            g->e[epos] = incidence_list[i][j];
            epos += 1;
        }
    }

    g->v_size = nv;
    g->d_size = nv;
    g->e_size = 2 * ne;

    g->max_degree = maxd;

    std::cout << "Vertices: \t" << g->v_size << std::endl;
    std::cout << "Edges: \t\t" << g->e_size << std::endl;
    std::cout << "Degrees: \t";
    for(auto it = degrees.begin(); it != degrees.end(); ++it)
        std::cout << *it << ", ";
    std::cout << std::endl;

    assert(nv == g->v_size);
    assert(nv == g->d_size);
    assert(2 * ne == g->e_size);
}
*/


void parser::parse_dimacs_file(std::string filename, sgraph* g) {
    std::cout << "Graph: \t\t" << filename << std::endl;
    std::ifstream infile(filename);
    vector<vector<int>> incidence_list;
    vector<std::set<int>>    incidence_set;
    std::set<int> degrees;
    string line;
    int nv, ne;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        //std::cout << line << std::endl;
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

    int maxd = 0;

    for(int i = 0; i < incidence_list.size(); ++i) {
        g->v[vpos] = epos;
        g->d[vpos] = incidence_list[i].size();
        degrees.insert(g->d[vpos]);
        if(g->d[vpos] > maxd)
            maxd = g->d[vpos];
        vpos += 1;
        for(int j = 0; j < incidence_list[i].size(); ++j) {
            g->e[epos] = incidence_list[i][j];
            epos += 1;
        }
    }

    g->v_size = nv;
    g->d_size = nv;
    g->e_size = 2 * ne;

    g->max_degree = maxd;

    std::cout << "Vertices: \t" << g->v_size << std::endl;
    std::cout << "Edges: \t\t" << g->e_size << std::endl;
    std::cout << "Degrees: \t";
    for(auto it = degrees.begin(); it != degrees.end(); ++it)
        std::cout << *it << ", ";
    std::cout << std::endl;

    assert(nv == g->v_size);
    assert(nv == g->d_size);
    assert(2 * ne == g->e_size);
}

void parser::parse_dimacs_file_digraph(std::string filename, sgraph* g) {
    std::cout << "Graph: \t\t" << filename << std::endl;
    std::ifstream infile(filename);
    vector<vector<int>>      incidence_list;
    vector<std::set<int>>    incidence_set;
    string line;
    int nv, ne;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        //std::cout << line << std::endl;
        char m;
        if (!(iss >> m)) break;
        switch (m) {
            case 'p':
                iss.ignore(6);
                iss >> nv >> ne;
                g->v = new int[nv];
                g->d = new int[nv];
                g->e = new int[ne];
                for(int i = 0; i < nv; ++i) {
                    incidence_list.emplace_back(vector<int>());
                    //incidence_set.emplace_back(std::set<int>());
                }
                break;
            case 'e':
                int nv1, nv2;
                iss >> nv1 >> nv2;
                incidence_list[nv1 - 1].push_back(nv2 - 1);
                //incidence_list[nv2 - 1].push_back(nv1 - 1);
                break;
            default:
                break;
        }
    }

    int epos = 0;
    int vpos = 0;

    int maxd = 0;

    for(int i = 0; i < incidence_list.size(); ++i) {
        g->v[vpos] = epos;
        g->d[vpos] = incidence_list[i].size();
        if(g->d[vpos] > maxd)
            maxd = g->d[vpos];
        vpos += 1;
        for(int j = 0; j < incidence_list[i].size(); ++j) {
            g->e[epos] = incidence_list[i][j];
            epos += 1;
        }
    }

    g->v_size = nv;
    g->d_size = nv;
    g->e_size = ne;

    g->max_degree = maxd;

    std::cout << "Vertices: \t" << g->v_size << std::endl;
    std::cout << "Edges: \t\t" << g->e_size << std::endl;

    assert(nv == g->v_size);
    assert(nv == g->d_size);
    assert(ne == g->e_size);
}
