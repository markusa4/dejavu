#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <assert.h>
#include <set>
#include <chrono>
#include "parser.h"

using std::string;
using std::vector;

void parser::parse_dimacs_file_g(std::string filename, sgraph_t<int, int, int>* g) {
    std::cout << "Graph: \t\t" << filename << std::endl;
    std::ifstream infile(filename);
    vector<vector<int>>   incidence_list;
    vector<std::set<int>> incidence_set;
    std::set<int> degrees;
    string line;
    int nv, ne;
    bool first = true;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if(line.empty())
            break;
        switch(first) {
            case true:
                iss >> nv >> ne;
                g->initialize(nv, ne * 2);
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
        }
    }

    int epos = 0;
    int vpos = 0;

    int maxd = 0;

    for(size_t i = 0; i < incidence_list.size(); ++i) {
        g->v[vpos] = epos;
        g->d[vpos] = incidence_list[i].size();
        degrees.insert(g->d[vpos]);
        if(g->d[vpos] > maxd)
            maxd = g->d[vpos];
        vpos += 1;
        for(size_t j = 0; j < incidence_list[i].size(); ++j) {
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

    assert(nv == g->v_size);
    assert(nv == g->d_size);
    assert(2 * ne == g->e_size);
}

void parser::parse_dimacs_file_fast(std::string filename, sgraph_t<int, int, int>* g, int** colmap) {
    //const size_t bufsize = 4*1024;
    //char buf[bufsize];
    std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
    std::cout << "Graph: \t\t" << filename << std::endl;
    std::ifstream infile(filename);
    //infile.rdbuf()->pubsetbuf(buf, bufsize);

    vector<vector<int>> incidence_list;
    std::set<int> degrees;
    std::set<int> colors;
    string line;
    string nv_str, ne_str;
    string nv1_string, nv2_string;
    int nv1, nv2;
    int i;
    int nv, ne;
    while (std::getline(infile, line)) {
        char m = line[0];
        int average_d;
        switch (m) {
            case 'p':
                for(i = 7; i < line.size() && line[i] != ' '; ++i) {
                    nv_str += line[i];
                }

                ++i;
                for(; i < line.size() && line[i] != ' '; ++i) {
                    ne_str += line[i];
                }
                nv = std::stoi(nv_str);
                ne = std::stoi(ne_str);
                average_d = (ne / nv) + 3;
                g->initialize(nv, ne * 2);
                incidence_list.reserve(nv);
                for(int i = 0; i < nv; ++i) {
                    incidence_list.emplace_back(vector<int>());
                    incidence_list[incidence_list.size() - 1].reserve(average_d);
                }
                break;
            case 'e':
                nv1_string = "";
                nv2_string = "";
                for(i = 2; i < line.size() && line[i] != ' '; ++i) {
                    nv1_string += line[i];
                }
                ++i;
                for(; i < line.size() && line[i] != ' '; ++i) {
                    nv2_string += line[i];
                }

                nv1 = std::stoi(nv1_string);
                nv2 = std::stoi(nv2_string);

                incidence_list[nv1 - 1].push_back(nv2 - 1);
                incidence_list[nv2 - 1].push_back(nv1 - 1);
                break;
            case 'n':
                if(*colmap == nullptr)
                    *colmap = new int[nv];
                nv1_string = "";
                nv2_string = "";
                for(i = 2; i < line.size() && line[i] != ' '; ++i) {
                    nv1_string += line[i];
                }
                ++i;
                for(; i < line.size() && line[i] != ' '; ++i) {
                    nv2_string += line[i];
                }

                nv1 = std::stoi(nv1_string);
                nv2 = std::stoi(nv2_string);
                (*colmap)[nv1 - 1] = nv2;
                break;
            default:
                break;
        }
    }

    int epos = 0;
    int vpos = 0;

    int maxd = 0;

    for(size_t i = 0; i < incidence_list.size(); ++i) {
        g->v[vpos] = epos;
        g->d[vpos] = incidence_list[i].size();
        //degrees.insert(g->d[vpos]);
        if(g->d[vpos] > maxd)
            maxd = g->d[vpos];
        vpos += 1;
        for(size_t j = 0; j < incidence_list[i].size(); ++j) {
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
    //std::cout << "Degrees: \t";
    //for(auto it = degrees.begin(); it != degrees.end(); ++it)
    //    std::cout << *it << ", ";
    //std::cout << std::endl;

    assert(nv == g->v_size);
    assert(nv == g->d_size);
    assert(2 * ne == g->e_size);
    const double parse_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
    std::cout << "Parse time: " << parse_time / 1000000.0 << "ms" << std::endl;
}

void parser::parse_dimacs_file(std::string filename, sgraph_t<int, int, int>* g, int** colmap) {
    //const size_t bufsize = 4*1024;
    //char buf[bufsize];
    std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
    std::cout << "Graph: \t\t" << filename << std::endl;
    std::ifstream infile(filename);
    //infile.rdbuf()->pubsetbuf(buf, bufsize);

    vector<vector<int>> incidence_list;
    std::set<int> degrees;
    std::set<int> colors;
    string line;
    int nv, ne;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        //std::cout << line << std::endl;
        char m;
        int average_d;
        if (!(iss >> m)) break;
        switch (m) {
            case 'p':
                iss.ignore(6);
                iss >> nv >> ne;
                average_d = (ne / nv) + 2;
                g->initialize(nv, ne * 2);
                incidence_list.reserve(nv);
                for(int i = 0; i < nv; ++i) {
                    incidence_list.emplace_back(vector<int>());
                    incidence_list[incidence_list.size() - 1].reserve(average_d);
                }
                break;
            case 'e':
                int nv1, nv2;
                iss >> nv1 >> nv2;
                incidence_list[nv1 - 1].push_back(nv2 - 1);
                incidence_list[nv2 - 1].push_back(nv1 - 1);
                break;
            case 'n':
                if(*colmap == nullptr)
                    *colmap = new int[nv];
                int v, col;
                iss >> v >> col;
                (*colmap)[v - 1] = col;
                break;
            default:
                break;
        }
    }

    int epos = 0;
    int vpos = 0;

    int maxd = 0;

    for(size_t i = 0; i < incidence_list.size(); ++i) {
        g->v[vpos] = epos;
        g->d[vpos] = incidence_list[i].size();
        //degrees.insert(g->d[vpos]);
        if(g->d[vpos] > maxd)
            maxd = g->d[vpos];
        vpos += 1;
        for(size_t j = 0; j < incidence_list[i].size(); ++j) {
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
    //std::cout << "Degrees: \t";
    //for(auto it = degrees.begin(); it != degrees.end(); ++it)
    //    std::cout << *it << ", ";
    //std::cout << std::endl;

    assert(nv == g->v_size);
    assert(nv == g->d_size);
    assert(2 * ne == g->e_size);
    const double parse_time = (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
    std::cout << "Parse time: " << parse_time / 1000000.0 << "ms" << std::endl;
}

void parser::parse_dimacs_file_digraph(std::string filename, sgraph_t<int, int, int>* g) {
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
                g->initialize(nv, ne);
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

    for(size_t i = 0; i < incidence_list.size(); ++i) {
        g->v[vpos] = epos;
        g->d[vpos] = incidence_list[i].size();
        if(g->d[vpos] > maxd)
            maxd = g->d[vpos];
        vpos += 1;
        for(size_t j = 0; j < incidence_list[i].size(); ++j) {
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
