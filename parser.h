#ifndef DEJAVU_PARSER_H
#define DEJAVU_PARSER_H
#include "sgraph.h"
#include <string>

static void parse_dimacs_file_fast(const std::string& filename, dejavu::sgraph* g, int** colmap) {
    //const size_t bufsize = 4*1024;
    //char buf[bufsize];
    std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();
    std::cout << "Graph: \t\t" << filename << std::endl;
    std::ifstream infile(filename);
    //infile.rdbuf()->pubsetbuf(buf, bufsize);

    std::vector<std::vector<int>> incidence_list;
    std::set<int> degrees;
    std::set<int> colors;
    std::string line;
    std::string nv_str, ne_str;
    std::string nv1_string, nv2_string;
    int nv1, nv2;
    size_t i;
    int nv = 0;
    int ne = 0;
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
                for(int j = 0; j < nv; ++j) {
                    incidence_list.emplace_back();
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

    for(auto & i : incidence_list) {
        g->v[vpos] = epos;
        g->d[vpos] = (int) i.size();
        //degrees.insert(g->d[vpos]);
        if(g->d[vpos] > maxd)
            maxd = g->d[vpos];
        vpos += 1;
        for(int j : i) {
            g->e[epos] = j;
            epos += 1;
        }
    }

    g->v_size = nv;
    g->e_size = 2 * ne;

    std::cout << "Vertices: \t" << g->v_size << std::endl;
    std::cout << "Edges: \t\t" << g->e_size << std::endl;
    //std::cout << "Degrees: \t";
    //for(auto it = degrees.begin(); it != degrees.end(); ++it)
    //    std::cout << *it << ", ";
    //std::cout << std::endl;

    assert(nv == g->v_size);
    assert(2 * ne == g->e_size);
    const double parse_time = (double) (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
    std::cout << "Parse time: " << parse_time / 1000000.0 << "ms" << std::endl;
}



#endif //DEJAVU_PARSER_H
