// Copyright 2025 Markus Anders
// This file is part of dejavu 2.1.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_PARSER_H
#define DEJAVU_PARSER_H

#include <cstdio>
#include <fstream>
#include <charconv>
#include "graph.h"
#ifdef OS_LINUX
#include <fcntl.h>
#endif

namespace dejavu {

    static inline bool read_unsigned(FILE* file, uint32_t& n) {
        int m = dejavu_getc(file);

        const unsigned d = (unsigned)(m - '0');
        if (d >= 10) {                 
            return false;
        }

        n = d;
        for (;;) {
            m = dejavu_getc(file);
            const unsigned e = (unsigned)(m - '0');

            if (e >= 10) {
                return (m == ' ' || m == '\t' || m == '\n' || m == EOF);
            }

            n = n * 10 + e;
        }
    }

    static inline bool parse_dimacs(const std::string& filename, sgraph* g, int** colmap, bool silent=true,
                                    int seed_permute=0) {
        std::chrono::high_resolution_clock::time_point timer = std::chrono::high_resolution_clock::now();

        FILE* file = fopen(filename.c_str(), "r");
#ifdef OS_LINUX
        posix_fadvise(fileno(file), 0, 0, POSIX_FADV_SEQUENTIAL);
#endif

        constexpr int line_buf_sz = 1024*64;
        char line_buffer[line_buf_sz];
        setvbuf(file, line_buffer, _IOFBF, line_buf_sz);
        dejavu_flockfile(file);

        std::vector<uint32_t> reshuffle;
        //std::vector<std::vector<int>> incidence_list;
        std::vector<std::pair<uint32_t,uint32_t>> edge;
        ds::workspace degree;
        const bool shuffle = seed_permute != 0;
        std::string line;
        uint32_t nv1, nv2;
        uint32_t nv = 0;
        uint32_t ne = 0;
        int line_number = 0;

        bool fail = false;
        bool initialized = false;

        char   m;

        while ((m = dejavu_getc(file)) != EOF) {
            ++line_number;
            switch (m) {
                case 'p': {
                              // eat characters
                              m = dejavu_getc(file); // == ' ' // should not be unsafe since getc will keep returning EOF
                              bool valid = (m == ' ' || m == '\t');
                              m = dejavu_getc(file); // == 'e'
                              valid = valid && (m == 'e' || m == 'E');
                              m = dejavu_getc(file); // == 'd'
                              valid = valid && (m == 'd' || m == 'D');
                              m = dejavu_getc(file); // == 'g'
                              valid = valid && (m == 'g' || m == 'G');
                              m = dejavu_getc(file); // == 'e'
                              valid = valid && (m == 'e' || m == 'E');
                              m = dejavu_getc(file); // == ' '
                              valid = valid && (m == ' ' || m == '\t');
                              if (!valid) {
                                  fail = true;
                                  break;
                              }

                              if (!read_unsigned(file, nv)) {
                                  line = "p unable to read number of vertices";
                                  fail = true;
                                  break;
                              }

                              if (!read_unsigned(file, ne)) {
                                  line = "p unable to read number of edges";
                                  fail = true;
                                  break;
                              }

                              initialized = true;

                              if(shuffle) {
                                  reshuffle.reserve(nv);
                                  for(size_t j = 0; j < nv; ++j) reshuffle.push_back(j + 1);
                                  if(seed_permute != 0) {
                                      std::mt19937 eng(seed_permute);
                                      std::shuffle(std::begin(reshuffle), std::end(reshuffle), eng);
                                  }
                              }

                              edge.reserve(ne);
                              degree.resize(nv);
                              break;
                          }

                case 'e': {

                              if(!initialized || (((m = dejavu_getc(file)) != ' ') && (m != '\t'))) {
                                  fail = true;
                                  line = !initialized?"not initialized":("invalid character " + std::to_string(m));
                                  break;
                              }

                              if (!read_unsigned(file, nv1) || nv1 > nv) {
                                  fail = true;
                                  line = "e v1 missing or invalid";
                                  break;
                              }

                              if (!read_unsigned(file, nv2) || nv2 > nv) {
                                  fail = true;
                                  line = "e v2 missing or invalid";
                                  break;
                              }

                              // potentially unsafe if shuffle active, and nv out of range
                              nv1 = shuffle?reshuffle[nv1-1]:nv1;
                              nv2 = shuffle?reshuffle[nv2-1]:nv2;

                              edge.emplace_back(nv2-1,nv1-1);
                              ++degree[nv2-1];
                              ++degree[nv1-1];
                              break;
                          }

                case 'n': {
                              if(!initialized) {
                                  fail = true;
                                  break;
                              }

                              if(((m = dejavu_getc(file)) != ' ') && (m != '\t')) {
                                  line = "n invalid character";
                                  fail = true;
                                  break;
                              }

                              if(*colmap == nullptr) *colmap = (int *) calloc(nv, sizeof(int));

                              if (!read_unsigned(file, nv1)) {
                                  fail = true;
                                  line = "n vertex missing or invalid";
                                  break;
                              }

                              if (!read_unsigned(file, nv2)) {
                                  fail = true;
                                  line = "n color missing or invalid";
                                  break;
                              }


                              nv1 = shuffle?reshuffle[nv1-1]:nv1;

                              if(nv1 < 0 || nv1 > nv) {

                                  line = "n vertex ID exceeds number of vertices";
                                  fail = true;
                                  break;
                              }

                              (*colmap)[nv1 - 1] = nv2;
                              break;
                          }
                case 'c': {
                              while(((m = dejavu_getc(file)) != '\n') && (m != EOF));
                              break;
                          }
                case ' ':
                case '\t':
                case '\n':
                          break;
                default: {
                             fail = true;
                             line = "invalid character " + std::to_string(m);
                             break;
                         }
            }
            if(fail) break;
        }

        if(fail) {
            std::clog << "parsing failed while reading file in line " << line_number << std::endl;
            std::clog << "fail: " << line << std::endl;
            return false;
        }

        if(!initialized) {
            std::clog << "fail: file does not contain a graph " << std::endl;
            return false;
        }

        if(ne != edge.size()) {
            std::clog << "fail: header claimed " << ne << " edges, but found " << edge.size() << " edges" << std::endl;
            return false;
        }

        g->initialize(nv, ne * 2);
        int epos = 0;
        for(size_t i = 0; i < nv; ++i) {
            g->v[i] = epos;
            g->d[i] = 0; 
            epos += degree[i];
        }

        for(const auto& e : edge) {
            const int v1 = (int) e.first;
            const int v2 = (int) e.second;
            g->e[g->v[v1] + g->d[v1]] = v2;
            ++g->d[v1];
            g->e[g->v[v2] + g->d[v2]] = v1;
            ++g->d[v2];
        }

        // sanity checks
        ds::markset neighbors(nv);
        for(int i = 0; i < static_cast<int64_t>(nv); ++i) {
            dej_assert(g->d[i] == degree[i]);
            int* start = g->e + g->v[i];
            int* end = start + g->d[i];
            neighbors.reset();
            for(int* l = start; l < end; ++l) {
                const int neighbor = *l;
                if(neighbor == i) {
                    std::clog << "fail: vertex " << i+1 << " has a self-loop" << std::endl;
                    return false;
                }
                if(neighbors.get(neighbor)) {
                    std::clog << "fail: vertex " << i+1 << " has multiple edges to " << neighbor+1 << std::endl;
                    return false;
                }
                neighbors.set(neighbor);
            }
        }


        const double read_mb = ftell(file)/1000000.0;

        dejavu_funlockfile(file);
        fclose(file);

        g->v_size = nv;
        g->e_size = 2 * ne;

        dej_assert(static_cast<int>(nv) == g->v_size);
        dej_assert(2 * static_cast<int>(ne) == g->e_size);
        const double parse_time = (double) (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - timer).count());
        if(!silent) std::clog << std::fixed     << std::setprecision(2) 
            << "c parse_time="  << parse_time / 1000000.0 << "ms"
                << ", file_size=" << read_mb << "mb";

        return true;
    }

}


#endif //DEJAVU_PARSER_H
