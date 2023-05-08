// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#include <atomic>
#include <iostream>
#include <mutex>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <set>
#include <cstring>
#include <queue>
#include <memory>
#include <chrono>
#include <iomanip>

#ifndef DEJAVU_UTILITY_H
#define DEJAVU_UTILITY_H

#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__))
    #define OS_WINDOWS
#endif

#define INV_MARK_ENDREF    (INT32_MAX - 5)
#define INV_MARK_STARTCELL (INT32_MAX - 6)
#define INV_MARK_ENDCELL   (INT32_MAX - 7)

#define PRINT(str) std::cout << str << std::endl;
//#define PRINT(str) (void)0;

// TODO maybe can be done faster
unsigned int hash(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

inline bool file_exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

typedef const std::function<void(int, const int *, int, const int *)> dejavu_hook;
extern volatile int dejavu_kill_request;

namespace dejavu {
    thread_local  std::chrono::high_resolution_clock::time_point first;
    thread_local  std::chrono::high_resolution_clock::time_point previous;
    thread_local bool set_first = false;

    static void progress_print_split() {
        PRINT("______________________________________________________________");
    }

    static void progress_print_header() {
        progress_print_split();
        PRINT(std::setw(9) << std::left <<"T (ms)" << std::setw(9) << "δ (ms)" << std::setw(14) << "proc"  << std::setw(16) << "p1"        << std::setw(16)        << "p2");
        progress_print_split();
        PRINT(std::setw(9) << std::left << 0 << std::setw(9) << 0 << std::setw(12) << "start" << std::setw(16) << "_" << std::setw(16) << "_" );
    }

    static void progress_print(const std::string
                               proc, const std::string p1, const std::string p2) {
        if(!set_first) {
            first     = std::chrono::high_resolution_clock::now();
            previous  = first;
            set_first = true;
        }
        auto now = std::chrono::high_resolution_clock::now();
        PRINT(std::fixed << std::setprecision(2) << std::setw(9) << std::left << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - first).count()) / 1000000.0  << std::setw(9) << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - previous).count()) / 1000000.0  << std::setw(12) << proc << std::setw(16) << p1 << std::setw(16) << p2);
        previous = now;
    }

    static void progress_print(const std::string
                               proc, const double p1, const double p2) {
        if(!set_first) {
            first     = std::chrono::high_resolution_clock::now();
            previous  = first;
            set_first = true;
        }

        auto now = std::chrono::high_resolution_clock::now();
        PRINT(std::fixed << std::setprecision(2) << std::setw(9) << std::left << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - first).count()) / 1000000.0  << std::setw(9) << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - previous).count()) / 1000000.0  << std::setw(12) << proc << std::setw(16) << p1 << std::setw(16) << p2);
        previous = now;
    }
}

#endif //DEJAVU_UTILITY_H
