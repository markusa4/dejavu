#include <atomic>
#include <iostream>
#include <mutex>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include "configuration.h"
#include <fstream>
#include <set>
#include <cstring>
#include <queue>
#include <memory>

#ifndef DEJAVU_UTILITY_H
#define DEJAVU_UTILITY_H

#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__))
    #define OS_WINDOWS
#endif

#define INV_MARK_ENDREF    (INT32_MAX - 5)
#define INV_MARK_STARTCELL (INT32_MAX - 6)
#define INV_MARK_ENDCELL   (INT32_MAX - 7)

#define MASH0(i) (i * (35235237 - i * 5))
#define MASH1(i) (i * (352355 - i * 3))
#define MASH2(i) ((i + 1) * (423733 - (i + 1)))
#define MASH3(i) ((i + 1) * (423233 - (i + 1)))
#define MASH4(i) ((i + 1) * (23524361 - i * 3))
#define MASH5(i) ((i + 1) * (23524361 - i * 3))

#define PRINT(str) std::cout << str << std::endl;
//#define PRINT(str) (void)0;

inline bool file_exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

#endif //DEJAVU_UTILITY_H
