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

#ifndef SASSY_UTILITY_H
#define SASSY_UTILITY_H

namespace sassy {

#define INV_MARK_ENDREF    (INT32_MAX - 5)
#define INV_MARK_STARTCELL (INT32_MAX - 6)
#define INV_MARK_ENDCELL   (INT32_MAX - 7)

#define MASH0(i) ((unsigned long) i * (35235237 - (unsigned long) i * 5))
#define MASH1(i) ((unsigned long) i * (352355 - (unsigned long) i * 3))
#define MASH2(i) (((unsigned long) i + 1) * (423733 - ((unsigned long) i + 1)))
#define MASH3(i) (((unsigned long) i + 1) * (423233 - ((unsigned long) i + 1)))
#define MASH4(i) (((unsigned long) i + 1) * (23524361 - (unsigned long) i * 3))
#define MASH5(i) (((unsigned long) i + 1) * (23524361 - (unsigned long) i * 3))

//#define PRINT(str) {if(config->CONFIG_PRINT) {std::cout << str << std::endl;}}
//#define PRINT(str) {}

// metrics used to compare strategies
    struct strategy_metrics {
        int color_refinement_cost = 0;
    };
}

#endif //SASSY_UTILITY_H
