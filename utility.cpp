#include <iostream>
#include "utility.h"
#include "configuration.h"

thread_local bool bulk_domain_reset = false;

int intRand(const int & min, const int & max, int seed) {
    static thread_local std::mt19937 generator(seed);
    std::uniform_int_distribution<int> distribution(min,max);
    return distribution(generator);
}

double doubleRand(const double & min, const double & max, int seed) {
    static thread_local std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distribution(min,max);
    return floor(distribution(generator));
}

// thread_local NFAllocBuf n_buffer;