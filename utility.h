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

#if __APPLE__
    #define OS_MAC
#endif

#if __linux__
    #define OS_LINUX
#endif

#define INV_MARK_ENDREF    (INT32_MAX - 5)
#define INV_MARK_STARTCELL (INT32_MAX - 6)
#define INV_MARK_ENDCELL   (INT32_MAX - 7)

#define PRINT(str) std::cout << str << std::endl;

// TODO maybe can be done faster
unsigned int hash(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

unsigned long add_to_hash(unsigned long hash, const int d) {
    const unsigned long ho = hash & 0xff00000000000000; // extract high-order 8 bits from hash
    hash    = hash << 8;                    // shift hash left by 5 bits
    hash    = hash ^ (ho >> 56);            // move the highorder 5 bits to the low-order
    hash    = hash ^ d;                     // XOR into hash

    return hash;
}

inline bool file_exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

typedef const std::function<void(int, const int *, int, const int *)> dejavu_hook;
extern volatile int dejavu_kill_request;

namespace dejavu {

    /**
     * A simple class to store big, positive numbers. Consists of a \a mantissa and a \a exponent, where the value of
     * the number is `mantissa^exponent`.
     *
     * Used to store group sizes.
     */
    class big_number {
    public:
        long double mantissa = 1.0; /**< mantissa, number is `mantissa^exponent`
                                          * \sa exponent */
        int         exponent = 0;   /**< exponent, number is `mantissa^exponent`
                                          * \sa mantissa  */

        friend bool operator<(const big_number& l, const big_number& r)
        {
            return (l.exponent < r.exponent) || (l.exponent == r.exponent && l.mantissa < r.mantissa);
        }

        /**
         * Multiply a \p number to this big_number.
         *
         * @param number The number to multiply.
         */
        void multiply(int number) {
            multiply(number, 0);
        }

        /**
         * Multiply a \p number to this big_number.
         *
         * @param number The number to multiply.
         */
        void multiply(big_number number) {
            multiply(number.mantissa, number.exponent);
        }


        /**
         * Multiply a number consisting of a mantissa (\p other_mantissa) and exponent (\p other_exponent) to this
         * big_number.
         *
         * @param other_mantissa Mantissa of number to multiply.
         * @param other_exponent Exponent of number to multiply.
         */
        void multiply(long double other_mantissa, int other_exponent) {
            if(std::fpclassify(other_mantissa) == FP_INFINITE ||  std::fpclassify(other_mantissa) == FP_NAN) {
                return;
            }
            while (other_mantissa >= 10.0) {
                exponent += 1;
                other_mantissa = other_mantissa / 10;
            }
            exponent += other_exponent;
            mantissa *= other_mantissa;
            man_to_exp();
        }

    private:
        void man_to_exp() {
            if(std::fpclassify(mantissa) == FP_INFINITE ||  std::fpclassify(mantissa) == FP_NAN) {
                return;
            }
            while(mantissa >= 10.0) {
                exponent += 1;
                mantissa = mantissa / 10;
            }
        }
    };

    inline std::ostream& operator<<(std::ostream& out, big_number& number) {
        return out << number.mantissa << "*10^" << number.exponent;
    }

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
