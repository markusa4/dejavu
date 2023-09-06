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
#include "ds.h"

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

#define PRINT_NO_NEWLINE(str) std::cout << str << std::flush;
#define PRINT(str) std::cout << str << std::endl;
//#define PRINT(str) {};

/**
 * Hash function for unsigned integers.
 *
 * @param x the unsigned integer
 * @return hashed integer
 */
static unsigned int hash(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

/**
 * Accumulate a hash, for example to be used to hash strings of integers.
 *
 * @param hash hash computed so far
 * @param d integer to accumulate to \p hash
 * @return the new hash
 */
static unsigned long add_to_hash(unsigned long hash, const int d) {
    const unsigned long ho = hash & 0xff00000000000000; // extract high-order 8 bits from hash
    hash    = hash << 8;                    // shift hash left by 5 bits
    hash    = hash ^ (ho >> 56);            // move the highorder 5 bits to the low-order
    hash    = hash ^ d;                     // XOR into hash

    return hash;
}

/**
 * Does the file with filename \p name exist?
 *
 * @param name filename to look for
 * @return whether file \p name exists
 */
static inline bool file_exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

typedef void type_dejavu_hook(int, const int*, int, const int*);
typedef std::function<void(int, const int*, int, const int*)> dejavu_hook;
extern volatile int dejavu_kill_request;

namespace dejavu {

    class multi_hook {
    private:
        std::vector<dejavu_hook*> hooks;
        dejavu_hook my_hook;
        void hook_func(int n, const int *p, int nsupp, const int *supp) {
            for(dejavu_hook* hook : hooks) {
                (*hook)(n, p, nsupp, supp);
            }
        }
    public:
        void add_hook(dejavu_hook* hook) {
            hooks.push_back(hook);
        }

        void clear() {
            hooks.clear();
        }

        [[nodiscard]] size_t size() const {
            return hooks.size();
        }

        dejavu_hook* get_hook() {
            my_hook = [this](auto && PH1, auto && PH2, auto && PH3, auto && PH4)
                    { return hook_func(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2),
                               std::forward<decltype(PH3)>(PH3), std::forward<decltype(PH4)>(PH4));
            };
            return &my_hook;
        }
    };

    class ostream_hook {
    private:
        dejavu_hook   my_hook;
        std::ostream& my_ostream;
        dejavu::ds::markset  test_set;

        void hook_func(int n, const int *p, int nsupp, const int *supp) {
            test_set.initialize(n);
            test_set.reset();
            for(int i = 0; i < nsupp; ++i) {
                const int v_from = supp[i];
                if(test_set.get(v_from)) continue;
                int v_next = p[v_from];
                if(v_from == v_next) continue;
                test_set.set(v_from);
                my_ostream << "(" << v_from;
                while(!test_set.get(v_next)) {
                    test_set.set(v_next);
                    my_ostream << " " << v_next;
                    v_next = p[v_next];
                }
                my_ostream << ")";
            }
            my_ostream << std::endl;
        }
    public:
        explicit ostream_hook(std::ostream& ostream) : my_ostream(ostream) {}

        dejavu_hook* get_hook() {
            my_hook = [this](auto && PH1, auto && PH2, auto && PH3, auto && PH4)
            { return hook_func(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2),
                               std::forward<decltype(PH3)>(PH3), std::forward<decltype(PH4)>(PH4));
            };
            return &my_hook;
        }
    };

    class random_source {
        bool true_random = false;
        std::mt19937 pseudo_random_device;
        //std::default_random_engine pseudo_random_device;
        std::random_device true_random_device;
    public:
        random_source(bool set_true_random, int set_seed) {
            true_random = set_true_random;
            pseudo_random_device.seed(set_seed);
        }
        int operator()() {
            return true_random?static_cast<int>(true_random_device()&INT32_MAX):
                               static_cast<int>(pseudo_random_device()&INT32_MAX);
        }
    };

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
            return (l.exponent < r.exponent) || (l.exponent == r.exponent && l.mantissa+0.01 < r.mantissa);
        }

        friend bool operator==(const big_number& l, const big_number& r)
        {
            return (l.exponent == r.exponent) && (l.mantissa > r.mantissa-0.01) && (l.mantissa < r.mantissa+0.01);
        }

        void set(long double set_mantissa, int set_exponent) {
            mantissa = set_mantissa;
            exponent = set_exponent;
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

    inline std::ostream& operator<<(std::ostream& out, big_number number) {
        return out << number.mantissa << "*10^" << number.exponent;
    }

    static void progress_print_split() {
        PRINT("\r______________________________________________________________");
    }

    static void progress_print_header() {
        progress_print_split();
        PRINT(std::setw(11) << std::left <<"T (ms)" << std::setw(11) << "delta(ms)" << std::setw(12) << "proc"
              << std::setw(16) << "p1"        << std::setw(16)        << "p2");
        progress_print_split();
    }


    class timed_print {
        std::chrono::high_resolution_clock::time_point first;
        std::chrono::high_resolution_clock::time_point previous;
    public:

        bool h_silent = false;

        timed_print() {
            first     = std::chrono::high_resolution_clock::now();
            previous  = first;
        }

        void print_header() const {
            if(h_silent) return;
            progress_print_header();
        }

        void print_split() {
            if(h_silent) return;
            PRINT("\r______________________________________________________________");
        }

        void print(const std::string& str) const {
            if(h_silent) return;
            PRINT("\r" << str);
        }

        void timer_print(const std::string& proc, const std::string& p1, const std::string& p2) {
            if(h_silent) return;
            auto now = std::chrono::high_resolution_clock::now();
            PRINT("\r" << std::fixed << std::setprecision(2) << std::setw(11) << std::left
                       << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - first).count()) / 1000000.0
                       << std::setw(11)
                       << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - previous).count()) / 1000000.0
                       << std::setw(12) << proc << std::setw(16) << p1 << std::setw(16) << p2);
            previous = now;
        }

        void timer_split() {
            previous = std::chrono::high_resolution_clock::now();
        }

        void timer_print(const std::string& proc, const int p1, const int p2) {
            if(h_silent) return;
            auto now = std::chrono::high_resolution_clock::now();
            PRINT("\r" << std::fixed << std::setprecision(2) << std::setw(11) << std::left
                       << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - first).count()) / 1000000.0
                       << std::setw(11)
                       << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - previous).count()) / 1000000.0
                       << std::setw(12) << proc << std::setw(16) << p1 << std::setw(16) << p2);
            previous = now;
        }

        void timer_print(const std::string& proc, const int p1, const double p2) {
            if(h_silent) return;
            auto now = std::chrono::high_resolution_clock::now();
            PRINT("\r" << std::fixed << std::setprecision(2) << std::setw(11) << std::left
                       << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - first).count()) / 1000000.0
                       << std::setw(11)
                       << (std::chrono::duration_cast<std::chrono::nanoseconds>(now - previous).count()) / 1000000.0
                       << std::setw(12) << proc << std::setw(16) << p1 << std::setw(16) << p2);
            previous = now;
        }

        void progress_current_method(const std::string& print) const  {
            if(h_silent) return;
            PRINT_NO_NEWLINE("\r>" << print);
        }
        void progress_current_method(const std::string& method_name, const std::string& var1, double var1_val,
                                            const std::string& var2, double var2_val) const  {
            if(h_silent) return;
            PRINT_NO_NEWLINE("\r>" << method_name << " " << var1 << "=" << var1_val << ", " << var2 << "=" << var2_val);
        }
        void progress_current_method(const std::string& method_name, const std::string& var1, int var1_val,
                                            const std::string& var2, int var2_val,
                                            const std::string& var3, double var3_val) const {
            if(h_silent) return;
            PRINT_NO_NEWLINE("\r>" << method_name << " "  << var1 << "=" << var1_val << ", " << var2 << "=" << var2_val
                                   << ", " << var3 << "=" << var3_val);
        }

        void progress_current_method(const std::string& method_name, const std::string& var1, double var1_val,
                                            const std::string& var2, int var2_val,
                                            const std::string& var3, int var3_val,
                                            const std::string& var4, int var4_val) const  {
            if(h_silent) return;
            PRINT_NO_NEWLINE("\r>" << method_name << " "  << var1 << "=" << var1_val << ", " << var2 << "=" << var2_val
                                   << ", " << var3 << "=" << var3_val << ", " << var4 << "=" << var4_val);
        }
    };
}

#endif //DEJAVU_UTILITY_H
