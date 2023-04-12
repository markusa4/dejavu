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

// set specialized for quick resets
class mark_set {
    int mark = 0;
    int *s = nullptr;
    int sz = -1;
    bool init = false;
public:
    mark_set() {};
    mark_set(int size) {
        initialize(size);
    }

    void initialize(int size) {
        //s = new int[size];
        s = (int*) calloc(size, sizeof(int));
        sz = size;
        init = true;
        //memset(s, mark, sz * sizeof(int));
        reset();
    }
    /*void initialize_from_array(int* arr, int size) {
        s  = arr;
        sz = size;
        init = false;
        memset(s, mark, sz * sizeof(int));
        reset();
    }*/
    bool get(int pos) {
        return s[pos] == mark;
    }
    void set(int pos) {
        s[pos] = mark;
    }
    void unset(int pos) {
        s[pos] = mark - 1;
    }
    void reset() {
        if(mark == -1) {
            memset(s, mark, sz * sizeof(int));
        }
        ++mark;
    }
    ~mark_set() {
        /*if(init)
            delete[] s;*/
        if(init)
            free(s);
    }
};

template<class T>
class concurrent_queue {
    std::unique_ptr<std::mutex> lock;
    std::deque<T> content;
public:
    concurrent_queue() {
        lock =  std::make_unique<std::mutex>();
    }

    void enqueue(T item) {
        lock->lock();
        content.emplace_back(item);
        lock->unlock();
    }

    void enqueue_bulk(T* items, int size) {
        lock->lock();
        for(int i = 0; i < size; ++i) {
            content.emplace_back(items[i]);
        }
        lock->unlock();
    }

    T dequeue() {
        lock->lock();
        T item = content.front();
        content.pop_front();
        lock->unlock();
        return item;
    }

    bool try_dequeue(T& item)  {
        int i = 0;
        lock->lock();
        if(content.size() > 0) {
            item = content.front();
            content.pop_front();
            ++i;
        }
        lock->unlock();
        return (i > 0);
    }

    int try_dequeue_bulk(T* arr, int chunk_size)  {
        int i = 0;
        lock->lock();
        while(chunk_size - i > 0 && content.size() > 0) {
            arr[i] = content.front();
            content.pop_front();
            ++i;
        }
        lock->unlock();
        return i;
    }

    void clear() {
        lock->lock();
        content.clear();
        lock->unlock();
    }

    int size() {
        int sz = 0;
        lock->lock();
        sz = content.size();
        lock->unlock();
        return sz;
    }
};

#endif //DEJAVU_UTILITY_H
