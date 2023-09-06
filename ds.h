// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_DS_H
#define DEJAVU_DS_H

#include <list>
#include <iostream>
#include <cstring>
#include <functional>
#include <mutex>
#include <cassert>
#include "coloring.h"

namespace dejavu {

    /**
     * \brief General-purpose datastructures.
     *
     */
    namespace ds {

        /**
         * Sorting utilizing minimal sorting networks for arrays of size `sz <= 6`, falling back to `std::sort` for
         * `sz > 0`.
         *
         * @tparam T Template parameter for the type of array elements.
         * @param arr Array of elements of type \p T.
         * @param sz Length of the array \p arr.
         */
        template<class T>
        void inline sort_t(T *arr, int sz) {
            std::sort(arr, arr + sz);
        }

        /**
         * \brief Shared queue datastructure
         *
         * A queue which can be accessed across multiple threads.
         *
         * @tparam T Type of elements stored in the queue.
         */
        template<class T>
        class shared_queue_t {
        private:
            std::mutex lock;
            std::vector<T> queue;
        public:
            /**
             * Add an element \p t to the queue.
             * @param t Element to be added.
             */
            void add(T &t) {
                lock.lock();
                queue.emplace_back(t);
                lock.unlock();
            }

            /**
             * Reserve space for at least \n elements.
             * @param n Space to be reserved.
             */
            void reserve(int n) {
                lock.lock();
                queue.reserve(n);
                lock.unlock();
            }

            /**
             * @return Whether the queue is empty.
             */
            bool empty() {
                return queue.empty();
            }

            /**
             * Pop an element from the queue and return it.
             *
             * @return The element popped from the queue.
             */
            T pop() {
                lock.lock();
                auto element = queue.back();
                queue.pop_back();
                lock.unlock();
                return element;
            }
        };

        /**
         * \brief Fixed-size array, uninitialized
         *
         * An array of fixed size, with some further convenience functions.
         *
         * @tparam T The type of array elements.
         */
        template<class T>
        class worklist_t {
        private:
            /**
             * Allocate an array of size \p size.
             * @param size Space to allocate.
             */
            void alloc(const int size) {
                dealloc();
                //arr = (T *) malloc(sizeof(T) * size);
                arr = new T[size];
                arr_sz = size;
            }

            void dealloc() {
                //if(arr) free(arr);
                if(arr) delete[] arr;
                arr = nullptr;
            }

        public:
            /**
             * Default constructor, does not allocate any memory.
             */
            worklist_t() = default;

            /**
             * Constructor that allocates the internal array with size \p size. The allocated memory is not
             * initialized.
             *
             * @param size Size to allocate.
             */
            explicit worklist_t(int size) {
                assert(size >= 0);
                allocate(size);
            }

            void copy(worklist_t<T>* other) {
                alloc(other->arr_sz);
                for(int i = 0; i < other->arr_sz; ++i) {
                    arr[i] = other->arr[i];
                }
                arr_sz  = other->arr_sz;
                cur_pos = other->cur_pos;
            }

            /**
             * Allocates the internal array with size \p size. The allocated memory is not
             * initialized. Initializes the internal position \a cur_pos of the array at 0.
             *
             * @param size Size to allocate.
             */
            void allocate(int size) {
                assert(size >= 0);
                alloc(size);
                cur_pos = 0;
            }

            /**
             * Push back an element at position \a cur_pos. Increments the internal position \a cur_pos.
             *
             * @param value Element to push back.
             */
            inline void push_back(T value) {
                assert(cur_pos >= 0 && cur_pos < arr_sz);
                arr[cur_pos] = value;
                cur_pos += 1;
            }

            /**
             * Pop an element at position \a cur_pos. Decreases the internal position \a cur_pos. There is no safeguard
             * in place if `cur_pos <= 0`.
             *
             * @return Element popped at \a cur_pos.
             *
             * \sa The function empty() tests whether `cur_pos == 0`.
             */
            inline T pop_back() {
                assert(cur_pos > 0);
                return arr[--cur_pos];
            }

            /**
             * @return Element at \a cur_pos.
             */
            inline T *last() const {
                return &arr[cur_pos - 1];
            }

            /**
             * @return Whether `cur_pos == 0`.
             */
            [[nodiscard]] bool empty() const {
                return cur_pos == 0;
            }

            /**
             * Sets \a cur_pos to \p size.
             *
             * @param size Value to set \a cur_pos to.
             */
            void set_size(const int size) {
                cur_pos = size;
            }

            /**
             * @return The current position \a cur_pos.
             */
            [[nodiscard]] int size() const {
                return cur_pos;
            }

            /**
             * Sets \a cur_pos to `0`.
             */
            inline void reset() {
                cur_pos = 0;
            }

            /**
             * Resizes the internal array to `size`. Copies the contents of the old array into the new one. If the new
             * size is larger than the old one, the new space is only allocated and not initialized.
             *
             * @param size New size to allocate the array to.
             */
            void resize(const int size) {
                assert(size >= 0);
                if (arr && size <= arr_sz) return;
                T *old_arr = nullptr;
                int old_arr_sz = arr_sz;
                if (arr) old_arr = arr;
                alloc(size);
                if (old_arr != nullptr) {
                    int cp_pt = std::min(old_arr_sz, arr_sz);
                    memcpy(arr, old_arr, cp_pt * sizeof(T));
                    delete[] old_arr;
                }
            }

            /**
             * Deallocates the internal array.
             */
            ~worklist_t() {
                dealloc();
            }

            /**
             * Sort the internal array up to position \a cur_pos.
             */
            void sort() {
                sort_t<T>(arr, cur_pos);
            }

            /**
             * @return A pointer to the internal memory.
             */
            inline T *get_array() const {
                return arr;
            }

            /**
             * Access element \p index in the internal array \p arr.
             *
             * @param index Index of the internal array.
             * @return The element `arr[index]`.
             */
            inline T &operator[](int index) const {
                assert(index >= 0);
                assert(index < arr_sz);
                return arr[index];
            }

            void sort_after_map(T *map) {
                struct comparator_map {
                    T *map;

                    explicit comparator_map(T *map) {
                        this->map = map;
                    }

                    bool operator()(const T &a, const T &b) {
                        return map[a] < map[b];
                    }
                };
                comparator_map c = comparator_map(map);
                std::sort(arr, arr + cur_pos, c);
            }

            int cur_pos = 0; /**< current position */
        private:
            int arr_sz = 0;       /**< size to which \a arr is currently allocated*/
            T *arr     = nullptr; /**< internal array */
        };

        typedef worklist_t<int> worklist;

        /**
         * \brief Fixed-size array, 0-initialized
         *
         * An array of fixed size, with some further convenience functions.
         *
         */
        class workspace {
        private:
            /**
             * Allocate an array of size \p size.
             * @param size Space to allocate.
             */
            void alloc(const int size) {
                assert(size >= 0);
                dealloc();
                arr = (int*) calloc(size, sizeof(int));
                arr_sz = size;
            }

            void dealloc() {
                if(arr) free(arr);
                arr = nullptr;
            }

        public:
            /**
             * Default constructor, does not allocate any memory.
             */
            workspace() = default;

            /**
             * Constructor that allocates the internal array with size \p size. The allocated memory is not
             * initialized.
             *
             * @param size Size to allocate.
             */
            explicit workspace(int size) {
                allocate(size);
            }

            void copy(workspace* other) {
                alloc(other->arr_sz);
                for(int i = 0; i < other->arr_sz; ++i) {
                    arr[i] = other->arr[i]; // TODO use memcpy
                }
                arr_sz  = other->arr_sz;
            }

            /**
             * Allocates the internal array with size \p size. The allocated memory is not
             * initialized. Initializes the internal position \a cur_pos of the array at 0.
             *
             * @param size Size to allocate.
             */
            void allocate(int size) {
                alloc(size);
            }

            /**
             * Sets the entire array to 0.
             */
            inline void reset() {
                memset(arr, 0, arr_sz * sizeof(int));
            }

            /**
             * Resizes the internal array to `size`. Copies the contents of the old array into the new one. If the new
             * size is larger than the old one, the new space is only allocated and not initialized.
             *
             * @param size New size to allocate the array to.
             */
            void resize(const int size) {
                assert(size >= 0);
                if (arr && size <= arr_sz) return;
                int *old_arr = arr;
                arr = nullptr;
                int old_arr_sz = arr_sz;
                alloc(size);
                if (old_arr) {
                    int cp_pt = std::min(old_arr_sz, arr_sz);
                    memcpy(arr, old_arr, cp_pt * sizeof(int));
                    free(old_arr);
                }
            }

            /**
             * Deallocates the internal array.
             */
            ~workspace() {
                dealloc();
            }

            /**
             * @return A pointer to the internal memory.
             */
            [[nodiscard]] inline int* get_array() const {
                return arr;
            }

            /**
             * Access element \p index in the internal array \p arr.
             *
             * @param index Index of the internal array.
             * @return The element `arr[index]`.
             */
            inline int &operator[](int index) const {
                assert(index >= 0);
                assert(index < arr_sz);
                return arr[index];
            }

        private:
            int arr_sz = 0;     /**< size to which \a arr is currently allocated*/
            int *arr = nullptr; /**< internal array */
        };

        // queue with fixed size limitation
        class work_queue {
        public:
            void initialize(int size) {
                if(init) delete[] queue;
                sz = size;
                pos = 0;
                queue = new int[size];
                init = true;
            }

            void push(int val) {
                assert(pos != sz);
                queue[pos] = val;
                pos++;
            }

            int pop() {
                assert(pos > 0);
                pos--;
                return queue[pos];
            }

            [[nodiscard]] bool empty() const {
                return (pos == 0);
            }

            ~work_queue() {
                if (init) {
                    delete[] queue;
                }
            }

            void reset() {
                pos = 0;
            }

        private:
            int *queue = nullptr;
            int pos    = 0;
            bool init  = false;
            int sz     = 0;
        };

        // work set with arbitrary type
        template<class T>
        class workset_t {
        public:
            void initialize(int size) {
                s = new T[size];
                reset_queue.initialize(size);

                memset(s, -1, size * sizeof(T)); // TODO should use calloc

                init = true;
                sz = size;
            }

            void set(int index, T value) {
                assert(index >= 0);
                assert(index < sz);
                s[index] = value;
            }

            T get(int index) {
                assert(index >= 0);
                assert(index < sz);
                return s[index];
            }

            void reset() {
                while (!reset_queue.empty()) s[reset_queue.pop()] = -1;
            }

            void reset_hard() {
                memset(s, -1, sz * sizeof(T));
                reset_queue.reset();
            }

            T inc(int index) {
                assert(index >= 0);
                assert(index < sz);
                if (s[index]++ == -1)
                    reset_queue.push(index);
                return s[index];
            }

            void inline inc_nr(int index) {
                assert(index >= 0 && index < sz);
                ++s[index];
            }

            ~workset_t() {
                if (init)
                    delete[] s;
            }

            work_queue reset_queue;
        private:
            bool init = false;
            T   *s = nullptr;
            int sz = 0;
        };

        typedef workset_t<int> work_set_int;

        /**
         * \brief Sets specialized for quick resets
         *
         * Set on a statically specified domain of elements 1, ..., size, with O(1) \a set and \a get.
         */
        class markset {
            int *s   = nullptr;
            int mark = 0;
            int sz = 0;

            void full_reset() {
                memset(s, mark, sz * sizeof(int));
            }

        public:
            /**
             * Initializes a set of size 0
             */
            markset() = default;

            /**
             * Initialize this set with the given \p size.
             * @param size size to initialize this set to
             */
            explicit markset(int size) {
                initialize(size);
            }

            /**
             * Resizes this set to size \p size, resets the set
             *
             * @param size new size of this set
             */
            void initialize(int size) {
                assert(size >= 0);
                if(s && sz == size) {
                    reset();
                    return;
                }
                if(s) free(s);
                s = (int*) calloc(size, sizeof(int));
                sz   = size;
                mark = 0;
                reset();
            }

            /**
             * @param pos element to check
             * @return Is element \p pos in set?
             */
            inline bool get(int pos) {
                assert(pos >= 0);
                assert(pos < sz);
                return s[pos] == mark;
            }

            /**
             * Adds element \p pos to set
             *
             * @param pos element to set
             */
            inline void set(int pos) {
                assert(pos >= 0);
                assert(pos < sz);
                s[pos] = mark;
            }

            /**
             * Removes element \p pos from set
             * @param pos element to remove
             */
            inline void unset(int pos) {
                assert(pos >= 0);
                assert(pos < sz);
                s[pos] = mark - 1;
            }
            /**
             * Resets this set to the empty set
             */
            void reset() {
                if(mark == -1) full_reset();
                ++mark;
            }

            void copy(markset* other) {
                initialize(other->sz);
                for(int i = 0; i < other->sz; ++i) {
                    s[i] = other->s[i];
                }
                mark = other->mark;
                sz   = other->sz;
            }

            ~markset() {
                if(s) free(s);
            }
        };
    }
}

#endif //DEJAVU_DS_H
