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
#define min(x, y) (x<y?x:y)
#define max(x, y) (x<y?y:x)
#define SWAP(x, y) { const T a = min(arr[x], arr[y]); \
                    const T b = max(arr[x], arr[y]); \
                    arr[x] = a; arr[y] = b; }
            switch (sz) {
                case 0:
                case 1:
                    break;
                case 2: SWAP(0, 1);
                    break;
                case 3: SWAP(0, 1);
                    SWAP(0, 2);
                    SWAP(1, 2);
                    break;
                case 4: SWAP(0, 1);
                    SWAP(2, 3);
                    SWAP(0, 2);
                    SWAP(1, 3);
                    SWAP(1, 2);
                    break;
                case 5: SWAP(0, 1);
                    SWAP(2, 3);
                    SWAP(0, 2);
                    SWAP(1, 4);
                    SWAP(0, 1);
                    SWAP(2, 3);
                    SWAP(1, 2);
                    SWAP(3, 4);
                    SWAP(2, 3);
                    break;
                case 6: SWAP(1, 2);
                    SWAP(4, 5);
                    SWAP(0, 2);
                    SWAP(3, 5);
                    SWAP(0, 1);
                    SWAP(3, 4);
                    SWAP(1, 4);
                    SWAP(0, 3);
                    SWAP(2, 5);
                    SWAP(1, 3);
                    SWAP(2, 4);
                    SWAP(2, 3);
                    break;
                default:
                    std::sort(arr, arr + sz);
            }
#undef SWAP
#undef min
#undef max
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
         * \brief Fixed-size array
         *
         * An array of fixed size, with some further convenience functions.
         *
         * @tparam T The type of array elements.
         */
        template<class T>
        class work_list_t {
        private:
            /**
             * Allocate an array of size \p size.
             * @param size Space to allocate.
             */
            void alloc(const unsigned int size) {
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
            work_list_t() = default;

            /**
             * Constructor that allocates the internal array with size \p size. The allocated memory is not
             * initialized.
             *
             * @param size Size to allocate.
             */
            explicit work_list_t(unsigned int size) {
                allocate(size);
            }

            /**
             * Allocates the internal array with size \p size. The allocated memory is not
             * initialized. Initializes the internal position \a cur_pos of the array at 0.
             *
             * @param size Size to allocate.
             */
            void allocate(unsigned int size) {
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
            void resize(const unsigned int size) {
                if (arr && size <= arr_sz) return;
                T *old_arr = nullptr;
                unsigned int old_arr_sz = arr_sz;
                if (arr) old_arr = arr;
                alloc(size);
                if (old_arr != nullptr) {
                    int cp_pt = std::min(old_arr_sz, arr_sz);
                    memcpy(arr, old_arr, cp_pt * sizeof(T));
                    free(old_arr);
                }
            }

            /**
             * Deallocates the internal array.
             */
            ~work_list_t() {
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
            unsigned int arr_sz = 0;      /**< size to which \a arr is currently allocated*/
            T *arr     = nullptr; /**< internal array */
        };

        // frequently used types
        typedef work_list_t<int> work_list;
        typedef work_list_t<std::pair<std::pair<int, int>, bool>> work_list_pair_bool;

        // queue with fixed size limitation
        class work_queue {
        public:
            void initialize(int size) {
                assert(!init);
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
        class work_set_t {
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
                while (!reset_queue.empty())
                    s[reset_queue.pop()] = -1;
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

            void inc_nr(int index) {
                assert(index >= 0 && index < sz);
                ++s[index];
            }

            ~work_set_t() {
                if (init)
                    delete[] s;
            }

            work_queue reset_queue;
        private:
            bool init = false;
            T   *s = nullptr;
            int sz = 0;
        };

        typedef work_set_t<int> work_set_int;

        /**
         * \brief Sets specialized for quick resets
         *
         * Set on a statically specified domain of elements 1, ..., size, with O(1) \a set and \a get.
         */
        class mark_set {
            int *s   = nullptr;
            int mark = 0;
            unsigned int sz = 0;
        public:
            mark_set() = default;
            explicit mark_set(int size) {
                initialize(size);
            }

            void initialize(unsigned int size) {
                if(s && sz == size) return;
                if(s) free(s);
                s = (int*) calloc(size, sizeof(int));
                sz   = size;
                mark = 0;
                reset();
            }

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
                if(s) free(s);
            }
        };
    }
}

#endif //DEJAVU_DS_H
