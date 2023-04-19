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
        // sorting utilizing minimal sorting networks for n <= 6
        template<class T>
        void sort_t(T *arr, int sz) {
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

        template<class T>
        class shared_queue_t {
        private:
            std::mutex lock;
            std::vector<T> queue;
        public:
            void add(T &t) {
                lock.lock();
                queue.emplace_back(t);
                lock.unlock();
            }

            void reserve(int n) {
                lock.lock();
                queue.reserve(n);
                lock.unlock();
            }

            bool empty() {
                return queue.empty();
            }

            T pop() {
                lock.lock();
                auto element = queue.back();
                queue.pop_back();
                lock.unlock();
                return element;
            }
        };

        // work list / stack with fixed size limitation
        template<class T>
        class work_list_t {
        private:
            void alloc(const int size) {
                arr = (T *) malloc(sizeof(T) * size);
                arr_sz = size;
            }

        public:
            work_list_t() {};

            work_list_t(int size) {
                initialize(size);
            }

            void initialize(int size) {
                assert(!init);
                //arr     = new T[size];
                alloc(size);
                init = true;
                cur_pos = 0;
            }

            void push_back(T value) {
                assert(cur_pos >= 0 && cur_pos < arr_sz);
                arr[cur_pos] = value;
                cur_pos += 1;
            }

            T pop_back() {
                assert(cur_pos > 0);
                return arr[--cur_pos];
            }

            T *last() {
                return &arr[cur_pos - 1];
            }

            bool empty() {
                return cur_pos == 0;
            }

            void set_size(const int size) {
                cur_pos = size;
            }

            int size() {
                return cur_pos;
            }

            void reset() {
                cur_pos = 0;
            }

            void resize(const int size) {
                if (init && size <= arr_sz) return;
                T *old_arr = nullptr;
                int old_arr_sz = arr_sz;
                if (init) old_arr = arr;
                alloc(size);
                init = true;
                if (old_arr != nullptr) {
                    int cp_pt = std::min(old_arr_sz, arr_sz);
                    memcpy(arr, old_arr, cp_pt * sizeof(T));
                    free(old_arr);
                }
            }

            ~work_list_t() {
                if (init) {
                    //delete[] arr;
                    free(arr);
                }
            }

            void sort() {
                sort_t<T>(arr, cur_pos);
            }

            T *get_array() {
                return arr;
            }

            T &operator[](int index) {
                assert(index >= 0);
                assert(index < arr_sz);
                return arr[index];
            }

            void sort_after_map(T *map) {
                struct comparator_map {
                    T *map;

                    comparator_map(T *map) {
                        this->map = map;
                    }

                    bool operator()(const T &a, const T &b) {
                        return map[a] < map[b];
                    }
                };
                comparator_map c = comparator_map(map);
                std::sort(arr, arr + cur_pos, c);
            }

            int cur_pos = 0;
            bool init = false;
        private:
            T *arr = nullptr;
            int arr_sz = -1;
        };

        typedef work_list_t<std::pair<std::pair<int, int>, bool>> work_list_pair_bool;
        typedef work_list_t<int> work_list;

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
                //assert(init);
                assert(pos != sz);
                queue[pos] = val;
                pos++;
            }

            int pop() {
                //assert(init);
                assert(pos > 0);
                pos--;
                return queue[pos];
            }

            bool empty() {
                //assert(init);
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

            void initialize_from_array(int *arr, int size) {
                assert(!init);
                sz = size;
                pos = 0;
                queue = arr;
                init = false;
            }

        private:
            int *queue;
            int pos;
            bool init = false;
            int sz;
        };

        // work set with arbitrary type
        template<class T>
        class alignas(16) work_set_t {
        public:
            void initialize(int size) {
                s = new T[size];
                reset_queue.initialize(size);

                memset(s, -1, size * sizeof(T)); // TODO should use calloc

                init = true;
                sz = size;
            }

            void initialize_from_array(T *arr, int size) {
                s = arr;
                reset_queue.initialize_from_array(arr + size, size);

                memset(s, -1, size * sizeof(T));

                init = false;
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
            T *s;
            int sz;
        };

        typedef work_set_t<int> work_set_int;
        typedef work_set_t<char> work_set_char;

        // ring queue for pairs of integers
        class ring_pair {
        public:
            void initialize(int size) {
                arr = new std::pair<int, int>[size];
                arr_sz = size;
                back_pos = 0;
                front_pos = 0;
                init = true;
            }

            void push_back(std::pair<int, int> value) {
                arr[back_pos] = value;
                back_pos = (back_pos + 1) % arr_sz;
            }

            std::pair<int, int> *front() {
                return &arr[front_pos];
            }

            void pop() {
                front_pos = (front_pos + 1) % arr_sz;
            }

            bool empty() {
                return (front_pos == back_pos);
            }

            ~ring_pair() {
                if (init)
                    delete[] arr;
            }

            void reset() {
                front_pos = back_pos;
            }

        private:
            std::pair<int, int> *arr;
            bool init = false;
            int arr_sz = -1;
            int front_pos = -1;
            int back_pos = -1;
        };

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
    }
}

#endif //DEJAVU_DS_H
