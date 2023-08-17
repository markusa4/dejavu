// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_GROUPS_H
#define DEJAVU_GROUPS_H

#include "coloring.h"
#include "sgraph.h"
#include "trace.h"

namespace dejavu {

    /**
     * \brief Data structures and algorithms to deal with groups.
     *
     * Contains basic data structures to construct and deal with automorphisms, as well as a Schreier structure.
     */
    namespace groups {

        /**
         * Reset an automorphism given in \p automorphism with support \p support to the identity.
         *
         * @param automorphism Dense notation of the given automorphism, i.e., vertex `i` is mapped to `automorphism[i]`.
         * @param support Support of the automorphism, contains all vertices where `i != automorphism[i]`.
         */
        static void reset_automorphism(int *automorphism, work_list *support) {
            for (int i = 0; i < support->cur_pos; ++i) {
                automorphism[(*support)[i]] = (*support)[i];
            }
            support->reset();
        }

        /**
         * Create an automorphism from two discrete vertex colorings.
         *
         * @param domain_size Size of the underlying domain (i.e., number of vertices of the graph).
         * @param vertex_to_col Vertex-to-color mapping of the first coloring, i.e., vertex `i` is mapped to color
         * `vertex_to_col[i]`.
         * @param col_to_vertex Color-to-vertex mapping of the second coloring, i.e., the color i contains the vertex
         * `col_to_vertex[i]` (since the colorings are assumed to be discrete, every `i` must be a distinct color).
         * @param automorphism Dense notation for the automorphism to be written. Assumed to be the identity upon
         * calling the method.
         * @param support Support for the automorphism \p automorphism.
         */
        static void color_diff_automorphism(int domain_size, const int *vertex_to_col, const int *col_to_vertex,
                                            int *automorphism, work_list *support) {
            for (int v1 = 0; v1 < domain_size; ++v1) {
                const int col = vertex_to_col[v1];
                const int v2  = col_to_vertex[col];
                if (v1 != v2) {
                    automorphism[v1] = v2;
                    support->push_back(v1);
                }
            }
        }

        /**
         * \brief Workspace for sparse automorphisms
         *
         * Enables O(1) lookup on a sparse automorphism by using an O(n) workspace.
         */
        class automorphism_workspace {
            work_list automorphism;
            work_list automorphism_supp;
            int domain_size = 0;

            bool support01 = false;
        public:
            /**
             * Initializes the stored automorphism to the identity.
             *
             * @param domain_size Size of the domain on which automorphisms operate
             */
            explicit automorphism_workspace(int domain_size) {
                automorphism.allocate(domain_size);
                for (int i = 0; i < domain_size; ++i)
                    automorphism[i] = i;
                automorphism_supp.allocate(domain_size);
                this->domain_size = domain_size;
            }

            /**
             * If 0/1 support is activated, this automorphism_workspace will only track whether the support of the
             * underlying automorphism is trivial (=0) or non-trivial (>=1). Some operations are more efficient, if
             * 0/1 support is used (e.g., \a apply)
             *
             * @param new_support01 Flag of whether to use 0/1 support, or not.
             */
            void set_support01(const bool new_support01) {
                this->support01 = new_support01;
            }

            /**
             * Create automorphism from two discrete vertex colorings.
             *
             * @param vertex_to_col Vertex-to-color mapping of the first coloring, i.e., vertex `i` is mapped to color
             * `vertex_to_col[i]`.
             * @param col_to_vertex Color-to-vertex mapping of the second coloring, i.e., the color i contains the vertex
             * `col_to_vertex[i]` (since the colorings are assumed to be discrete, every `i` must be a distinct color).
             */
            void write_color_diff(const int *vertex_to_col, const int *col_to_vertex) {
                color_diff_automorphism(domain_size, vertex_to_col, col_to_vertex, automorphism.get_array(),
                                        &automorphism_supp);
            }

            /**
             * Updates the support using the internal dense notation.
             */
            void update_support() {
                // rewrite support
                if (!support01) {
                    automorphism_supp.reset();
                    for (int i = 0; i < domain_size; ++i) {
                        if (i != automorphism[i])
                            automorphism_supp.push_back(i);
                    }
                } else {
                    automorphism_supp.reset();
                    int i;
                    for (i = 0; i < domain_size; ++i) {
                        if (i != automorphism[i]) break;
                    }
                    automorphism_supp.cur_pos = (i != domain_size);
                }
            }

            /**
             * Apply another automorphism to the stored automorphism. To be more precise, `other^pwr` is applied to
             * the automorphism stored in this object.
             *
             * Closely follows the implementation in nauty / Traces.
             *
             * @param scratch_apply1 Auxiliary workspace used for the operation.
             * @param scratch_apply2 Auxiliary workspace used for the operation.
             * @param scratch_apply3 Auxiliary workspace used for the operation.
             * @param other Automorphism in sparse notation that is applied to this automorphism in.
             * @param pwr Power with which \p other is applied to this automorphism.
             */
            [[maybe_unused]] void apply(work_list &scratch_apply1, work_list &scratch_apply2, mark_set &scratch_apply3,
                       automorphism_workspace *other, int pwr = 1) {
                apply(scratch_apply1, scratch_apply2, scratch_apply3, other->perm(), pwr);
            }

            /**
             * Apply another automorphism to the stored automorphism. To be more precise, `other^pwr` is applied to
             * the automorphism stored in this object.
             *
             * Closely follows the implementation in nauty / Traces.
             *
             * @param scratch_apply1 Auxiliary workspace used for the operation.
             * @param scratch_apply2 Auxiliary workspace used for the operation.
             * @param scratch_apply3 Auxiliary workspace used for the operation.
             * @param other Automorphism in dense notation that is applied to this automorphism.
             * @param pwr Power with which \p other is applied to this automorphism.
             */
            void apply(work_list &scratch_apply1, work_list &scratch_apply2, mark_set &scratch_apply3,
                  const int *p, int pwr = 1) {
                if (pwr == 0) return;
                if (pwr <= 5) {
                    if (pwr == 1)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[automorphism[i]];
                    else if (pwr == 2)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[automorphism[i]]];
                    else if (pwr == 3)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[p[automorphism[i]]]];
                    else if (pwr == 4)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[p[p[automorphism[i]]]]];
                    else if (pwr == 5)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[p[p[p[automorphism[i]]]]]];
                } else if (pwr <= 19) {
                    // apply other automorphism
                    for (int j = 0; j < domain_size; ++j) {
                        scratch_apply1[j] = p[p[p[j]]];
                    }
                    for (; pwr >= 6; pwr -= 6)
                        for (int j = 0; j < domain_size; ++j)
                            automorphism[j] = scratch_apply1[scratch_apply1[automorphism[j]]];

                    if (pwr == 1)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[automorphism[i]];
                    else if (pwr == 2)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[automorphism[i]]];
                    else if (pwr == 3)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = scratch_apply1[automorphism[i]];
                    else if (pwr == 4)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[scratch_apply1[automorphism[i]]];
                    else if (pwr == 5)
                        for (int i = 0; i < domain_size; ++i) automorphism[i] = p[p[scratch_apply1[automorphism[i]]]];
                } else {
                    // 1 cycle at a time

                    scratch_apply3.reset();
                    for (int i = 0; i < domain_size; ++i) {
                        if (scratch_apply3.get(i)) continue;
                        if (p[i] == i)
                            scratch_apply2[i] = i;
                        else {
                            int cyclen = 1;
                            scratch_apply1[0] = i;
                            for (int j = p[i]; j != i; j = p[j]) {
                                scratch_apply1[cyclen++] = j;
                                scratch_apply3.set(j);
                            }
                            int kk = pwr % cyclen;
                            for (int j = 0; j < cyclen; ++j) {
                                scratch_apply2[scratch_apply1[j]] = scratch_apply1[kk];
                                if (++kk == cyclen) kk = 0;
                            }
                        }
                    }
                    for (int i = 0; i < domain_size; ++i) automorphism[i] = scratch_apply2[automorphism[i]];
                    scratch_apply3.reset();
                }

                // rewrite support
                update_support();
            }


            /**
             * Create automorphism from two canonically-ordered vectors of singletons. The resulting automorphism maps
             * `singletons1[i]` to `singletons2[i]` for `i` in `pos_start, ..., pos_end`.
             *
             * @param singletons1 first vector of singletons
             * @param singletons2 second vector of singletons
             * @param pos_start start reading the vectors at this position
             * @param pos_end stop reading the vecvtors at this position.
             */
            void __attribute__((noinline))  write_singleton(const std::vector<int> *singletons1, const std::vector<int> *singletons2,
                                 const int pos_start, const int pos_end) {
                for (int i = pos_start; i < pos_end; ++i) {
                    const int from = (*singletons1)[i];
                    const int to   = (*singletons2)[i];
                    assert(automorphism[from] == from);
                    if (from != to) {
                        automorphism_supp.push_back(from);
                        automorphism[from] = to;
                    }
                }
            }

            void write_single_map(const int from, const int to) {
                assert(automorphism[from] == from);
                if (from != to) {
                    automorphism_supp.push_back(from);
                    automorphism[from] = to;
                }
            }

            /**
             * Reset the contained automorphism back to the identity.
             */
            void reset() {
                reset_automorphism(automorphism.get_array(), &automorphism_supp);
            }

            /**
             * @return Integer array \p p describing the stored automorphism, where point v is mapped to \p p[v].
             */
            [[nodiscard]] int *perm() const {
                return automorphism.get_array();
            }

            /**
             * @return Integer array which contains all vertices in the support of the contained automorphism.
             */
            [[nodiscard]] int* support() const {
                return automorphism_supp.get_array();
            }

            /**
             * @return Size of the support.
             */
            [[nodiscard]] int nsupport() const {
                return automorphism_supp.cur_pos;
            }

            /** Set support (unsafe).
             */
            void set_nsupport(int nsupp) {
                automorphism_supp.set_size(nsupp);
            }
        };

        /**
         * \brief Orbit partition
         *
         * Keeps track of an orbit partition, and provides tools to manipulate the orbits within.
         */
        class orbit {
            int        sz = 0;
            mark_set   touched;
            work_list  reset_arr;
            work_list  map_arr;
            work_list  orb_sz;
        public:

            /**
             * Retrieve the orbit of the given vertex.
             *
             * @param v Vertex of the specified domain.
             * @return The orbit of \p v.
             */
            int find_orbit(const int v) {
                assert(v >= 0);
                assert(v < sz);
                int last_map;
                int next_map = v;
                do {
                    last_map = next_map;
                    next_map = map_arr[last_map];
                } while(next_map != last_map);
                map_arr[v] = next_map;
                return next_map;
            }

            /**
             * Returns the size of an orbit.
             *
             * @param v Vertex of the specified domain.
             * @return Size of the orbit of \p v.
             */
            int orbit_size(const int v) {
                assert(v >= 0);
                assert(v < sz);
                return orb_sz[find_orbit(v)];
            }

            /**
             * Every orbit has precisely one representative. This function enables to test this.
             *
             * @param v Vertex of the specified domain.
             * @return Whether \p v is the representative of the orbit of \p v.
             */
            bool represents_orbit(const int v) {
                return v == map_arr[v];
            }

            /**
             * Combines the orbits of two given vertices.
             *
             * @param v1 The first vertex.
             * @param v2  The second vertex.
             */
            void combine_orbits(const int v1, const int v2) {
                assert(v1 >= 0);
                assert(v2 >= 0);
                assert(v1 < sz);
                assert(v2 < sz);
                if(v1 != v2) {
                    if(!touched.get(v1)) reset_arr.push_back(v1);
                    if(!touched.get(v2)) reset_arr.push_back(v2);
                    touched.set(v1);
                    touched.set(v2);
                    int orbit1 = find_orbit(v1);
                    int orbit2 = find_orbit(v2);
                    if(orbit1 == orbit2) return;
                    if(orbit1 < orbit2) {
                        map_arr[orbit2] = orbit1;
                        orb_sz[orbit1] += orb_sz[orbit2];
                    } else {
                        map_arr[orbit1] = orbit2;
                        orb_sz[orbit2] += orb_sz[orbit1];
                    }
                }
            }


            /**
             * Checks whether two given vertices are in the same orbit
             *
             * @param v1 The first vertex.
             * @param v2  The second vertex.
             * @return Whether \p v1 and \p v2 are in the same orbit.
             */
            bool are_in_same_orbit(const int v1, const int v2) {
                assert(v1 >= 0);
                assert(v2 >= 0);
                assert(v1 < sz);
                assert(v2 < sz);
                if(v1 == v2) return true;
                const int orbit1 = find_orbit(v1);
                const int orbit2 = find_orbit(v2);
                return (orbit1 == orbit2);
            }

            void reset() {
                while(!reset_arr.empty()) {
                    const int v = reset_arr.pop_back();
                    map_arr[v] = v;
                    orb_sz[v]  = 1;
                }
                touched.reset();
            }

            /**
             * Applies an automorphism to the orbit structure.
             *
             * @param aut Automorphism workspace which is applied.
             */
            void __attribute__((noinline))  add_automorphism_to_orbit(groups::automorphism_workspace& aut) {
                const int  nsupp = aut.nsupport();
                const int* supp  = aut.support();
                const int* p     = aut.perm();
                for (int i = 0; i < nsupp; ++i) {
                    combine_orbits(p[supp[i]], supp[i]);
                }
            }

            /**
             * Initializes the orbit structure with the given size.
             *
             * @param domain_size size of the underlying domain
             */
            void initialize(int domain_size) {
                sz = domain_size;
                touched.initialize(domain_size);
                reset_arr.allocate(domain_size);
                map_arr.allocate(domain_size);
                orb_sz.allocate(domain_size);
                for(int i = 0; i < domain_size; ++i) {
                    map_arr.push_back(i);
                    orb_sz.push_back(1);
                }
            }

            /**
             * Initializes orbit with the given domain size.
             *
             * @param domain_size size of the underlying domain
             */
            explicit orbit(int domain_size) {
                initialize(domain_size);
            }

            orbit() = default;
        };

        /**
         * \brief Stores a link to an automorphism.
         *
         * Used to implement indiscriminate loading of dense and sparse automorphisms.
         */
        class dense_sparse_arbiter {
            int *loaded_automorphism = nullptr;
        public:
            void load(int *automorphism) {
                loaded_automorphism = automorphism;
            }

            int *p() {
                return loaded_automorphism;
            }
        };

        /**
         * \brief Stores an automorphism in a dense or sparse manner, dynamically.
         */
        class stored_automorphism {
        public:
            enum stored_automorphism_type { STORE_DENSE, ///< stored densely, in size O(\a domain_size)
                                            STORE_SPARSE ///< stored in minimal encoding size of automorphism
                                          };

        private:
            work_list data;
            int domain_size = 0;
            stored_automorphism_type store_type = STORE_SPARSE;

        public:
            /**
             * @return whether the automorphism is stored in a dense or sparse manner
             */
            stored_automorphism_type get_store_type() {
                return store_type;
            }

            /**
             * Load this stored automorphism into workspace.
             *
             * @param loader Store a pointer to permutation of automorphism in the arbiter.
             * @param space Auxiliary space that may or may not be used, depending on whether the loaded automorphism is
             *              sparse. Should only be reset once \p loader is not used anymore.
             */
            void load(dense_sparse_arbiter &loader, automorphism_workspace &space) {
                if (store_type == STORE_SPARSE) {
                    space.reset();
                    int first_of_cycle = 0;
                    for (int i = 0; i < data.size(); ++i) {
                        const int j = abs(data[i]) - 1;
                        const bool is_last = data[i] < 0;
                        assert(i == data.size() - 1 ? is_last : true);
                        if (is_last) {
                            space.write_single_map(j, abs(data[first_of_cycle]) - 1);
                            first_of_cycle = i + 1;
                        } else {
                            assert(i + 1 < data.size());
                            space.write_single_map(j, abs(data[i + 1]) - 1);
                        }
                    }
                    loader.load(space.perm());
                } else {
                    // store_type == STORE_DENSE
                    loader.load(data.get_array());
                }
            }

            /**
             * Load inverse of stored automorphism into workspace.
             *
             * @param loader Store a pointer to permutation of automorphism in the arbiter.
             * @param space Auxiliary space that may or may not be used, depending on whether the loaded automorphism is
             *              sparse. Should only be reset once \p loader is not used anymore.
             */
            void load_inverse(dense_sparse_arbiter &loader, automorphism_workspace &space) {
                if (store_type == STORE_SPARSE) {
                    space.reset();
                    int first_of_cycle = data.size() - 1;
                    for (int i = data.size() - 1; i >= 0; --i) {
                        if (data[i] < 0) first_of_cycle = i;

                        const int j = abs(data[i]) - 1;
                        const bool is_last = i == 0 || (data[i - 1] < 0);
                        if (is_last) {
                            space.write_single_map(j, abs(data[first_of_cycle]) - 1);
                        } else {
                            space.write_single_map(j, abs(data[i - 1]) - 1);
                        }
                    }
                    loader.load(space.perm());
                } else {
                    space.reset();
                    for (int i = 0; i < domain_size; ++i) {
                        if (i != data[i]) {
                            space.write_single_map(data[i], i);
                        }
                    }
                    loader.load(space.perm());
                }
            }

            /**
             * Store the given automorphism workspace.
             *
             * @param automorphism The automorphism to be stored.
             */
            void store(int new_domain_size, automorphism_workspace &automorphism, mark_set &helper) {
                domain_size = new_domain_size;
                assert(data.empty());

                int support = 0;
                for (int i = 0; i < domain_size; ++i) support += (automorphism.perm()[i] != i);

                // decide whether to store dense or sparse representation
                if (support < domain_size / 4) {
                    store_type = STORE_SPARSE;
                    helper.reset();

                    data.allocate(support);
                    for (int i = 0; i < domain_size; ++i) {
                        if (automorphism.perm()[i] == i) continue;
                        const int j = i;
                        if (helper.get(j)) continue;
                        helper.set(j);
                        int map_j = automorphism.perm()[j];
                        assert(map_j != j);
                        while (!helper.get(map_j)) {
                            data.push_back(map_j + 1);
                            helper.set(map_j);
                            map_j = automorphism.perm()[map_j];
                        }
                        assert(map_j == j);
                        data.push_back(-(j + 1));
                    }
                    helper.reset();
                    assert(data.size() == support);
                } else {
                    store_type = STORE_DENSE;
                    data.allocate(domain_size);
                    data.set_size(domain_size);
                    memcpy(data.get_array(), automorphism.perm(), domain_size * sizeof(int));
                    assert(data.size() == domain_size);
                }
            }
        };

        /**
         * \brief Auxiliary workspace used for Schreier computations
         *
         * A global (thread local) state used for computations in Schreier structures. Used such that auxiliary space
         * does not have to be allocated or re-allocated for every single operation, but only once. Also, the same space
         * is used across different operations.
         */
        class schreier_workspace {
        public:
            /**
             * Initialize this workspace.
             *
             * @param domain_size Size of the underlying domain (i.e., number of vertices of the graph).
             */
            explicit schreier_workspace(int new_domain_size) : scratch_auto(new_domain_size) {
                scratch1.initialize(new_domain_size);
                scratch2.initialize(new_domain_size);
                scratch_apply1.allocate(new_domain_size);
                scratch_apply2.allocate(new_domain_size);
                scratch_apply3.initialize(new_domain_size);
            }

            dense_sparse_arbiter loader; /**< used for indiscriminate loading of dense and sparse automorphisms */

            mark_set scratch1;        /**< auxiliary space */
            mark_set scratch2;        /**< auxiliary space */
            work_list scratch_apply1; /**< auxiliary space used for `apply` operations */
            work_list scratch_apply2; /**< auxiliary space used for `apply` operations */
            mark_set  scratch_apply3; /**< auxiliary space used for `apply` operations */
            automorphism_workspace scratch_auto; /**< used to store a sparse automorphism*/
        };

        /**
         * \brief Stores a generating set.
         *
         * Can be used across multiple threads.
         */
        class generating_set {
            std::mutex lock_generators;          /**< locks the generators */
            std::vector<stored_automorphism *> generators; /** list of generators */
            int domain_size = 0;
        public:
            int s_stored_sparse = 0; /**< how many generators are stored in a sparse manner */
            int s_stored_dense = 0;  /**< how many generators are stored in a dense manner */

            /**
             * Set up this generating set.
             * @param domain_size Size of the domain of the stored generators.
             */
            void initialize(int new_domain_size) {
                this->domain_size = new_domain_size;
            }

            /**
             * Add a generator to this generating set.
             *
             * @param w The Schreier workspace.
             * @param automorphism The automorphism to be stored as a generator
             * @return Identifier of the new generator in the generating set.
             */
            int add_generator(schreier_workspace &w, automorphism_workspace &automorphism) {
                lock_generators.lock();
                generators.emplace_back(new stored_automorphism);
                const auto num = generators.size() - 1;
                generators[num]->store(domain_size, automorphism, w.scratch2);
                lock_generators.unlock();

                s_stored_sparse += (generators[num]->get_store_type() ==
                                    stored_automorphism::stored_automorphism_type::STORE_SPARSE);
                s_stored_dense += (generators[num]->get_store_type() ==
                                   stored_automorphism::stored_automorphism_type::STORE_DENSE);

                return static_cast<int>(num);
            }

            void remove_generator(size_t num) {
                assert(num >= 0);
                assert(num < generators.size());
                delete generators[num];
                generators[num] = nullptr;
            }

            /**
             * Remove generators from generating set which do not fix points in \p fix_points.
             *
             * @param w A Schreier workspace.
             * @param fix_points Points to be fixed by generating set.
             */
            void filter(schreier_workspace& w, std::vector<int> &fix_points) {
                for(int test_pt : fix_points) {
                    for(size_t j = 0; j < generators.size(); ++j) {
                        const auto gen = generators[j];
                        if(gen == nullptr) continue;
                        gen->load(w.loader, w.scratch_auto);
                        if(w.loader.p()[test_pt] != test_pt) remove_generator(j);
                        w.scratch_auto.reset();
                    }
                }
                compact_generators();
            }

            /**
             * Remove generators marked for deletion.
             */
            void compact_generators() {
                std::vector<stored_automorphism *> new_gens;
                for(auto gen : generators) {
                    if(gen) new_gens.push_back(gen);
                }
                generators.swap(new_gens);
            }

            /**
             * Retrieve a generator from the stored generating set.
             *
             * @param num Identifier of a generator.
             * @return A pointer to the generator which corresponds to the identifier.
             */
            [[nodiscard]] stored_automorphism* get_generator(const size_t num) const {
                return generators[num];
            }

            stored_automorphism* operator[](const size_t num) const {
                return get_generator(num);
            }

            /**
             * @return number of stored generators
             */
            [[nodiscard]] int size() const {
                return static_cast<int>(generators.size());
            }

            void clear() {
                for(auto & generator : generators) {
                    delete generator;
                }
                generators.clear();
            }

            ~generating_set() {
                clear();
            }
        };

        /**
         * \brief A transversal in a Schreier structure.
         *
         * Can be used across multiple threads.
         */
        class shared_transversal {
            /*enum stored_transversal_type {
                STORE_DENSE, STORE_SPARSE
            };*/

            std::mutex lock_transversal;          /**< locks this transversal */

            int fixed = -1;                       /**< vertex fixed by this transversal */
            int sz_upb = INT32_MAX;               /**< upper bound for size of the transversal (e.g. color class size) */
            int level = -1;
            bool finished = false;

            std::vector<int> fixed_orbit;         /**< contains vertices of orbit at this schreier level */
            std::vector<int> fixed_orbit_to_perm; /**< maps fixed_orbit[i] to generators[i] in class \ref schreier. */
            std::vector<int> fixed_orbit_to_pwr;  /**< power that should to be applied to generators[i] in class \ref schreier. */

            //stored_transversal_type store_type = STORE_SPARSE; /**< whether above structures are stored dense or sparse */

            /**
             * We load the stored orbit to a given schreier_workspace to unlock O(1) access.
             *
             * @param w The schreier_workspace to which the orbit is loaded.
             */
            void load_orbit_to_scratch(schreier_workspace &w) const {
                w.scratch1.reset();
                for (int p : fixed_orbit) {
                    w.scratch1.set(p);
                }
            }

            /**
             * Add \p vertex to the orbit of vertex \a fixed. Store in the transversal that we begin mapping \p vertex
             * to \a fixed by applying `perm^pwr` (further applications of permutations might be necessary, the
             * transversal stores a Schreier vector).
             *
             * @param vertex Vertex added to orbit of \a fixed.
             * @param perm Corresponding permutation used in Schreier vector stored in the transversal.
             * @param pwr Corresponding power used in Schreier vector stored in the transversal.
             */
            void add_to_fixed_orbit(const int vertex, const int perm, const int pwr) {
                fixed_orbit.push_back(vertex);
                fixed_orbit_to_perm.push_back(perm);
                fixed_orbit_to_pwr.push_back(pwr);
            }

            /**
             * The method performs two operations: first, it loads the generator \p gen_num of the generating set
             * \p generators. We denote the generator with \p gen.
             *
             * Second, it applies \p gen with power \p pwr (i.e., `gen^pwr`) to the permutation stored in \p automorphism.
             *
             *
             * @param w Schreier workspace used as auxiliary space for computations.
             * @param automorphism The automorphism to which the permutation is applied.
             * @param generators A generating set.
             * @param gen_num The number of the generator in \p generators, which is applied to \p automorphism.
             * @param pwr A power that is applied to the generator first.
             */
            static void apply_perm(schreier_workspace &w, automorphism_workspace &automorphism,
                       generating_set &generators, const int gen_num, const int pwr) {
                // load perm into workspace
                auto generator = generators[gen_num];

                // apply generator
                if (pwr < 0) { // use inverse
                    generator->load_inverse(w.loader, w.scratch_auto);
                    // multiply
                    automorphism.apply(w.scratch_apply1, w.scratch_apply2, w.scratch_apply3, w.loader.p(), abs(pwr));
                    w.scratch_auto.reset();
                } else if (pwr > 0) { // use generator
                    generator->load(w.loader, w.scratch_auto);
                    // multiply
                    automorphism.apply(w.scratch_apply1, w.scratch_apply2, w.scratch_apply3, w.loader.p(), pwr);
                    w.scratch_auto.reset();
                }
            }

        public:

            /**
             * @return Size of the transversal.
             */
            [[nodiscard]] int size() const {
                return (int) fixed_orbit.size();
            }

            void set_size_upper_bound(const int new_sz_upb) {
                sz_upb = new_sz_upb;
            }

            [[nodiscard]] int get_size_upper_bound() {
                return sz_upb;
            }
            /**
             * Check whether a point \p p is contained in transversal.
             *
             * @param p Point to check.
             * @return Position of point \p in transversal, or -1 if not contained.
             */
            [[nodiscard]] int find_point(const int p) const {
                for (size_t i = 0; i < fixed_orbit.size(); ++i) {
                    if (p == fixed_orbit[i]) {
                        return (int) i;
                    }
                }
                return -1;
            }

            /**
             * Reduces a vector of vertices \p selection to contain only points not contained in this transversal
             *
             * @param w A Schreier workspace.
             * @param selection Vector to be reduced.
             */
            void reduce_to_unfinished(schreier_workspace &w, std::vector<int> &selection) {
                load_orbit_to_scratch(w);
                int back_swap = (int) selection.size() - 1;
                int front_pt;
                for (front_pt = 0; front_pt <= back_swap;) {
                    if (!w.scratch1.get(selection[front_pt])) {
                        ++front_pt;
                    } else {
                        selection[front_pt] = selection[back_swap];
                        --back_swap;
                    }
                }
                selection.resize(front_pt);
                w.scratch1.reset();
            }

            /**
             * @return Whether size of this transversal matches the given upper bound.
             */
            [[nodiscard]] bool is_finished() const {
                return finished;
            }

            /**
             * @return Point fixed by this transversal.
             */
            [[nodiscard]] int fixed_point() const {
                return fixed;
            }

            /**
             * Initialize this transversal.
             *
             * @param fixed_vertex Vertex fixed by the transversal.
             * @param level Position of the transversal in the base of Schreier structure.
             * @param sz_upb Upper bound for the size of transversal (e.g., color class size in combinatorial base).
             */
            void initialize(const int fixed_vertex, const int new_level, const int new_sz_upb) {
                assert(fixed_vertex >= 0);
                assert(new_level >= 0);
                fixed = fixed_vertex;
                this->level  = new_level;
                this->sz_upb = new_sz_upb;
                add_to_fixed_orbit(fixed_vertex, -1, 0);
            }

            /**
             * Extend this transversal using the given \p automorphism.
             *
             * @param generators Generating set of this Schreier structure.
             * @param automorphism Extend the transversal using this automorphism.
             * @return whether transversal was extended
             */
            bool extend_with_automorphism(schreier_workspace &w, generating_set &generators,
                                          automorphism_workspace &automorphism) {
                if (finished) return false; // Already finished transversal? Nothing to do, then!

                // TODO maybe attempt to not load the orbit sometimes?

                // load orbit of this transversal to our workspace, such that we have random access to points in O(1)
                load_orbit_to_scratch(w);
                bool changed = false; /*< we changed this transversal? */
                int  gen_num = -1;    /*< the generator we added to extend this transversal, -1 means no generator */

                // watch out, we may enlarge fixed_orbit in the loop below
                for (int i = 0; i < static_cast<int>(fixed_orbit.size()); ++i) {
                    const int p = fixed_orbit[i];
                    int j = automorphism.perm()[p];
                    if (j == p || w.scratch1.get(j)) continue;

                    int pwr = 0; // power we save in Schreier vector
                    for (int jj = j; !w.scratch1.get(jj); jj = automorphism.perm()[jj]) ++pwr;

                    while (!w.scratch1.get(j)) {
                        // we change this traversal
                        changed = true;

                        // add generator to generating set (once)
                        if (gen_num == -1) gen_num = generators.add_generator(w, automorphism);

                        // add entry to traversal
                        add_to_fixed_orbit(j, gen_num, pwr);
                        w.scratch1.set(j);

                        // we check out the entire cycle of j now
                        j = automorphism.perm()[j];
                        --pwr;
                    }
                }

                // We reached upper bound for the size of this transversal? Then mark transversal as "finished"!
                if (sz_upb == (int) fixed_orbit.size() && !finished) {
                    finished = true;
                }

                assert(sz_upb >= (int) fixed_orbit.size());

                // reset our workspace
                w.scratch1.reset();
                return changed;
            }

            /**
             * Fix \a fixed in \p automorphism.
             *
             * @param generators The underlying generating set.
             * @param automorphism Automorphism where \a fixed will be fixed.
             * @return whether \p automorphism is now the identity
             */
            bool fix_automorphism(schreier_workspace &w, generating_set &generators,
                                  automorphism_workspace &automorphism) const {
                int fixed_map = automorphism.perm()[fixed]; // Where is fixed mapped to?

                // as long as `fixed` is not yet fixed, we apply automorphisms from the transversal as prescribed by
                // the Schreier vector
                while (fixed != fixed_map) {
                    const int pos = find_point(fixed_map); // Where are we storing the information for `fixed_map`?
                    assert(pos >= 0);
                    const int perm = fixed_orbit_to_perm[pos]; // generator to apply for `fixed_map`
                    const int pwr  = fixed_orbit_to_pwr[pos];  // power to use for `fixed_map`
                    apply_perm(w, automorphism, generators, perm, pwr);
                    fixed_map = automorphism.perm()[fixed]; // Fixed now? Or we need to go again?
                }
                assert(automorphism.perm()[fixed] == fixed);
                return automorphism.nsupport() == 0;
            }


            // TODO: flip to dense over certain threshold -- could also include "semi-dense" with sorted list or hash
        };

        /**
         * \brief Schreier structure with fixed base.
         *
         * Enables sifting of automorphisms into a Schreier structure with fixed base. Can be used across multiple threads
         * in a safe manner, i.e., the structure can lock appropriate parts of itself.
         *
         */
        class shared_schreier {
        private:
            int domain_size    = -1;
            int finished_up_to = -1;

            generating_set generators;
            std::vector<shared_transversal *> transversals;

            bool init = false;

        public:
            int s_consecutive_success = 0;  /**< track repeated sifts for probabilistic abort criterion */
            int h_error_bound         = 10; /**< determines error probability                           */
            big_number s_grp_sz; /**< size of the automorphism group computed */

            /**
             * @return Number of stored generators using a sparse data structure.
             */
            [[nodiscard]] int s_sparsegen() const {
                return generators.s_stored_sparse;
            }

            /**
             * @return Number of stored generators using a dense data structure.
             */
            [[nodiscard]] int s_densegen() const {
                return generators.s_stored_dense;
            }

            /**
             * Set up this Schreier structure using the given base. The base is then fixed and can not be adjusted
             * later on.
             *
             * @param base the base
             * @param base_sizes upper bounds for the size of transversals
             * @param stop integer which indicates to stop reading the base at this position
             */
            void initialize(const int new_domain_size, std::vector<int> &base, std::vector<int> &base_sizes,
                            const int stop) {
                assert(base.size() >= stop);
                domain_size = new_domain_size;
                assert(this->domain_size > 0);
                generators.initialize(domain_size);
                generators.clear();
                transversals.reserve(stop);
                transversals.clear();

                finished_up_to = -1;

                //transversals.set_size(stop);
                for (int i = 0; i < stop; ++i) {
                    transversals.push_back(new shared_transversal());
                    //transversals[i] = new shared_transversal();
                    transversals[i]->initialize(base[i], i, base_sizes[i]);
                }
                init = true;
            }

            /**
             * Reset up this Schreier structure with a new base.
             *
             * @param new_base the new base
             * @param new_base_sizes upper bounds for the size of transversals
             * @param stop integer which indicates to stop reading the base at this position
             * @param keep_old If true, attempt to keep parts of the base that is already stored.
             */
            bool reset(int new_domain_size, schreier_workspace& w, std::vector<int> &new_base,
                       std::vector<int> &new_base_sizes, const int stop, bool keep_old,
                       std::vector<int> &global_fixed_points) {
                const int old_size = static_cast<int>(transversals.size());
                if(!init || new_domain_size != domain_size) {
                    initialize(new_domain_size, new_base, new_base_sizes, stop);
                    return false;
                }
                const int new_size = stop;

                // compare with stored base, keep whatever is possible
                int keep_until = 0;
                if(keep_old) {
                    for (; keep_until < old_size && keep_until < new_size; ++keep_until) {
                        if (transversals[keep_until]->fixed_point() != new_base[keep_until]) break;
                    }
                } else {
                    //generators.clear();
                    generators.filter(w, global_fixed_points);
                }

                if(keep_until == new_size && new_size == old_size) return false;

                finished_up_to = -1;

                transversals.resize(new_size);
                //transversals.set_size(new_size);


                for(int i = 0; i < keep_until; ++i) {
                    transversals[i]->set_size_upper_bound(new_base_sizes[i]);
                }
                for (int i = keep_until; i < stop; ++i) {
                    if(i < old_size) delete transversals[i];
                    transversals[i] = new shared_transversal();
                    assert(new_base[i] >= 0);
                    transversals[i]->initialize(new_base[i], i, new_base_sizes[i]);
                }

                // TODO resift old generators if desired

                return true;
            }

            /**
             * Returns a vertex to individualize for each color of \p root_coloring that matches in size a corresponding
             * transversal.
             *
             * @param save_to_individualize Vector in which vertices deemed save to individualize are pushed.
             * @param root_coloring The coloring with which the stored transversals are compared.
             */
            void determine_potential_individualization(std::vector<std::pair<int, int>>* save_to_individualize,
                                                       coloring* root_coloring) {
                for (int i = base_size()-1; i >= 0; --i) {
                    const int corresponding_root_color_sz = root_coloring->ptn[root_coloring->vertex_to_col[transversals[i]->fixed_point()]] + 1;
                    if(transversals[i]->size() >= corresponding_root_color_sz && corresponding_root_color_sz > 1) {
                        save_to_individualize->emplace_back(transversals[i]->fixed_point(), corresponding_root_color_sz);
                    }
                }
            }

            /**
             * @param pos Position in base.
             * @return Vertex fixed at position \p pos in base.
             */
            [[nodiscard]] int base_point(int pos) const {
                return transversals[pos]->fixed_point();
            }

            /**
             * @return Size of base of this Schreier structure.
             */
            [[nodiscard]] int base_size() const {
                return transversals.size();
            }

            /**
             * Checks whether a vertex \p v is contained in the traversal at position \p s_base_pos.
             *
             * @param base_pos Position in base.
             * @param v Vertex to check.
             * @return Bool indicating whether \p v is contained in the traversal at position \p s_base_pos.
             */
            bool is_in_base_orbit(const int base_pos, const int v) {
                if (base_pos >= static_cast<int>(transversals.size())) return false;
                if(v < 0) return false;
                assert(base_pos >= 0);
                assert(base_pos < transversals.size());
                const int search = transversals[base_pos]->find_point(v);
                return search != -1;
            }

            [[nodiscard]]int get_base_orbit_size(const int base_pos) {
                if (base_pos >= static_cast<int>(transversals.size())) return 0;
                return transversals[base_pos]->size();
            }

            /**
             * Reduces a vector of vertices \p selection to contain only points not contained in transversal at position
             * \p s_base_pos in Schreier structure.
             *
             * @param w A Schreier workspace.
             * @param selection Vector to be reduced.
             * @param base_pos Position in base.
             */
            void reduce_to_unfinished(schreier_workspace &w, std::vector<int> &selection, int base_pos) {
                transversals[base_pos]->reduce_to_unfinished(w, selection);
            }

            /**
             * Checks whether the traversal at position \p s_base_pos matches its size upper bound.
             *
             * @param base_pos Position in base.
             * @return Bool that indicates whether the traversal at position \p s_base_pos matches its size upper bound.
             */
            bool is_finished(const int base_pos) {
                return transversals[base_pos]->is_finished();
            }

            /**
             * Sift automorphism into the Schreier structure.
             *
             * @param w Auxiliary workspace used for procedures.
             * @param automorphism Automorphism to be sifted. Will be manipulated by the method.
             * @return Whether automorphism was added to the Schreier structure or not.
             */
            bool sift(schreier_workspace &w, automorphism_workspace &automorphism, bool uniform = false) {
                bool changed = false; /*< keeps track of whether we changed the Schreier structure while sifting */

                automorphism.set_support01(true); // we don't need to track full support
                for (int level = 0; level < static_cast<int>(transversals.size()); ++level) { // sift level-by-level
                    // first, we try to extend the traversal using the new automorphism
                    changed = transversals[level]->extend_with_automorphism(w, generators, automorphism)
                              || changed;

                    // secondly, we fix the point of this transversal in automorphism
                    const bool identity = transversals[level]->fix_automorphism(w, generators, automorphism);
                    if (identity) {
                        //std::cout << "sift until " << level << ", c: " << changed << std::endl;
                        break; // if automorphism is the identity now, no need to sift further
                    }
                }

                // keep track of how far this Schreier structure is finished (matches given upper bounds)
                for (int level = finished_up_to + 1; level < static_cast<int>(transversals.size()); ++level) {
                    if (finished_up_to == level - 1 && transversals[level]->is_finished()) {
                        ++finished_up_to;
                    } else {
                        break;
                    }
                }

                    // reset the automorphism_workspace
                automorphism.set_support01(false);
                automorphism.update_support();
                automorphism.reset();

                if(uniform) record_sift_result(changed); // uniform automorphisms count towards abort criterion

                return changed;
            }

            /**
             * Generate a (semi-)random element from the generators.
             *
             * @param w Auxiliary workspace used for procedures.
             * @param automorphism Random element is saved in this automorphism_workspace.
             * @param rn_generator Random number generator used for the generation.
             */
            void generate_random(schreier_workspace& w, automorphism_workspace& automorphism, std::default_random_engine& rn_generator) {
                automorphism.reset();

                const int random_mult = static_cast<int>(rn_generator() % 7); // 7
                const int num_mult = 1 + (random_mult);
                for(int i = 0; i < num_mult; ++i) {
                    // load generator
                    const int next_gen_num = static_cast<int>(rn_generator() % generators.size());
                    assert(next_gen_num >= 0);
                    assert(next_gen_num < generators.size());
                    auto next_gen = generators[next_gen_num];
                    assert(next_gen != nullptr);
                    next_gen->load(w.loader, w.scratch_auto);

                    // multiply
                    automorphism.apply(w.scratch_apply1, w.scratch_apply2, w.scratch_apply3, w.loader.p(), 1);
                    w.scratch_auto.reset();
                }
            }

            /**
             * Sift a (semi-)randomly generated element into the Schreier structure.
             *
             * @param w Auxiliary workspace used for procedures.
             * @param automorphism An automorphism_workspace used to store the randomly generated element.
             * @return Whether the generated automorphism was added to the Schreier structure or not.
             */
            bool sift_random(schreier_workspace &w, automorphism_workspace& automorphism,
                             std::default_random_engine& rn_generator) {
                if(generators.size() <= 0) {
                    return false;
                }
                automorphism.reset();
                automorphism.set_support01(true);
                generate_random(w, automorphism, rn_generator);
                return sift(w, automorphism, false);
            }

            /**
             * Records a sift result for the probabilistic abort criterion.
             *
             * @param changed Whether the sift was successful or not.
             */
            void record_sift_result(const bool changed) {
                if(!changed) {
                    ++s_consecutive_success;
                } else {
                    if(s_consecutive_success > 0) {
                        ++h_error_bound;
                        s_consecutive_success = 0;
                    }
                }
            }

            /**
             * Reset the probabilistic abort criterion.
             */
            void reset_probabilistic_criterion() {
                s_consecutive_success = 0;
            }

            /**
             * @return Whether the probabilistic abort criterion allows termination or not.
             */
            [[nodiscard]] bool probabilistic_abort_criterion() const {
                return (s_consecutive_success > h_error_bound);
            }

            /**
             * @return Whether the deterministic abort criterion allows termination or not.
             */
            [[nodiscard]] bool deterministic_abort_criterion() const {
                /*int num_finished = 0;
                int num_not_finished = 0;
                for (int level = 0; level < static_cast<int>(transversals.size()); ++level) {
                    num_finished += transversals[level]->is_finished();
                    num_not_finished += !transversals[level]->is_finished();
                    if(!transversals[level]->is_finished()) std::cout << level << " bp: " << transversals[level]->fixed_point()  << "sz: " << transversals[level]->size() << "/" << transversals[level]->get_size_upper_bound()  << std::endl;
                }
                    std::cout << finished_up_to_level() + 1 << "/" << num_finished << "-" << num_not_finished << "/" << base_size() << std::endl;*/
                return (finished_up_to_level() + 1 == base_size());
            }

            /**
             * @return Level up to which Schreier structure is guaranteed to be complete according to given upper bounds.
             *         -1 indicates no level has been finished.
             */
            [[nodiscard]] int finished_up_to_level() const {
                return finished_up_to;
            }

            /**
             * @return Size of group described by this Schreier structure.
             */
            void compute_group_size() {
                s_grp_sz.mantissa = 1.0;
                s_grp_sz.exponent = 0;
                // multiply the sizes of the individual levels in the Schreier table
                for (int level = 0; level < static_cast<int>(transversals.size()); ++level) {
                    s_grp_sz.multiply(transversals[level]->size());
                }
            }
        };



        /**
         * \brief Schreier structure with fixed base. Implements precisely the interface as used by the dejavu
         * automorphism solver.
         *
         * Internally, stores a shared_schreier structure in a compressed form using a domain_compressor.
         *
         */
        class compressed_schreier {
        private:
            shared_schreier         internal_schreier;
            domain_compressor*      compressor;
            automorphism_workspace* compressed_automorphism;
            std::vector<int> original_base;
            std::vector<int> original_base_sizes;
            int              original_stop = 0;

            int              s_compressed_until = 0;

        public:
            double h_min_compression_ratio = 0.4; /*< only use compression if at least this ratio can be achieved */
            double s_compression_ratio = 1.0;

            /**
             * @return Number of stored generators using a sparse data structure.
             */
            [[nodiscard]] int s_sparsegen() const {
                return internal_schreier.s_sparsegen();
            }

            /**
             * @return Number of stored generators using a dense data structure.
             */
            [[nodiscard]] int s_densegen() const {
                return internal_schreier.s_densegen();
            }

            /**
             * Reset up this Schreier structure with a new base.
             *
             * @param new_base the new base
             * @param new_base_sizes upper bounds for the size of transversals
             * @param stop integer which indicates to stop reading the base at this position
             * @param keep_old If true, attempt to keep parts of the base that is already stored.
             */
            bool reset(ds::domain_compressor* new_compressor, int new_domain_size, schreier_workspace& w, std::vector<int> &new_base,
                       std::vector<int> &new_base_sizes, const int stop, bool keep_old,
                       std::vector<int> &global_fixed_points) {
                s_compressed_until = 0;
                compressor = new_compressor;
                // do not use compression at all if ration too bad (to save on translation cost, and such that keep_old
                // works for difficult graphs)
                if(compressor != nullptr &&
                   compressor->s_compression_ratio > h_min_compression_ratio) compressor = nullptr;

                if(compressor != nullptr) {
                    // need more management if we are using compression
                    original_base       = new_base;
                    original_base_sizes = new_base_sizes;
                    original_stop = stop;
                    delete compressed_automorphism;
                    compressed_automorphism = new automorphism_workspace(compressor->compressed_domain_size());
                    keep_old = false;
                    std::vector<int> new_basec;
                    new_basec.reserve(new_base.size());
                    for (int i = 0; i < new_base.size(); ++i) {
                        new_basec.push_back(compressor->compress(new_base[i]));
                    }
                    s_compression_ratio = compressor->s_compression_ratio;
                    return internal_schreier.reset(new_compressor->compressed_domain_size(), w, new_basec,
                                                    new_base_sizes, stop, keep_old,
                                                    global_fixed_points);
                } else {
                    // we are not using compression, so simply use the Schreier structure normally
                    s_compression_ratio = 1.0;
                    return internal_schreier.reset(new_domain_size, w, new_base,
                                                    new_base_sizes, stop, keep_old,
                                                    global_fixed_points);
                }
            }

            void consider_compress_more_from_start(schreier_workspace& w, int compress_until, coloring& c) {
                // TODO breaks some large-cfi
                if(compressor == nullptr || s_compression_ratio < 0.1) return;
                assert(compressor != nullptr);
                std::vector<int> new_base;
                std::vector<int> new_base_sizes;
                new_base.resize(original_base.size() - compress_until);
                new_base_sizes.resize(original_base.size() - compress_until);
                int j = 0;
                for(int i = compress_until; i < original_base.size(); ++i) {
                    new_base[j]       = original_base[i];
                    new_base_sizes[j] = original_base_sizes[i];
                    ++j;
                }
                const int new_stop = original_stop - compress_until;

                if(compressor->determine_compression_ratio(c, new_base, new_stop) >= 0.1) return;
                s_compressed_until = compress_until;
                compressor->determine_compression_map(c, new_base, new_stop);
                s_compression_ratio = compressor->s_compression_ratio;
                std::vector<int> new_basec;
                for (int &b: new_base) {
                    new_basec.push_back(compressor->compress(b));
                }
                delete compressed_automorphism;
                compressed_automorphism = new automorphism_workspace(compressor->compressed_domain_size());
                std::vector<int> empty_vec;
                internal_schreier.reset(compressor->compressed_domain_size(), w, new_basec,
                                               new_base_sizes, new_stop, false,empty_vec);
            }

            /**
             * Returns a vertex to individualize for each color of \p root_coloring that matches in size a corresponding
             * transversal.
             *
             * @param save_to_individualize Vector in which vertices deemed save to individualize are pushed.
             * @param root_coloring The coloring with which the stored transversals are compared.
             */
            void determine_potential_individualization(std::vector<std::pair<int, int>>* save_to_individualize,
                                                       coloring* root_coloring) {
                if(s_compressed_until > 0) return; // let's not bother with this case...

                if(compressor == nullptr) {
                    internal_schreier.determine_potential_individualization(save_to_individualize, root_coloring);
                    return;
                }

                for (int i = base_size()-1; i >= 0; --i) {
                    const int corresponding_root_color_sz = root_coloring->ptn[root_coloring->vertex_to_col[original_base[i]]] + 1;
                    if(internal_schreier.get_base_orbit_size(i) >= corresponding_root_color_sz && corresponding_root_color_sz > 1) {
                        save_to_individualize->emplace_back(original_base[i], corresponding_root_color_sz);
                    }
                }
            }

            /**
             * @param pos Position in base.
             * @return Vertex fixed at position \p pos in base.
             */
            [[nodiscard]] int base_point(int pos) const {
                if(compressor != nullptr) {
                    return original_base[pos];
                } else {
                    return internal_schreier.base_point(pos);
                }
            }

            /**
             * @return Size of base of this Schreier structure.
             */
            [[nodiscard]] int base_size() const {
                return internal_schreier.base_size();
            }

            /**
             * Checks whether a vertex \p v is contained in the traversal at position \p s_base_pos.
             *
             * @param base_pos Position in base.
             * @param v Vertex to check.
             * @return Bool indicating whether \p v is contained in the traversal at position \p s_base_pos.
             */
            bool is_in_base_orbit(const int base_pos, const int v) {
                int v_translate = compressor != nullptr? compressor->compress(v) : v;
                return internal_schreier.is_in_base_orbit(base_pos-s_compressed_until, v_translate);
            }

            /**
             * Reduces a vector of vertices \p selection to contain only points not contained in transversal at position
             * \p s_base_pos in Schreier structure.
             *
             * @param w A Schreier workspace.
             * @param selection Vector to be reduced.
             * @param base_pos Position in base.
             */
            void reduce_to_unfinished(schreier_workspace &w, std::vector<int> &selection, int base_pos) {
                if(compressor != nullptr) {
                    for(int i = 0; i < selection.size(); ++i) {
                        assert(compressor->compress(selection[i]) >= 0);
                        selection[i] = compressor->compress(selection[i]);
                    }
                }
                internal_schreier.reduce_to_unfinished(w, selection, base_pos-s_compressed_until);
                if(compressor != nullptr) {
                    for(int i = 0; i < selection.size(); ++i) {
                        selection[i] = compressor->decompress(selection[i]);
                        assert(selection[i] >= 0);
                    }
                }
            }

            void compress_automorphism(automorphism_workspace &automorphism, automorphism_workspace &automorphism_compress) {
                automorphism_compress.reset();

                const int support = automorphism.nsupport();
                for(int i = 0; i < support; ++i) {
                    const int vsupport = automorphism.support()[i];
                    const int vsupportc = compressor->compress(vsupport);
                    const int vmapsto = automorphism.perm()[vsupport];
                    const int vmapstoc = compressor->compress(vmapsto);
                    if(vsupportc >= 0) {
                        assert(vmapstoc >= 0);
                        automorphism_compress.write_single_map(vsupportc, vmapstoc);
                    }
                }
            }

            /**
             * Sift automorphism into the Schreier structure.
             *
             * @param w Auxiliary workspace used for procedures.
             * @param automorphism Automorphism to be sifted. Will be manipulated by the method.
             * @return Whether automorphism was added to the Schreier structure or not.
             */
            bool sift(schreier_workspace &w, automorphism_workspace &automorphism, bool uniform = false) {
                if(compressor != nullptr) {
                    compress_automorphism(automorphism, *compressed_automorphism);
                    return internal_schreier.sift(w, *compressed_automorphism, uniform);
                } else {
                    return internal_schreier.sift(w, automorphism, uniform);
                }
            }

            /**
             * Sift a (semi-)randomly generated element into the Schreier structure.
             *
             * @param w Auxiliary workspace used for procedures.
             * @param automorphism An automorphism_workspace used to store the randomly generated element.
             * @return Whether the generated automorphism was added to the Schreier structure or not.
             */
            bool sift_random(schreier_workspace &w, automorphism_workspace& automorphism,
                             std::default_random_engine& rn_generator) {
                if(compressor != nullptr) {
                    return internal_schreier.sift_random(w, *compressed_automorphism, rn_generator);
                } else {
                    return internal_schreier.sift_random(w, automorphism, rn_generator);
                }
            }

            /**
             * Records a sift result for the probabilistic abort criterion.
             *
             * @param changed Whether the sift was successful or not.
             */
            void record_sift_result(const bool changed) {
                internal_schreier.record_sift_result(changed);
            }

            /**
             * Reset the probabilistic abort criterion.
             */
            void reset_probabilistic_criterion() {
                internal_schreier.reset_probabilistic_criterion();
            }

            /**
             * @return Whether the probabilistic abort criterion allows termination or not.
             */
            [[nodiscard]] bool probabilistic_abort_criterion() const {
                return internal_schreier.probabilistic_abort_criterion();
            }

            /**
             * @return Whether the deterministic abort criterion allows termination or not.
             */
            [[nodiscard]] bool deterministic_abort_criterion() const {
                //std::cout << internal_schreier.base_point(14) << "->" << compressor->decompress(internal_schreier.base_point(14)) << std::endl;
                return internal_schreier.deterministic_abort_criterion();
            }

            void set_error_bound(int error_bound) {
                internal_schreier.h_error_bound = error_bound;
            }

            int get_consecutive_success() {
                return internal_schreier.s_consecutive_success;
            }

            big_number get_group_size() {
                return internal_schreier.s_grp_sz;
            }

            /**
             * @return Level up to which Schreier structure is guaranteed to be complete according to given upper bounds.
             *         -1 indicates no level has been finished.
             */
            [[nodiscard]] int finished_up_to_level() const {
                return s_compressed_until + internal_schreier.finished_up_to_level();
            }

            /**
             * @return Size of group described by this Schreier structure.
             */
            void compute_group_size() {
                internal_schreier.compute_group_size();
            }
        };
    }
}

#endif //DEJAVU_GROUPS_H
