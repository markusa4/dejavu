#ifndef DEJAVU_GROUPS_H
#define DEJAVU_GROUPS_H

#include "refinement.h"
#include "bijection.h"
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

        // reset internal automorphism structure to the identity
        static void reset_automorphism(int *rautomorphism, work_list *automorphism_supp) {
            for (int i = 0; i < automorphism_supp->cur_pos; ++i) {
                rautomorphism[(*automorphism_supp)[i]] = (*automorphism_supp)[i];
            }
            automorphism_supp->reset();
        };

        // make automorphism from colorings
        static void color_diff_automorphism(int n, const int *vertex_to_col, const int *col_to_vertex,
                                            int *automorphism_map, work_list *automorphism_supp) {
            for (int v1 = 0; v1 < n; ++v1) {
                const int col = vertex_to_col[v1];
                const int v2 = col_to_vertex[col];
                if (v1 != v2) {
                    automorphism_map[v1] = v2;
                    automorphism_supp->push_back(v1);
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
            int domain_size;

            bool support01 = false;
        public:
            /**
             * Initializes the stored automorphism to the identity.
             *
             * @param domain_size Size of the domain on which automorphisms operate
             */
            automorphism_workspace(int domain_size) {
                automorphism.initialize(domain_size);
                for (int i = 0; i < domain_size; ++i)
                    automorphism[i] = i;
                automorphism_supp.initialize(domain_size);
                this->domain_size = domain_size;
            }

            void set_support01(bool support01) {
                this->support01 = support01;
            }

            /**
             * Create automorphism from two given vertex colorings.
             *
             * @param vertex_to_col a vertex-to-color mapping which describes the first coloring
             * @param col_to_vertex a color-to-vertex mapping which describes the second coloring
             */
            void write_color_diff(const int *vertex_to_col, const int *col_to_vertex) {
                color_diff_automorphism(domain_size, vertex_to_col, col_to_vertex, automorphism.get_array(),
                                        &automorphism_supp);
            }

            /**
             * Apply another automorphism to the stored automorphism. Closely follows the implementation in nauty / Traces.
             *
             * @param other
             * @param pwr
             */
            void apply(work_list &scratch_apply1, work_list &scratch_apply2, mark_set &scratch_apply3,
                       automorphism_workspace *other, int pwr = 1) {
                apply(scratch_apply1, scratch_apply2, scratch_apply3, other->perm(), pwr);
            }

            void __attribute__ ((noinline)) update_support() {
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
             * Apply another automorphism to the stored automorphism. Closely follows the implementation in nauty / Traces.
             *
             * @param other
             * @param pwr
             */
            void apply(work_list &scratch_apply1, work_list &scratch_apply2, mark_set &scratch_apply3,
                  const int *p, int pwr = 1) {
                if (pwr == 0)
                    return;
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
                    if (pwr >= 6) {
                        for (int j = 0; j < domain_size; ++j) {
                            scratch_apply1[j] = p[p[p[j]]];
                        }
                        for (; pwr >= 6; pwr -= 6)
                            for (int j = 0; j < domain_size; ++j)
                                automorphism[j] = scratch_apply1[scratch_apply1[automorphism[j]]];
                    }

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
                }

                // rewrite support
                update_support();
            }


            /**
             * Create automorphism from two canonically-ordered vectors of singletons. The resulting automorphism maps
             * \p singletons1[i] to \p singletons2[i] for i in \p pos_start, ..., \p pos_end.
             *
             * @param singletons1 first vector of singletons
             * @param singletons2 second vector of singletons
             * @param pos_start start reading the vectors at this position
             * @param pos_end stop reading the vecvtors at this position.
             */
            void write_singleton(std::vector<int> *singletons1, std::vector<int> *singletons2, int pos_start, int pos_end) {
                for (int i = pos_start; i < pos_end; ++i) {
                    const int from = (*singletons1)[i];
                    const int to = (*singletons2)[i];
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
            int *perm() {
                return automorphism.get_array();
            }

            /**
             * @return Integer array which contains all vertices in the support of the contained automorphism.
             */
            int *support() {
                return automorphism_supp.get_array();
            }

            /**
             * @return Size of the support.
             */
            int nsupport() {
                return automorphism_supp.cur_pos;
            }
        };

        /**
         * \brief Orbit partition
         *
         * Keeps track of an orbit partition, and provides tools to manipulate the orbits within.
         */
        class orbit {
            int               sz;
            mark_set          touched;
            work_list_t<int>  reset_arr;
            work_list_t<int>  map_arr;
            work_list_t<int>  orb_sz;
        public:

            /**
             * Retrieve the orbit of the given vertex.
             *
             * @param v Vertex of the specified domain.
             * @return The orbit of \p v.
             */
            int find_and_cut_orbit(const int v) {
                assert(v >= 0);
                assert(v < sz);
                int orbit1 = map_arr[v];
                while(orbit1 != map_arr[orbit1])
                    orbit1 = map_arr[orbit1];
                map_arr[v] = orbit1;
                return orbit1;
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
                return orb_sz[find_and_cut_orbit(v)];
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
                    if(!touched.get(v1))
                        reset_arr.push_back(v1);
                    if(!touched.get(v2))
                        reset_arr.push_back(v2);
                    touched.set(v1);
                    touched.set(v2);
                    int orbit1 = find_and_cut_orbit(v1);
                    int orbit2 = find_and_cut_orbit(v2);
                    if(orbit1 == orbit2)
                        return;
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
                if(v1 == v2)
                    return true;
                const int orbit1 = find_and_cut_orbit(v1);
                const int orbit2 = find_and_cut_orbit(v2);
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
            void add_automorphism_to_orbit(groups::automorphism_workspace& aut) {
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
             * @param domain_size The size of the underlying domain.
             */
            void initialize(int domain_size) {
                sz = domain_size;
                touched.initialize(domain_size);
                reset_arr.initialize(domain_size);
                map_arr.initialize(domain_size);
                orb_sz.initialize(domain_size);
                for(int i = 0; i < domain_size; ++i) {
                    map_arr.push_back(i);
                    orb_sz.push_back(1);
                }
            }
        };

        typedef void type_unload_hook();

        /**
         * \brief Stores a link to an automorphism.
         *
         * Used to implement indiscriminate loading of dense and sparse automorphisms.
         */
        class dense_sparse_arbiter {
            int *loaded_automorphism;
            std::function<type_unload_hook> *unload_hook;
        public:
            void load(int *automorphism, std::function<type_unload_hook> *unload_hook) {
                loaded_automorphism = automorphism;
                this->unload_hook = unload_hook;
            }

            void unload() {
                if (unload_hook) (*unload_hook)();
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
            int domain_size;
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
             * @param automorphism Auxiliary space used if loaded automorphism is sparse.
             */
            void load(dense_sparse_arbiter &loader, automorphism_workspace &automorphism) {
                if (store_type == STORE_SPARSE) {
                    automorphism.reset();
                    int first_of_cycle = 0;
                    for (int i = 0; i < data.size(); ++i) {
                        const int j = abs(data[i]) - 1;
                        const bool is_last = data[i] < 0;
                        assert(i == data.size() - 1 ? is_last : true);
                        if (is_last) {
                            automorphism.write_single_map(j, abs(data[first_of_cycle]) - 1);
                            first_of_cycle = i + 1;
                        } else {
                            assert(i + 1 < data.size());
                            automorphism.write_single_map(j, abs(data[i + 1]) - 1);
                        }
                    }
                    loader.load(automorphism.perm(), nullptr);
                } else {
                    loader.load(data.get_array(), nullptr);
                }
            }

            /**
             * Load inverse of stored automorphism into workspace.
             *
             * @param loader Store a pointer to permutation of automorphism in the arbiter.
             * @param automorphism Auxiliary space used if loaded automorphism is sparse.
             */
            void load_inverse(dense_sparse_arbiter &loader, automorphism_workspace &automorphism) {
                if (store_type == STORE_SPARSE) {
                    automorphism.reset();
                    int first_of_cycle = data.size() - 1;
                    for (int i = data.size() - 1; i >= 0; --i) {
                        if (data[i] < 0)
                            first_of_cycle = i;

                        const int j = abs(data[i]) - 1;
                        const bool is_last = i == 0 || (data[i - 1] < 0);
                        if (is_last) {
                            automorphism.write_single_map(j, abs(data[first_of_cycle]) - 1);
                        } else {
                            automorphism.write_single_map(j, abs(data[i - 1]) - 1);
                        }
                    }
                    loader.load(automorphism.perm(), nullptr);
                } else {
                    automorphism.reset();
                    for (int i = 0; i < domain_size; ++i) {
                        if (i != data[i]) {
                            automorphism.write_single_map(data[i], i);
                        }
                    }
                    loader.load(automorphism.perm(), nullptr);
                }
            }

            /**
             * Store the given automorphism workspace.
             *
             * @param automorphism The automorphism to be stored.
             */
            void store(int domain_size, automorphism_workspace &automorphism, mark_set &helper) {
                this->domain_size = domain_size;
                assert(data.empty());

                int support = 0;
                for (int i = 0; i < domain_size; ++i) support += (automorphism.perm()[i] != i);

                // decide whether to store dense or sparse representation
                if (support < domain_size / 4) {
                    store_type = STORE_SPARSE;
                    helper.reset();

                    data.initialize(support);
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
                    assert(data.size() == support);
                } else {
                    store_type = STORE_DENSE;
                    data.initialize(domain_size);
                    data.set_size(domain_size);
                    memcpy(data.get_array(), automorphism.perm(), domain_size * sizeof(int));
                    assert(data.size() == domain_size);
                }
            }
        };

        /**
         * \brief Workspace used for Schreier computations
         *
         * A global (thread local) state used for computations in Schreier structures.
         */
        class schreier_workspace {
        public:
            /**
             * Initialize this workspace.
             *
             * @param domain_size
             * @param R only used for debug purposes
             * @param g only used for debug purposes
             */
            schreier_workspace(int domain_size, refinement *R, sgraph *g) : scratch_auto(domain_size) {
                scratch1.initialize(domain_size);
                scratch2.initialize(domain_size);
                scratch_apply1.initialize(domain_size);
                scratch_apply2.initialize(domain_size);
                scratch_apply3.initialize(domain_size);
                this->R = R;
                this->g = g;
            }

            dense_sparse_arbiter loader;

            mark_set scratch1;
            mark_set scratch2;
            work_list scratch_apply1;
            work_list scratch_apply2;
            mark_set scratch_apply3;
            automorphism_workspace scratch_auto;

            // TODO debug
            refinement *R;
            sgraph *g;
        };

        /**
         * \brief Stores a generating set.
         *
         * Can be used across multiple threads.
         */
        class generating_set {
            std::mutex lock_generators;          /**< locks the generators */
            std::vector<stored_automorphism *> generators; /** list of generators */
            int domain_size;
        public:
            int s_stored_sparse = 0; /**< how many generators are stored in a sparse manner */
            int s_stored_dense = 0;  /**< how many generators are stored in a dense manner */

            /**
             * Set up this generating set.
             * @param domain_size Size of the domain of the stored generators.
             */
            void initialize(int domain_size) {
                this->domain_size = domain_size;
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
                const int num = generators.size() - 1;
                assert(w.R->certify_automorphism_sparse(w.g, automorphism.perm(), automorphism.nsupport(),
                                                        automorphism.support()));
                generators[num]->store(domain_size, automorphism, w.scratch2);
                lock_generators.unlock();

                s_stored_sparse += (generators[num]->get_store_type() ==
                                    stored_automorphism::stored_automorphism_type::STORE_SPARSE);
                s_stored_dense += (generators[num]->get_store_type() ==
                                   stored_automorphism::stored_automorphism_type::STORE_DENSE);

                return num;
            }

            /**
             * Retrieve a generator from the stored generating set.
             *
             * @param num Identifier of a generator.
             * @return A pointer to the generator which corresponds to the identifier.
             */
            stored_automorphism *get_generator(const int num) {
                // TODO locks...
                return generators[num];
            }

            /**
             * @return number of stored generators
             */
            int size() {
                return generators.size();
            }
        };

        /**
         * \brief A transversal in a Schreier structure.
         *
         * Can be used across multiple threads.
         */
        class transversal {
            enum stored_transversal_type {
                STORE_DENSE, STORE_SPARSE
            };

            std::mutex lock_transversal;          /**< locks this transversal */

            int fixed;                            /**< vertex fixed by this transversal */
            int sz_upb = INT32_MAX;               /**< upper bound for size of the transversal (e.g. color class size) */
            int level;
            bool finished = false;

            std::vector<int> fixed_orbit;         /**< contains vertices of orbit at this schreier level */
            std::vector<int> fixed_orbit_to_perm; /**< maps fixed_orbit[i] to generators[i] in class \ref schreier. */
            std::vector<int> fixed_orbit_to_pwr;  /**< power that should to be applied to generators[i] in class \ref schreier. */

            stored_transversal_type store_type = STORE_SPARSE; /**< whether above structures are stored dense or sparse */

            void __attribute__ ((noinline)) load_orbit_to_scratch(schreier_workspace &w) {
                for (int i = 0; i < fixed_orbit.size(); ++i) {
                    w.scratch1.set(fixed_orbit[i]);
                }
            }

            void add_to_fixed_orbit(const int vertex, const int perm, const int pwr) {
                //std::cout << "add " << vertex << ", " << perm << ", " << pwr << std::endl;
                fixed_orbit.push_back(vertex);
                fixed_orbit_to_perm.push_back(perm);
                fixed_orbit_to_pwr.push_back(pwr);
            }


            void apply_perm(schreier_workspace &w, automorphism_workspace &automorphism,
                       generating_set &generators, const int perm, const int pwr) {
                // load perm into workspace
                auto generator = generators.get_generator(perm);

                // apply generator
                if (pwr < 0) { // use inverse
                    generator->load_inverse(w.loader, w.scratch_auto);
                    // multiply
                    automorphism.apply(w.scratch_apply1, w.scratch_apply2, w.scratch_apply3, w.loader.p(), abs(pwr));
                } else if (pwr > 0) { // use generator
                    generator->load(w.loader, w.scratch_auto);
                    // multiply
                    automorphism.apply(w.scratch_apply1, w.scratch_apply2, w.scratch_apply3, w.loader.p(), pwr);
                }
            }

        public:

            /**
             * @return Size of the transversal.
             */
            int size() {
                return fixed_orbit.size();
            }

            /**
             * Check whether a point \p p is contained in transversal.
             *
             * @param p Point to check.
             * @return Position of point \p in transversal, or -1 if not contained.
             */
            int find_point(const int p) {
                for (int i = 0; i < fixed_orbit.size(); ++i) {
                    if (p == fixed_orbit[i]) {
                        return i;
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
                int back_swap = selection.size() - 1;
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
            bool is_finished() {
                return finished;
            }

            /**
             * @return Point fixed by this transversal.
             */
            int fixed_point() {
                return fixed;
            }

            /**
             * Initialize this transversal.
             *
             * @param fixed_vertex Vertex fixed by the transversal.
             * @param level Position of the transversal in the base of Schreier structure.
             * @param sz_upb Upper bound for the size of transversal (e.g., color class size in combinatorial base).
             */
            void initialize(const int fixed_vertex, const int level, const int sz_upb) {
                fixed = fixed_vertex;
                this->level = level;
                this->sz_upb = sz_upb;
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
                if (finished)
                    return false;

                load_orbit_to_scratch(w);
                bool changed = false;
                // TODO probably acquire lock already here
                // TODO how this is performed precisely should depend on fixed_orbit size versus support of generator
                int gen_num = -1;

                int self_apply_after = 0;

                for (int i = 0; i < fixed_orbit.size(); ++i) {
                    int j = automorphism.perm()[fixed_orbit[i]];
                    if (w.scratch1.get(j))
                        continue;

                    int ipwr = 0; // inverted power
                    int npwr = 1; //non-inverted power
                    int jj;
                    for (jj = j; !w.scratch1.get(jj); jj = automorphism.perm()[jj]) ++ipwr;
                    //bool reversible = (jj == fixed_orbit[i]) && (ipwr > 100);

                    /*if(!w.scratch1.get(j) && reversible && i == 0) {
                        assert(fixed_orbit[i] == fixed);
                        self_apply_after = ipwr;
                    }*/

                    //int pwr = 1;
                    while (!w.scratch1.get(j)) {
                        // we change this traversal
                        changed = true;

                        // add generator to generating set (once)
                        if (gen_num == -1) gen_num = generators.add_generator(w, automorphism);

                        // add entry to traversal
                        //if(ipwr <= npwr || !reversible) {
                        add_to_fixed_orbit(j, gen_num, ipwr);
                        // } else {
                        //    add_to_fixed_orbit(j, gen_num, -npwr);
                        // }
                        w.scratch1.set(j);

                        // we check out the entire cycle of j now
                        j = automorphism.perm()[j];
                        --ipwr;
                        ++npwr;
                    }
                }

                /*if(self_apply_after > 0) {
                    assert(changed);
                    assert(gen_num >= 0);
                    apply_perm(w, automorphism, generators, gen_num, self_apply_after);
                    assert(automorphism.perm()[fixed] == fixed);
                }*/

                if (sz_upb == fixed_orbit.size() && !finished) {
                    finished = true;
                }

                w.scratch1.reset();
                return changed;
            }

            /**
             * @param generators
             * @param automorphism
             * @return whether transversal was extended
             */
            bool fix_automorphism(schreier_workspace &w, generating_set &generators,
                                  automorphism_workspace &automorphism) {
                int fixed_map = automorphism.perm()[fixed];
                while (fixed != fixed_map) {
                    const int pos = find_point(fixed_map);
                    const int perm = fixed_orbit_to_perm[pos];
                    const int pwr = fixed_orbit_to_pwr[pos];
                    apply_perm(w, automorphism, generators, perm, pwr);
                    fixed_map = automorphism.perm()[fixed];
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
        class schreier {
            // TODO could use color sizes for memory allocation, predicting dense/sparse?
            // TODO using inverses to reduce pwr seems like low-hanging fruit

            int domain_size;
            int finished_up_to = -1;

            generating_set generators;
            work_list_t<transversal *> transversals;

        public:

            /**
             * @return Number of stored generators using a sparse data structure.
             */
            int stat_sparsegen() {
                return generators.s_stored_sparse;
            }

            /**
             * @return Number of stored generators using a dense data structure.
             */
            int stat_densegen() {
                return generators.s_stored_dense;
            }

            /**
             * Set up this Schreier structure using the given base. The base is then fixed and can not be adjusted
             * later on.
             *
             * @param base the base
             * @param stop integer which indicates to stop reading the base at this position
             */
            void initialize(const int domain_size, std::vector<int> &base, std::vector<int> &base_sizes, const int stop) {
                assert(base.size() >= stop);
                this->domain_size = domain_size;
                generators.initialize(domain_size);
                transversals.initialize(stop);
                transversals.set_size(stop);
                for (int i = 0; i < stop; ++i) {
                    transversals[i] = new transversal();
                    transversals[i]->initialize(base[i], i, base_sizes[i]);
                }
            }

            /**
             * @param pos Position in base.
             * @return Vertex fixed at position \p pos in base.
             */
            int base_point(int pos) {
                return transversals[pos]->fixed_point();
            }

            /**
             * @return Size of base of this Schreier structure.
             */
            int base_size() {
                return transversals.size();
            }

            /**
             * Checks whether a vertex \p v is contained in the traversal at position \p base_pos.
             *
             * @param base_pos Position in base.
             * @param v Vertex to check.
             * @return Bool indicating whether \p v is contained in the traversal at position \p base_pos.
             */
            bool is_in_base_orbit(const int base_pos, const int v) {
                if (base_pos >= transversals.size()) return false;
                assert(base_pos >= 0);
                assert(base_pos < transversals.size());
                const int search = transversals[base_pos]->find_point(v);
                return search != -1;
            }

            /**
             * Reduces a vector of vertices \p selection to contain only points not contained in transversal at position
             * \p base_pos in Schreier structure.
             *
             * @param w A Schreier workspace.
             * @param selection Vector to be reduced.
             * @param base_pos Position in base.
             */
            void reduce_to_unfinished(schreier_workspace &w, std::vector<int> &selection, int base_pos) {
                transversals[base_pos]->reduce_to_unfinished(w, selection);
            }

            /**
             * Checks whether the traversal at position \p base_pos matches its size upper bound.
             *
             * @param base_pos Position in base.
             * @return Bool that indicates whether the traversal at position \p base_pos matches its size upper bound.
             */
            bool is_finished(const int base_pos) {
                return transversals[base_pos]->is_finished();
            }

            /**
             * Sift automorphism into the Schreier structure.
             *
             * @param automorphism Automorphism to be sifted. Will be manipulated by the method.
             * @return Whether automorphism was added to the Schreier structure or not.
             */
            bool sift(schreier_workspace &w, sgraph *g, refinement *R, automorphism_workspace &automorphism) {
                bool changed = false;

                automorphism.set_support01(true);
                for (int level = 0; level < transversals.size(); ++level) {
                    //std::cout << "extending transversal..." << std::endl;
                    changed = transversals[level]->extend_with_automorphism(w, generators, automorphism) || changed;

                    if (finished_up_to == level - 1 && transversals[level]->is_finished()) {
                        ++finished_up_to;
                    }

                    //std::cout << "fixing automorphism..." << std::endl;
                    const bool is_identity = transversals[level]->fix_automorphism(w, generators, automorphism);
                    if (is_identity) break;
                }
                automorphism.set_support01(false);
                automorphism.update_support();
                automorphism.reset();

                //std::cout << generators.size() << ", " << finished_up_to << "/" << transversals.size() << std::endl;
                return changed;
            }

            /**
             * @return Level up to which Schreier structure is guaranteed to be complete according to given upper bounds.
             *         -1 indicates no level has been finished.
             */
            int finished_up_to_level() {
                return finished_up_to;
            }

            /**
             * @return Size of group described by this Schreier structure.
             */
            std::pair<long double, int> compute_group_size() {
                long double grp_sz_man = 1;
                int         grp_sz_exp = 0;
                for (int level = 0; level < transversals.size(); ++level) {
                    grp_sz_man *= transversals[level]->size();
                    while (grp_sz_man > 10) {
                        grp_sz_exp += 1;
                        grp_sz_man = grp_sz_man / 10;
                    }
                }
                return {grp_sz_man, grp_sz_exp};
            }
        };
    }
}

#endif //DEJAVU_GROUPS_H