// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_COMPONENTS_H
#define DEJAVU_COMPONENTS_H

#include "sgraph.h"
#include "coloring.h"
#include "ds.h"

namespace dejavu::ir {

    /**
     * Compute the quotient components of graph \p g colored with vertex coloring \p c.
     *
     * @param g the graph
     * @param c the vertex coloring
     * @param vertex_to_component map from vertex to component, where components are enumerated from 0 to number of
     * components - 1
     * @returns number of components
     */
    static int quotient_components(sgraph *g, int* colmap, ds::worklist *vertex_to_component) {
        coloring c;
        g->initialize_coloring(&c, colmap);

        ds::mark_set  handled(g->v_size);
        ds::mark_set  col_handled(g->v_size);
        ds::worklist  wl(g->v_size);
        int current_component      = 0;
        int total_size             = 0;
        int non_trivial_components = 0;

        for(int i = 0; i < g->v_size; ++i) {
            if(handled.get(i)) continue;
            handled.set(i);

            wl.push_back(i);
            int current_component_size = 0;

            while(!wl.empty()) {
                const int k = wl.pop_back();

                const int col    = c.vertex_to_col[k];
                const int col_sz = c.ptn[col] + 1;

                if(col_sz == 1) {
                    (*vertex_to_component)[k] = -1;
                    continue;
                }

                (*vertex_to_component)[k] = current_component;
                ++current_component_size;

                const int deg = g->d[k];
                const int vpt = g->v[k];
                for (int j = vpt; j < vpt + deg; ++j) {
                    const int neighbour = g->e[j];
                    if(!handled.get(neighbour)) {
                        handled.set(neighbour);
                        wl.push_back(neighbour);
                    }
                }
                if(!col_handled.get(col)) {
                    col_handled.set(col);
                    for (int j = col; j < col + col_sz; ++j) {
                        const int neighbour = c.lab[j];
                        if (!handled.get(neighbour)) {
                            handled.set(neighbour);
                            wl.push_back(neighbour);
                        }
                    }
                }
            }

            total_size += current_component_size;

            if(current_component_size > 0) {
                ++current_component;
            }
            if(current_component_size > 1) ++non_trivial_components;

            if(total_size == g->v_size) break;
        }

        return current_component;
    }

    /**
     * \brief Decomposes graphs and manages decomposition information
     *
     */
    class graph_decomposer {
        int num_components = 0;
        std::vector<sgraph>  components_graph;
        std::vector<int*>    components_coloring;
        std::vector<int>     backward_translation;
        std::vector<int>     component_to_backward_translation;

    public:
        /**
         * Maps back vertex \p vertex of component \p component
         * @param component
         * @param vertex
         * @return vertex of the original graph
         */
        int map_back(int component, int vertex) {
            assert(component >= 0);
            assert(component < num_components);
            assert(component_to_backward_translation[component] + vertex < backward_translation.size());
            return (num_components <= 1? vertex :
                    backward_translation[component_to_backward_translation[component] + vertex]);
        }

        /**
         * Decompose the given graph into components, as defined by \p vertex_to_component. Rearranges \p g and stores
         * decomposition information internally.
         *
         * @param g the graph
         * @param c vertex coloring of \p g
         * @param vertex_to_component maps vertices of \p g to their components
         * @param new_num_components how many components
         */
        void decompose(sgraph *g, int* colmap, ds::worklist& vertex_to_component, int new_num_components) {
            // set up forward / backward maps
            num_components = new_num_components; // new_num_components
            if(num_components <= 1) return;

            std::vector<int> vertices_in_component;
            vertices_in_component.resize(num_components);

            std::vector<int> forward_translation;
            forward_translation.resize(g->v_size);

            int backward_size = 0;

            for(int v = 0; v < g->v_size; ++v) {
                const int component = vertex_to_component[v];
                if(component  < 0) {
                    forward_translation[v] = -1;
                    //assert(false);
                    continue;
                }
                const int v_in_component = vertices_in_component[component];
                ++vertices_in_component[component];
                ++backward_size;
                forward_translation[v] = v_in_component;
            }

            backward_translation.resize(backward_size);
            component_to_backward_translation.reserve(num_components);

            int current_pos = 0;
            for(int i = 0; i < num_components; ++i) {
                const int component_size = vertices_in_component[i];
                component_to_backward_translation.push_back(current_pos);
                current_pos += component_size;
            }

            for(int v = 0; v < g->v_size; ++v) {
                const int component = vertex_to_component[v];
                if(component  < 0) continue;
                const int v_in_component = forward_translation[v];
                assert(v_in_component >= 0);
                assert(v_in_component < vertices_in_component[component]);
                backward_translation[component_to_backward_translation[component] + v_in_component] = v;
            }

            for(int v = 0; v < g->v_size; ++v) {
                const int deg = g->d[v];
                const int vpt = g->v[v];
                int ept = vpt;
                for (int j = vpt; j < vpt + deg; ++j) {
                    const int neighbour = g->e[j];
                    const int fwd_neighbour = forward_translation[neighbour];
                    if(fwd_neighbour >= 0) {
                        g->e[ept] = fwd_neighbour;
                        assert(vertex_to_component[v] == vertex_to_component[neighbour]);
                        //assert(g->e[ept] >= 0 && g->e[ept] < vertices_in_component[vertex_to_component[ept]]);
                        ++ept;
                    }
                }
            }

            std::vector<int> original_v;
            std::vector<int> original_d;
            original_v.reserve(g->v_size);
            original_d.reserve(g->v_size);
            for(int v = 0; v < g->v_size; ++v) original_v.push_back(g->v[v]);
            for(int v = 0; v < g->v_size; ++v) original_d.push_back(g->d[v]);

            for(int v = 0; v < g->v_size; ++v) {
                const int component = vertex_to_component[v];
                if(component  < 0) continue;
                const int v_in_component = forward_translation[v];
                const int bw_translate = component_to_backward_translation[component];
                g->v[bw_translate + v_in_component] = original_v[v];
                g->d[bw_translate + v_in_component] = original_d[v];
            }

            for(int v = 0; v < g->v_size; ++v) original_v[v] = colmap[v];

            components_graph.reserve(num_components);
            for(int i = 0; i < num_components; ++i) {
                const int bw_translate = component_to_backward_translation[i];
                components_graph.emplace_back();
                sgraph* component_i = &components_graph[components_graph.size() - 1];
                component_i->v_size = vertices_in_component[i];
                component_i->e_size = g->e_size;
                component_i->v = g->v + bw_translate;
                component_i->d = g->d + bw_translate;
                component_i->e = g->e;
                component_i->initialized = false;

                for(int j = 0; j < vertices_in_component[i]; ++j) {
                    colmap[bw_translate + j] = original_v[backward_translation[bw_translate + j]];
                }
                components_coloring.push_back(colmap + bw_translate);
            }

        }

        /**
         * Returns the graph of component \p i.
         *
         * @param i number of component
         * @return graph of component \p i
         */
        sgraph* get_component(int i) {
            return &components_graph[i];
        }

        /**
         * Returns the vertex coloring of component \p i.
         *
         * @param i number of component
         * @return vertex coloring of component \p i
         */
        int* get_colmap(int i) {
            return components_coloring[i];
        }
    };
}

#endif //DEJAVU_COMPONENTS_H
