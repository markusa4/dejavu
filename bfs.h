// Copyright 2023 Markus Anders
// This file is part of dejavu 2.0.
// See LICENSE for extended copyright information.

#ifndef DEJAVU_BFS_H
#define DEJAVU_BFS_H

#include "ir.h"

namespace dejavu {
    namespace search_strategy {

        /**
         * \brief Breadth-first search.
         */
        class bfs_ir {
            groups::automorphism_workspace* gws_automorphism = nullptr;

        public:
            bool h_use_deviation_pruning = true; /**< use pruning using deviation maps */

            // TODO some of this should go into shared_tree
            // statistics
            int s_total_prune              = 0; /**< how many nodes were pruned on last level */
            int s_total_kept               = 0; /**< how many nodes were kept, i.e., not pruned */
            int s_total_automorphism_prune = 0; /**< how many nodes were pruned using automorphism pruning */
            int s_total_leaves             = 0; /**< how many of the computed nodes were leaves */
            int s_deviation_prune          = 0; /**< how many nodes were pruned using deviation maps */

            void link_to_workspace(groups::automorphism_workspace* automorphism) {
                gws_automorphism = automorphism;
            }

            void do_a_level(sgraph* g, dejavu_hook* hook, ir::shared_tree& ir_tree, ir::controller& local_state, std::function<ir::type_selector_hook> *selector) {
                int current_level = ir_tree.get_finished_up_to();

                s_deviation_prune = 0;
                s_total_prune     = 0;
                s_total_kept      = 0;
                s_total_automorphism_prune = 0;
                s_total_leaves = 0;

                assert(ir_tree.get_current_level_size() > 0);

                queue_up_level(selector, ir_tree, current_level);
                work_on_todo(g, hook, &ir_tree, local_state);
                ir_tree.set_finished_up_to(current_level + 1);
            }

            static int next_level_estimate(ir::shared_tree& ir_tree, std::function<ir::type_selector_hook> *selector) {
                const int base_pos = ir_tree.get_finished_up_to();
                const auto start_node = ir_tree.get_level(base_pos);
                assert(start_node != nullptr);
                const auto level_size = ir_tree.get_level_size(base_pos);
                auto next_node_save = start_node->get_save();
                auto c = next_node_save->get_coloring();
                auto base_pos_   = next_node_save->get_base_position();
                int col = (*selector)(c, base_pos_);
                return level_size * (c->ptn[col] + 1);
            }

            static void queue_up_level(std::function<ir::type_selector_hook> *selector, ir::shared_tree& ir_tree,
                                       int base_pos) {
                auto start_node = ir_tree.get_level(base_pos);
                assert(start_node != nullptr);
                while(!start_node->get_base()) {
                    start_node = start_node->get_next();
                }
                start_node = start_node->get_next();

                const auto level_size = ir_tree.get_level_size(base_pos);
                auto next_node = start_node;
                bool reserve = false;

                do {
                    auto next_node_save = next_node->get_save();
                    auto c = next_node_save->get_coloring();
                    auto this_base_pos    = next_node_save->get_base_position();
                    int col = (*selector)(c, this_base_pos);
                    if(!reserve && col >= 0) {
                        //expected_for_base = c->ptn[col] + 1;
                        ir_tree.stored_deviation.start(c->ptn[col] + 1);
                        ir_tree.queue_reserve((c->ptn[col] + 1) * level_size);
                        reserve = true;
                    }

                    if(col >= 0) {
                        for (int i = 0; i < c->ptn[col] + 1; ++i) {
                            ir_tree.queue_missing_node({next_node, c->lab[col + i]});
                        }
                    }
                    next_node = next_node->get_next();
                } while(next_node != start_node);
            }

            void compute_node(sgraph* g, dejavu_hook* hook, ir::shared_tree* ir_tree, ir::controller& local_state,
                              ir::tree_node* node, const int v, ir::limited_save* last_load) {
                auto next_node_save = node->get_save();

                // TODO consider base size 1 and top-level automorphisms

                // node is already pruned
                const bool is_pruned = node->get_prune();
                if(is_pruned && h_use_deviation_pruning) {
                    ++s_total_prune;
                    ++s_deviation_prune;
                    assert(!node->get_base());
                    return;
                }

                // special code for automorphism pruning on base size 2
                const int parent_node_base_pos  = node->get_save()->get_base_position()-1;
                const int parent_node_base_vert = parent_node_base_pos>=0?node->get_save()->get_base()[parent_node_base_pos]:-1;
                const int vert_on_base          = parent_node_base_pos>=0?local_state.compare_base[parent_node_base_pos]:-1;
                const int vert_on_base_sl       = parent_node_base_pos==-1?local_state.compare_base[0]:-1;
                if(parent_node_base_pos == 0 && !ir_tree->h_bfs_top_level_orbit.represents_orbit(parent_node_base_vert)) {
                    ++s_total_automorphism_prune;
                    return;
                }

                if(parent_node_base_pos == -1 && v != vert_on_base_sl &&
                   ir_tree->h_bfs_top_level_orbit.are_in_same_orbit(v, vert_on_base_sl)) {
                    ++s_total_automorphism_prune;
                    return;
                }

                // do efficient loading if parent is the same as previous load
                if(next_node_save != last_load || g->v_size < 2000) { // TODO heuristic to check how much has changed
                    local_state.use_reversible(false); // potentially loads more efficiently
                    local_state.load_reduced_state(*next_node_save);
                } else {
                    local_state.move_to_parent();
                    local_state.load_reduced_state_without_coloring(*next_node_save); // TODO <- this should be unecessary, right?
                }

                if(local_state.s_base_pos > 0) local_state.use_increase_deviation(true);

                // do computation
                local_state.reset_trace_equal();
                local_state.use_reversible(g->v_size >= 2000);
                local_state.use_trace_early_out(true);

                local_state.move_to_child(g, v);

                // we want to keep track of whether we are on the base or not
                const bool parent_is_base = node->get_base();
                const bool is_base = parent_is_base && (v == local_state.compare_base[local_state.s_base_pos - 1]);

                bool cert = true;

                if(g->v_size == local_state.c->cells && local_state.T->trace_equal()) {
                    gws_automorphism->write_color_diff(local_state.c->vertex_to_col, local_state.leaf_color.lab);
                    cert = local_state.certify(g, *gws_automorphism);
                    if(cert) {
                        ir_tree->h_bfs_top_level_orbit.add_automorphism_to_orbit(*gws_automorphism);

                        // Output automorphism
                        if (hook)
                            (*hook)(g->v_size, gws_automorphism->perm(), gws_automorphism->nsupport(),
                                    gws_automorphism->support());
                    }
                    ++s_total_leaves;
                    gws_automorphism->reset();

                    if(parent_node_base_pos == 0 && vert_on_base == parent_node_base_vert) ++ir_tree->h_bfs_automorphism_pw;
                }

                // TODO work on control flow below
                if(local_state.T->trace_equal() && cert) {
                    ++s_total_kept;
                    auto new_save = new ir::limited_save();
                    local_state.save_reduced_state(*new_save);
                    ir_tree->add_node(local_state.s_base_pos, new_save, node, is_base);
                    if(local_state.s_base_pos > 1) ir_tree->record_add_invariant(v, local_state.T->get_hash());
                } else {
                    assert(!is_base);
                    // deviation map
                    if(local_state.s_base_pos > 1) {
                        /*if(!h_use_deviation_pruning) {
                            const int first_level_v = node->get_save()->get_base()[0];
                            ir_tree->record_add_invariant(first_level_v, local_state.T->get_hash());
                            ir_tree->record_add_invariant(first_level_v, local_state.T->get_hash());
                        }*/
                        ++s_total_prune;
                        if (parent_is_base) ir_tree->stored_deviation.record_deviation(local_state.T->get_hash());
                        else {
                            if (!ir_tree->stored_deviation.check_deviation(local_state.T->get_hash())) {
                                assert(!parent_is_base);
                                node->prune();
                            }
                        }
                    } else {
                        ir_tree->record_invariant(v, local_state.T->get_hash());
                    }
                }

                // keep track how many we computed for deviation map
                if(parent_is_base && local_state.s_base_pos > 1 && local_state.T->trace_equal()) {
                    ir_tree->stored_deviation.record_no_deviation();
                }
            }

            void work_on_todo(sgraph* g, dejavu_hook* hook, ir::shared_tree* ir_tree, ir::controller& local_state) {
                ir::limited_save* last_load = nullptr;
                while(!ir_tree->queue_missing_node_empty()) {
                    const auto todo = ir_tree->queue_missing_node_pop();
                    compute_node(g, hook, ir_tree, local_state, todo.first, todo.second, last_load);
                    last_load = todo.first->get_save();
                }
            }
        };
    }
}

#endif //DEJAVU_BFS_H
