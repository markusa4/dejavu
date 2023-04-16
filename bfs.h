#ifndef DEJAVU_BFS_H
#define DEJAVU_BFS_H

#include "ir.h"

namespace dejavu {
    namespace search_strategy {

        /**
         * \brief Breadth-first search.
         */
        class bfs_ir {
            // TODO deviation map should go into shared_tree as well!
            std::unordered_set<long> deviation_map;
            int computed_for_base = 0;
            int expected_for_base = 0;
            bool deviation_done = false;
            // TODO: bfs_ir should be workspace for local thread, and not a shared structure!

            int s_deviation_prune = 0;
            int s_total_prune     = 0;
            int s_total_kept      = 0;

        public:
            void do_a_level(refinement* R, sgraph* g, ir::shared_tree& ir_tree, ir::controller& local_state, std::function<ir::type_selector_hook> *selector) {
                int current_level = ir_tree.get_finished_up_to();

                s_deviation_prune = 0;
                s_total_prune     = 0;
                s_total_kept      = 0;

                queue_up_level(selector, ir_tree, current_level);
                work_on_todo(R, g, &ir_tree, local_state);
                ir_tree.set_finished_up_to(current_level + 1);
                //std::cout << s_deviation_prune << "/" << s_total_prune << " - " << s_total_kept << std::endl;
            }

            int next_level_estimate(ir::shared_tree& ir_tree, std::function<ir::type_selector_hook> *selector) {
                const int base_pos = ir_tree.get_finished_up_to();
                const auto start_node = ir_tree.get_level(base_pos);
                const auto level_size = ir_tree.get_level_size(base_pos);
                auto next_node_save = start_node->get_save();
                auto c = next_node_save->get_coloring();
                auto base_pos_   = next_node_save->get_base_position();
                int col = (*selector)(c, base_pos_);
                return level_size * (c->ptn[col] + 1);
            }

            void queue_up_level(std::function<ir::type_selector_hook> *selector, ir::shared_tree& ir_tree, int base_pos) {
                auto start_node = ir_tree.get_level(base_pos);
                while(!start_node->get_base()) {
                    start_node = start_node->get_next();
                }
                start_node = start_node->get_next();

                const auto level_size = ir_tree.get_level_size(base_pos);
                auto next_node = start_node;
                bool reserve = false;
                reset_deviation_map();

                do {
                    auto next_node_save = next_node->get_save();
                    auto c = next_node_save->get_coloring();
                    auto base_pos    = next_node_save->get_base_position();
                    int col = (*selector)(c, base_pos);
                    if(!reserve && col >= 0) {
                        expected_for_base = c->ptn[col] + 1;
                        ir_tree.queue_reserve((c->ptn[col] + 1) * level_size);
                        reserve = true;
                    }

                    if(col >= 0) {
                        for (int i = 0; i < c->ptn[col] + 1; ++i) ir_tree.queue_missing_node({next_node, c->lab[col + i]});
                    }
                    next_node = next_node->get_next();
                } while(next_node != start_node);
            }

            void reset_deviation_map() {
                deviation_map.clear();
                computed_for_base = 0;
                deviation_done = false;
            }

            void add_deviation(long hash) {
                deviation_map.insert(hash);
            }

            void finish_deviation() {
                deviation_done = true;
            }

            bool check_deviation(long hash) {
                return !deviation_done || deviation_map.contains(hash);
            }

            void compute_node(refinement* R, sgraph* g, ir::shared_tree* ir_tree, ir::controller& local_state, ir::tree_node* node, const int v, ir::reduced_save* last_load) {
                auto next_node_save = node->get_save();

                // node is already pruned
                const bool is_pruned = node->get_prune();
                if(is_pruned) {
                    ++s_total_prune;
                    ++s_deviation_prune;
                    return;
                }

                // do efficient loading if parent is the same as previous load
                if(next_node_save != last_load || g->v_size < 500) {
                    local_state.load_reduced_state(*next_node_save);
                } else {
                    local_state.move_to_parent();
                    local_state.load_reduced_state_without_coloring(*next_node_save);
                }

                if(local_state.base_pos > 0) local_state.use_increase_deviation_hash(true);

                // do computation
                local_state.reset_trace_equal();
                if(g->v_size >= 500) local_state.use_limited_reversible_for_next();
                local_state.use_trace_early_out(true);
                local_state.move_to_child(R, g, v);

                // we want to keep track of whether we are on the base or not
                const bool parent_is_base = node->get_base();
                const bool is_base = parent_is_base && (v == local_state.compare_base[local_state.base_pos-1]);

                if(local_state.T->trace_equal()) { // TODO: what if leaf?
                    ++s_total_kept;
                    auto new_save = new ir::reduced_save();
                    local_state.save_reduced_state(*new_save);
                    ir_tree->add_node(local_state.base_pos, new_save, is_base);
                } else {
                    // deviation map
                    if(local_state.base_pos > 1) {
                        ++s_total_prune;
                        if (parent_is_base) add_deviation(local_state.T->get_hash());
                        else {
                            if (!check_deviation(local_state.T->get_hash())) {
                                node->prune();
                            }
                        }
                    }
                }

                // keep track how many we computed for deviation map
                if(parent_is_base && local_state.base_pos > 1) {
                    computed_for_base += 1;
                    if(computed_for_base == expected_for_base) {
                        finish_deviation();
                    }
                }
            }

            void work_on_todo(refinement* R, sgraph* g, ir::shared_tree* ir_tree, ir::controller& local_state) {
                ir::reduced_save* last_load = nullptr;
                while(!ir_tree->queue_missing_node_empty()) {
                    const auto todo = ir_tree->queue_missing_node_pop();
                    compute_node(R, g, ir_tree, local_state, todo.first, todo.second, last_load);
                    last_load = todo.first->get_save();
                }
            }
        };
    }
}

#endif //DEJAVU_BFS_H
