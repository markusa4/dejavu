#ifndef DEJAVU_TRACE_H
#define DEJAVU_TRACE_H

#include <vector>
#include <cassert>

namespace dejavu {
    namespace ir {
        #define TRACE_MARKER_INDIVIDUALIZE     (-1)
        #define TRACE_MARKER_REFINE_START      (-2)
        #define TRACE_MARKER_REFINE_END        (-3)
        #define TRACE_MARKER_REFINE_CELL_START (-4)
        #define TRACE_MARKER_REFINE_CELL_END   (-5)

        /**
         * \brief The trace invariant.
         *
         * Class that serves to store and compare the trace of a walk in an individualization-refinement tree. The class
         * provides several different modes in which information is recorded and/or compared.
         *
         * Specifically, it is possible to (1) record a full trace, (2) compare to a full trace, or (3) compare to a full
         * trace and recording a hash invariant as soon as the new computation deviates from the stored trace (see also
         * \ref ir_mode ).
         *
         * While comparing to a stored trace (2, 3), the class facilitates the use of the blueprint heuristic, which enables
         * skipping of non-splitting cells in the stored trace.
        */
        class trace {
            // the trace
        protected:
            std::vector<int> data; /**< keeps all the data of the trace */
        private:
            trace *compare_trace = nullptr; /**< link to a stored trace to compare to */
            long hash = 0; /**< hash value to summarize all operations performed on this trace */

            // mode
            bool compare = false; /**< whether to compare operations to a stored trace*/
            bool freeze = false;
            bool record = false; /**< whether to record a trace */

            // housekeeping
            int cell_act_spot = -1;
            int cell_old_color = -1;
            bool assert_cell_act = false;
            bool assert_refine_act = false;

            // comparison variables
            int position = 0;
            bool comp = true;

            void add_to_hash(int d) {
                hash += (d * (352355 - d * 3));
            }

            void write_compare(int d) {
                add_to_hash(d);
                if (record)
                    data.push_back(d);
                if (compare) {
                    if(position < compare_trace->data.size()) {
                        assert(compare_trace->data.size() > position);
                        comp = comp && compare_trace->data[position] == d;
                    } else {
                        comp = false;
                    }
                }
                ++position;
            }

            void write_skip_compare(int d) {
                if (record)
                    data.push_back(d);
                ++position;
            }

        public:

            void update_blueprint_hash() {
                bool skipping_cell = false;
                bool write_false_next = false;
                hash = 0;
                for (int i = 0; i < data.size(); ++i) {
                    const int dt = data[i];
                    if(dt == TRACE_MARKER_REFINE_END)
                        continue;
                    if (dt == TRACE_MARKER_REFINE_CELL_START) {
                        if (data[i + 1] == false) {
                            skipping_cell = true;
                        } else {
                            write_false_next = true;
                        }
                    }
                    if (!skipping_cell) {
                        if (write_false_next) {
                            // uses write_skip_compare
                            write_false_next = false;
                        } else {
                            add_to_hash(dt);
                        }
                    }
                    if (dt == TRACE_MARKER_REFINE_CELL_END) {
                        skipping_cell = false;
                        write_false_next = false;
                    }
                }
            }

            /**
             * Records an individualization.
             * @param color The color being individualized.
             */
            void op_individualize(int old_color, int ind_color) {
                assert(ind_color >= 0);
                assert(old_color >= 0);
                assert(ind_color != old_color);
                write_compare(TRACE_MARKER_INDIVIDUALIZE);
                write_compare(old_color);
                write_compare(ind_color);
            }

            /**
             * Records the start of a refinement.
             */
            void op_refine_start() {
                assert(!comp || !assert_refine_act);
                write_compare(TRACE_MARKER_REFINE_START);
                assert_refine_act = true;
            }

            /**
             * Records the start of a refinement with respect to a color.
             * @param color The color in respect to which the coloring is refined.
             */
            void op_refine_cell_start(int color) {
                assert(!comp || !assert_cell_act);
                write_compare(TRACE_MARKER_REFINE_CELL_START);
                cell_old_color = color;
                cell_act_spot = data.size();
                write_skip_compare(false);
                assert_cell_act = true;
            }

            /**
             * Records a that a new color appeared while refining with respect to a color.
             * @param new_color The new color that was refined.
             * @param new_color_size The size of the new color.
             * @param new_color_deg The color degree which lead to the new color class.
             */
            void op_refine_cell_record(int new_color, int new_color_size, int new_color_deg) {
                assert(!comp || assert_cell_act);
                write_compare(new_color);
                //write_compare(new_color_size);
                //write_compare(new_color_deg);
                if (new_color != cell_old_color && record)
                    data[cell_act_spot] = true;
            }

            void op_additional_info(long d) {
                write_compare(d);
            }

            /**
             * Records the end of a refinement with respect to a color.
             */
            void op_refine_cell_end() {
                assert(!comp || assert_cell_act);
                assert_cell_act = false;
                write_compare(TRACE_MARKER_REFINE_CELL_END);
            }

            /**
             * Records the end of a refinement.
             */
            void op_refine_end() {
                assert(!comp || !assert_cell_act);
                assert(assert_refine_act);
                assert_refine_act = false;
                write_skip_compare(TRACE_MARKER_REFINE_END);
            }

            /**
             * Only applicable if we are currently comparing to a stored trace (if \a compare is set).
             *
             * @return Determines whether in the stored trace, the next color in respect to which was refined created new colors
             * (i.e., whether the next color is splitting).
             */
            bool blueprint_is_next_cell_active() {
                if (!compare)
                    return true;

                assert(compare_trace);
                size_t read_pt = position;
                assert(compare_trace->data.size() > read_pt);
                assert(compare_trace->data[read_pt] == TRACE_MARKER_REFINE_CELL_START);
                ++read_pt;
                assert(compare_trace->data.size() > read_pt);
                assert(compare_trace->data[read_pt] == false || compare_trace->data[read_pt] == true);
                return compare_trace->data[read_pt];
            }

            /**
             * Only applicable if we are currently comparing to a stored trace (if \a compare is set).
             *
             * Skips the \a position to the start of the next refinement with respect to a color. To be used after
             * \a blueprint_is_next_cell_active() determined the current color to be non-splitting.
             */
            void blueprint_skip_to_next_cell() {
                assert(compare_trace->data[position] == TRACE_MARKER_REFINE_CELL_START);
                while (compare_trace->data[position] != TRACE_MARKER_REFINE_CELL_END) {
                    assert(compare_trace->data.size() > (size_t) position);
                    ++position;
                }
                ++position;
            }

            /**
             * Rewinds the \a position to the previous individualization.
             */
            void rewind_to_individualization() {
                if (record) {
                    int read_pt = data.size() - 1;
                    while (read_pt >= 0 && data[read_pt] != TRACE_MARKER_INDIVIDUALIZE) {
                        --read_pt;
                    }
                    data.resize(read_pt);
                    position = read_pt;
                }
                if (compare) {
                    int read_pt = position - 1;
                    while (read_pt >= 0 && compare_trace->data[read_pt] != TRACE_MARKER_INDIVIDUALIZE) {
                        --read_pt;
                    }
                    position = read_pt;
                }
            }

            /**
             * Only applicable if we are currently comparing to a stored trace (if \a compare is set).
             *
             * Skips the \a position to the next individualization in the \a compare_trace.
             */
            void skip_to_individualization() {
                assert_cell_act = false;
                assert_refine_act = false;
                if (compare) {
                    int read_pt = position - 1;
                    while ((size_t) read_pt < compare_trace->data.size() &&
                           compare_trace->data[read_pt] != TRACE_MARKER_INDIVIDUALIZE) {
                        ++read_pt;
                    }
                    position = read_pt;
                }
            }

            /**
             * Stores a new trace to compare to.
             *
             * @param compare_trace
             */
            void set_compare_trace(trace *compare_trace) {
                this->compare_trace = compare_trace;
                compare_trace->freeze = true;
            }

            /**
             * Determines whether we want to compare the following operations with those of the trace stored in
             * \a compare_trace.
             *
             * @param compare
             */
            void set_compare(bool compare) {
                this->compare = compare;
            }

            /**
             * @return A hash value summarizing the operations recorded in this trace.
             */
            long get_hash() {
                return hash;
            }

            /**
             * Sets the hash value to a pre-determined value.
             * @param hash The hash value.
             */
            void set_hash(long hash) {
                this->hash = hash;
            }

            /**
             * @return Whether the recorded operations deviated from the stored trace in \a compare_trace.
             */
            bool trace_equal() {
                return comp;
            }

            /**
             * Resets the trace to find new deviations from the stored trace.
             */
            void reset_trace_equal() {
                comp = true;
            }

            void reset() {
                data.clear();
                compare_trace = nullptr;
                position = 0;
                hash = 0;
                reset_trace_equal();
            }

            // record to trace
            void set_record(bool record) {
                this->record = record;
            }

            // position the trace
            void set_position(int position) {
                assert_cell_act   = false;
                assert_refine_act = false;
                this->position = position;
            }

            int get_position() {
                return position;
            }
        };
    }
}


#endif //DEJAVU_TRACE_H