#ifndef DEJAVU_TRACE_H
#define DEJAVU_TRACE_H

#include<vector>
#include "assert.h"

#define TRACE_MARKER_INDIVIDUALIZE     -1
#define TRACE_MARKER_REFINE_START      -2
#define TRACE_MARKER_REFINE_END        -3
#define TRACE_MARKER_REFINE_CELL_START -4
#define TRACE_MARKER_REFINE_CELL_END   -5

namespace dejavu {
    class trace {
    protected:
        std::vector<int> data;

    private:
        // the trace
        trace *compare_trace = nullptr;
        std::vector<long> hash;

        // mode
        bool compare = false;
        bool freeze = false;
        bool record = false;

        // housekeeping
        int cell_act_spot = -1;
        int cell_old_color = -1;
        bool assert_cell_act = false;
        bool assert_refine_act = false;

        // comparison variables
        int  position = 0;
        bool comp = true;

        void write_compare(int d) {
            if (record)
                data.push_back(d);
            if (compare) {
                assert(compare_trace->data.size() > position);
                comp = comp && compare_trace->data[position] == d;
            }
            ++position;
        }
    public:
        // recording & comparing
        void op_individualize(int color) {
            assert(color >= 0);
            write_compare(TRACE_MARKER_INDIVIDUALIZE);
            write_compare(color);
        }

        void op_refine_start() {
            assert(!assert_refine_act);
            write_compare(TRACE_MARKER_REFINE_START);
            assert_refine_act = true;
        }

        void op_refine_cell_start(int old_color) {
            assert(!assert_cell_act);
            write_compare(TRACE_MARKER_REFINE_CELL_END);
            cell_old_color = old_color;
            cell_act_spot = data.size();
            write_compare(false);
            assert_cell_act = true;
        }

        void op_refine_cell_record(int new_color, int new_color_size, int new_color_deg) {
            assert(assert_cell_act);
            write_compare(new_color);
            write_compare(new_color_size);
            //write_compare(new_color_deg);
            if (new_color != cell_old_color && record)
                data[cell_act_spot] = true;
        }

        void op_refine_cell_end() {
            assert(assert_cell_act);
            assert_cell_act = false;
            write_compare(TRACE_MARKER_REFINE_CELL_END);
        }

        void op_refine_end() {
            assert(assert_refine_act);
            assert_refine_act = false;
            write_compare(TRACE_MARKER_REFINE_END);
        }

        // blueprint usage
        bool blueprint_is_next_cell_active() {
            assert(compare_trace);
            size_t read_pt = position;
            assert(compare_trace->data.size() > read_pt);
            assert(compare_trace->data[read_pt] == TRACE_MARKER_REFINE_CELL_START);
            ++read_pt;
            assert(compare_trace->data.size() > read_pt);
            assert(compare_trace->data[read_pt] == false || compare_trace->data[read_pt] == true);
            return compare_trace->data[read_pt];
        }

        void blueprint_skip_to_next_cell() {
            assert(compare_trace->data[position] == TRACE_MARKER_REFINE_CELL_START);
            while(compare_trace->data[position] != TRACE_MARKER_REFINE_CELL_END) {
                assert(compare_trace->data.size() > (size_t) position);
                ++position;
            }
            ++position;
        }

        // rewind the invariant to last individualization
        void rewind_to_individualization() {
            if(record) {
                int read_pt = data.size() - 1;
                while(read_pt >= 0 && data[read_pt] != TRACE_MARKER_INDIVIDUALIZE) {
                    --read_pt;
                }
                data.resize(read_pt);
            }
            if(compare) {
                int read_pt = position - 1;
                while(read_pt >= 0 && compare_trace->data[read_pt] != TRACE_MARKER_INDIVIDUALIZE) {
                    --read_pt;
                }
                position = read_pt;
            }
        }

        // compare to other traces
        void set_compare_trace(trace *compare_trace) {
            this->compare_trace = compare_trace;
            compare_trace->freeze = true;
        }

        void set_compare(bool compare) {
            this->compare = compare;
        }

        // does this deviate from compare trace yet?
        bool trace_equal() {
            return comp;
        }

        // reset comparison failure
        void reset_trace_equal() {
            comp = true;
        }

        // record to trace
        void set_record(bool record) {
            this->record = record;
        }

        // position the trace
        void set_position(int position) {
            this->position = position;
        }

        int get_position() {
            return position;
        }
    };
}


#endif //DEJAVU_TRACE_H