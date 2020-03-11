#ifndef DEJAVU_INVARIANT_H
#define DEJAVU_INVARIANT_H


#include <stack>
#include <vector>
#include <iostream>
#include "assert.h"
#include "utility.h"

enum cell_state {CELL_ACTIVE, CELL_IDLE, CELL_END};

class alignas(16) invariant {
public:
    std::vector<int>*        vec_cells     = nullptr;
    std::vector<int>*        vec_selections= nullptr;
    std::vector<cell_state>* vec_protocol  = nullptr;
    std::vector<int>*        vec_invariant = nullptr;

    invariant*        compareI;
    std::vector<int>* compare_vec;
    bool has_compare = false;
    bool no_write    = false;
    bool never_fail  = false;
    int  comp_fail_pos = -2;
    int  comp_fail_val = -1;
    long comp_fail_acc = -1;
    int  cur_pos = -1;
    int  protocol_pos = -1;
    long acc     = 0;

    // currently a bit convoluted, really should be split into 2 functions...
    inline bool write_top_and_compare(int i) {
        acc += MASH0(i) ;
        if(no_write) {
            const bool comp = (comp_fail_pos == -2) && ((*compare_vec)[++cur_pos] == i);
            if(!comp) {
                if(comp_fail_pos == -2) {
                    comp_fail_pos = cur_pos;
                    comp_fail_val = i;
                    comp_fail_acc = i;
                }
                comp_fail_acc += MASH1(i); // could just use acc instead
            }
            return (comp || never_fail);
        } else {
            vec_invariant->push_back(i);
            cur_pos += 1;
            if (has_compare) {
                if ((compareI->vec_invariant)->size() < vec_invariant->size())
                    return false;
                return (*vec_invariant)[cur_pos] == (*compareI->vec_invariant)[cur_pos];
            } else {
                return true;
            }
        }
    }

    inline void write_cells(int cells) {
        if(!no_write) {
            vec_cells->push_back(cells);
        }
    }

    inline void selection_write(int cell) {
        if(!no_write) {
            vec_selections->push_back(cell);
        }
    }

    inline int selection_read(int base_point) {
        return (*vec_selections)[base_point];
    }

    inline void protocol_write(bool active_cell, int c) {
        if(!no_write) {
            ++protocol_pos;
            vec_protocol->push_back(active_cell?cell_state::CELL_ACTIVE:cell_state::CELL_IDLE);
        }
    }

    inline void protocol_mark() {
        if(!no_write) {
            ++protocol_pos;
            vec_protocol->push_back(CELL_END);
        }
    }

    inline bool protocol_read(int c) {
        ++protocol_pos;
        return (*compareI->vec_protocol)[protocol_pos] == cell_state::CELL_ACTIVE;
    }

    inline void protocol_skip_to_mark() {
        while((*compareI->vec_protocol)[protocol_pos] != cell_state::CELL_END) {
            ++protocol_pos;
        }
    }

    inline void fast_forward(int abort_marker) {
        while((*compare_vec)[cur_pos] != abort_marker) {
            write_top_and_compare((*compare_vec)[cur_pos + 1]);
        }
    }

    void reset_deviation() {
        comp_fail_pos = -2;
        comp_fail_val = -1;
        comp_fail_acc = -1;
    }

    void set_compare_invariant(invariant *I);

    void create_vector(int prealloc) {
        vec_cells     = new std::vector<int>();
        vec_protocol  = new std::vector<cell_state>();
        vec_invariant = new std::vector<int>();
        vec_selections= new std::vector<int>();
        vec_cells->reserve(prealloc + 16);
        vec_protocol->reserve(prealloc + 16);
        vec_invariant->reserve(prealloc * 20);
        vec_selections->reserve(prealloc + 16);
    }

    /*~invariant() {
        if(vec_invariant != nullptr)
            delete vec_invariant;
    }*/
};


#endif //DEJAVU_INVARIANT_H
