#ifndef BRUTUS_INVARIANT_H
#define BRUTUS_INVARIANT_H


#include <stack>
#include <vector>

class invariant {
    std::vector<std::vector<int>> vec_invariant;
    invariant* compareI;
    bool has_compare = false;
    bool no_write = false;
    int cur_pos = -1;
    int fake_sz = 0;
    std::vector<int>* compare_level;
public:
    virtual void push_level();
    virtual void pop_level();
    virtual int current_level();
    virtual std::vector<int>* get_level(int i);
    virtual std::vector<int> top_level();
    virtual int top_is_geq(std::vector<int> *other);
    virtual void write_top(int i);
    virtual void print();
    virtual bool top_is_eq(std::vector<int> *other);

    bool level_is_eq(invariant *other, int level);

    inline bool write_top_and_compare(int i) {
        if(no_write) {
            int pos2 = fake_sz;
            fake_sz += 1;
            return (compare_level->size() > pos2) && (i == (*compare_level)[pos2]);
        } else {
            vec_invariant[cur_pos].push_back(i);
            int pos2 = vec_invariant[cur_pos].size() - 1;
            if (has_compare) {
                if ((*compareI->get_level(cur_pos)).size() <= pos2)
                    return false;
                return vec_invariant[cur_pos][pos2] == (*compareI->get_level(cur_pos))[pos2];
            } else {
                return true;
            }
        }
    }

    bool compare_sizes();

    void set_compare_invariant(invariant *I);
};


#endif //BRUTUS_INVARIANT_H
