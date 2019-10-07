#ifndef BRUTUS_INVARIANT_H
#define BRUTUS_INVARIANT_H


#include <stack>
#include <vector>

class invariant {
    std::vector<std::vector<int>> vec_invariant;
    invariant* compareI;
    bool has_compare = false;
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

    bool write_top_and_compare(int i);
    bool compare_sizes();

    void set_compare_invariant(invariant *I);
};


#endif //BRUTUS_INVARIANT_H
