#ifndef BRUTUS_IR_TOOLS_H
#define BRUTUS_IR_TOOLS_H


#include <set>
#include <stack>
#include "bijection.h"
#include "sgraph.h"

class ir_tools {
public:
    void label_graph(sgraph* g, bijection* canon_p);
};

#endif //BRUTUS_IR_TOOLS_H
