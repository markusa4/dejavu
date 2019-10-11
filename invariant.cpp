#include <assert.h>
#include <iostream>
#include "invariant.h"


void invariant::set_compare_invariant(invariant* I) {
    no_write = true;
    has_compare = true;
    compareI = I;
}