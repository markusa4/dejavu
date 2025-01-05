# Compilation, Library, Tests
dejavu is a solver and C++ library for the fast detection and manipulation of [combinatorial symmetry](https://automorphisms.org/quick_start/symmetry/). 
Below, you can find some information on how to build the solver and include the library.
More detailed information can be found in our [get started guide](https://automorphisms.org/) or in the [full documentation](https://automorphisms.org/documentation/), which can also be built from the source code using [doxygen](https://www.doxygen.nl/).


## Compilation
Using *cmake*, the project should compile without any further dependencies:
```text
cmake .
make
```
Compilation produces a binary *dejavu*. It accepts a DIMACS graph as input, and computes the automorphism group of the graph. For available options and more descriptions, please refer to our [guide](https://automorphisms.org/quick_start/standalone/).

## Use dejavu as a library
dejavu is a header-only library. You can simply add dejavu to your C++ project by including the respective header file: 
```cpp
#include "dejavu.h"
```

Note that currently, dejavu requires to be *compiled with C++ version 14*. For a more thorough description, please refer to our [guide](https://automorphisms.org/quick_start/cpp_api/).

By default, dejavu is compiled without assertions. We recommend activating assertions for debugging purposes (by adding the definition `DEJDEBUG`). Assertions do however slow the code considerably.

## Running the tests
Using *cmake*, a test target `dejavu_test` can be produced by setting the following flag:
```text
cmake . -DCOMPILE_TEST_SUITE=1
```

In order to run all the tests, the [test graphs](https://automorphisms.org/graphs/graphs.zip) are required to be placed into `tests/graphs/`.
