# Documentation

This is the source code documentation for dejavu 2.0. It contains detailed descriptions for all classes and methods.
If you are looking for a guide to get started, please instead refer to our [get started guide](https://automorphisms.org) on our main page.

## Potentially Interesting Pages

Below, we list some potentially interesting pages of this documentation.

* [solver](@ref dejavu::solver): the dejavu solver, used to compute automorphisms of a graph 
* [static_graph](@ref dejavu::static_graph): the graph interface
* [refinement](@ref dejavu::ir::refinement): the color refinement algorithm
* [orbit](@ref dejavu::groups::orbit): can be used to compute orbits
* [hooks](@ref dejavu::hooks): hooks are used to interact with the computed symmetries
* [random_schreier](@ref dejavu::groups::random_schreier): a implementation of the random Schreier algorithm

## Bug Reports & Feedback

If you come across any bugs or have any feedback to share, please always feel free to reach out to me at `markus (at) automorphisms.org`.

## How does it work?

An up-to-date description of the algorithms will be available in Markus Anders' PhD thesis.

The underlying algorithms are based on a series of papers by Markus Anders and Pascal Schweitzer, in particular

* Search Problems in Trees with Symmetries, ICALP 2021
* Engineering a Fast Probabilistic Isomorphism Test, ALENEX 2021
* Parallel Computation of Combinatorial Symmetries, ESA 2021

The preprocessing routines are described in

* Engineering a Preprocessor for Symmetry Detection, SEA 2023, additionally co-authored by Julian Stieß

## Copyright & License

The solver is released under the MIT license. For more information, see the [license](LICENSE_source.html).

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
