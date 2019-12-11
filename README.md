# Compilation
The implementation has binary dependencies to nauty. The easiest way to compile dejavu is to place a compiled version of nauty / Traces (available a thttp://pallini.di.uniroma1.it/) in "nauty/".

The solver itself technically only depends on "nauty.h" and "naurng.h", utilizing some of the included data structures and randomization functions. The distribution does however also include the benchmark code used for the benchmarks of the paper, which call nauty and Traces.

# Usage
