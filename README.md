# Compilation
The implementation has binary dependencies to nauty. The easiest way to compile dejavu is to place a compiled version of nauty / Traces (available at http://pallini.di.uniroma1.it/) in "nauty/".

The solver itself technically only depends on "nauty.h" and "naurng.h", utilizing some of the included data structures and randomization functions. The distribution does however also include the benchmark code used for the benchmarks of the paper, which call nauty and Traces.

# Usage
Compilation produces two binaries (dejavu and bench), which are however similar in usage. dejavu is the binary of the actual solver. bench is basically just a frontend which can call nauty, Traces and dejavu, while recording and tracking time measurements, as well as manage timeouts. 

## dejavu
The solver only accepts files in the DIMACS graph format. 

## bench

