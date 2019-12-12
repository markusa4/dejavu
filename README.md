# Compilation
The implementation has binary dependencies to nauty. The easiest way to compile dejavu is to place a configured and compiled distribution of nauty / Traces (available at http://pallini.di.uniroma1.it/) in "nauty/".

The solver itself technically only depends on "nauty.h" and "naurng.h", utilizing some of the included data structures and randomization functions. The distribution does however also include the benchmark code used for the benchmarks of the paper, which call nauty and Traces.

After placing nauty and Traces into the respective folder, the project can be compiled by simply running cmake.

# Usage
Compilation produces two binaries (dejavu and bench), which are however similar in usage. dejavu is the binary of the actual solver. bench is basically just a frontend which can call nauty, Traces and dejavu, while recording and tracking time measurements, as well as manage timeouts. 

## dejavu
The solver only accepts files in the DIMACS graph format. The solver accepts the following command line arguments:

Command Line Argument | Effect
--- | ---
`--file` | specifies the graph input file in DIMACS format
`--threads` | specifies the amount of additional threads used (total is threads + 1, default is 0)
`--write_auto` | writes all of the generators of the found automorphism group
`--timeout` | specifies a timeout (in seconds)
`--force_selector` | forces a specific cell selector (0 is first, 1 is smallest, 2 is largest)
`--permute` | permutes the input graph randomly
`--permute_seed` | specify a fixed seed for the previous argument

## bench
`bench` understands the same command line arguments, with the addition of the following: 

Command Line Argument | Effect
--- | ---
`--no_nauty` | does not run nauty
`--no_traces` | does not run Traces
`--no_dejavu` | does not run dejavu
`--stat_file` | specify a file to which measurements shall be written