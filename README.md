# Compilation
The benchmark utility has binary dependencies to nauty. The easiest way to compile it is to place a configured and compiled distribution of nauty / Traces (available at http://pallini.di.uniroma1.it/) in "nauty/". After placing nauty and Traces into the respective folder, the project can be compiled with cmake:
```
cmake -DBENCH=1 --build "."
make 
```
If only the dejavu solver is to be compiled, the `bench` target can be ignored and the project will compile without any further dependencies:
```
cmake .
make 
```

# Usage
Compilation produces two binaries (`dejavu` and `bench`), which are however similar in usage. dejavu is the binary of the actual solver. bench is basically just a frontend which can call nauty, Traces and dejavu, while recording and tracking time measurements, as well as manage timeouts. 

## dejavu
The solver only accepts files in the DIMACS graph format. The solver accepts the following command line arguments:

Command Line Argument | Effect
--- | ---
`--file` | specifies the graph input file in DIMACS format
`--threads` | specifies the amount of additional threads used (total is threads + 1, default is 0)
`--write_auto` | outputs a generating set of the found automorphism group
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

# API
By including "dejavu.h" you can call the automorphism computation directly. The API currently only consists of one function, which can be used as follows. If you have an `sgraph` g, you can compute its automorphism group with a call to `automorphisms`:
```cpp
shared_permnode* gens;
dejavu d;
d.automorphisms(&g, &gens);
```
Configurations can be made using the global struct `config`, in which things such as the thread count can be defined (similar to the commandline arguments).

# Remarks
The source code contains modified source code of the [nauty / Traces](http://pallini.di.uniroma1.it) distribution, as well as the lock-free queue implementation described [here](http://moodycamel.com/blog/2014/a-fast-general-purpose-lock-free-queue-for-c++).