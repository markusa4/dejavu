# Compilation
If only the dejavu automorphism and isomorphism solvers are to be compiled the project should compile without any further dependencies:
```
cmake .
make 
```

The benchmark utility has binary dependencies to nauty. The easiest way to compile it is to place a configured and compiled distribution of nauty / Traces (available at http://pallini.di.uniroma1.it/) in "nauty/". After placing nauty and Traces into the respective folder, the project can be compiled with cmake:
```
cmake -DBENCH=1 --build "."
make 
```

# Usage
Compilation produces (up to) three binaries (`dejavu`, `dejavu-iso` and `bench`), which are however similar in usage. dejavu is the binary of the automorphism solver, dejavu-iso of the isomorphism solver. bench is a frontend which can call nauty, Traces and dejavu, while recording and tracking time measurements, as well as manage timeouts. 

## dejavu
The automorphism solver only accepts files in the DIMACS graph format. Only undirected graphs can be handled at this point -- but the graphs may be colored. The solver accepts the following command line arguments:

Command Line Argument | Effect
--- | ---
`--file` | specifies the graph input file in DIMACS format
`--threads` | specifies the amount of additional threads used (total is threads + 1, default is 0)
`--write-auto` | outputs a generating set of the found automorphism group
`--timeout` | specifies a timeout (in seconds)
`--permute` | permutes the input graph randomly
`--permute-seed` | specify a fixed seed for the previous argument
`--force-selector` | forces a specific cell selector (0 is first, 1 is smallest, 2 is largest)
`--no-idleskip` | suppresses skipping of idle cells, i.e., cells that cause splits
`--kdeviation` | determine the number of additional refinement steps for deviation values

## dejavu-iso
The isomorphism solver accepts the same arguments, however, it requires two input files.

Command Line Argument | Effect
--- | ---
`--file1` | specifies the first graph input file in DIMACS format
`--file2` | specifies the second graph input file in DIMACS format

## bench
`bench` understands the same command line arguments, with the addition of the following: 

Command Line Argument | Effect
--- | ---
`--no-nauty` | does not run nauty
`--no-traces` | does not run Traces
`--no-dejavu` | does not run dejavu
`--stat-file` | specify a file to which measurements shall be written

# API
By including "dejavu.h" you can call the automorphism computation directly. The API currently only consists of one function, which can be used as follows. If you have an `sgraph` g, you can compute its automorphism group with a call to `dejavu_automorphisms`:
```cpp
shared_permnode* gens;
dejavu_automorphisms(&g, nullptr, &gens);
```
Configurations can be made using the global struct `config`, in which things such as the thread count can be defined (analogous to the commandline arguments). The second argument defines the initial colors of vertices. If the graph is uncolored, a null pointer can be passed. Otherwise, an integer array mapping vertices to colors is expected.

# License
The source code contains modified source code of the [nauty / Traces](http://pallini.di.uniroma1.it) distribution, as well as the lock-free queue implementation described [here](http://moodycamel.com/blog/2014/a-fast-general-purpose-lock-free-queue-for-c++). 

1. All files not listed below, copyright Markus Anders.
2. `naudefs.h`, `naurng.cpp`, `naurng.h`, `schreier_sequential.cpp`, `schreier_sequential.h`, `schreier_shared.cpp`, `schreier_shared.h` are (derivative) work, originally copyright Brendan McKay. Licensed under the Apache License, Version 2.0 (see `LICENSE`).
3. `concurrentqueue.h`, copyright Cameron Desrochers. Licensed under the Simplified BSD license (see `concurrentqueue.h`)
