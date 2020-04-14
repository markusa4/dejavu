cmake_minimum_required(VERSION 3.5)
project(dejavu)

set(CMAKE_CXX_STANDARD 14)
add_compile_options("-O3")
add_definitions(-DNDEBUG)
set(CMAKE_CXX_STANDARD_LIBRARIES -lpthread)

set(source_list "nauty/gtnauty.o;nauty/gtools.o;nauty/gutil1.o;nauty/gutil2.o;nauty/naugraph.o;nauty/naugroup.o;nauty/naurng.o;nauty/nausparse.o;nauty/nautil.o;nauty/nautinv.o;nauty/naututil.o;nauty/nauty.o;nauty/traces.o;nauty/schreier.o")
add_custom_target(nauty_target
        COMMAND make all
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/nauty
        )

add_executable(dejavu schreier_sequential.cpp utility.cpp main.cpp sgraph.h naurng.h naurng.cpp invariant.cpp invariant.h refinement.cpp refinement.h selector.h coloring.h parser.cpp parser.h bijection.h dejavu.cpp dejavu.h configuration.h schreier_shared.h schreier_shared.cpp bfs.h)
add_executable(vujade schreier_sequential.cpp utility.cpp vujade.h sgraph.h naurng.h naurng.cpp invariant.cpp invariant.h refinement.cpp refinement.h selector.h coloring.h parser.cpp parser.h bijection.h configuration.h schreier_shared.h schreier_shared.cpp bfs.h vujade.cpp)

set(conauto_source_list "conauto/conauto.o;conauto/factor.o;conauto/graph.o;conauto/groups.o;conauto/mismatches.o;conauto/orbits.o;conauto/refinements.o;conauto/seqpart.o;conauto/sort.o")

add_executable(bench_iso ${conauto_source_list} schreier_sequential.cpp utility.cpp bench_iso.cpp sgraph.h naurng.h naurng.cpp invariant.cpp invariant.h refinement.cpp refinement.h selector.h coloring.h parser.cpp parser.h bijection.h vujade.h configuration.h schreier_shared.h schreier_shared.cpp bfs.h ${source_list} naudefs.h)
add_dependencies(bench_iso nauty_target)

if (DEFINED BENCH)
	add_executable(bench schreier_sequential.cpp utility.cpp bench.cpp sgraph.h naurng.h naurng.cpp invariant.cpp invariant.h refinement.cpp refinement.h selector.h coloring.h parser.cpp parser.h bijection.h dejavu.cpp dejavu.h configuration.h schreier_shared.h schreier_shared.cpp bfs.h ${source_list} naudefs.h conauto/conauto.h)
	add_dependencies(bench nauty_target)
endif ()
