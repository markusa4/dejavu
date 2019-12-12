cmake_minimum_required(VERSION 3.5)
project(dejavu)

set(CMAKE_CXX_STANDARD 14)
add_compile_options("-O3")
add_definitions(-DNDEBUG)
set(CMAKE_CXX_STANDARD_LIBRARIES -lpthread)

set(source_list "nauty/gtnauty.o;nauty/gtools.o;nauty/gutil1.o;nauty/gutil2.o;nauty/naugraph.o;nauty/naugroup.o;nauty/naurng.o;nauty/nausparse.o;nauty/nautil.o;nauty/nautinv.o;nauty/naututil.o;nauty/nauty.o;nauty/traces.o;nauty/schreier.o" configuration.h schreier_shared.h schreier_shared.cpp bfs.cpp bfs.h lowdeg.cpp lowdeg.h)
add_custom_target(nauty_target
        COMMAND make all
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/nauty
        )

add_executable(dejavu schreier_sequential.cpp group_shared.cpp utility.cpp main.cpp sgraph.cpp sgraph.h invariant.cpp invariant.h refinement.cpp refinement.h selector.cpp selector.h coloring.cpp coloring.h parser.cpp parser.h bijection.cpp bijection.h dejavu.cpp dejavu.h ${source_list})
add_executable(bench schreier_sequential.cpp group_shared.cpp utility.cpp bench.cpp sgraph.cpp sgraph.h invariant.cpp invariant.h refinement.cpp refinement.h selector.cpp selector.h coloring.cpp coloring.h parser.cpp parser.h bijection.cpp bijection.h dejavu.cpp dejavu.h ${source_list})
add_dependencies(dejavu nauty_target)
add_dependencies(bench nauty_target)
