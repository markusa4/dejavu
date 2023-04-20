cmake_minimum_required(VERSION 3.5)
project(dejavu)

set(CMAKE_CXX_STANDARD 20)
add_compile_options("-O3")
add_compile_options("-W")
add_compile_options("-Wall")
add_definitions(-DNDEBUG)
#add_definitions(-COMPILE_BENCH_AUTO)
add_definitions(-DDEJAVU_STANDALONE)
add_definitions(-g)
set(CMAKE_CXX_STANDARD_LIBRARIES -lpthread)

add_executable(dejavu-auto dejavu_auto.cpp sgraph.h refinement.h coloring.h parser.h)