cmake_minimum_required(VERSION 3.14)
project(brutus)

set(CMAKE_CXX_STANDARD 14)
add_compile_options("-D_GLIBCXX_ASSERTIONS" "-D_FORTIFY_SOURCE=2" "-O0")
add_executable(brutus main.cpp graph.cpp graph.h invariant.cpp invariant.h refinement.cpp refinement.h selector.cpp selector.h coloring.cpp coloring.h parser.cpp parser.h bijection.cpp bijection.h)