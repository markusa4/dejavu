cmake_minimum_required(VERSION 3.5)
project(dejavu)

set(CMAKE_CXX_STANDARD 14)
add_compile_options("-D_GLIBCXX_ASSERTIONS" "-D_FORTIFY_SOURCE=2" "-O3")
# add_definitions(-DNDEBUG)
set(CMAKE_CXX_STANDARD_LIBRARIES -lpthread)

set(source_list "nauty/gtnauty.o;nauty/gtools.o;nauty/gutil1.o;nauty/gutil2.o;nauty/naugraph.o;nauty/naugroup.o;nauty/naurng.o;nauty/nausparse.o;nauty/nautil.o;nauty/nautinv.o;nauty/naututil.o;nauty/nauty.o;nauty/traces.o;nauty/schreier.o" configuration.h invariant_acc.cpp invariant_acc.h pipeline_schreier.h pipeline_schreier.cpp pipeline_group.cpp pipeline_group.h)
add_custom_target(nauty_target
        COMMAND make all
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/nauty
        )

add_executable(dejavu main.cpp sgraph.cpp sgraph.h invariant.cpp invariant.h refinement.cpp refinement.h selector.cpp selector.h coloring.cpp coloring.h parser.cpp parser.h bijection.cpp bijection.h ir_tools.cpp ir_tools.h sequential_group.cpp sequential_group.h auto_blaster.cpp auto_blaster.h trail.cpp trail.h ${source_list})
add_dependencies(dejavu nauty_target)