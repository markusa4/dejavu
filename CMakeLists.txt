cmake_minimum_required(VERSION 3.14)
project(dejavu)

set(CMAKE_CXX_STANDARD 14)
add_compile_options("-D_GLIBCXX_ASSERTIONS" "-D_FORTIFY_SOURCE=2" "-O3")
set(CMAKE_CXX_STANDARD_LIBRARIES -lpthread)

set(source_list "nauty/gtnauty.o;nauty/gtools.o;nauty/gutil1.o;nauty/gutil2.o;nauty/naugraph.o;nauty/naugroup.o;nauty/naurng.o;nauty/nausparse.o;nauty/nautil.o;nauty/nautinv.o;nauty/naututil.o;nauty/nauty.o;nauty/traces.o;nauty/schreier.o" badtree.cpp badtree.h configuration.h invariant_acc.cpp invariant_acc.h)
add_custom_target(nauty_target
        COMMAND make all
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/nauty
        )

add_executable(dejavu main.cpp sgraph.cpp sgraph.h invariant.cpp invariant.h refinement.cpp refinement.h selector.cpp selector.h coloring.cpp coloring.h parser.cpp parser.h bijection.cpp bijection.h ir_tools.cpp ir_tools.h group.cpp group.h auto_blaster.cpp auto_blaster.h trail.cpp trail.h ${source_list})
add_dependencies(dejavu nauty_target)