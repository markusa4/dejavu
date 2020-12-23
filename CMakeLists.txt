cmake_minimum_required(VERSION 3.5)
project(dejavu)

set(CMAKE_CXX_STANDARD 14)
add_compile_options("-O3")
add_definitions(-DNDEBUG)
set(CMAKE_CXX_STANDARD_LIBRARIES -lpthread)

add_executable(dejavu-auto schreier_sequential.cpp utility.cpp sgraph.h naurng.h naurng.cpp invariant.cpp invariant.h refinement.cpp refinement.h selector.h coloring.h parser.cpp parser.h bijection.h dejavu_auto.cpp dejavu_auto.h configuration.h schreier_shared.h schreier_shared.cpp bfs.h)
add_executable(dejavu-iso  utility.cpp dejavu_iso.h sgraph.h invariant.cpp invariant.h refinement.cpp refinement.h selector.h coloring.h parser.cpp parser.h bijection.h configuration.h bfs.h dejavu_iso.cpp)
# add_executable(dejavu-api-test utility.cpp dejavu_iso.h dejavu_api.h sgraph.h invariant.cpp invariant.h refinement.cpp refinement.h selector.h coloring.h parser.cpp parser.h bijection.h configuration.h bfs.h dejavu_api.cpp)
add_library(dejavu-api utility.cpp dejavu_iso.h dejavu_api.h sgraph.h invariant.cpp invariant.h refinement.cpp refinement.h selector.h coloring.h parser.cpp parser.h bijection.h configuration.h bfs.h)

if (DEFINED COMPILE_PYTHON2API)
    # might need to replace 2.7 with appropriate version
    find_package(Python2 2.7 COMPONENTS Interpreter Development REQUIRED)
    message(STATUS "PYTHON_LIBRARIES = ${Python2_LIBRARIES}")
    message(STATUS "PYTHON_EXECUTABLE = ${Python2_EXECUTABLE}")
    message(STATUS "PYTHON_INCLUDE_DIRS = ${Python2_INCLUDE_DIRS}")

    # might need to replace "python" with "python2" or "python27"
    find_package(Boost COMPONENTS python2 REQUIRED)
    set(CMAKE_SHARED_MODULE_PREFIX "")
    add_library(dejavupy2 MODULE dejavu_api_python.cpp utility.cpp dejavu_iso.h dejavu_api.h sgraph.h invariant.cpp invariant.h refinement.cpp refinement.h selector.h coloring.h parser.cpp parser.h bijection.h configuration.h bfs.h)
    set_target_properties(dejavupy2 PROPERTIES SUFFIX ".pyd")
    target_link_libraries(dejavupy2 ${Boost_LIBRARIES} ${Python2_LIBRARIES})
    target_include_directories(dejavupy2 PRIVATE ${Python2_INCLUDE_DIRS})
endif ()

if (DEFINED COMPILE_PYTHON3API)
    # might need to replace 3.6 with appropriate version
    find_package(Python3 3.6 COMPONENTS Interpreter Development REQUIRED)
    message(STATUS "PYTHON_LIBRARIES = ${Python3_LIBRARIES}")
    message(STATUS "PYTHON_EXECUTABLE = ${Python3_EXECUTABLE}")
    message(STATUS "PYTHON_INCLUDE_DIRS = ${Python3_INCLUDE_DIRS}")

    # might need to replace "python" with "python3" or "python36"
    find_package(Boost COMPONENTS python3 REQUIRED)
    set(CMAKE_SHARED_MODULE_PREFIX "")
    add_library(dejavupy3 MODULE dejavu_api_python.cpp utility.cpp dejavu_iso.h dejavu_api.h sgraph.h invariant.cpp invariant.h refinement.cpp refinement.h selector.h coloring.h parser.cpp parser.h bijection.h configuration.h bfs.h)
    set_target_properties(dejavupy3 PROPERTIES SUFFIX ".pyd")
    target_link_libraries(dejavupy3 ${Boost_LIBRARIES} ${Python3_LIBRARIES})
    target_include_directories(dejavupy3 PRIVATE ${Python3_INCLUDE_DIRS})
endif ()

####################################################################
# If cmake can not properly detect boost python, we suggest        #
# trying to link the library manually (see below).                 #
####################################################################
if (DEFINED COMPILE_PYTHON3API_LINK_BOOST_MANUAL)
    find_package(Python3 3.6 COMPONENTS Interpreter Development REQUIRED)
    message(STATUS "PYTHON_LIBRARIES = ${Python3_LIBRARIES}")
    message(STATUS "PYTHON_EXECUTABLE = ${Python3_EXECUTABLE}")
    message(STATUS "PYTHON_INCLUDE_DIRS = ${Python3_INCLUDE_DIRS}")

    ####################################################################
    # Insert approriate paths here                                     #
    ####################################################################
    set( BOOST_ROOT "C:/boost" )
    set (BOOST_INCLUDEDIR "C:/boost/include/boost-1_74")
    set (BOOST_LIBRARYDIR "C:/boost/lib")
    ####################################################################


    set(CMAKE_SHARED_MODULE_PREFIX "")
    add_library(dejavupy3 MODULE dejavu_api_python.cpp utility.cpp dejavu_iso.h dejavu_api.h sgraph.h invariant.cpp invariant.h refinement.cpp refinement.h selector.h coloring.h parser.cpp parser.h bijection.h configuration.h bfs.h)
    set_target_properties(dejavupy3 PROPERTIES SUFFIX ".pyd")
    # uncomment for static link
    #set_target_properties(libdejavu_python PROPERTIES LINK_FLAGS "-static")
    #add_definitions(/DBOOST_PYTHON_STATIC_LIB)

    target_link_directories(dejavupy3 PRIVATE ${BOOST_LIBRARYDIR})
    ####################################################################
    # Insert approriate library names here                             #
    ####################################################################
    target_link_libraries(dejavupy3 libboost_python36-mgw9-mt-x64-1_74.dll.a libboost_python36-mgw9-mt-x64-1_74.dll  ${Python3_LIBRARIES})
    ####################################################################
    target_include_directories(dejavupy3 PRIVATE ${Python3_INCLUDE_DIRS} ${BOOST_INCLUDEDIR})
endif ()

set(source_list "nauty/gtnauty.o;nauty/gtools.o;nauty/gutil1.o;nauty/gutil2.o;nauty/naugraph.o;nauty/naugroup.o;nauty/naurng.o;nauty/nausparse.o;nauty/nautil.o;nauty/nautinv.o;nauty/naututil.o;nauty/nauty.o;nauty/traces.o;nauty/schreier.o")
add_custom_target(nauty_target
        COMMAND make all
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/nauty
        )
set(conauto_source_list "conauto/conauto.o;conauto/factor.o;conauto/graph.o;conauto/groups.o;conauto/mismatches.o;conauto/orbits.o;conauto/refinements.o;conauto/seqpart.o;conauto/sort.o")

if (DEFINED COMPILE_BENCH_ISO)
        add_executable(bench-iso ${conauto_source_list} schreier_sequential.cpp utility.cpp bench_iso.cpp sgraph.h naurng.h naurng.cpp invariant.cpp invariant.h refinement.cpp refinement.h selector.h coloring.h parser.cpp parser.h bijection.h dejavu_iso.h configuration.h schreier_shared.h schreier_shared.cpp bfs.h ${source_list} naudefs.h)
        add_dependencies(bench-iso nauty_target)
endif ()

if (DEFINED COMPILE_BENCH_AUTO)
	add_executable(bench schreier_sequential.cpp utility.cpp bench.cpp sgraph.h naurng.h naurng.cpp invariant.cpp invariant.h refinement.cpp refinement.h selector.h coloring.h parser.cpp parser.h bijection.h dejavu_auto.cpp dejavu_auto.h configuration.h schreier_shared.h schreier_shared.cpp bfs.h ${source_list} naudefs.h conauto/conauto.h)
	add_dependencies(bench nauty_target)
endif ()