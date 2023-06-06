cmake_minimum_required(VERSION 3.5)
project(dejavu)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options("-O0")
add_compile_options("-W")
add_compile_options("-Wall")
#add_compile_options("-funroll-loops")
#add_compile_options("-march=native")
add_definitions(-DNDEBUG)
#add_definitions(-COMPILE_BENCH_AUTO)
add_definitions(-DDEJAVU_STANDALONE)
add_definitions(-g)
set(CMAKE_CXX_STANDARD_LIBRARIES -lpthread)

add_executable(dejavu-auto dejavu_auto.cpp sgraph.h refinement.h coloring.h parser.h)

set(CMAKE_C_COMPILER /usr/bin/cc)
set(CMAKE_CXX_COMPILER /usr/bin/c++)
set(COMPILE_FLAGS "--coverage")
set(CMAKE_EXE_LINKER_FLAGS "--coverage")

remove_definitions(-DNDEBUG)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
add_executable(
        dejavu_test
        tests/markset_test.cpp
        tests/simple_graphs_test.cpp
        tests/main_test.cpp
        tests/refinement_test.cpp
        tests/orbit_test.cpp
        tests/graphs_test.cpp
)
target_link_libraries(
        dejavu_test
        GTest::gtest
)

target_compile_definitions(dejavu_test PUBLIC TEST_RESOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/graphs/")

include(GoogleTest)
gtest_discover_tests(dejavu_test)