cmake_minimum_required(VERSION 3.5)
project(dejavu)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options("-O3")
add_compile_options("-Wall" "-Wpointer-arith" "-Wcast-align" "-Wwrite-strings" "-Wshadow" "-Wredundant-decls"
                    "-Wdisabled-optimization" "-Wnon-virtual-dtor" "-Wreorder" "-Woverloaded-virtual" "-Wsign-promo"
                    "-Wsynth" "-Wcast-qual" "-Wno-long-long" "-Wno-unknown-pragmas" "-Wno-unused-parameter"
                    "-Wno-strict-overflow")
add_compile_options("-march=native")
add_definitions(-DNDEBUG)
set(COMPILE_TEST_SUITE FALSE CACHE BOOL "Whether to compile the test suite")
#set(COMPILE_TEST_SUITE TRUE)

add_executable(dejavu dejavu.cpp)

if (${COMPILE_TEST_SUITE})
    message("Tests active...")

    add_compile_definitions("-D_GLIBCXX_ASSERTIONS")
    remove_definitions(-DNDEBUG)
    add_definitions(-DDEJDEBUG)
    add_definitions(-g)

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
            tests/worklist_test.cpp
            tests/workspace_test.cpp
            tests/automorphism_test.cpp
            tests/simple_graphs_test.cpp
            tests/main_test.cpp
            tests/refinement_test.cpp
            tests/orbit_test.cpp
            tests/static_graph_test.cpp
            tests/schreier_test.cpp
            tests/graphs_test.cpp
    )
    target_link_libraries(
            dejavu_test
            GTest::gtest
    )

    target_compile_definitions(dejavu_test PUBLIC TEST_RESOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/graphs/")

    include(GoogleTest)
    gtest_discover_tests(dejavu_test)
endif()