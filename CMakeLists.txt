cmake_minimum_required(VERSION 3.5)
project(dejavu)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options("-O3")
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options("-Wall" "-Wpointer-arith" "-Wcast-align" "-Wwrite-strings" "-Wshadow" "-Wredundant-decls"
        "-Wdisabled-optimization" "-Wnon-virtual-dtor" "-Wreorder" "-Woverloaded-virtual" "-Wsign-promo"
        "-Wsynth" "-Wcast-qual" "-Wno-long-long" "-Wno-unknown-pragmas" "-Wno-unused-parameter"
        "-Wno-strict-overflow")
    add_compile_options("-march=native")
endif()
add_definitions(-DNDEBUG)

add_executable(dejavu dejavu.cpp)

if (${TEST})
    message("Tests target active...")
    message("WARNING: dejavu will be compiled with assertions ON.")

    add_compile_definitions("-D_GLIBCXX_ASSERTIONS")
    remove_definitions(-DNDEBUG)
    add_definitions(-DDEJDEBUG)
    add_definitions(-g)

    # Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24:
	if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
		cmake_policy(SET CMP0135 NEW)
    endif()
    include(FetchContent)

    # Download googletest
    FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    # Download graphs and groups used for testing
    set(TEST_GRAPH_LOCATION ${CMAKE_SOURCE_DIR}/tests/graphs)
    FetchContent_Declare(
            testgraphs
            URL https://automorphisms.org/graphs/graphs.zip
            SOURCE_DIR ${TEST_GRAPH_LOCATION})

    FetchContent_MakeAvailable(testgraphs)

    enable_testing()
    add_executable(
            dejavu_test
            tests/markset_test.cpp
            tests/workset_test.cpp
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
