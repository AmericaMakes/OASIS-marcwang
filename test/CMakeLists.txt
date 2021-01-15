cmake_minimum_required(VERSION 3.16)

project(oasis_test LANGUAGES CXX)

enable_testing()

include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        release-1.10.0
)

FetchContent_GetProperties(googletest)

if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)

    # Prevent GoogleTest from overriding our compiler/linker options
    # when building with Visual Studio
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    # Prevent GoogleTest from using PThreads
    set(gtest_disable_pthreads ON CACHE BOOL "" FORCE)

    # adds the targers: gtest, gtest_main, gmock, gmock_main
    add_subdirectory(
        ${googletest_SOURCE_DIR}
        ${googletest_BINARY_DIR}
    )

endif()


set(TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/hello_world.cpp
)

add_executable(oasis_test ${TEST_SOURCES})

add_test(NAME oasis_test COMMAND generic_test)

target_link_libraries(oasis_test PUBLIC gtest_main)