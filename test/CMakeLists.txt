cmake_minimum_required(VERSION 3.16)

project(oasis_test LANGUAGES CXX)

set(TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/hello_world.cpp
)

add_executable(oasis_test ${TEST_SOURCES})

add_test(NAME oasis_test COMMAND oasis_test)

target_link_libraries(oasis_test PUBLIC gtest)