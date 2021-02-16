cmake_minimum_required(VERSION 3.13.4)
project(oasis_lib)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/geogram/cmake)
find_package(Catch2 REQUIRED)
find_package(Geogram REQUIRED)
find_package(pybind11 REQUIRED)

set(LIB_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/OasisLib/src/GeoGramVoronoi.cpp
)

set(LIB_SRC_WRAP
    ${CMAKE_CURRENT_SOURCE_DIR}/OasisLib/src/wrap_geogram.cpp
)

pybind11_add_module(OasisLib ${LIB_SRC} ${LIB_SRC_WRAP})

target_link_libraries(OasisLib PUBLIC Geogram::geogram)
target_include_directories(OasisLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/OasisLib/include) 

######### TEST #############
SET(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)
SET(TEST_SRC ${TEST_DIR}/test_voronoi.cpp
             ${TEST_DIR}/test_main.cpp)

add_executable(OasisLibtest ${TEST_SRC} ${LIB_SRC})
target_include_directories(OasisLibtest PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/include
)
target_link_libraries(OasisLibtest Catch2::Catch2WithMain Geogram::geogram)

enable_testing()

add_test(NAME OasisLibtest COMMAND RunUnitTest)
