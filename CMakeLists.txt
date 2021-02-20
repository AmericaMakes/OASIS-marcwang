cmake_minimum_required(VERSION 3.13.4)
project(oasis_lib)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/geogram/cmake)
find_package(Catch2 REQUIRED)
find_package(Geogram REQUIRED)
find_package(pybind11 REQUIRED)
find_package(Boost REQUIRED)

set(LIB_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/OasisLib/src/GeogramVoronoi.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OasisLib/src/MeshGraph.cpp
)

set(LIB_HEADER
    ${CMAKE_CURRENT_SOURCE_DIR}/OasisLib/include/GeogramVoronoi.h
    ${CMAKE_CURRENT_SOURCE_DIR}/OasisLib/include/GeogramBase.h
    ${CMAKE_CURRENT_SOURCE_DIR}/OasisLib/include/MeshGraph.h
)

set(LIB_SRC_WRAP
    ${CMAKE_CURRENT_SOURCE_DIR}/OasisLib/src/wrap_geogram.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OasisLib/include/GeogramCaster.h
)

pybind11_add_module(OasisLib ${LIB_SRC} ${LIB_SRC_WRAP} ${LIB_HEADER})

target_link_libraries(OasisLib PUBLIC Geogram::geogram)
target_include_directories(OasisLib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/OasisLib/include) 

######### TEST #############

SET(TEST_ARTIFACT 
    ${CMAKE_CURRENT_SOURCE_DIR}/test/test_artifact/3DBenchy.stl
    ${CMAKE_CURRENT_SOURCE_DIR}/test/test_artifact/Parameter_quality_nut_1.stl
    ${CMAKE_CURRENT_SOURCE_DIR}/test/test_artifact/Spiral_lattice_cone_1.stl
)

file(COPY 
    ${CMAKE_CURRENT_SOURCE_DIR}/test/test_artifact/
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test_artifact/
)

SET(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)
SET(TEST_SRC ${TEST_DIR}/test_voronoi.cpp
             ${TEST_DIR}/test_main.cpp)

add_executable(OasisLibtest ${TEST_SRC} ${LIB_SRC} ${LIB_HEADER})
target_include_directories(OasisLibtest PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/OasisLib/include
    ${CMAKE_CURRENT_SOURCE_DIR}/test/
)
target_link_libraries(OasisLibtest Catch2::Catch2WithMain Geogram::geogram)

include(CTest)
include(Catch)

catch_discover_tests(OasisLibtest)
