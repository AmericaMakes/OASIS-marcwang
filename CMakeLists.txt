cmake_minimum_required(VERSION 3.19)

project(OASIS_AmericaMake)

include(ExternalProject)
set(CMAKE_CXX_STANDARD 14)

set (VTK_LIB "VTK")
ExternalProject_Add(
    ${VTK_LIB} 
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/vtk
    BUILD_ALWAYS      OFF
    INSTALL_DIR       ${CMAKE_CURRENT_BINARY_DIR}/ext/${VTK_LIB}
    CMAKE_CACHE_ARGS
		-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    BUILD_COMMAND     ${CMAKE_COMMAND} --build <BINARY_DIR> --config Release --target INSTALL
)

add_subdirectory(test)
add_subdirectory(lib/googletest)