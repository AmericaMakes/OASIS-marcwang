SET(OASIS_SRC 
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkVoxelizer.cpp
)

SET(OASIS_HEADER
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkVoxelizer.h
)

add_library(oasis_lib ${OASIS_HEADER} ${OASIS_SRC})