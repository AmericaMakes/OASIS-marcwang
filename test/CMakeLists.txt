set(TEST_EXEC ${CMAKE_PROJECT_NAME}_test)

set(TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/hello_world.cpp
)

add_executable(${TEST_EXEC} ${TEST_SOURCES})

add_test(NAME ${TEST_EXEC} COMMAND ${TEST_EXEC})

target_link_libraries(${TEST_EXEC} PUBLIC gtest)