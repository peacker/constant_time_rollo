cmake_minimum_required(VERSION 3.7)
project(${PROJECT_NAME})

set(HDR
        hexhelper.h
        rollo-i-128-tests/performance_comparisons.h
        rollo-i-128-tests/polynomial_arithmetic_tests.h
        rollo-i-128-tests/keygen_tests.h
        rollo-i-128-tests/encapsulation_tests.h
        rollo-i-128-tests/decapsulation_tests.h
        rollo-i-128-tests/vector_space_tests.h
        rollo-i-128-tests/bf_tests.h
        )

set(SRC
        hexhelper.c
        test_all.c
        rollo-i-128-tests/all_rollo-i-128_tests.c
        rollo-i-128-tests/performance_comparisons.c
        rollo-i-128-tests/polynomial_arithmetic_tests.c
        rollo-i-128-tests/keygen_tests.c
        rollo-i-128-tests/encapsulation_tests.c
        rollo-i-128-tests/decapsulation_tests.c
        rollo-i-128-tests/vector_space_tests.c
        rollo-i-128-tests/bf_tests.c
        )


include_directories(
        ${CMAKE_SOURCE_DIR}
)

set(EXECUTABLE_NAME rollo_i_128_tests)
add_executable(${EXECUTABLE_NAME} ${SRC})
add_dependencies(${EXECUTABLE_NAME} ${PROJECT_NAME})
target_link_libraries(${EXECUTABLE_NAME} ${PROJECT_NAME})

if (UNIX AND NOT DEFINED ARM)
    target_link_libraries(${EXECUTABLE_NAME} pthread)
endif ()
