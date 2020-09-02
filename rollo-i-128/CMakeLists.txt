cmake_minimum_required(VERSION 3.7)
project(${PROJECT_NAME})

set(HDR
        rollo_i_polynomial_arithmetic.h
        rollo_i_polynomial_arithmetic_precomputed_matrices.h
        rollo-i-128_parameters.h
        rollo-i-keygen.h
        rollo-i-encapsulation.h
        rollo-i-decapsulation.h
        vector_space.h
        bf.h
        )

set(SRC
        rollo_i_polynomial_arithmetic.c
        rollo-i-keygen.c
        rollo-i-encapsulation.c
        rollo-i-decapsulation.c
        vector_space.c
        bf.c
        bf_uint128_t.c
        main.c
        )

include_directories(
        ${CMAKE_SOURCE_DIR})

set(EXECUTABLE_NAME rollo_i_128_simple_test)
add_executable(${EXECUTABLE_NAME} ${SRC})
#
add_dependencies(${EXECUTABLE_NAME} ${PROJECT_NAME})
target_link_libraries(${EXECUTABLE_NAME} ${PROJECT_NAME})

if (UNIX AND NOT DEFINED ARM)
    target_link_libraries(${EXECUTABLE_NAME} pthread)
endif ()

install(TARGETS ${EXECUTABLE_NAME} DESTINATION bin)


set(LIBRARY_NAME rollo-i-128)
add_library(${LIBRARY_NAME} OBJECT ${SRC})
install(FILES ${HDR} DESTINATION include/${LIBRARY_NAME})