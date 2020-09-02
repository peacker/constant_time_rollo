cmake_minimum_required(VERSION 3.7)
project(tii-pqclib)

set(HDR
        avx2_conversion.h
        )

set(SRC
        avx2_conversion.c)

include_directories(
        ${CMAKE_SOURCE_DIR})

add_library(avx2_conversion OBJECT ${SRC})

install(FILES ${HDR} DESTINATION include/avx2_conversion)
