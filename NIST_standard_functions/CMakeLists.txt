cmake_minimum_required(VERSION 3.7)
project(tii-pqclib)

set(HDR
        rng/nist-rng.h
        hash/hash.h
        )

set(SRC
        rng/nist-rng.c
        hash/hash.c)

include_directories(
        ${CMAKE_SOURCE_DIR}
        ${CRYPTOLIB_INCLUDE_DIR}
        /usr/local/include)

add_library(nist OBJECT ${SRC})

install(FILES ${HDR} DESTINATION include/nist)
