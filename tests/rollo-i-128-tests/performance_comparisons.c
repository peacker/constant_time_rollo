#include <stdio.h>
#include <time.h>
#include <stdint.h>
//#include <rollo-i-128/bf_3149.h>
#include <rollo-i-128/rollo_i_polynomial_arithmetic.h>
#include <rollo-i-128/rollo-i-keygen.h>
#include <rollo-i-128/rollo-i-encapsulation.h>
#include <rollo-i-128/rollo-i-decapsulation.h>
#include <rollo-i-128/vector_space.h>
#include "performance_comparisons.h"

#ifndef TESTRUN_LEN
#define TESTRUN_LEN 2
#endif

void performance_of_unreduced_bf_multiplications() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;
    bf_element_t p1, p2;
    bf_double_element_t p3;

    do {
        bf_get_random_element(&prng, &p1);
        bf_get_random_element(&prng, &p2);

        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            bf_unreduced_multiplication(&p3, &p2, &p1);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("bf_unreduced_multiplication", total_time, total_iterations, total_cycles);
}

#if AVX2

void performance_of_unreduced_bf_multiplications_karatsuba() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;
    bf_element_t p1, p2;
    bf_double_element_t p3;
    do {
        bf_get_random_element(&prng, &p1);
        bf_get_random_element(&prng, &p2);

        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            bf_unreduced_multiplication_karatsuba(&p3, &p2, &p1);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("bf_unreduced_multiplication_karatsuba", total_time, total_iterations, total_cycles);
}

#endif


void performance_of_bf_reductions() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;
    bf_element_t p1[50], p2[50];
    bf_double_element_t p3[50];
    do {
        bf_get_random_element(&prng, p1);
        bf_get_random_element(&prng, p2);
        bf_unreduced_multiplication(p3, p2, p1);

        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            bf_reduction(p1 , p3);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("bf_reduction", total_time, total_iterations, total_cycles);
}

void performance_of_bf_reductions_uint64() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;
    uint64_t output[2];
    uint64_t input[3];
    do {
        randombytes((unsigned char*) input, 24);
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            bf_reduction_uint64(output, input);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("bf_reduction_uint64", total_time, total_iterations, total_cycles);
}

void performance_of_unreduced_squaring() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;
    bf_element_t p1;
    bf_double_element_t p3;
    bf_get_random_element(&prng, &p1);
    do {
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            bf_unreduced_square(&p3, &p1);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("bf_unreduced_square", total_time, total_iterations, total_cycles);
}

void performance_of_bf_inversion() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;
    bf_element_t p1, p2;
    bf_get_random_element(&prng, &p1);
    do {
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            bf_inversion(&p2, &p1);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("bf_inversion", total_time, total_iterations, total_cycles);
}

void performance_of_to_row_echelon_form_10_rows() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t NROWS = 10;
    bf_element_t mat[NROWS];
    size_t dim = 0;

    do {
        //    // fill matrix
        for (i = 0; i < NROWS; i++) {
            bf_get_random_element(&prng, &mat[i]);
        }
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = to_row_echelon_form(mat, NROWS);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("to_row_echelon_form_10_rows", total_time, total_iterations, total_cycles);
}

void performance_of_to_row_echelon_form_20_rows() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t NROWS = 20;
    bf_element_t mat[NROWS];
    size_t dim = 0;

    do {
        //    // fill matrix
        for (i = 0; i < NROWS; i++) {
            bf_get_random_element(&prng, &mat[i]);
        }
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = to_row_echelon_form(mat, NROWS);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("to_row_echelon_form_20_rows", total_time, total_iterations, total_cycles);
}

void performance_of_to_row_echelon_form_30_rows() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t NROWS = 30;
    bf_element_t mat[NROWS];
    size_t dim = 0;

    do {
        //    // fill matrix
        for (i = 0; i < NROWS; i++) {
            bf_get_random_element(&prng, &mat[i]);
        }
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = to_row_echelon_form(mat, NROWS);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("to_row_echelon_form_30_rows", total_time, total_iterations, total_cycles);
}

void performance_of_to_row_echelon_form_100_rows() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t NROWS = 100;
    bf_element_t mat[NROWS];
    size_t dim = 0;

    do {
        //    // fill matrix
        for (i = 0; i < NROWS; i++) {
            bf_get_random_element(&prng, &mat[i]);
        }
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = to_row_echelon_form(mat, NROWS);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("to_row_echelon_form_100_rows", total_time, total_iterations, total_cycles);
}

void performance_of_to_reduced_row_echelon_form_10_rows() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t NROWS = 10;
    bf_element_t mat[NROWS];
    size_t dim = 0;

    do {
        //    // fill matrix
        for (i = 0; i < NROWS; i++) {
            bf_get_random_element(&prng, &mat[i]);
        }
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = to_reduced_row_echelon_form(mat, NROWS);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("to_reduced_row_echelon_form_10_rows", total_time, total_iterations, total_cycles);
}

void performance_of_to_reduced_row_echelon_form_20_rows() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t NROWS = 20;
    bf_element_t mat[NROWS];
    size_t dim = 0;

    do {
        //    // fill matrix
        for (i = 0; i < NROWS; i++) {
            bf_get_random_element(&prng, &mat[i]);
        }
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = to_reduced_row_echelon_form(mat, NROWS);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("to_reduced_row_echelon_form_20_rows", total_time, total_iterations, total_cycles);
}

void performance_of_to_reduced_row_echelon_form_30_rows() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t NROWS = 30;
    bf_element_t mat[NROWS];
    size_t dim = 0;

    do {
        //    // fill matrix
        for (i = 0; i < NROWS; i++) {
            bf_get_random_element(&prng, &mat[i]);
        }
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = to_reduced_row_echelon_form(mat, NROWS);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("to_reduced_row_echelon_form_30_rows", total_time, total_iterations, total_cycles);
}

void performance_of_to_reduced_row_echelon_form_100_rows() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t NROWS = 100;
    bf_element_t mat[NROWS];
    size_t dim = 0;

    do {
        //    // fill matrix
        for (i = 0; i < NROWS; i++) {
            bf_get_random_element(&prng, &mat[i]);
        }
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = to_reduced_row_echelon_form(mat, NROWS);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("to_reduced_row_echelon_form_100_rows", total_time, total_iterations, total_cycles);
}

void performance_of_to_systematic_form_10_rows() {
    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t NROWS = 10;
    bf_element_t mat[NROWS];
    size_t dim = 0;

    do {
        //    // fill matrix
        for (i = 0; i < NROWS; i++) {
            bf_get_random_element(&prng, &mat[i]);
        }
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = to_systematic_form(mat, NROWS);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("to_systematic_form_10_rows", total_time, total_iterations, total_cycles);
}

void performance_of_to_systematic_form_20_rows() {
    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t NROWS = 20;
    bf_element_t mat[NROWS];
    size_t dim = 0;

    do {
        //    // fill matrix
        for (i = 0; i < NROWS; i++) {
            bf_get_random_element(&prng, &mat[i]);
        }
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = to_systematic_form(mat, NROWS);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("to_systematic_form_20_rows", total_time, total_iterations, total_cycles);
}

void performance_of_to_systematic_form_30_rows() {
    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t NROWS = 30;
    bf_element_t mat[NROWS];
    size_t dim = 0;

    do {
        //    // fill matrix
        for (i = 0; i < NROWS; i++) {
            bf_get_random_element(&prng, &mat[i]);
        }
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = to_systematic_form(mat, NROWS);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("to_systematic_form_30_rows", total_time, total_iterations, total_cycles);
}

void performance_of_to_systematic_form_100_rows() {
    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t NROWS = 100;
    bf_element_t mat[NROWS];
    size_t dim = 0;

    do {
        //    // fill matrix
        for (i = 0; i < NROWS; i++) {
            bf_get_random_element(&prng, &mat[i]);
        }
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = to_systematic_form(mat, NROWS);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("to_systematic_form_100_rows", total_time, total_iterations, total_cycles);
}

void performance_generate_list_of_vectors_with_given_rank_RANK_7_LIST_LEN_10() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;


    uint16_t rank = 7;
    bf_element_t list_basis[rank];
    uint16_t list_length = 10;
    bf_element_t list[list_length];


    do {
        // init

        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            generate_list_of_vectors_with_given_rank(list, list_length, rank,
                                                     list_basis, &prng);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("gen_vectors_given_RANK_7_LIST_LEN_10", total_time,
                 total_iterations, total_cycles);
}

void performance_generate_list_of_vectors_with_given_rank_RANK_7_LIST_LEN_100() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;


    uint16_t rank = 7;
    bf_element_t list_basis[rank];
    uint16_t list_length = 100;
    bf_element_t list[list_length];


    do {
        // init

        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            generate_list_of_vectors_with_given_rank(list, list_length, rank,
                                                     list_basis, &prng);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("gen_vectors_given_RANK_7_LIST_LEN_100", total_time,
                 total_iterations, total_cycles);
}

void performance_generate_two_lists_of_vectors_with_given_rank_RANK_7_LIST_LEN_10() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;


    uint16_t rank = 7;
    bf_element_t list_basis[rank];
    uint16_t list_length = 10;
    bf_element_t list1[list_length];
    bf_element_t list2[list_length];


    do {
        // init

        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            generate_two_list_of_vectors_with_given_rank_and_their_basis(list1, list_length,
                                                                         list2, list_length, rank, list_basis, &prng);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("gen_2_vectors_lists_RANK_7_LIST_LEN_10", total_time,
                 total_iterations, total_cycles);
}

void performance_generate_two_lists_of_vectors_with_given_rank_RANK_7_LIST_LEN_100() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;


    uint16_t rank = 7;
    bf_element_t list_basis[rank];
    uint16_t list_length = 100;
    bf_element_t list1[list_length];
    bf_element_t list2[list_length];


    do {
        // init

        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            generate_two_list_of_vectors_with_given_rank_and_their_basis(list1, list_length,
                                                                         list2, list_length, rank, list_basis, &prng);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("gen_2_vectors_lists_RANK_7_LIST_LEN_100", total_time,
                 total_iterations, total_cycles);
}

void performance_zassenhaus_algorithm_sum_reduction_ROWS_83() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t s1_size = ROLLO_I_CODE_LENGTH;
    uint16_t s2_size = ROLLO_I_CODE_LENGTH;
    uint16_t out_size = 2 * ROLLO_I_CODE_LENGTH;
    bf_element_t s1[s1_size];
    bf_element_t s2[s2_size];
    bf_element_t out[out_size];

    size_t dim = 0;

    do {
        // fill matrix
        for (uint16_t i = 0; i < s1_size; i++) {
            bf_get_random_element(&prng, &s1[i]);
        }
        for (uint16_t i = 0; i < s2_size; i++) {
            bf_get_random_element(&prng, &s2[i]);
        }

        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = zassenhaus_algorithm_sum(out, out_size, s1, s1_size, s2, s2_size);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("zassenhaus_sum_ROWS_83", total_time, total_iterations, total_cycles);
}

void performance_zassenhaus_algorithm_intersection_reduction_ROWS_83() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t s1_size = ROLLO_I_CODE_LENGTH;
    uint16_t s2_size = ROLLO_I_CODE_LENGTH;
    uint16_t out_size = 2 * ROLLO_I_CODE_LENGTH;
    bf_element_t s1[s1_size];
    bf_element_t s2[s2_size];
    bf_element_t out[out_size];

    size_t dim = 0;

    do {
        // fill matrix
        for (uint16_t i = 0; i < s1_size; i++) {
            bf_get_random_element(&prng, &s1[i]);
        }
        for (uint16_t i = 0; i < s2_size; i++) {
            bf_get_random_element(&prng, &s2[i]);
        }

        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = zassenhaus_algorithm_sum(out, out_size, s1, s1_size, s2, s2_size);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("zassenhaus_intersection_ROWS_83", total_time, total_iterations, total_cycles);
}

void performance_zassenhaus_algorithm_sum_reduction_ROWS_100() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t s1_size = 100;
    uint16_t s2_size = 100;
    uint16_t out_size = 2 * s1_size;
    bf_element_t s1[s1_size];
    bf_element_t s2[s2_size];
    bf_element_t out[out_size];

    size_t dim = 0;

    do {
        // fill matrix
        for (uint16_t i = 0; i < s1_size; i++) {
            bf_get_random_element(&prng, &s1[i]);
        }
        for (uint16_t i = 0; i < s2_size; i++) {
            bf_get_random_element(&prng, &s2[i]);
        }

        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = zassenhaus_algorithm_sum(out, out_size, s1, s1_size, s2, s2_size);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("zassenhaus_sum_ROWS_100", total_time, total_iterations, total_cycles);
}

void performance_zassenhaus_algorithm_intersection_reduction_ROWS_100() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint16_t s1_size = 100;
    uint16_t s2_size = 100;
    uint16_t out_size = 2 * s1_size;
    bf_element_t s1[s1_size];
    bf_element_t s2[s2_size];
    bf_element_t out[out_size];

    size_t dim = 0;

    do {
        // fill matrix
        for (uint16_t i = 0; i < s1_size; i++) {
            bf_get_random_element(&prng, &s1[i]);
        }
        for (uint16_t i = 0; i < s2_size; i++) {
            bf_get_random_element(&prng, &s2[i]);
        }

        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            dim = zassenhaus_algorithm_sum(out, out_size, s1, s1_size, s2, s2_size);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("zassenhaus_intersection_ROWS_100", total_time, total_iterations, total_cycles);
}

void performance_rank_support_recover() {
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    // F rank(F) = d
    bf_element_t basis_of_subspaceF[ROLLO_I_KEY_VECTORS_RANK_WEIGHT];

    // E rank(E) = r
    uint16_t subspaceE_length = 1u << ROLLO_I_ERROR_VECTORS_RANK_WEIGHT;
    bf_element_t subspaceE[subspaceE_length];
    bf_element_t basis_of_subspaceE[ROLLO_I_ERROR_VECTORS_RANK_WEIGHT];
    bf_element_t recovered_subspaceE[subspaceE_length];

    polynomial_t x;
    polynomial_t y;
    polynomial_t h;
    polynomial_t e1, e2;
    polynomial_t c;
    polynomial_t s;

    do {
        // initialize subspaces
        // generate vector subspace F of GF(2^m) of rank d and its basis
        generate_basis_and_sample_two_polynomials_from_basis(
                &x,
                &y,
                basis_of_subspaceF, ROLLO_I_KEY_VECTORS_RANK_WEIGHT,
                &prng);
        polynomial_inversion_mod_p(&h, &x);
        polynomial_multiplication_mod_p(&h, &h, &y);

        // generate c = e1 + H * e2
        generate_basis_and_sample_two_polynomials_from_basis(
                &e1,
                &e2,
                basis_of_subspaceE, ROLLO_I_ERROR_VECTORS_RANK_WEIGHT,
                &prng);

        // construct E
        generate_support_from_basis(subspaceE, basis_of_subspaceE, ROLLO_I_ERROR_VECTORS_RANK_WEIGHT);
        polynomial_multiplication_mod_p(&c, &h, &e2);
        polynomial_addition(&c, &c, &e1);
        // compute syndrome
        polynomial_multiplication_mod_p(&s, &x, &c);

        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            rank_support_recover(recovered_subspaceE, basis_of_subspaceF, s);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("rank_support_recover", total_time, total_iterations, total_cycles);
}

void performance_polynomial_multiplication() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;
    polynomial_t p1, p2, p3;
    for (i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
        bf_get_random_element(&prng, &p1.coefficients[i]);
        bf_get_random_element(&prng, &p2.coefficients[i]);
    }
    do {
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            polynomial_multiplication_mod_p(&p3, &p2, &p1);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("polynomial_multiplication_mod_p", total_time, total_iterations, total_cycles);
}

void performance_polynomial_inversion() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;
    polynomial_t p1, p2;
    for (i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
        bf_get_random_element(&prng, &p1.coefficients[i]);
    }
    do {
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            polynomial_inversion_mod_p(&p2, &p1);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("polynomial_inversion_mod_p", total_time, total_iterations, total_cycles);
}


void performance_polynomial_inversion_with_precomputed_matrices() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);
    uint64_t iterations = 0x100;
    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    uint64_t total_iterations = iterations;
    int i;
    polynomial_t p1, p2;
    for (i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
        bf_get_random_element(&prng, &p1.coefficients[i]);
    }
    do {
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            polynomial_inversion_mod_p_with_precomputes_matrices(&p2, &p1);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("polynomial_inversion_mod_p [precomputed matrices]", total_time, total_iterations, total_cycles);
}

void performance_keygen() {

    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    int i;
    uint64_t iterations = 0x100;
    uint64_t total_iterations = iterations;
    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    rollo_I_key_pair_t key_pair;
    do {
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            rollo_I_keygen(&key_pair, &prng);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("rollo_I_128_keygen", total_time, total_iterations, total_cycles);
}


void performance_encapsulation() {

    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    int i;
    uint64_t iterations = 0x100;
    uint64_t total_iterations = iterations;
    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    rollo_I_key_pair_t key_pair;
    rollo_I_keygen(&key_pair, &prng);

    uint8_t shared_secret[ROLLO_I_SHARED_SECRET_BYTE_SIZE];
    polynomial_t cipher;
    do {
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            rollo_I_encapsulation(shared_secret, &cipher, &prng, &key_pair.public_key);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("rollo_I_128_encryption", total_time, total_iterations, total_cycles);
}

void performance_decapsulation() {

    uint64_t start, lap, total_cycles;
    clock_t clock_start, clock_lap, total_time;
    total_time = 0;
    total_cycles = 0;
    int i;
    uint64_t iterations = 0x100;
    uint64_t total_iterations = iterations;
    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    rollo_I_key_pair_t key_pair;
    rollo_I_keygen(&key_pair, &prng);

    uint8_t shared_secret[ROLLO_I_SHARED_SECRET_BYTE_SIZE];
    polynomial_t cipher;
    do {
        rollo_I_keygen(&key_pair, &prng);
        rollo_I_encapsulation(shared_secret, &cipher, &prng, &key_pair.public_key);
        start = rdtsc();
        clock_start = clock();
        for (i = 0; i < iterations; i++) {
            rollo_I_decapsulation(shared_secret, &cipher, &key_pair.private_key);
        }
        lap = rdtsc() - start;
        clock_lap = clock() - clock_start;
        total_time += clock_lap;
        total_cycles += lap;
        total_iterations += iterations;

        iterations <<= 1u;
        printf(".");
        fflush(stdout);

    } while (clock_lap < TESTRUN_LEN * CLOCKS_PER_SEC);

    printresults("rollo_I_128_decryption", total_time, total_iterations, total_cycles);
}
