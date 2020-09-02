#include <stdio.h>
#include "vector_space_tests.h"
#include "rollo-i-128/vector_space.h"
#include "rollo-i-128/rollo-i-128_parameters.h"
#include "all_rollo-i-128_tests.h"

int test_to_reduced_row_echelon_form() {
    size_t NROWS = 6;
    bf_element_t mat[NROWS];

    uint8_t seed[32] = {10,};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    size_t dim = 0;
    size_t expected_dim = NROWS;

    uint8_t verb = 0;

//    // fill matrix
    for (uint8_t i = 0; i < NROWS; i++) {
//        bf113_weak_rand(&mat[i]);
        bf_get_random_element(&prng, &mat[i]);
    }

    // fill matrix
//    for (uint8_t i = 0; i < NROWS; i++) {
//        bf_set_to_ones(&mat[i]);
//    }

    if (verb) {
        printf("\nbinary: A = \n");
        bf_print_matrix_bin(mat, NROWS);
    }

    dim = to_reduced_row_echelon_form(mat, NROWS);

    if (verb) {
        printf("binary: A = \n");
        bf_print_matrix_bin(mat, NROWS);
    }

    if (dim != expected_dim) {
        printf("computed dim = %ld\n", dim);
        printf("expected dim = %ld\n", expected_dim);
        printf("\n\n");
        return 1;
    }

    return 0;
}

int test_to_systematic_form() {
    size_t NROWS = 6;
    bf_element_t mat[NROWS];

    uint8_t seed[32] = {10,};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    size_t is_systematic = 0;
    size_t expected_dim = NROWS;

    uint8_t verb = 0;

//    // fill matrix
    for (uint16_t i = 0; i < NROWS; i++) {
        bf_get_random_element(&prng, &mat[i]);
    }

    if (verb) {
        printf("\nbinary: A = \n");
        bf_print_matrix_bin(mat, NROWS);
    }

    is_systematic = to_systematic_form(mat, NROWS);

    if (verb) {
        printf("binary: A = \n");
        bf_print_matrix_bin(mat, NROWS);
    }

    // TODO: implement a test that check if the matrix is systematic (by checking if it contains the identity matrix)
    if (is_systematic != 1) {
        printf("The matrix is not systematic\n");
        printf("\n\n");
        return 1;
    }

    return 0;
}

int test_uint16_t_is_greater_than() {
    uint16_t a = 0;
    uint16_t b = 0;

    for (uint16_t i = 0; i < 100; ++i) {
        for (uint16_t j = 0; j < 100; ++j) {
            randombytes((uint8_t *) &a, 2);
            randombytes((uint8_t *) &b, 2);
            if (!((a > b & uint16_t_is_greater_than(a, b) == 1) || (a <= b & uint16_t_is_greater_than(a, b) == 0))) {
                printf("\na = %d\n", a);
                printf("b = %d\n", b);
                printf("output = %d\n", uint16_t_is_greater_than(a, b));

                return EXIT_FAILURE;
            }
        }
    }


    return EXIT_SUCCESS;
}

int test_to_row_echelon_form_for_a_matrix_with_identical_entries() {
    uint16_t NROWS = ROLLO_I_CODE_LENGTH;
    bf_element_t mat[NROWS];

    size_t dim = 0;

    //    // fill matrix
    //    for (i = 0; i < NROWS; i++) {
    //        bf113_weak_rand(&mat[i]);
    //    }

    // fill matrix
    for (uint16_t i = 0; i < NROWS; ++i)
        bf_set_to_ones(&mat[i]);
    //        bf_set_from_hex_string(&mat[i], "0xF");

//        printf("binary: A = \n");
//        bf_print_matrix_bin(mat, NROWS);

    dim = to_row_echelon_form(mat, NROWS);

//        printf("binary: A = \n");
//        bf_print_matrix_bin(mat, NROWS);

    if (dim != 1) {
        printf("computed dim = %ld\n", dim);
        printf("expected dim = 1\n");
        printf("\n\n");
        return 1;
    }

    return 0;
}

//int test_to_row_echelon_form_with_five_independent_linear_rows() {
int test_to_row_echelon_form_with_five_independent_linear_rows() {
    uint16_t NROWS = ROLLO_I_CODE_LENGTH;
    bf_element_t mat[NROWS];
    bf_element_t tmp;
    size_t dim = 0;

//    // fill matrix
//    for (i = 0; i < NROWS; i++) {
//        bf113_weak_rand(&mat[i]);
//    }

    // fill matrix
    bf_set_from_hex_string(&mat[0], "7FFFFFFFFFFFFFFF1");
    bf_set_from_hex_string(&mat[1], "7FFFFFFFFFFFFFFF2");
    bf_set_from_hex_string(&mat[2], "7FFFFFFFFFFFFFFF4");
    bf_set_from_hex_string(&mat[3], "7FFFFFFFFFFFFFFF8");
    bf_set_from_hex_string(&mat[4], "7FFFFFFFFFFFFFFF0");
    for (uint16_t i = 5; i < NROWS; ++i)
        bf_set_from_hex_string(&mat[i], "7FFFFFFFFFFFFFFF0");

//    printf("binary: A = \n");
//    bf_print_matrix_bin(mat, NROWS);

    dim = to_row_echelon_form(mat, NROWS);

//    printf("binary: A = \n");
//    bf_print_matrix_bin(mat, NROWS);

    if (dim != 5) {
        printf("computed dim = %ld\n", dim);
        printf("expected dim = 5\n");
        printf("\n\n");
        return 1;
    }

    return 0;
}


int test_zassenhaus_algorithm_reduction_with_random_values() {
    uint16_t NROWS = 2 * ROLLO_I_FINITE_FIELD_DEGREE;
//    uint16_t NROWS = 2 * 2;
//    uint16_t NROWS = 5;
    uint16_t expected_dimension = NROWS;

    if (NROWS >= 2 * ROLLO_I_FINITE_FIELD_DEGREE) {
        expected_dimension = 2 * ROLLO_I_FINITE_FIELD_DEGREE;
    }

//    uint16_t NROWS = 128;
    bf_element_t matrix_right[NROWS];
    bf_element_t matrix_left[NROWS];

    uint8_t seed[32] = {10,};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    size_t dim = 0;

    uint8_t verb = 0;

    // fill matrix
    // fill first half
    bf_element_t x;
    bf_set_from_hex_string(&x, "2");
    bf_set_from_hex_string(&matrix_left[0], "1");
    bf_get_random_element(&prng, &matrix_right[0]);
    for (uint16_t i = 1; i < NROWS / 2; i++) {
        bf_multiplication(&matrix_left[i], &matrix_left[i - 1], &x);
        bf_get_random_element(&prng, &matrix_right[i]);
    }
    // fill second half
    bf_set_to_zero(&matrix_left[NROWS / 2]);
    bf_set_from_hex_string(&matrix_right[NROWS / 2], "1");
    for (uint16_t i = NROWS / 2 + 1; i < NROWS; i++) {
        bf_set_to_zero(&matrix_left[i]);
        bf_multiplication(&matrix_right[i], &matrix_right[i - 1], &x);
    }

//    for (uint16_t i = 0; i < NROWS; i++) {
//        bf_get_random_element(&prng, &matrix_left[i]);
//        bf_get_random_element(&prng, &matrix_right[i]);
//    }

    if (verb) {
        printf("\nbinary: left || right \n");
        for (uint16_t i = 0; i < NROWS; ++i) {
            printf("row: %3d = ", i);
            bf_print_binary(&matrix_left[i]);
            printf(" || ");
            bf_print_binary(&matrix_right[i]);
            printf("\n");
        }
    }


    dim = zassenhaus_algorithm_reduction(matrix_left, matrix_right, NROWS);

    if (dim != expected_dimension) {
        printf("\nbinary: left || right \n");
        for (uint16_t i = 0; i < NROWS; ++i) {
            printf("row: %3d = ", i);
            bf_print_binary(&matrix_left[i]);
            printf(" || ");
            bf_print_binary(&matrix_right[i]);
            printf("\n");
        }

        printf("computed dim = %ld\n", dim);
        printf("expected dim = %d\n", expected_dimension);
        printf("\n\n");
        return 1;
    }

    return 0;
}

int test_zassenhaus_algorithm_returns_correct_dimension() {
    uint16_t s1_size = ROLLO_I_CODE_LENGTH / 2;
    uint16_t s2_size = ROLLO_I_CODE_LENGTH / 2;
    uint16_t out_size = 2 * ROLLO_I_CODE_LENGTH / 2;
    bf_element_t s1[s1_size];
    bf_element_t s2[s2_size];
    bf_element_t out[out_size];
    size_t out_dim_sum, out_dim_intersection;
    size_t dim_s1, dim_s2;

    uint8_t seed[32] = {10,};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    uint8_t verb = 0;

    // fill matrix
    for (uint16_t i = 0; i < s1_size; i++) {
        bf_get_random_element(&prng, &s1[i]);
    }
    for (uint16_t i = 0; i < s2_size; i++) {
        bf_get_random_element(&prng, &s2[i]);
    }

    if (verb) {
        printf("\nbinary: s1 || s2 \n");
        for (uint16_t i = 0; i < s1_size; ++i) {
            printf("row: %3d = ", i);
            bf_print_binary(&s1[i]);
            printf(" || ");
            bf_print_binary(&s2[i]);
            printf("\n");
        }
    }

    // s1 + s2 (sum)

    out_dim_sum = zassenhaus_algorithm_sum(out, out_size, s1, s1_size, s2, s2_size);

    if (verb) {
        printf("\nbinary: s1+s2 = \n");
        for (uint16_t i = 0; i < out_size; ++i) {
            printf("row: %3d = ", i);
            bf_print_binary(&out[i]);
            printf("\n");
        }
        printf("dim(s1+s2) = %lu\n", out_dim_sum);
    }

    // s1 ^ s2 (intersection)
    out_dim_intersection = zassenhaus_algorithm_intersection(out, out_size, s1, s1_size, s2, s2_size);

    if (verb) {
        printf("binary: s1*s2 = \n");
        for (uint16_t i = 0; i < out_size; ++i) {
            printf("row: %3d = ", i);
            bf_print_binary(&out[i]);
            printf("\n");
        }
        printf("dim(s1*s2) = %lu\n", out_dim_intersection);
    }

    dim_s1 = to_row_echelon_form(s1, s1_size);
    dim_s2 = to_row_echelon_form(s2, s2_size);
    if (out_dim_sum + out_dim_intersection != dim_s1 + dim_s2) {
        printf("\n");
        printf("dim(S1+S2) = %lu\n", out_dim_sum);
        printf("dim(S1^S2) = %lu\n", out_dim_intersection);
        printf("dim(S1) = %lu\n", dim_s1);
        printf("dim(S2) = %lu\n", dim_s2);
        return 1;
    }
    return 0;
}

int test_generate_random_support_basis() {
    uint32_t basis_length = 10;
    bf_element_t basis[basis_length];

    uint8_t seed[32] = {10,};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    uint8_t rank;
    uint8_t verb = 0;

    generate_random_support_basis(basis, basis_length, &prng);

    if (verb) {
        printf("\nError list:\n");
        bf_print_matrix_bin(basis, basis_length);
        printf("\n");
    }

    rank = to_row_echelon_form(basis, basis_length);

    if (verb) {
        printf("\nReduced error list:\n");
        bf_print_matrix_bin(basis, basis_length);
        printf("\n");
    }

    if (rank != basis_length) {
        printf("\nError! Found rank = %d, expecting rank = %d\n", rank, basis_length);
        return TEST_FAILURE;
    }

    return TEST_SUCCESS;
}

int test_generate_random_linear_combination() {
    uint8_t basis_length = 10;
    bf_element_t basis[basis_length];
    bf_element_t linear_combination;
    bf_element_t matrix[basis_length + 1];

    uint8_t seed[32] = {10,};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    uint8_t rank;
    uint8_t verb = 0;

    generate_random_support_basis(basis, basis_length, &prng);
    generate_random_linear_combination(&linear_combination, basis, basis_length, &prng);

    for (uint8_t i = 0; i < basis_length; ++i) {
        bf_copy(&matrix[i], &basis[i]);
    }
    bf_copy(&matrix[basis_length], &linear_combination);

    if (verb) {
        printf("\nLinear combination = \n");
        bf_print_binary(&linear_combination);
        printf("\n");
        printf("\nBasis list:\n");
        bf_print_matrix_bin(matrix, basis_length + 1);
        printf("\n");
    }

    rank = to_row_echelon_form(matrix, basis_length + 1);

    if (verb) {
        printf("\nReduced basis list:\n");
        bf_print_matrix_bin(matrix, basis_length + 1);
        printf("\n");
    }

    if (rank != basis_length) {
        printf("\nError! Found rank = %d, expecting rank = %d\n", rank, basis_length);
        return TEST_FAILURE;
    }

    return TEST_SUCCESS;
}

int test_that_generate_list_of_vectors_with_given_rank_returns_correct_rank() {
    // NOTE: this test might fail with probability 1/2^ROLLO_I_ERROR_VECTORS_RANK_WEIGHT
    // thus approximately every 1/2^ROLLO_I_ERROR_VECTORS_RANK_WEIGHT tests
    uint8_t log2_vector_list_length = ROLLO_I_ERROR_VECTORS_RANK_WEIGHT;
    uint16_t vector_list_length = 1u << ROLLO_I_ERROR_VECTORS_RANK_WEIGHT;
    bf_element_t vector_list[vector_list_length];
    bf_element_t vector_list_basis[log2_vector_list_length];

    uint8_t seed[32] = {10,};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    uint8_t rank;


    uint8_t verb = 0;

//    for (uint8_t i = 0; i < log2_vector_list_length; ++i) {
    uint8_t i = log2_vector_list_length;
    for (uint32_t j = 0; j < (1u << log2_vector_list_length); ++j) {
        bf_set_to_zero(&vector_list[j]);
    }
    generate_list_of_vectors_with_given_rank(vector_list,
                                             1u << i,
                                             i,
                                             vector_list_basis,
                                             &prng);

    if (verb) {
        printf("\ni = %d\n", i);
        printf("\nVector list:\n");
        bf_print_matrix_bin(vector_list, 1u << i);
        printf("\n");

        printf("\nVector list basis:\n");
        bf_print_matrix_bin(vector_list_basis, i);
        printf("\n");
    }

    rank = to_row_echelon_form(vector_list, 1u << i);

    if (verb) {
        printf("\nReduced vector list:\n");
        bf_print_matrix_bin(vector_list, 1u << i);
        printf("\n");
    }

    if (rank != i) {
        printf("\nError! Found rank = %d, expecting rank = %d\n", rank, i);
        return TEST_FAILURE;
    }
//    }
    return TEST_SUCCESS;
}

int test_rank_support_recover_returns_the_same_reduced_basis() {

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    // F rank(F) = d
    bf_element_t basis_of_subspaceF[ROLLO_I_KEY_VECTORS_RANK_WEIGHT];

    // E rank(E) = r
    uint16_t subspaceE_basis_length = ROLLO_I_ERROR_VECTORS_RANK_WEIGHT;
    bf_element_t subspaceE_basis[subspaceE_basis_length];
    bf_element_t recovered_subspaceE_basis[subspaceE_basis_length];

    polynomial_t x;
    polynomial_t y;
    polynomial_t h;

    int16_t i, j;
    polynomial_t e1, e2;
    polynomial_t c;

    polynomial_t s;

    uint8_t verb = 0;

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
            subspaceE_basis, ROLLO_I_ERROR_VECTORS_RANK_WEIGHT,
            &prng);
    polynomial_multiplication_mod_p(&c, &h, &e2);
    polynomial_addition(&c, &c, &e1);

    // compute reduced row echelon form of E basis
    to_reduced_row_echelon_form(subspaceE_basis, subspaceE_basis_length);

    // compute syndrome
    polynomial_multiplication_mod_p(&s, &x, &c);

    if (verb) {
        printf("\nRank(x) = %d\n", rank((bf_element_t *) &x, ROLLO_I_CODE_LENGTH));
        printf("Rank(y) = %d\n", rank((bf_element_t *) &y, ROLLO_I_CODE_LENGTH));
        printf("Rank(h) = %d\n", rank((bf_element_t *) &h, ROLLO_I_CODE_LENGTH));
        printf("Rank(e1) = %d\n", rank((bf_element_t *) &e1, ROLLO_I_CODE_LENGTH));
        printf("Rank(e2) = %d\n", rank((bf_element_t *) &e2, ROLLO_I_CODE_LENGTH));
        printf("Rank(c) = %d\n", rank((bf_element_t *) &c, ROLLO_I_CODE_LENGTH));
    }

    rank_support_recover(recovered_subspaceE_basis, basis_of_subspaceF, s);

    // checks
    uint8_t are_equal = 1;

    // recovered E contained in E
    for (i = 0; i < (subspaceE_basis_length); ++i) {
        if (!bf_is_equal(&recovered_subspaceE_basis[i], &subspaceE_basis[i])) {
            are_equal = 0;
        }
    }

    if (!are_equal) {
        printf("\nERROR! The sorted recovered subspace E is NOT equal to the initial sorted subspace E!\n");
        printf("\nrecovered subspace E ?= subspace E\n");
        for (i = 0; i < (subspaceE_basis_length); ++i) {
            bf_print_binary(&recovered_subspaceE_basis[i]);
            printf("  =?  ");
            bf_print_binary(&subspaceE_basis[i]);
            printf(" equals? %d\n", bf_is_equal(&recovered_subspaceE_basis[i], &subspaceE_basis[i]));
        }
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}