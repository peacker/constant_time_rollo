#include <stdio.h>
#include <rollo-i-128/bf.h>
#include "bf_tests.h"
#include "all_rollo-i-128_tests.h"
#include "../../rollo-i-128/rollo-i-128_parameters.h"

// test utils


int test_uint8_t_ith_coefficient() {
    uint8_t coefficient;
    uint8_t expected_coefficient;
    bf_element_t a;

    bf_set_to_ones(&a);
    expected_coefficient = 1;

    for (uint8_t i = 0; i < ROLLO_I_FINITE_FIELD_DEGREE; ++i) {
        coefficient = uint8_t_ith_coefficient(&a, i);
        if ( coefficient!= expected_coefficient) {
            printf("\ni = %d\n", i);
            printf("input vector:  ");
            bf_print_binary(&a);
            printf(" = ");
            bf_print(&a);
            printf("\n");
            printf("computed coefficient = %d\n", coefficient);
            printf("expected coefficient = %d\n", expected_coefficient);
            return EXIT_FAILURE;
        }
    }

    bf_set_to_zero(&a);
    expected_coefficient = 0;

    for (uint8_t i = 0; i < ROLLO_I_FINITE_FIELD_DEGREE; ++i) {
        coefficient = uint8_t_ith_coefficient(&a, i);
        if ( coefficient!= expected_coefficient) {
            printf("\ni = %d\n", i);
            printf("input vector:  ");
            bf_print_binary(&a);
            printf(" = ");
            bf_print(&a);
            printf("\n");
            printf("computed coefficient = %d\n", coefficient);
            printf("expected coefficient = %d\n", expected_coefficient);
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

int test_bf_set_mask() {
    bf_element_t mask;
    bf_element_t expected_mask;

    bf_set_to_ones(&expected_mask);

    bf_set_mask(&mask, 1);
    if (!bf_is_equal(&mask, &expected_mask)) {
        printf("\nbit = 1\n");
        printf("computed mask: ");
        bf_print_binary(&mask);
        printf(" = ");
        bf_print(&mask);
        printf("\n");
        printf("expected mask: ");
        bf_print_binary(&expected_mask);
        printf(" = ");
        bf_print(&expected_mask);
        printf("\n");
        return EXIT_FAILURE;
    }

    bf_set_to_zero(&expected_mask);

    bf_set_mask(&mask, 0);
    if (!bf_is_equal(&mask, &expected_mask)) {
        printf("\nbit = 1\n");
        printf("computed mask: ");
        bf_print_binary(&mask);
        printf(" = ");
        bf_print(&mask);
        printf("\n");
        printf("expected mask: ");
        bf_print_binary(&expected_mask);
        printf(" = ");
        bf_print(&expected_mask);
        printf("\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int test_bf_compute_mask() {
    bf_element_t mask;
    bf_element_t expected_mask;
    bf_element_t a;

    bf_set_to_ones(&a);
    bf_set_to_ones(&expected_mask);

    for (uint8_t i = 0; i < ROLLO_I_FINITE_FIELD_DEGREE; ++i) {
        bf_compute_mask(&mask, &a, i);
        if (!bf_is_equal(&mask, &expected_mask)) {
            printf("\ni = %d\n", i);
            printf("input vector:  ");
            bf_print_binary(&a);
            printf(" = ");
            bf_print(&a);
            printf("\n");
            printf("computed mask: ");
            bf_print_binary(&mask);
            printf(" = ");
            bf_print(&mask);
            printf("\n");
            printf("expected mask: ");
            bf_print_binary(&expected_mask);
            printf(" = ");
            bf_print(&expected_mask);
            printf("\n");
            return EXIT_FAILURE;
        }
    }

    bf_set_to_zero(&a);
    bf_set_to_zero(&expected_mask);

    for (uint8_t i = 0; i < ROLLO_I_FINITE_FIELD_DEGREE; ++i) {
        bf_compute_mask(&mask, &a, i);
        if (!bf_is_equal(&mask, &expected_mask)) {
            printf("\ni = %d\n", i);
            printf("input vector:  ");
            bf_print_binary(&a);
            printf(" = ");
            bf_print(&a);
            printf("\n");
            printf("computed mask: ");
            bf_print_binary(&mask);
            printf(" = ");
            bf_print(&mask);
            printf("\n");
            printf("expected mask: ");
            bf_print_binary(&expected_mask);
            printf(" = ");
            bf_print(&expected_mask);
            printf("\n");
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}


int test_bf_compute_mask_inverse() {
    bf_element_t mask;
    bf_element_t expected_mask;
    bf_element_t a;

    bf_set_to_ones(&a);
    bf_set_to_zero(&expected_mask);

    for (uint8_t i = 0; i < ROLLO_I_FINITE_FIELD_DEGREE; ++i) {
        bf_compute_mask_inverse(&mask, &a, i);
        if (!bf_is_equal(&mask, &expected_mask)) {
            printf("\ni = %d\n", i);
            printf("input vector:  ");
            bf_print_binary(&a);
            printf(" = ");
            bf_print(&a);
            printf("\n");
            printf("computed mask: ");
            bf_print_binary(&mask);
            printf(" = ");
            bf_print(&mask);
            printf("\n");
            printf("expected mask: ");
            bf_print_binary(&expected_mask);
            printf(" = ");
            bf_print(&expected_mask);
            printf("\n");
            return EXIT_FAILURE;
        }
    }

    bf_set_to_zero(&a);
    bf_set_to_ones(&expected_mask);

    for (uint8_t i = 0; i < ROLLO_I_FINITE_FIELD_DEGREE; ++i) {
        bf_compute_mask_inverse(&mask, &a, i);
        if (!bf_is_equal(&mask, &expected_mask)) {
            printf("\ni = %d\n", i);
            printf("input vector:  ");
            bf_print_binary(&a);
            printf(" = ");
            bf_print(&a);
            printf("\n");
            printf("computed mask: ");
            bf_print_binary(&mask);
            printf(" = ");
            bf_print(&mask);
            printf("\n");
            printf("expected mask: ");
            bf_print_binary(&expected_mask);
            printf(" = ");
            bf_print(&expected_mask);
            printf("\n");
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

// test vectors for bf 67

int test_bf_add_test_vector() {

    bf_element_t a, b, c, exp_c;

    bf_set_from_hex_string(&a, "74830087672108887");
    bf_set_from_hex_string(&b, "3aade08767f10a7c3");
    bf_set_from_hex_string(&exp_c, "4E2EE00000D002F44");

    bf_addition(&c, &a, &b);

    if (!bf_is_equal(&c, &exp_c)){
        return TEST_FAILURE;
    }

    return TEST_SUCCESS;

}

int test_bf_double_addition_test_vector() {

    bf_double_element_t a, b, c, test;

    bf_double_set_from_string(&a, "183c687647f884d7d0a8d62fb6f0edae9c");
    bf_double_set_from_string(&b, "183c687647f884d7d0a8d62fb6f0edae9c");
    bf_double_set_from_string(&test, "0000000000000000000000000000000000");

    bf_double_addition(&c, &a, &b);

    if (!bf_double_is_equal(&test, &c)) {
        printf("\na          = ");
        bf_double_element_print(&a);
        printf("\nb          = ");
        bf_double_element_print(&b);
        printf("\nc = a+b    = ");
        bf_double_element_print(&c);
        printf("\nexpected c = ");
        bf_double_element_print(&test);
        printf("\n");
        return TEST_FAILURE;
    }

    return TEST_SUCCESS;

}


int test_bf_multiplication_without_reduction_test_vector() {

    bf_element_t a, b;
    bf_double_element_t c, test;

    bf_set_from_hex_string(&a, "583c687645c872c00");
    bf_set_from_hex_string(&b, "7f884d7d0a8d62fb6");
    bf_double_set_from_string(&test, "1BECCDF4F046965D4A33164487B137A800");


    bf_unreduced_multiplication(&c, &a, &b);
//    bf_unreduced_multiplication_karatsuba(&c, &a, &b);

    if (!bf_double_is_equal(&test, &c)) {
        printf("\na          = ");
        bf_print(&a);
        printf("\nb          = ");
        bf_print(&b);
        printf("\nc = a*b    = ");
        bf_double_element_print(&c);
        printf("\nexpected c = ");
        bf_double_element_print(&test);
        printf("\n");
        return TEST_FAILURE;
    }
    return TEST_SUCCESS;

}



int test_bf_reduction_test_vector() {
    bf_element_t a, b, axb, test_vector;
    bf_double_element_t c;

    bf_set_from_hex_string(&a, "583c687645c872c00");
    bf_set_from_hex_string(&b, "7f884d7d0a8d62fb6");

    bf_unreduced_multiplication(&c, &a, &b);
    bf_reduction(&axb, &c);

    bf_set_from_hex_string(&test_vector, "5f0edae9c9152aedb");

    if (!bf_is_equal(&axb, &test_vector)) {
        printf("\n");
        printf("result          = ");
        bf_print(&axb);
        printf("\n");
        printf("expected result = ");
        bf_print(&test_vector);
        printf("\n");

        return TEST_FAILURE;
    }
    return TEST_SUCCESS;
}

int test_bf_multiplication_test_vector() {
    bf_element_t a, b, axb, test_vector;

    bf_set_from_hex_string(&a, "583c687645c872c00");
    bf_set_from_hex_string(&b, "7f884d7d0a8d62fb6");
    bf_set_from_hex_string(&test_vector, "5f0edae9c9152aedb");

    bf_multiplication(&axb, &a, &b);

    if (!bf_is_equal(&axb, &test_vector)) {
        printf("\n");
        printf("result          = ");
        bf_print(&axb);
        printf("\n");
        printf("expected result = ");
        bf_print(&test_vector);
        printf("\n");

        return TEST_FAILURE;
    }
    return TEST_SUCCESS;
}

int test_bf_inversion_test_vector() {
    bf_element_t a, b, axb, test_vector;

    bf_set_from_hex_string(&a, "583c687645c872c00");
    bf_set_from_hex_string(&test_vector, "1");

    bf_inversion(&b, &a);

    bf_multiplication(&axb, &a, &b);

    if (!bf_is_equal(&axb, &test_vector)) {
        printf("\n");
        printf("result          = ");
        bf_print(&axb);
        printf("\n");
        printf("expected result = ");
        bf_print(&test_vector);
        printf("\n");

        return TEST_FAILURE;
    }
    return TEST_SUCCESS;
}

// unit tests

int test_bf_is_greater_than() {
    bf_element_t a, b;

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    for (size_t i = 0; i < 100; ++i) {
        bf_get_random_element(&prng, &a);
        bf_get_random_element(&prng, &b);

        if (bf_greater_than(&a,&b) != bf_greater_than_not_constant_time(&a,&b)) {
            printf("\ni = %lu\n", i);
            printf("a = ");
            bf_print(&a);
            printf("\n");
            printf("b = ");
            bf_print(&b);
            printf("\n");
            printf("bf_greater_than(&a,&b)       = %d\n", bf_greater_than(&a,&b));
            printf("bf_greater_than_NO_CT(&a,&b) = %d\n", bf_greater_than_not_constant_time(&a,&b));
            printf("\n\n");
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

int test_bf_add_is_commutative() {

    bf_element_t a, b, axorb, bxora;

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    for (size_t i = 0; i < 100; ++i) {

        bf_get_random_element(&prng, &a);
        bf_get_random_element(&prng, &b);

        bf_addition(&axorb, &a, &b);
        bf_addition(&bxora, &b, &a);
        if (!bf_is_equal(&axorb, &bxora)) {
            printf("a+b = ");
            bf_print(&axorb);
            printf("\n\n");
            printf("b+a = ");
            bf_print(&bxora);
            printf("\n\n");
            return TEST_FAILURE;
        }
    }

    return TEST_SUCCESS;
}


int test_bf_mul_is_commutative() {

    bf_element_t a, b, axb, bxa;

    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    for (size_t i = 0; i < 100; ++i) {

        bf_get_random_element(&prng, &a);
        bf_get_random_element(&prng, &b);

        bf_multiplication(&axb, &a, &b);
        bf_multiplication(&bxa, &b, &a);
        if (!bf_is_equal(&axb, &bxa)) {
            printf("a*b = ");
            bf_print(&axb);
            printf("\n\n");
            printf("b*a = ");
            bf_print(&bxa);
            printf("\n\n");
            return TEST_FAILURE;
        }
    }

    return TEST_SUCCESS;
}

int test_bf_inv() {
    bf_element_t a, b, c;
    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    for (size_t i = 0; i < 100; ++i) {
        bf_get_random_element(&prng, &a);

        bf_inversion(&b, &a);
        bf_inversion(&c, &b);

        if (!bf_is_equal(&c, &a)) {
            printf("a = ");
            bf_print(&a);
            printf("\n\n");
            printf("(a^-1)^-1 = ");
            bf_print(&c);
            printf("\n\n");
            return TEST_FAILURE;
        }
    }

    return TEST_SUCCESS;
}