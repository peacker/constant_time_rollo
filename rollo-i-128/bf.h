#ifndef BF_H
#define BF_H

#define BF_UINT128_T
//#define BF_UINT64_UINT8_T

#ifdef AVX2

#include <immintrin.h>

#endif

#include <NIST_standard_functions/rng/nist-rng.h>
#include <stdint.h>

typedef __uint128_t uint128_t;

#ifdef BF_UINT64_UINT8_T
typedef struct {
    uint8_t high;
    uint64_t low;
} bf_element_t;
#endif

#ifdef BF_UINT128_T
typedef struct {
    uint128_t value;
} bf_element_t;
#endif

typedef struct {
    __m256i value;
} bf_double_element_t;

// print

void bf_print(const bf_element_t *a);

void bf_double_element_print(const bf_double_element_t *a);

void bf_print_binary(const bf_element_t *a);

void bf_print_matrix_bin(bf_element_t *mat, size_t NROWS);

// set

void bf_set_to_zero(bf_element_t *c);

void bf_set_to_ones(bf_element_t *c);

uint8_t bf_set_from_hex_string(bf_element_t *result, const char *h);

uint8_t bf_double_set_from_string(bf_double_element_t *result, const char *h);

// basic

void bf_copy(bf_element_t *a, bf_element_t *b);

void bf_swap(bf_element_t *a, bf_element_t *b);

void bf_get_random_element(AES_XOF_struct *prng, bf_element_t *c);
void bf_get_list_of_random_elements(AES_XOF_struct *prng, bf_element_t *, uint16_t number_of_elements);

int bf_is_zero(bf_element_t *a);

int bf_is_one(bf_element_t *a);

int bf_is_equal(bf_element_t *a, bf_element_t *b);
int bf_double_is_equal(bf_double_element_t *a, bf_double_element_t *b);

int bf_greater_than(bf_element_t *a, bf_element_t *b);
int bf_greater_than_not_constant_time(bf_element_t *a, bf_element_t *b);

// constant time branching

int bf_set_mask(bf_element_t *mask, uint8_t bit);
int bf_compute_mask(bf_element_t *mask, bf_element_t *a, uint8_t bit_position);
int bf_compute_mask_inverse(bf_element_t *mask, bf_element_t *a, uint8_t bit_position);

// arithmetic operations

void bf_addition(bf_element_t *c, bf_element_t *a, bf_element_t *b);

void bf_double_addition(bf_double_element_t *c, bf_double_element_t *a, bf_double_element_t *b);

void bf_unreduced_multiplication(bf_double_element_t *c, const bf_element_t *a, const bf_element_t *b);

void bf_unreduced_multiplication_karatsuba(bf_double_element_t *c, const bf_element_t *a, const bf_element_t *b);

void bf_reduction(bf_element_t *reduced_a, bf_double_element_t *a);
void bf_reduction_uint64(uint64_t *o, const uint64_t *e);

void bf_multiplication(bf_element_t *c, bf_element_t *a, bf_element_t *b);

void bf_unreduced_square(bf_double_element_t *c, bf_element_t *a);

void bf_square(bf_element_t *c, bf_element_t *a);

void bf_inversion(bf_element_t *x_inverse, bf_element_t *x);

uint8_t uint8_t_ith_coefficient(bf_element_t *a, uint16_t i);

void bf_ith_coefficient(bf_element_t *b, bf_element_t *a, int i);

void bf_and(bf_element_t *c, bf_element_t *b, bf_element_t *a); // c = a & b

void bf_serialize_list(uint8_t * serialized_list, uint16_t serialized_list_len, bf_element_t *list, uint16_t number_of_elements);

#endif //BF_H
