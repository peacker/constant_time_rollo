#ifndef BF_TESTS_H
#define BF_TESTS_H

// test utils

int test_uint8_t_ith_coefficient();

int test_bf_set_mask();

int test_bf_compute_mask();

int test_bf_compute_mask_inverse();

// Test vectors bf 67

int test_bf_add_test_vector();

int test_bf_double_addition_test_vector();

int test_bf_multiplication_without_reduction_test_vector();

int test_bf_reduction_test_vector();

int test_bf_multiplication_test_vector();

int test_bf_inversion_test_vector();

// Unit tests bf

int test_bf_is_greater_than();

int test_bf_add_is_commutative();

int test_bf_mul_is_commutative();

int test_bf_sqr_x_equals_mul_x_x();
int test_bf_inv();


#endif // BF_TESTS_H
