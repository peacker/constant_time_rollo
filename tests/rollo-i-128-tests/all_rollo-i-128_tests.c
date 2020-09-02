#include <stdio.h>
#include <stdlib.h>
#include <rollo-i-128/rollo-i-128_parameters.h>
#include "all_rollo-i-128_tests.h"
#include "bf_tests.h"
#include "vector_space_tests.h"
#include "polynomial_arithmetic_tests.h"
#include "encapsulation_tests.h"
#include "keygen_tests.h"
#include "decapsulation_tests.h"
#include "performance_comparisons.h"

static void print_test_result(int result) {
    if (result == 0) {
        printf("PASS\n");
    } else {
        printf("*** FAILED ***\n");
    }
}

int rollo_i_128_unit_tests() {

    int result;
    int total_failed = 0;

    printf("===== ROLLO-I-128 Tests ===== \n");

    printf("===== ROLLO-I-128 Binary field GF(2^%d) ===== \n", ROLLO_I_FINITE_FIELD_DEGREE);

    printf("Executing test_uint8_t_ith_coefficient()... ");
    result = test_uint8_t_ith_coefficient();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_bf_set_mask()... ");
    result = test_bf_set_mask();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_bf_compute_mask()... ");
    result = test_bf_compute_mask();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_bf_compute_mask_inverse()... ");
    result = test_bf_compute_mask_inverse();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_bf_add_test_vector()... ");
    result = test_bf_add_test_vector();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_bf_double_addition_test_vector()... ");
    result = test_bf_double_addition_test_vector();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_bf_multiplication_without_reduction_test_vector()... ");
    result = test_bf_multiplication_without_reduction_test_vector();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_bf_reduction_test_vector()... ");
    result = test_bf_reduction_test_vector();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_bf_multiplication_test_vector()... ");
    result = test_bf_multiplication_test_vector();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_bf_inversion_test_vector()... ");
    result = test_bf_inversion_test_vector();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_bf_is_greater_than()... ");
    result = test_bf_is_greater_than();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_bf_add_is_commutative()... ");
    result = test_bf_add_is_commutative();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_bf_mul_is_commutative()... ");
    result = test_bf_mul_is_commutative();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_bf_inv()... ");
    result = test_bf_inv();
    print_test_result(result);
    total_failed += result;

    printf("===== ROLLO-I-128 Binary Vector Space of length %d over GF(2) ===== \n", ROLLO_I_FINITE_FIELD_DEGREE);

    printf("Executing test_uint16_t_is_greater_than()... ");
    result = test_uint16_t_is_greater_than();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_to_reduced_row_echelon_form()... ");
    result = test_to_reduced_row_echelon_form();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_to_systematic_form()... ");
    result = test_to_systematic_form();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_to_row_echelon_form_for_a_matrix_with_identical_entries()... ");
    result = test_to_row_echelon_form_for_a_matrix_with_identical_entries();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_to_row_echelon_form_with_five_independent_linear_rows()... ");
    result = test_to_row_echelon_form_with_five_independent_linear_rows();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_zassenhaus_algorithm_reduction_with_random_values()... ");
    result = test_zassenhaus_algorithm_reduction_with_random_values();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_zassenhaus_algorithm_returns_correct_dimension()... ");
    result = test_zassenhaus_algorithm_returns_correct_dimension();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_generate_random_support_basis()... ");
    result = test_generate_random_support_basis();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_generate_random_linear_combination()... ");
    result = test_generate_random_linear_combination();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_that_generate_list_of_vectors_with_given_rank_returns_correct_rank()... ");
    result = test_that_generate_list_of_vectors_with_given_rank_returns_correct_rank();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_rank_support_recover_returns_the_same_reduced_basis()... ");
    result = test_rank_support_recover_returns_the_same_reduced_basis();
    print_test_result(result);
    total_failed += result;

    printf("===== ROLLO-I-128 Polynomial Arithmetic ===== \n");

    printf("Executing test_polynomial_multiplication_mod_p... ");
    result = test_polynomial_multiplication_mod_p();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_polynomial_inversion... ");
    result = test_polynomial_inversion();
    print_test_result(result);
    total_failed += result;

    printf("===== ROLLO-I-128 KEM ===== \n");

    printf("Executing test_keygen... ");
    result = test_keygen();
    print_test_result(result);
    total_failed += result;

    printf("Executing test_encapsulation()... ");
    result = test_encapsulation();
    total_failed += result;
    print_test_result(result);

    printf("Executing test_decapsulation()... ");
    result = test_decapsulation();
    print_test_result(result);
    total_failed += result;

    return total_failed;
}

int rollo_i_128_performance_tests(){

    printf("===== ROLLO-I-128 Performance ===== \n");

    printf("\nBinary field arithmetic\n\n");

    printf("Executing performance_of_unreduced_bf_multiplications:");
    performance_of_unreduced_bf_multiplications();

    #if AVX2
    printf("Executing performance_of_unreduced_bf_multiplications_karatsuba:");
    performance_of_unreduced_bf_multiplications_karatsuba();
    #endif

    printf("Executing performance_of_unreduced_squaring:");
    performance_of_unreduced_squaring();

    printf("Executing performance_bf_reductions:");
    performance_of_bf_reductions();

    printf("Executing performance_bf_reductions:");
    performance_of_bf_reductions_uint64();

    printf("Executing performance_of_bf_inversion:");
    performance_of_bf_inversion();

    printf("\nBinary vector space arithmetic\n\n");

    printf("Executing performance_of_to_row_echelon_form_10_rows:");
    performance_of_to_row_echelon_form_10_rows();

    printf("Executing performance_of_to_row_echelon_form_20_rows:");
    performance_of_to_row_echelon_form_20_rows();

    printf("Executing performance_of_to_row_echelon_form_30_rows:");
    performance_of_to_row_echelon_form_30_rows();

    printf("Executing performance_of_to_row_echelon_form_100_rows:");
    performance_of_to_row_echelon_form_100_rows();

    printf("Executing performance_of_to_reduced_row_echelon_form_10_rows:");
    performance_of_to_reduced_row_echelon_form_10_rows();

    printf("Executing performance_of_to_reduced_row_echelon_form_20_rows:");
    performance_of_to_reduced_row_echelon_form_20_rows();

    printf("Executing performance_of_to_reduced_row_echelon_form_30_rows:");
    performance_of_to_reduced_row_echelon_form_30_rows();

    printf("Executing performance_of_to_reduced_row_echelon_form_100_rows:");
    performance_of_to_reduced_row_echelon_form_100_rows();

    printf("Executing performance_of_to_systematic_form_10_rows:");
    performance_of_to_systematic_form_10_rows();

    printf("Executing performance_of_to_systematic_form_20_rows:");
    performance_of_to_systematic_form_20_rows();

    printf("Executing performance_of_to_systematic_form_30_rows:");
    performance_of_to_systematic_form_30_rows();

    printf("Executing performance_of_to_systematic_form_100_rows:");
    performance_of_to_systematic_form_100_rows();

    printf("Executing performance_generate_list_of_vectors_with_given_rank_RANK_7_LIST_LEN_10:");
    performance_generate_list_of_vectors_with_given_rank_RANK_7_LIST_LEN_10();

    printf("Executing performance_generate_list_of_vectors_with_given_rank_RANK_7_LIST_LEN_100:");
    performance_generate_list_of_vectors_with_given_rank_RANK_7_LIST_LEN_100();

    printf("Executing performance_generate_two_lists_of_vectors_with_given_rank_RANK_7_LIST_LEN_10:");
    performance_generate_two_lists_of_vectors_with_given_rank_RANK_7_LIST_LEN_10();

    printf("Executing performance_generate_two_lists_of_vectors_with_given_rank_RANK_7_LIST_LEN_100:");
    performance_generate_two_lists_of_vectors_with_given_rank_RANK_7_LIST_LEN_100();

    printf("Executing performance_zassenhaus_algorithm_sum_reduction_ROWS_83:");
    performance_zassenhaus_algorithm_sum_reduction_ROWS_83();

    printf("Executing performance_zassenhaus_algorithm_intersection_reduction_ROWS_83:");
    performance_zassenhaus_algorithm_intersection_reduction_ROWS_83();

    printf("Executing performance_zassenhaus_algorithm_sum_reduction_ROWS_100:");
    performance_zassenhaus_algorithm_sum_reduction_ROWS_100();

    printf("Executing performance_zassenhaus_algorithm_intersection_reduction_ROWS_100:");
    performance_zassenhaus_algorithm_intersection_reduction_ROWS_100();

    printf("Executing performance_rank_support_recover:");
    performance_rank_support_recover();

    printf("\nComposite Galois field arithmetic\n\n");

    printf("Executing performance_polynomial_multiplication:");
    performance_polynomial_multiplication();

    printf("Executing performance_polynomial_inversion:");
    performance_polynomial_inversion();

    printf("Executing performance_polynomial_inversion_with_precomputed_matrices:");
    performance_polynomial_inversion_with_precomputed_matrices();

    printf("\nScheme routines\n\n");

    printf("Executing performance_keygen:");
    performance_keygen();

    printf("Executing performance_encapsulation:");
    performance_encapsulation();

    printf("Executing performance_decapsulation:");
    performance_decapsulation();

    return EXIT_SUCCESS;
}