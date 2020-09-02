#ifndef VECTOR_SPACE_TESTS_H
#define VECTOR_SPACE_TESTS_H

// Unit test vector space

int test_uint16_t_is_greater_than();

int test_to_reduced_row_echelon_form();

int test_to_systematic_form();

int test_to_row_echelon_form_for_a_matrix_with_identical_entries();

int test_to_row_echelon_form_with_five_independent_linear_rows();

int test_zassenhaus_algorithm_reduction_with_random_values();

int test_zassenhaus_algorithm_returns_correct_dimension();

int test_generate_random_support_basis();

int test_generate_random_linear_combination();

int test_that_generate_list_of_vectors_with_given_rank_returns_correct_rank();

int test_rank_support_recover_returns_the_same_reduced_basis();

#endif //VECTOR_SPACE_TESTS_H
