#ifndef PERFORMANCE_COMPARISONS_H
#define PERFORMANCE_COMPARISONS_H

#include <stdint.h>
#include <time.h>

uint64_t rdtsc();

void printresults(const char* measured_function, clock_t total_time, uint64_t total_iterations, uint64_t total_cycles);

void performance_of_unreduced_bf_multiplications();

#if AVX2
void performance_of_unreduced_bf_multiplications_karatsuba();
#endif

void performance_of_unreduced_squaring();
void performance_of_bf_reductions();
void performance_of_bf_reductions_uint64();
void performance_of_bf_inversion();

void performance_of_to_row_echelon_form_10_rows();
void performance_of_to_row_echelon_form_20_rows();
void performance_of_to_row_echelon_form_30_rows();
void performance_of_to_row_echelon_form_100_rows();

void performance_of_to_reduced_row_echelon_form_10_rows();
void performance_of_to_reduced_row_echelon_form_20_rows();
void performance_of_to_reduced_row_echelon_form_30_rows();
void performance_of_to_reduced_row_echelon_form_100_rows();

void performance_of_to_systematic_form_10_rows();
void performance_of_to_systematic_form_20_rows();
void performance_of_to_systematic_form_30_rows();
void performance_of_to_systematic_form_100_rows();

void performance_generate_list_of_vectors_with_given_rank_RANK_7_LIST_LEN_10();
void performance_generate_list_of_vectors_with_given_rank_RANK_7_LIST_LEN_100();
void performance_generate_two_lists_of_vectors_with_given_rank_RANK_7_LIST_LEN_10();
void performance_generate_two_lists_of_vectors_with_given_rank_RANK_7_LIST_LEN_100();
void performance_zassenhaus_algorithm_sum_reduction_ROWS_83();
void performance_zassenhaus_algorithm_intersection_reduction_ROWS_83();
void performance_zassenhaus_algorithm_sum_reduction_ROWS_100();
void performance_zassenhaus_algorithm_intersection_reduction_ROWS_100();

void performance_rank_support_recover();

void performance_polynomial_multiplication();
void performance_polynomial_inversion();
void performance_polynomial_inversion_with_precomputed_matrices();

void performance_keygen();
void performance_encapsulation();
void performance_decapsulation();

#endif //PERFORMANCE_COMPARISONS_H
