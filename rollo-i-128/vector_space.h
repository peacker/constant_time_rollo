#ifndef VECTOR_SPACE_H
#define VECTOR_SPACE_H

#include "bf.h"
#include "rollo_i_polynomial_arithmetic.h"


int8_t uint16_t_is_greater_than(uint16_t a, uint16_t b);

int to_reduced_row_echelon_form(bf_element_t *matrix, uint16_t number_of_rows);
int to_systematic_form(bf_element_t *matrix, uint16_t number_of_rows);

int to_row_echelon_form(bf_element_t *matrix, uint16_t number_of_rows);
int rank(bf_element_t *input_matrix, uint16_t number_of_rows);

int generate_random_support_basis(bf_element_t *support_basis, uint16_t rank, AES_XOF_struct *prng);

int generate_support_from_basis(bf_element_t *support, bf_element_t *support_basis, uint16_t support_basis_size);

int generate_random_linear_combination(bf_element_t *output, bf_element_t *basis, u_int16_t basis_length, AES_XOF_struct *prng);

int generate_list_of_vectors_with_given_rank(bf_element_t *list, uint16_t list_length, u_int16_t rank, bf_element_t *list_basis, AES_XOF_struct *prng);
int generate_two_list_of_vectors_with_given_rank_and_their_basis(
        bf_element_t *list1, uint16_t list1_length,
        bf_element_t *list2, uint16_t list2_length,
        u_int16_t desired_rank, bf_element_t *list_basis, // TODO invert position of list basis and rank
        AES_XOF_struct *prng);

int zassenhaus_algorithm_reduction(bf_element_t *matrix_left, bf_element_t *matrix_right, uint16_t number_of_lines);

int zassenhaus_algorithm_intersection(bf_element_t *output, size_t output_length,
                                      bf_element_t *S1, size_t S1_length,
                                      bf_element_t *S2, size_t S2_length);

int zassenhaus_algorithm_sum(bf_element_t *output, size_t output_length,
                             bf_element_t *S1, size_t S1_length,
                             bf_element_t *S2, size_t S2_length);

int sort_list_of_vectors(bf_element_t *list, uint16_t list_length);

int rank_support_recover(bf_element_t *subspaceE, bf_element_t *subspaceF_basis, polynomial_t s);

#endif //VECTOR_SPACE_H
