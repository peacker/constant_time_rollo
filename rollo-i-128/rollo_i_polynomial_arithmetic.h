#ifndef ROLLO_I_128_POLYNOMIAL_ARITHMETIC_H
#define ROLLO_I_128_POLYNOMIAL_ARITHMETIC_H

#include "bf.h"
#include "rollo-i-128_parameters.h"

typedef struct polynomial {
    // coefficients[0] is for degree 0 ???
    bf_element_t coefficients[ROLLO_I_CODE_LENGTH];
} polynomial_t;

void polynomial_set_to_zero(polynomial_t *p);

// returns 1 if the two polynomials have the same coefficients
int rollo_i_polynomial_is_equal(polynomial_t *a, polynomial_t *b);

int polynomial_is_one(polynomial_t *a);

void polynomial_addition(polynomial_t *c, polynomial_t *a, polynomial_t *b);

void polynomial_multiplication_mod_p(polynomial_t *c, polynomial_t *a, polynomial_t *b);

void polynomial_to_the_2_to_the_n_times_m(polynomial_t *a_to_the_2_to_the_n, polynomial_t *a,
                                          const uint128_t *matrix_s_to_the_n);

void polynomial_inversion_mod_p(polynomial_t *a_inverse, polynomial_t *a);
void polynomial_inversion_mod_p_with_precomputes_matrices(polynomial_t *a_inverse, polynomial_t *a);

int sample_two_polynomials_from_vector_space(polynomial_t *p1, polynomial_t *p2, bf_element_t *list, uint16_t list_length, uint8_t mask, AES_XOF_struct *prng);
int generate_basis_and_sample_two_polynomials_from_basis(
        polynomial_t *p1, polynomial_t *p2,
        bf_element_t *basis, uint16_t basis_length,
        AES_XOF_struct *prng);

#endif //ROLLO_I_128_POLYNOMIAL_ARITHMETIC_H
