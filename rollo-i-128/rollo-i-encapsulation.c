#include <NIST_standard_functions/hash/hash.h>
#include "rollo-i-encapsulation.h"
#include "vector_space.h"

int rollo_I_encapsulation(uint8_t *shared_secret, polynomial_t *cipher, AES_XOF_struct *prng,
                          rollo_I_128_public_key_t *public_key) {

    uint16_t error_basis_list_len = ROLLO_I_ERROR_VECTORS_RANK_WEIGHT;
    bf_element_t error_basis_list[error_basis_list_len];

    polynomial_t e1, e2;

    //    initialize and error_basis_list
    for (uint16_t i = 0; i < ROLLO_I_ERROR_VECTORS_RANK_WEIGHT; ++ i) {
        bf_set_to_zero(&error_basis_list[i]);
    }

    generate_basis_and_sample_two_polynomials_from_basis(
            &e1,
            &e2,
            error_basis_list, ROLLO_I_ERROR_VECTORS_RANK_WEIGHT,
            prng);

    polynomial_multiplication_mod_p(cipher, &public_key->h, &e2);

    polynomial_addition(cipher, cipher, &e1);

    to_reduced_row_echelon_form(error_basis_list, ROLLO_I_ERROR_VECTORS_RANK_WEIGHT);

    uint16_t serialized_error_list_basis_len = ROLLO_I_BF_ELEMENT_BYTE_SIZE * ROLLO_I_ERROR_VECTORS_RANK_WEIGHT;
    uint8_t serialized_error_list[serialized_error_list_basis_len];
    bf_serialize_list(serialized_error_list, serialized_error_list_basis_len, error_basis_list, error_basis_list_len);

    sha256(shared_secret, (uint8_t *) serialized_error_list, serialized_error_list_basis_len);

    return 0;
}