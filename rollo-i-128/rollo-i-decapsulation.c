#include <string.h>
#include <NIST_standard_functions/hash/hash.h>
#include "rollo-i-decapsulation.h"
#include "vector_space.h"


void rollo_I_decapsulation(uint8_t *shared_secret, polynomial_t *cipher, rollo_I_128_private_key_t *private_key) {

    polynomial_t syndrome;

    polynomial_multiplication_mod_p(&syndrome, cipher, &private_key->x);

    size_t error_basis_list_length = ROLLO_I_ERROR_VECTORS_RANK_WEIGHT;
    bf_element_t error_basis_list[ROLLO_I_ERROR_VECTORS_RANK_WEIGHT];

    //    initialize error_basis_list
    for (uint16_t i = 0; i < ROLLO_I_ERROR_VECTORS_RANK_WEIGHT; ++ i) {
        bf_set_to_zero(&error_basis_list[i]);
    }

    rank_support_recover(error_basis_list, private_key->basis_of_support, syndrome);

    uint16_t serialized_error_list_basis_len = ROLLO_I_BF_ELEMENT_BYTE_SIZE * ROLLO_I_ERROR_VECTORS_RANK_WEIGHT;
    uint8_t serialized_error_list[serialized_error_list_basis_len];
    bf_serialize_list(serialized_error_list, serialized_error_list_basis_len, error_basis_list, error_basis_list_length);

    sha256(shared_secret, (uint8_t *) serialized_error_list, serialized_error_list_basis_len);
}

