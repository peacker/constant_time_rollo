#include "rollo-i-keygen.h"

void rollo_I_keypair_init(rollo_I_key_pair_t *key_pair) {
    polynomial_set_to_zero(&key_pair->private_key.x);
    polynomial_set_to_zero(&key_pair->private_key.y);
    for (uint8_t i = 0; i < ROLLO_I_KEY_VECTORS_RANK_WEIGHT; ++i) {
        bf_set_to_zero(&key_pair->private_key.basis_of_support[i]);
    }
    polynomial_set_to_zero(&key_pair->public_key.h);
}

int rollo_I_keygen(rollo_I_key_pair_t *key_pair, AES_XOF_struct *prng) {

    generate_basis_and_sample_two_polynomials_from_basis(
            &key_pair->private_key.x,
            &key_pair->private_key.y,
            key_pair->private_key.basis_of_support, ROLLO_I_KEY_VECTORS_RANK_WEIGHT,
            prng);

//    polynomial_inversion_mod_p(&key_pair->public_key.h, &key_pair->private_key.x);
    polynomial_inversion_mod_p_with_precomputes_matrices(&key_pair->public_key.h, &key_pair->private_key.x);

    polynomial_multiplication_mod_p(&key_pair->public_key.h, &key_pair->public_key.h, &key_pair->private_key.y);

    return 0;
}