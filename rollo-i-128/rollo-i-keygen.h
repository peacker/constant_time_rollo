#ifndef ROLLO_I_KEYGEN_H
#define ROLLO_I_KEYGEN_H

#include "rollo_i_polynomial_arithmetic.h"

typedef struct rollo_I_128_private_key {
    polynomial_t x;
    polynomial_t y;
    bf_element_t basis_of_support[ROLLO_I_KEY_VECTORS_RANK_WEIGHT];

} rollo_I_128_private_key_t;

typedef struct rollo_I_128_public_key {
    polynomial_t h;
} rollo_I_128_public_key_t;

typedef struct rollo_I_128_key_pair {
    rollo_I_128_private_key_t private_key;
    rollo_I_128_public_key_t public_key;
} rollo_I_key_pair_t;

void rollo_I_keypair_init(rollo_I_key_pair_t *key_pair);
int rollo_I_keygen(rollo_I_key_pair_t *key_pair, AES_XOF_struct *prng);

#endif //ROLLO_I_KEYGEN_H
