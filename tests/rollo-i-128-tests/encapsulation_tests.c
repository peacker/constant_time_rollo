#include <rollo-i-128/rollo-i-keygen.h>
#include <rollo-i-128/rollo-i-encapsulation.h>
#include <rollo-i-128/vector_space.h>
#include "encapsulation_tests.h"

int test_encapsulation() {

    int i;
    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    rollo_I_key_pair_t key_pair;
    rollo_I_keygen(&key_pair, &prng);

    uint8_t shared_secret[ROLLO_I_SHARED_SECRET_BYTE_SIZE] = {0,};
    polynomial_t cipher, cipher_prime;
    rollo_I_encapsulation(shared_secret, &cipher, &prng, &key_pair.public_key);

    uint8_t accumulator = 0;
    for (i = 0; i < ROLLO_I_SHARED_SECRET_BYTE_SIZE; i++) {
        accumulator |= shared_secret[i];
    }
    if (accumulator == 0) {
        return 1;
    }

    polynomial_multiplication_mod_p(&cipher_prime, &cipher, &key_pair.private_key.x);

    bf_element_t matrix[ROLLO_I_CODE_LENGTH];
    for (i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
        bf_copy(&matrix[i], &cipher_prime.coefficients[i]);
    }

    int rank = to_row_echelon_form(matrix, ROLLO_I_CODE_LENGTH);

    if (rank > ROLLO_I_ERROR_VECTORS_RANK_WEIGHT * ROLLO_I_KEY_VECTORS_RANK_WEIGHT) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}