#include <rollo-i-128/rollo-i-keygen.h>
#include "keygen_tests.h"

int test_keygen() {
    // check that h*x = y

    rollo_I_key_pair_t key_pair;

    uint8_t seed[32] = {10,};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    rollo_I_keygen(&key_pair, &prng);

    polynomial_t x_test;
    polynomial_multiplication_mod_p(&x_test, &key_pair.public_key.h, &key_pair.private_key.x);

    if (!rollo_i_polynomial_is_equal(&x_test, &key_pair.private_key.y)) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}