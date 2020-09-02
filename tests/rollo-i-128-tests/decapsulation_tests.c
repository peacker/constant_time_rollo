#include <stdio.h>
#include <rollo-i-128/rollo-i-keygen.h>
#include <rollo-i-128/rollo-i-encapsulation.h>
#include <rollo-i-128/rollo-i-decapsulation.h>
#include "decapsulation_tests.h"

int test_decapsulation() {

    size_t i;
    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, 0x100000000 - 1);

    rollo_I_key_pair_t key_pair;
    rollo_I_keygen(&key_pair, &prng);

    uint8_t shared_secret[ROLLO_I_SHARED_SECRET_BYTE_SIZE] = {0,};
    uint8_t recovered_shared_secret[ROLLO_I_SHARED_SECRET_BYTE_SIZE] = {0,};
    polynomial_t cipher;

    uint8_t verb = 0;

    rollo_I_encapsulation(shared_secret, &cipher, &prng, &key_pair.public_key);

    rollo_I_decapsulation(recovered_shared_secret, &cipher, &key_pair.private_key);

    uint8_t tmp = 0;
    for (i = 0; i < ROLLO_I_SHARED_SECRET_BYTE_SIZE; i++) {
        tmp = tmp | (uint8_t) (shared_secret[i] ^ recovered_shared_secret[i]);
    }

    if (tmp == 0) {
        if (verb) {
            printf("\nshared secret = \n");
            for (i = 0; i < ROLLO_I_SHARED_SECRET_BYTE_SIZE; ++i) {
                printf("%02x", shared_secret[i]);
            }
            printf("\nrecovered shared secret = \n");
            for (i = 0; i < ROLLO_I_SHARED_SECRET_BYTE_SIZE; ++i) {
                printf("%02x", recovered_shared_secret[i]);
            }
            printf("\n");
        }
        return EXIT_SUCCESS;
    } else {
        printf("\nshared secret = \n");
        for (i = 0; i < ROLLO_I_SHARED_SECRET_BYTE_SIZE; ++i) {
            printf("%02x", shared_secret[i]);
        }
        printf("\nrecovered shared secret = \n");
        for (i = 0; i < ROLLO_I_SHARED_SECRET_BYTE_SIZE; ++i) {
            printf("%02x", recovered_shared_secret[i]);
        }
        printf("\n");
        return EXIT_FAILURE;
    }
}