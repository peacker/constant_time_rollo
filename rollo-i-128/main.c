#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "rollo-i-encapsulation.h"
#include "rollo-i-keygen.h"
#include "rollo-i-decapsulation.h"

#define NBYTES 32
#define SEEDEXPANDER_MAX_LENGTH 4294967295

int main()
{
    printf("Generate keypair: ");
    AES_XOF_struct prng;
    uint8_t seed[32];
    seedexpander_init(&prng, seed, seed + 32, SEEDEXPANDER_MAX_LENGTH);

    rollo_I_key_pair_t keypair;
    rollo_I_keygen(&keypair, &prng);
    printf("DONE\n");


    printf("Encapsulation: ");
    uint8_t e_ss[NBYTES];
    polynomial_t ciphertext;
    rollo_I_encapsulation(e_ss, &ciphertext, &prng, &keypair.public_key);
    printf("DONE\n");


    printf("Decapsulation: ");
    uint8_t d_ss[NBYTES];
    rollo_I_decapsulation(d_ss, &ciphertext, &keypair.private_key);
    assert(memcmp((uint8_t *)e_ss, (uint8_t *)d_ss, NBYTES) == 0);
    printf("DONE\n");

    return 0;
}
