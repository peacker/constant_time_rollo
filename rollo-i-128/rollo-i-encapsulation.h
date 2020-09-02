#ifndef ROLLO_I_ENCAPSULATION_H
#define ROLLO_I_ENCAPSULATION_H

#include <NIST_standard_functions/rng/nist-rng.h>
#include "rollo_i_polynomial_arithmetic.h"
#include "rollo-i-keygen.h"

int rollo_I_encapsulation(uint8_t *shared_secret, polynomial_t *cipher, AES_XOF_struct *prng,
                          rollo_I_128_public_key_t *public_key);

#endif //ROLLO_I_ENCAPSULATION_H
