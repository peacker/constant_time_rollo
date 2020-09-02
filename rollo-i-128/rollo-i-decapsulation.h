#ifndef ROLLO_I_DECRYPTION_H
#define ROLLO_I_DECRYPTION_H

#include "rollo-i-keygen.h"

void rollo_I_decapsulation(uint8_t *shared_secret, polynomial_t *cipher, rollo_I_128_private_key_t *private_key);

#endif //ROLLO_I_DECRYPTION_H
