/**
 * \brief A wrapper for OpenSSL SHA512
 */

#include "hash.h"
#include <openssl/sha.h>

void sha512(unsigned char* output, unsigned char* input, size_t size) {
    SHA512_CTX sha512;
    SHA512_Init(&sha512);
    SHA512_Update(&sha512, input, size);
    SHA512_Final(output, &sha512);
}

void sha256(unsigned char* output, unsigned char* input, size_t size) {
    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    SHA256_Update(&sha256, input, size);
    SHA256_Final(output, &sha256);
}