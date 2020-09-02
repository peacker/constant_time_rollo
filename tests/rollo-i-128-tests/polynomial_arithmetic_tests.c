#include <rollo-i-128/rollo_i_polynomial_arithmetic.h>
#include <memory.h>
#include "polynomial_arithmetic_tests.h"


int test_polynomial_multiplication_mod_p() {

    polynomial_t p1, p2, p3, test_vector;
    // set to 0 all coefficients
    memset(&p1, 0x00, sizeof(polynomial_t));
    memset(&p2, 0x00, sizeof(polynomial_t));
    memset(&p3, 0x00, sizeof(polynomial_t));
    memset(&test_vector, 0x00, sizeof(polynomial_t));

//    // p1 = x^66
//    bf_set_from_hex_string(&p1.coefficients[66], "1");
//
//    // p2 = x^66
//    bf_set_from_hex_string(&p2.coefficients[66], "1");
//
//    //test_vector = x^66 + x^65 + x^8 + x^4 + x^3 + x^2 + x + 1
//    bf_set_from_hex_string(&test_vector.coefficients[66], "1");
//    bf_set_from_hex_string(&test_vector.coefficients[65], "1");
//    bf_set_from_hex_string(&test_vector.coefficients[8], "1");
//    bf_set_from_hex_string(&test_vector.coefficients[4], "1");
//    bf_set_from_hex_string(&test_vector.coefficients[3], "1");
//    bf_set_from_hex_string(&test_vector.coefficients[2], "1");
//    bf_set_from_hex_string(&test_vector.coefficients[1], "1");
//    bf_set_from_hex_string(&test_vector.coefficients[0], "1");

    // p1 = x^82
    bf_set_from_hex_string(&p1.coefficients[82], "1");

    // p2 = x
    bf_set_from_hex_string(&p2.coefficients[1], "1");

    //test_vector = x^7  + x^4 + x^2 + 1
    bf_set_from_hex_string(&test_vector.coefficients[7], "1");
    bf_set_from_hex_string(&test_vector.coefficients[4], "1");
    bf_set_from_hex_string(&test_vector.coefficients[2], "1");
    bf_set_from_hex_string(&test_vector.coefficients[0], "1");

    polynomial_multiplication_mod_p(&p3, &p1, &p2);

    if (!rollo_i_polynomial_is_equal(&test_vector, &p3)) {
        return 1;
    }

    return 0;
}

int test_polynomial_inversion() {

    int i;
    uint8_t seed[32] = {0};
    AES_XOF_struct prng;
    seedexpander_init(&prng, seed, seed, ((uint64_t) 1) << 32u);

    polynomial_t p1, p2, p3;
    for (i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
        bf_get_random_element(&prng, &p1.coefficients[i]);
    }

    polynomial_inversion_mod_p(&p2, &p1);
    polynomial_multiplication_mod_p(&p3, &p2, &p1);

    if (!polynomial_is_one(&p3)) {
        return 1;
    }

    return 0;
}