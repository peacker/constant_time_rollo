#include <stdlib.h>
#include "rollo_i_polynomial_arithmetic.h"
#include "rollo_i_polynomial_arithmetic_precomputed_matrices.h"
#include "vector_space.h"

void polynomial_set_to_zero(polynomial_t *p) {
    for (size_t i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
        bf_set_to_zero(&p->coefficients[i]);
    }
}

int rollo_i_polynomial_is_equal(polynomial_t *a, polynomial_t *b) {
    for (size_t i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
        if (!bf_is_equal(&a->coefficients[i], &b->coefficients[i])) {
            return 0;
        }
    }
    return 1;
}

int polynomial_is_one(polynomial_t *a) {
    if (!bf_is_one(&a->coefficients[0])) {
        return 0;
    }
    for (size_t i = 1; i < ROLLO_I_CODE_LENGTH; i++) {
        if (!bf_is_zero(&a->coefficients[i])) {
            return 0;
        }
    }

    return 1;
}

void polynomial_addition(polynomial_t *c, polynomial_t *a, polynomial_t *b) {
    for (uint8_t i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
        bf_addition(&c->coefficients[i], &a->coefficients[i], &b->coefficients[i]);
    }
}

//static inline void shift_polynomial_mod_p(polynomial_t *polynomial) {
//    // P =  X^67 + X^5 + X^2 + X + 1
//    int q;
//    bf_element_t tmp = polynomial->coefficients[ROLLO_I_CODE_LENGTH - 1];
//
//    // shift the polynomial by 1
//    for (q = ROLLO_I_CODE_LENGTH - 1; q > 0; q--) {
//        polynomial->coefficients[q] = polynomial->coefficients[q - 1];
//    }
//
//    // Do the reduction of X^67 by P, i.e.  X^67 = X^5 + X^2 + X + 1
//    polynomial->coefficients[0] = tmp;
//
//    bf_addition(&polynomial->coefficients[1], &polynomial->coefficients[1], &tmp);
//    bf_addition(&polynomial->coefficients[2], &polynomial->coefficients[2], &tmp);
//    bf_addition(&polynomial->coefficients[5], &polynomial->coefficients[5], &tmp);
//}

static inline void shift_polynomial_mod_p(polynomial_t *polynomial) {
    // P =  X^83 + X^7 + X^4 + X^2 + 1
    int q;
    bf_element_t tmp = polynomial->coefficients[ROLLO_I_CODE_LENGTH - 1];

    // shift the polynomial by 1
    for (q = ROLLO_I_CODE_LENGTH - 1; q > 0; q--) {
        polynomial->coefficients[q] = polynomial->coefficients[q - 1];
    }

    // Do the reduction of X^83 by P, i.e.  X^83 + X^7 + X^4 + X^2 + 1
    polynomial->coefficients[0] = tmp;

    bf_addition(&polynomial->coefficients[2], &polynomial->coefficients[2], &tmp);
    bf_addition(&polynomial->coefficients[4], &polynomial->coefficients[4], &tmp);
    bf_addition(&polynomial->coefficients[7], &polynomial->coefficients[7], &tmp);
}

void polynomial_multiplication_mod_p(polynomial_t *c, polynomial_t *a, polynomial_t *b) {

    int i, j;
    polynomial_t b_copy = *b;
    bf_double_element_t accumulators[ROLLO_I_CODE_LENGTH];
    bf_double_element_t tmp;

    for (i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
        bf_unreduced_multiplication(&accumulators[i], &a->coefficients[0], &b_copy.coefficients[i]);
    }
    shift_polynomial_mod_p(&b_copy);
    for (j = 1; j < ROLLO_I_CODE_LENGTH - 1; j++) {
        for (i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
            bf_unreduced_multiplication(&tmp, &a->coefficients[j], &b_copy.coefficients[i]);
            bf_double_addition(&accumulators[i], &accumulators[i], &tmp);
        }
        shift_polynomial_mod_p(&b_copy);
    }
    for (i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
        bf_unreduced_multiplication(&tmp, &a->coefficients[ROLLO_I_CODE_LENGTH - 1],
                                    &b_copy.coefficients[i]);
        bf_double_addition(&accumulators[i], &accumulators[i], &tmp);
        bf_reduction(&c->coefficients[i], &accumulators[i]);
    }
}

//using http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.46.2899&rep=rep1&type=pdf

static void rollo_i_binary_matrix_transpose(uint128_t *transpose, const uint128_t *matrix) {

    uint8_t i, j;
    uint128_t mask;

    for (i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
        mask = (uint128_t) 1 << i;
        transpose[i] = 0;
        for (j = 0; j < ROLLO_I_CODE_LENGTH; j++) {
            if (j - i >= 0) {
                transpose[i] |= ((matrix[j] & mask) << (j - i));
            } else {
                transpose[i] |= ((matrix[j] & mask) >> (i - j));
            }
        }
    }
}

static uint128_t bit_parity(uint128_t x) {

    uint128_t y = x ^(x >> 1u);
    y = y ^ (y >> 2u);
    y = y ^ (y >> 4u);
    y = y ^ (y >> 8u);
    y = y ^ (y >> 16u);
    y = y ^ (y >> 32u);
    y = y ^ (y >> 64u);

    return y & (uint128_t) 1;
}

static void rollo_i_binary_matrix_multiplication(uint128_t *c, uint128_t *a, uint128_t *b) {

    uint8_t i, j;

    uint128_t a_transposed[ROLLO_I_CODE_LENGTH];

    rollo_i_binary_matrix_transpose(a_transposed, a);

    for (i = 0; i < ROLLO_I_CODE_LENGTH; i++) {

        c[i] = 0;
        for (j = 0; j < ROLLO_I_CODE_LENGTH; j++) {
            c[i] |= bit_parity(b[i] & a_transposed[j]) << j;
        }
    }
}

void polynomial_to_the_2_to_the_n_times_m(polynomial_t *a_to_the_2_to_the_n, polynomial_t *a,
                                          const uint128_t *matrix_s_to_the_n) {

    size_t row, column;
    bf_element_t mask;

    for (row = 0; row < ROLLO_I_CODE_LENGTH; row++) {
        bf_set_to_zero(&a_to_the_2_to_the_n->coefficients[row]);
        for (column = 0; column < ROLLO_I_CODE_LENGTH; column++) {
            bf_set_mask(&mask, (matrix_s_to_the_n[column] >> row) & ((uint128_t) 1));
            bf_and(&mask, &a->coefficients[column], &mask);
            bf_addition(&a_to_the_2_to_the_n->coefficients[row], &a_to_the_2_to_the_n->coefficients[row], &mask);
        }
    }
}

void polynomial_inversion_mod_p(polynomial_t *a_inverse, polynomial_t *a) {

    polynomial_t a_to_the_r, a_to_r_minus_1, tmp;

    uint128_t matrix_s[ROLLO_I_CODE_LENGTH] = {
            (uint128_t) 0x1,
            (uint128_t) 0x683881de4774db7c  ^ (((uint128_t) 0x500f3 ) << 64u),
            (uint128_t) 0xca6c51380dc9ac1b  ^ (((uint128_t) 0x1f151 ) << 64u),
            (uint128_t) 0x756820de80373326  ^ (((uint128_t) 0x150b6 ) << 64u),
            (uint128_t) 0x98e676a2a15338cf  ^ (((uint128_t) 0x67f73 ) << 64u),
            (uint128_t) 0x572a8f874ca8cf2a  ^ (((uint128_t) 0x6ae2 ) << 64u),
            (uint128_t) 0xb141f60cc12d4fbe  ^ (((uint128_t) 0x3c9f0 ) << 64u),
            (uint128_t) 0x3617794327c765f0  ^ (((uint128_t) 0x2c28a ) << 64u),
            (uint128_t) 0x167d8f08cca1b703  ^ (((uint128_t) 0x5b3ff ) << 64u),
            (uint128_t) 0xb0abc065a53607fd  ^ (((uint128_t) 0xe838 ) << 64u),
            (uint128_t) 0x589bfd923e3329c6  ^ (((uint128_t) 0x42472 ) << 64u),
            (uint128_t) 0xcbc38f35e29df344  ^ (((uint128_t) 0x46ffd ) << 64u),
            (uint128_t) 0x2a25f4629a423b49  ^ (((uint128_t) 0x1a6ef ) << 64u),
            (uint128_t) 0x80d4538b23d8e7b8  ^ (((uint128_t) 0x406c9 ) << 64u),
            (uint128_t) 0x28f7895587c207b5  ^ (((uint128_t) 0x6124b ) << 64u),
            (uint128_t) 0xbd5f38f94a8f0c8f  ^ (((uint128_t) 0x4fe99 ) << 64u),
            (uint128_t) 0xafaa77c0b2a11ccc  ^ (((uint128_t) 0x1503d ) << 64u),
            (uint128_t) 0x68ed02a705334ade  ^ (((uint128_t) 0x7e5ee ) << 64u),
            (uint128_t) 0x189df530e660f751  ^ (((uint128_t) 0x5d839 ) << 64u),
            (uint128_t) 0xff835298d1dfd4d5  ^ (((uint128_t) 0x75307 ) << 64u),
            (uint128_t) 0x773fb909ba14bfd7  ^ (((uint128_t) 0x30aaf ) << 64u),
            (uint128_t) 0x3337886fe942b7b8  ^ (((uint128_t) 0x6ea03 ) << 64u),
            (uint128_t) 0xabb167a4f7598df3  ^ (((uint128_t) 0x6604e ) << 64u),
            (uint128_t) 0xe815033aec0ca36b  ^ (((uint128_t) 0x3281a ) << 64u),
            (uint128_t) 0x0c1e7e66cef86a61  ^ (((uint128_t) 0x657e1 ) << 64u),
            (uint128_t) 0xa19b21baf4ff525e  ^ (((uint128_t) 0x7d6e5 ) << 64u),
            (uint128_t) 0x02bb5945c59e6d81  ^ (((uint128_t) 0x74120 ) << 64u),
            (uint128_t) 0x16350260d849eec2  ^ (((uint128_t) 0x6d90 ) << 64u),
            (uint128_t) 0xeee53e29db43dde5  ^ (((uint128_t) 0x183fa ) << 64u),
            (uint128_t) 0x5cc41e2056050fb7  ^ (((uint128_t) 0x5e1e6 ) << 64u),
            (uint128_t) 0x2e5a15e49d2a2cbc  ^ (((uint128_t) 0x78aad ) << 64u),
            (uint128_t) 0x33e90dc07389e182  ^ (((uint128_t) 0x2c78f ) << 64u),
            (uint128_t) 0x18b0035fc733a270  ^ (((uint128_t) 0x6c8a0 ) << 64u),
            (uint128_t) 0x71c502f303a76d2b  ^ (((uint128_t) 0x82a3 ) << 64u),
            (uint128_t) 0x6d58e128dd54e422  ^ (((uint128_t) 0x3ada8 ) << 64u),
            (uint128_t) 0x75e6d02e0c93a905  ^ (((uint128_t) 0x4fd94 ) << 64u),
            (uint128_t) 0x088ea20cea207b40  ^ (((uint128_t) 0x43108 ) << 64u),
            (uint128_t) 0xddeac347fc86057c  ^ (((uint128_t) 0x79b82 ) << 64u),
            (uint128_t) 0x7076aeabd9495cc7  ^ (((uint128_t) 0x759b8 ) << 64u),
            (uint128_t) 0xf551c9f74b40a87c  ^ (((uint128_t) 0x165a4 ) << 64u),
            (uint128_t) 0x221f524238ab808a  ^ (((uint128_t) 0x124a7 ) << 64u),
            (uint128_t) 0xcc62b760ddb610fd  ^ (((uint128_t) 0xeeeb ) << 64u),
            (uint128_t) 0xd419687c3c43ddb4  ^ (((uint128_t) 0x55837 ) << 64u),
            (uint128_t) 0xf8dfd9c24df3865a  ^ (((uint128_t) 0x78839 ) << 64u),
            (uint128_t) 0x7ad97735906280f3  ^ (((uint128_t) 0x4a491 ) << 64u),
            (uint128_t) 0x9a7fa52852134cc9  ^ (((uint128_t) 0x376b8 ) << 64u),
            (uint128_t) 0x42174c78579db152  ^ (((uint128_t) 0x56b6c ) << 64u),
            (uint128_t) 0x073ff0822f434f47  ^ (((uint128_t) 0x34d20 ) << 64u),
            (uint128_t) 0xb8c6f0c2032df7ff  ^ (((uint128_t) 0x5adc3 ) << 64u),
            (uint128_t) 0x9c5f40ce2649e901  ^ (((uint128_t) 0x7320d ) << 64u),
            (uint128_t) 0x1cac17c30f546ea2  ^ (((uint128_t) 0x6f481 ) << 64u),
            (uint128_t) 0x0ed93a612c3b378a  ^ (((uint128_t) 0x569c2 ) << 64u),
            (uint128_t) 0xba91411a6080497f  ^ (((uint128_t) 0x21a03 ) << 64u),
            (uint128_t) 0xb6e0736943a97589  ^ (((uint128_t) 0x18c66 ) << 64u),
            (uint128_t) 0x4be003818d4cf004  ^ (((uint128_t) 0x57114 ) << 64u),
            (uint128_t) 0x65d6601bb1575c2f  ^ (((uint128_t) 0xadb6 ) << 64u),
            (uint128_t) 0xaf825d4818e9a899  ^ (((uint128_t) 0x20c1c ) << 64u),
            (uint128_t) 0xbcbdeafd4c24c94e  ^ (((uint128_t) 0x3db7e ) << 64u),
            (uint128_t) 0x78f5bd36837f437c  ^ (((uint128_t) 0x16c97 ) << 64u),
            (uint128_t) 0x7cc1558bc3bb3395  ^ (((uint128_t) 0x19636 ) << 64u),
            (uint128_t) 0x2641eb0db812c3ac  ^ (((uint128_t) 0x7d8f6 ) << 64u),
            (uint128_t) 0x9e38fe3a8febc7d1  ^ (((uint128_t) 0x1efa4 ) << 64u),
            (uint128_t) 0x98df383cb83be1bb  ^ (((uint128_t) 0x65376 ) << 64u),
            (uint128_t) 0xff27f071a6ade069  ^ (((uint128_t) 0x3705c ) << 64u),
            (uint128_t) 0xb2951309b424015c  ^ (((uint128_t) 0x17f9 ) << 64u),
            (uint128_t) 0x53e8602e385dae76  ^ (((uint128_t) 0x6a4f5 ) << 64u),
            (uint128_t) 0x62dc16223769a445  ^ (((uint128_t) 0xfd4b ) << 64u),
            (uint128_t) 0xdd8982a1e715d2fc  ^ (((uint128_t) 0x60b4 ) << 64u),
            (uint128_t) 0x17787401c81dc993  ^ (((uint128_t) 0x5c354 ) << 64u),
            (uint128_t) 0xc0cac207d0be06e2  ^ (((uint128_t) 0x515eb ) << 64u),
            (uint128_t) 0x1bdb32182da400d0  ^ (((uint128_t) 0x7dbe0 ) << 64u),
            (uint128_t) 0x51566dd42dfdafe8  ^ (((uint128_t) 0x10df7 ) << 64u),
            (uint128_t) 0xf0ec00acb88962e1  ^ (((uint128_t) 0x65862 ) << 64u),
            (uint128_t) 0xbc00ed9ebdfaea9c  ^ (((uint128_t) 0x4484e ) << 64u),
            (uint128_t) 0xdd1f55216c760226  ^ (((uint128_t) 0x38e41 ) << 64u),
            (uint128_t) 0x81c68a115acafadb  ^ (((uint128_t) 0x27e02 ) << 64u),
            (uint128_t) 0x05c8626062142261  ^ (((uint128_t) 0x7585b ) << 64u),
            (uint128_t) 0x6c1f84e8d9281a18  ^ (((uint128_t) 0x3f729 ) << 64u),
            (uint128_t) 0x53eacb19dcf732d2  ^ (((uint128_t) 0x294a9 ) << 64u),
            (uint128_t) 0xf0b5134b07316378  ^ (((uint128_t) 0xb003 ) << 64u),
            (uint128_t) 0xe6b22e8fd77b59ce  ^ (((uint128_t) 0x6fbaa ) << 64u),
            (uint128_t) 0x57df6a8df8786219  ^ (((uint128_t) 0x3d0f3 ) << 64u),
            (uint128_t) 0x1d24e29761ce245b  ^ (((uint128_t) 0xd94d ) << 64u)
    };
    int i;

    uint128_t matrix_s_to_the_n[ROLLO_I_CODE_LENGTH];

    for (i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
        matrix_s_to_the_n[i] = matrix_s[i];
    }

    polynomial_to_the_2_to_the_n_times_m(&a_to_r_minus_1, a, matrix_s_to_the_n);
    int n;

//    FILE *f = fopen("/Users/ema/Desktop/tmp_rollo.h", "w");
//    if (f == NULL)
//    {
//        printf("Error opening file!\n");
//        exit(1);
//    }

//    fprintf(f, "matrix_s_to_the_n_list[ROLLO_I_CODE_LENGTH-2][ROLLO_I_CODE_LENGTH] = {\n", n);
    for (n = 2; n < ROLLO_I_CODE_LENGTH; n++) {
        rollo_i_binary_matrix_multiplication(matrix_s_to_the_n, matrix_s_to_the_n, matrix_s);

//        fprintf(f, "{ // matrix %d\n", n);
//        for (i = 0; i < ROLLO_I_CODE_LENGTH; i++) {
//            fprintf(f, "(uint128_t) 0x%016llx  ^ (((uint128_t) 0x%05llx ) << 64u)",
//                    (uint64_t) (matrix_s_to_the_n[i] & 0xFFFFFFFFFFFFFFFFu), (uint64_t) (matrix_s_to_the_n[i] >> 64u));
//            if (i == ROLLO_I_CODE_LENGTH-1)
//                fprintf(f, "  // i = %d\n", i);
//            else
//                fprintf(f, ",  // i = %d\n", i);
//        }
//        fprintf(f, "},\n");

        polynomial_to_the_2_to_the_n_times_m(&tmp, a, matrix_s_to_the_n);
        polynomial_multiplication_mod_p(&a_to_r_minus_1, &tmp, &a_to_r_minus_1);
    }
//    fprintf(f, "};\n\n");
    polynomial_multiplication_mod_p(&a_to_the_r, a, &a_to_r_minus_1);

    bf_inversion(&a_to_the_r.coefficients[0], &a_to_the_r.coefficients[0]);
    for (n = 0; n < ROLLO_I_CODE_LENGTH; n++) {
        bf_multiplication(&a_inverse->coefficients[n], &a_to_the_r.coefficients[0], &a_to_r_minus_1.coefficients[n]);
    }
}

void polynomial_inversion_mod_p_with_precomputes_matrices(polynomial_t *a_inverse, polynomial_t *a) {

    polynomial_t a_to_the_r, a_to_r_minus_1, tmp;

    uint128_t matrix_s[ROLLO_I_CODE_LENGTH] = {
            (uint128_t) 0x1,
            (uint128_t) 0x683881de4774db7c  ^ (((uint128_t) 0x500f3 ) << 64u),
            (uint128_t) 0xca6c51380dc9ac1b  ^ (((uint128_t) 0x1f151 ) << 64u),
            (uint128_t) 0x756820de80373326  ^ (((uint128_t) 0x150b6 ) << 64u),
            (uint128_t) 0x98e676a2a15338cf  ^ (((uint128_t) 0x67f73 ) << 64u),
            (uint128_t) 0x572a8f874ca8cf2a  ^ (((uint128_t) 0x6ae2 ) << 64u),
            (uint128_t) 0xb141f60cc12d4fbe  ^ (((uint128_t) 0x3c9f0 ) << 64u),
            (uint128_t) 0x3617794327c765f0  ^ (((uint128_t) 0x2c28a ) << 64u),
            (uint128_t) 0x167d8f08cca1b703  ^ (((uint128_t) 0x5b3ff ) << 64u),
            (uint128_t) 0xb0abc065a53607fd  ^ (((uint128_t) 0xe838 ) << 64u),
            (uint128_t) 0x589bfd923e3329c6  ^ (((uint128_t) 0x42472 ) << 64u),
            (uint128_t) 0xcbc38f35e29df344  ^ (((uint128_t) 0x46ffd ) << 64u),
            (uint128_t) 0x2a25f4629a423b49  ^ (((uint128_t) 0x1a6ef ) << 64u),
            (uint128_t) 0x80d4538b23d8e7b8  ^ (((uint128_t) 0x406c9 ) << 64u),
            (uint128_t) 0x28f7895587c207b5  ^ (((uint128_t) 0x6124b ) << 64u),
            (uint128_t) 0xbd5f38f94a8f0c8f  ^ (((uint128_t) 0x4fe99 ) << 64u),
            (uint128_t) 0xafaa77c0b2a11ccc  ^ (((uint128_t) 0x1503d ) << 64u),
            (uint128_t) 0x68ed02a705334ade  ^ (((uint128_t) 0x7e5ee ) << 64u),
            (uint128_t) 0x189df530e660f751  ^ (((uint128_t) 0x5d839 ) << 64u),
            (uint128_t) 0xff835298d1dfd4d5  ^ (((uint128_t) 0x75307 ) << 64u),
            (uint128_t) 0x773fb909ba14bfd7  ^ (((uint128_t) 0x30aaf ) << 64u),
            (uint128_t) 0x3337886fe942b7b8  ^ (((uint128_t) 0x6ea03 ) << 64u),
            (uint128_t) 0xabb167a4f7598df3  ^ (((uint128_t) 0x6604e ) << 64u),
            (uint128_t) 0xe815033aec0ca36b  ^ (((uint128_t) 0x3281a ) << 64u),
            (uint128_t) 0x0c1e7e66cef86a61  ^ (((uint128_t) 0x657e1 ) << 64u),
            (uint128_t) 0xa19b21baf4ff525e  ^ (((uint128_t) 0x7d6e5 ) << 64u),
            (uint128_t) 0x02bb5945c59e6d81  ^ (((uint128_t) 0x74120 ) << 64u),
            (uint128_t) 0x16350260d849eec2  ^ (((uint128_t) 0x6d90 ) << 64u),
            (uint128_t) 0xeee53e29db43dde5  ^ (((uint128_t) 0x183fa ) << 64u),
            (uint128_t) 0x5cc41e2056050fb7  ^ (((uint128_t) 0x5e1e6 ) << 64u),
            (uint128_t) 0x2e5a15e49d2a2cbc  ^ (((uint128_t) 0x78aad ) << 64u),
            (uint128_t) 0x33e90dc07389e182  ^ (((uint128_t) 0x2c78f ) << 64u),
            (uint128_t) 0x18b0035fc733a270  ^ (((uint128_t) 0x6c8a0 ) << 64u),
            (uint128_t) 0x71c502f303a76d2b  ^ (((uint128_t) 0x82a3 ) << 64u),
            (uint128_t) 0x6d58e128dd54e422  ^ (((uint128_t) 0x3ada8 ) << 64u),
            (uint128_t) 0x75e6d02e0c93a905  ^ (((uint128_t) 0x4fd94 ) << 64u),
            (uint128_t) 0x088ea20cea207b40  ^ (((uint128_t) 0x43108 ) << 64u),
            (uint128_t) 0xddeac347fc86057c  ^ (((uint128_t) 0x79b82 ) << 64u),
            (uint128_t) 0x7076aeabd9495cc7  ^ (((uint128_t) 0x759b8 ) << 64u),
            (uint128_t) 0xf551c9f74b40a87c  ^ (((uint128_t) 0x165a4 ) << 64u),
            (uint128_t) 0x221f524238ab808a  ^ (((uint128_t) 0x124a7 ) << 64u),
            (uint128_t) 0xcc62b760ddb610fd  ^ (((uint128_t) 0xeeeb ) << 64u),
            (uint128_t) 0xd419687c3c43ddb4  ^ (((uint128_t) 0x55837 ) << 64u),
            (uint128_t) 0xf8dfd9c24df3865a  ^ (((uint128_t) 0x78839 ) << 64u),
            (uint128_t) 0x7ad97735906280f3  ^ (((uint128_t) 0x4a491 ) << 64u),
            (uint128_t) 0x9a7fa52852134cc9  ^ (((uint128_t) 0x376b8 ) << 64u),
            (uint128_t) 0x42174c78579db152  ^ (((uint128_t) 0x56b6c ) << 64u),
            (uint128_t) 0x073ff0822f434f47  ^ (((uint128_t) 0x34d20 ) << 64u),
            (uint128_t) 0xb8c6f0c2032df7ff  ^ (((uint128_t) 0x5adc3 ) << 64u),
            (uint128_t) 0x9c5f40ce2649e901  ^ (((uint128_t) 0x7320d ) << 64u),
            (uint128_t) 0x1cac17c30f546ea2  ^ (((uint128_t) 0x6f481 ) << 64u),
            (uint128_t) 0x0ed93a612c3b378a  ^ (((uint128_t) 0x569c2 ) << 64u),
            (uint128_t) 0xba91411a6080497f  ^ (((uint128_t) 0x21a03 ) << 64u),
            (uint128_t) 0xb6e0736943a97589  ^ (((uint128_t) 0x18c66 ) << 64u),
            (uint128_t) 0x4be003818d4cf004  ^ (((uint128_t) 0x57114 ) << 64u),
            (uint128_t) 0x65d6601bb1575c2f  ^ (((uint128_t) 0xadb6 ) << 64u),
            (uint128_t) 0xaf825d4818e9a899  ^ (((uint128_t) 0x20c1c ) << 64u),
            (uint128_t) 0xbcbdeafd4c24c94e  ^ (((uint128_t) 0x3db7e ) << 64u),
            (uint128_t) 0x78f5bd36837f437c  ^ (((uint128_t) 0x16c97 ) << 64u),
            (uint128_t) 0x7cc1558bc3bb3395  ^ (((uint128_t) 0x19636 ) << 64u),
            (uint128_t) 0x2641eb0db812c3ac  ^ (((uint128_t) 0x7d8f6 ) << 64u),
            (uint128_t) 0x9e38fe3a8febc7d1  ^ (((uint128_t) 0x1efa4 ) << 64u),
            (uint128_t) 0x98df383cb83be1bb  ^ (((uint128_t) 0x65376 ) << 64u),
            (uint128_t) 0xff27f071a6ade069  ^ (((uint128_t) 0x3705c ) << 64u),
            (uint128_t) 0xb2951309b424015c  ^ (((uint128_t) 0x17f9 ) << 64u),
            (uint128_t) 0x53e8602e385dae76  ^ (((uint128_t) 0x6a4f5 ) << 64u),
            (uint128_t) 0x62dc16223769a445  ^ (((uint128_t) 0xfd4b ) << 64u),
            (uint128_t) 0xdd8982a1e715d2fc  ^ (((uint128_t) 0x60b4 ) << 64u),
            (uint128_t) 0x17787401c81dc993  ^ (((uint128_t) 0x5c354 ) << 64u),
            (uint128_t) 0xc0cac207d0be06e2  ^ (((uint128_t) 0x515eb ) << 64u),
            (uint128_t) 0x1bdb32182da400d0  ^ (((uint128_t) 0x7dbe0 ) << 64u),
            (uint128_t) 0x51566dd42dfdafe8  ^ (((uint128_t) 0x10df7 ) << 64u),
            (uint128_t) 0xf0ec00acb88962e1  ^ (((uint128_t) 0x65862 ) << 64u),
            (uint128_t) 0xbc00ed9ebdfaea9c  ^ (((uint128_t) 0x4484e ) << 64u),
            (uint128_t) 0xdd1f55216c760226  ^ (((uint128_t) 0x38e41 ) << 64u),
            (uint128_t) 0x81c68a115acafadb  ^ (((uint128_t) 0x27e02 ) << 64u),
            (uint128_t) 0x05c8626062142261  ^ (((uint128_t) 0x7585b ) << 64u),
            (uint128_t) 0x6c1f84e8d9281a18  ^ (((uint128_t) 0x3f729 ) << 64u),
            (uint128_t) 0x53eacb19dcf732d2  ^ (((uint128_t) 0x294a9 ) << 64u),
            (uint128_t) 0xf0b5134b07316378  ^ (((uint128_t) 0xb003 ) << 64u),
            (uint128_t) 0xe6b22e8fd77b59ce  ^ (((uint128_t) 0x6fbaa ) << 64u),
            (uint128_t) 0x57df6a8df8786219  ^ (((uint128_t) 0x3d0f3 ) << 64u),
            (uint128_t) 0x1d24e29761ce245b  ^ (((uint128_t) 0xd94d ) << 64u)
    };
    int i;

    polynomial_to_the_2_to_the_n_times_m(&a_to_r_minus_1, a, matrix_s);
    int n;
    for (n = 2; n < ROLLO_I_CODE_LENGTH; n++) {
        polynomial_to_the_2_to_the_n_times_m(&tmp, a, matrix_s_to_the_n_list[n-2]);
        polynomial_multiplication_mod_p(&a_to_r_minus_1, &tmp, &a_to_r_minus_1);
    }
    polynomial_multiplication_mod_p(&a_to_the_r, a, &a_to_r_minus_1);

    bf_inversion(&a_to_the_r.coefficients[0], &a_to_the_r.coefficients[0]);
    for (n = 0; n < ROLLO_I_CODE_LENGTH; n++) {
        bf_multiplication(&a_inverse->coefficients[n], &a_to_the_r.coefficients[0], &a_to_r_minus_1.coefficients[n]);
    }
}

int generate_basis_and_sample_two_polynomials_from_basis(
        polynomial_t *p1, polynomial_t *p2,
        bf_element_t *basis, uint16_t basis_length,
        AES_XOF_struct *prng) {

    generate_two_list_of_vectors_with_given_rank_and_their_basis(
            (bf_element_t *) p1, ROLLO_I_CODE_LENGTH,
            (bf_element_t *) p2, ROLLO_I_CODE_LENGTH,
            basis_length, basis,
            prng);

    return EXIT_SUCCESS;
}

// TODO: vector space length could be removed
//       and also prng if using random bytes rather than seedexpander
int sample_two_polynomials_from_vector_space(polynomial_t *p1, polynomial_t *p2,
                                             bf_element_t *list, uint16_t list_length, uint8_t mask, AES_XOF_struct *prng) {
    uint8_t random_combination = 0;
    for (uint16_t i = 0; i < ROLLO_I_CODE_LENGTH; i++) {

        randombytes(&random_combination, 1);
//        seedexpander(prng, &random_combination, 1);
        bf_copy(&p1->coefficients[i],
                &list[random_combination & (mask)]);

        randombytes(&random_combination, 1);
//        seedexpander(prng, &random_combination, 1);
        bf_copy(&p2->coefficients[i],
                &list[random_combination & (mask)]);
    }
    return EXIT_SUCCESS;
}