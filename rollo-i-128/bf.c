
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "bf.h"
#include "rollo-i-128_parameters.h"
#include "avx2_conversion/avx2_conversion.h"

#ifdef BF_UINT64_UINT8_T

static int x_hexdigit(int ch) {
    if (ch >= '0' && ch <= '9')
        return ch - '0';
    if (ch >= 'A' && ch <= 'F')
        return ch - 'A' + 10;
    if (ch >= 'a' && ch <= 'f')
        return ch - 'a' + 10;
    return -1;
}

uint8_t bf_set_from_hex_string(bf_element_t *result, const char *h) {

    size_t l = strlen(h);
    if (l > ROLLO_I_BF_ELEMENT_BYTE_SIZE * 2) {
        return 0;
    }

    result->low = 0;
    result->high = 0;

    if (l > 16) {
        result->high = x_hexdigit(h[0]);
        for (size_t i = 0; i < l-1; ++i) {
            result->low = result->low + (((uint64_t) x_hexdigit(h[l - i - 1])) << (4 * i));
        }
    } else {
        for (size_t i = 0; i < l; ++i) {
            result->low = result->low + (((uint64_t) x_hexdigit(h[l - i - 1])) << (4 * i));
        }
    }

    result->high = result->high & ROLLO_I_BF_MASK_HIGH;

    return 1;
}


// the input hex string must have length 34
uint8_t bf_double_set_from_string(bf_double_element_t *result, const char *h) {
    size_t l = strlen(h);

    if (l != 34) {
        printf("\n Error! The input string should contain 34 hex chatacters!");
        return 0;
    }

    result->value[3] = 0;

    result->value[0] = 0;
    for (size_t i = 0; i < 16; ++i) {
        result->value[0] = result->value[0] + (((uint128_t) x_hexdigit(h[l - i - 1])) << (4 * i));
    }

    result->value[1] = 0;
    for (size_t i = 0; i < 16; ++i) {
        result->value[1] = result->value[1] + (((uint128_t) x_hexdigit(h[l - i - 1 - 16])) << (4 * i));
    }

    result->value[2] = 0;
    for (size_t i = 0; i < 16; ++i) {
        result->value[2] = result->value[2] + (((uint128_t) x_hexdigit(h[l - i - 1 - 32])) << (4 * i));
    }

    result->value[2] = result->value[2] & 0x3Fu;

    return 1;
}

void bf_set_to_zero(bf_element_t *c) {
    c->low = 0;
    c->high = 0;
}

void bf_set_to_ones(bf_element_t *c) {
    c->low = ROLLO_I_BF_MASK_LOW;
    c->high = ROLLO_I_BF_MASK_HIGH;
}

void bf_print(const bf_element_t *a) {
//    printf("0x%01llx %016llx", (uint64_t) (a->value >> 64u), (uint64_t) a->value);
    printf("%0x %016llx", (uint8_t) (a->high), a->low);
}

void bf_double_element_print(const bf_double_element_t *a) {
    printf("%016llx %016llx %016llx %016llx", a->value[3], a->value[2], a->value[1], a->value[0] );
}

void bf_print_binary(const bf_element_t *a) {
    for (int i = 5; i < 8; i++) {
        printf("%d", (int) ((uint64_t) ((uint64_t) a->high >> (8u - 1u - i)) & 0x1u));
    }
    for (int i = 0; i < 64; i++) {
        printf("%d", (int) ((a->low >> (64u - 1u - i)) & 0x1u));
    }
//    for (int i = 0; i < ROLLO_I_FINITE_FIELD_DEGREE; i++) {
//        printf("%d", (int) ((a->low >> (ROLLO_I_FINITE_FIELD_DEGREE - 1 - i)) & 0x1u));
//    }
}

void bf_print_matrix_bin(bf_element_t *mat, size_t NROWS) {
    for (size_t i = 0; i < NROWS; ++i) {
        bf_print_binary(&mat[i]);
        printf("\n");
    }
}

inline void bf_copy(bf_element_t *a, bf_element_t *b) {
    a->low = b->low;
    a->high = b->high;
//    *a = *b;
}

//inline void bf_swap(bf_element_t *a, bf_element_t *b) {
//    bf_element_t tmp;
//    tmp.value = b->value;
//    b->value = a->value;
//    a->value = tmp.value;
//}

void bf_get_random_element(AES_XOF_struct *prng, bf_element_t *c) {
//    seedexpander(prng, (uint8_t *) &c->value, ROLLO_I_BF_ELEMENT_BYTE_SIZE);
//    seedexpander(prng, (uint8_t *) &c->low, 8);
//    seedexpander(prng, (uint8_t *) &c->high, 1);

    //    randombytes((uint8_t *) &c->value, ROLLO_I_BF_ELEMENT_BYTE_SIZE);

//    randombytes((uint8_t *) &c->low, 8);
//    randombytes((uint8_t *) &c->high, 1);

    randombytes((uint8_t *) &c->low, 9);
    c->high &= ROLLO_I_BF_MASK_HIGH;
}

void bf_get_list_of_random_elements(AES_XOF_struct *prng, bf_element_t *c, uint16_t number_of_elements){
    //    randombytes((uint8_t *) support_basis, rank*ROLLO_I_BF_ELEMENT_BYTE_SIZE);
    uint8_t verb = 0;

    if (verb) {
        printf(("\nBEFORE RANDOMBYTES\n"));
        for (uint16_t i = 0; i < number_of_elements; ++i) {
            printf("i: %3d = ", i);
            bf_print_binary(&c[i]);
            printf("\n");
        }
    }
    // TODO: for some reason this does not work...
//    randombytes((uint8_t *) c, ROLLO_I_BF_ELEMENT_BYTE_SIZE*number_of_elements);

    for (uint16_t i = 0; i < number_of_elements; i++) {
//        bf_get_random_element(prng, &c[i]);

//        randombytes((uint8_t *) &c[i], 9);

        randombytes((uint8_t *) &c[i].low, 8);
        randombytes((uint8_t *) &c[i].high, 1);

        c[i].high = c[i].high & ROLLO_I_BF_MASK_HIGH;
    }

    if (verb) {
        printf(("\nAFTER RANDOMBYTES\n"));
        for (uint16_t i = 0; i < number_of_elements; ++i) {
            printf("i: %3d = ", i);
            bf_print_binary(&c[i]);
            printf("\n");
        }
    }
}

// returns 1 if a is 0, 0 otherwise
int bf_is_zero(bf_element_t *a) {
    uint8_t result = 0;
    for (uint8_t i = 0; i < 64; ++i)
        result |= (uint8_t) (a->low >> i) & 0x1u;
    for (uint8_t i = 0; i < 3; ++i)
        result |= (uint8_t) (a->high >> i) & 0x1u;
    return (uint8_t) (result ^ 0x1u); // flip bit
}

// returns 1 if a is 1, 0 otherwise
int bf_is_one(bf_element_t *a) {
    uint8_t otherbits = 0;
    uint8_t firstbit = (uint8_t) (a->low & 0x1u);
    for (uint8_t i = 1; i < 64; ++i)
        otherbits |= (uint8_t) (a->low >> i) & 0x1u;
    for (uint8_t i = 0; i < 3; ++i)
        otherbits |= (uint8_t) (a->high >> i) & 0x1u;
    return (uint8_t) (otherbits ^ 0x1u) & firstbit; // (0,1) is ok, (0,0),(1,0),(1,1) are not
}

// return 1 if a = b, 0 otherwise
int bf_is_equal(bf_element_t *a, bf_element_t *b) {
    uint8_t result = 0;
    for (uint8_t i = 0; i < 64; ++i)
        result |= (uint8_t) ((a->low >> i) ^ (b->low >> i)) & 0x1u;
    for (uint8_t i = 0; i < 3; ++i)
        result |= (uint8_t) (((uint64_t) a->high >> i) ^ ((uint64_t) b->high >> i)) & 0x1u;
    return (uint8_t) (result ^ 0x1u); // flip bit
}

// return 1 if a = b, 0 otherwise
int bf_double_is_equal(bf_double_element_t *a, bf_double_element_t *b) {
    uint8_t result = 0;
    for (uint8_t j = 0; j < 4; ++j) {
        for (uint8_t i = 0; i < 128; ++i)
            result |= (uint8_t) ((a->value[j] >> i) ^ (b->value[j] >> i)) & 0x1u;
    }
    return (uint8_t) (result ^ 0x1u); // flip bit
}


// returns 1 if a > b, 0 if a <= b
inline int bf_greater_than(bf_element_t *a, bf_element_t *b) {
    uint64_t aH, aL, bH, bL, u, c, t;

    uint64_t a0,a1,a2,a3,b0,b1,b2,b3;

    c = 0;
    t = 0;

    aH = a->high;
    aL = a->low;
    bH = b->high;
    bL = b->low;

    a3 = (uint64_t) (aH >> 32u);
    a2 = aH & 0xFFFFFFFFu;
    a1 = (uint64_t) (aL >> 32u);
    a0 = aL & 0xFFFFFFFFu;

    b3 = (uint64_t) (bH >> 32u);
    b2 = bH & 0xFFFFFFFFu;
    b1 = (uint64_t) (bL >> 32u);
    b0 = bL & 0xFFFFFFFFu;

    u = b0 - a0 - c;
    c = u >> 63u;
    t |= u;
    u = b1 - a1 - c;
    c = u >> 63u;
    t |= u;
    u = b2 - a2 - c;
    c = u >> 63u;
    t |= u;
    u = b3 - a3 - c;
    c = u >> 63u;
    t |= u;

//    u = bL - aL - c;
//    c = u >> 63u;
//    t |= u;
//    u = bH - aH - c;
//    c = u >> 63u;
//    t |= u;

    t = (t | -t) >> 63u;

    return t * c;
}

// used only for testing
inline int bf_greater_than_not_constant_time(bf_element_t *a, bf_element_t *b) {
    if (a->high > b->high || (a->high == b->high && a->low > b->low)) 
        return 1;
    else
        return 0;
}

// returns 1 if a > b, 0 if a = b, -1 if a < b
inline int bf_compare(bf_element_t *a, bf_element_t *b) {
    uint64_t aH, aL, bH, bL, u, c, t;

    c = 0;
    t = 0;

    aH = a->high;
    aL = a->low;
    bH = b->high;
    bL = b->low;

    u = aL - bL - c;
    c = (u >> 63u) & 1u;
    t |= u;
    u = aH - bH - c;
    c = (u >> 63u) & 1u;
    t |= u;

    c = c * -2 + 1;
    t = ((t | -t) >> 63u) & 1u;

    return t * c;
}

// constant time branching

// set the mask to all ones if bit is 1, to all zeroes otherwise
int bf_set_mask(bf_element_t *mask, uint8_t bit) {
    mask->low = -((uint64_t) bit);
//    mask->high = (uint64_t) -((uint8_t) bit) & ROLLO_I_BF_MASK_HIGH;
    mask->high = -((uint8_t) bit);

    mask->high &= ROLLO_I_BF_MASK_HIGH;

    return EXIT_SUCCESS;
}

// if j-th bit of a is 1 then mask is set to all ones, otherwise to all zeros
int bf_compute_mask(bf_element_t *mask, bf_element_t *a, uint8_t bit_position) {

    uint8_t pos = bit_position / 64u;
    uint8_t bit =  ((uint64_t) (pos * (a->high >> (bit_position - 64u))) ^ (1u-pos)*(a->low >> bit_position)) & 0x1u;
    mask->low = -((uint64_t) bit);
    mask->high = (uint64_t) -((uint8_t) bit) & ROLLO_I_BF_MASK_HIGH;

    return EXIT_SUCCESS;
}

// if j-th bit of a is 0 then mask is set to all ones, otherwise to all zeros
int bf_compute_mask_inverse(bf_element_t *mask, bf_element_t *a, uint8_t bit_position) {
//    uint128_t tmp = ((uint128_t) a->high << 64u) & (a->low);
//    uint16_t bit = 0x1u - (uint128_t) (tmp >> bit_position) & 0x1u;
//    mask->low = -((uint64_t) bit);
//    mask->high = -((uint8_t) bit);

    uint8_t pos = bit_position / 64u;
    uint8_t bit =  0x1u - (((uint64_t) (pos * (a->high >> (bit_position - 64u))) ^ (1u-pos)*(a->low >> bit_position)) & 0x1u);
    mask->low = -((uint64_t) bit);
    mask->high = (uint64_t) -((uint8_t) bit) & ROLLO_I_BF_MASK_HIGH;

    return EXIT_SUCCESS;
}

void bf_addition(bf_element_t *c, bf_element_t *a, bf_element_t *b) {
    c->high = a->high ^ b->high;
    c->low = a->low ^ b->low;
}

void bf_double_addition(bf_double_element_t *c, bf_double_element_t *a, bf_double_element_t *b) {
    c->value = _mm256_xor_si256(a->value, b->value);
}

void bf_unreduced_multiplication_karatsuba(bf_double_element_t *c, const bf_element_t *a,
                                           const bf_element_t *b) {
    __m128i a_low_times_b_low, a_low_times_b_high_xor_a_high_times_b_low, a_high_times_b_high;
    __m256i a_low_times_b_high_xor_a_high_times_b_low256;

    // 2 128x128 multiplications: AlBl, AlBh
    a_low_times_b_low = _mm_clmulepi64_si128((uint128_t) a->low, (uint128_t) b->low, 0x00);
    // karatsuba to compute AhBl, AhBh (1 128x128 multiplication + 2 128x128 xor + 2 256x256 xor)
    a_high_times_b_high = _mm_clmulepi64_si128((uint128_t) a->high, (uint128_t) b->high, 0x11);
    a_low_times_b_high_xor_a_high_times_b_low = _mm_clmulepi64_si128(
                                                        (__m128i) ((uint128_t) a->low ^
                                                                   a->high),
                                                        (__m128i) ((uint128_t) b->low ^
                                                                   b->high),
                                                        0x00) ^
                                                a_low_times_b_low ^ a_high_times_b_high;

    // convert AlBl from 128 to 256
    __m128i zero = _mm_setzero_si128();

    // convert AhBl ^ AlBh from 128 to 256
    a_low_times_b_high_xor_a_high_times_b_low256 = _mm256_set_m128i(zero, a_low_times_b_high_xor_a_high_times_b_low);
    a_low_times_b_high_xor_a_high_times_b_low256 = _mm256_permute4x64_epi64(
            a_low_times_b_high_xor_a_high_times_b_low256, 0xD2); //0x93 works too

    // concatenate AlBl || AhBh
    c->value = _mm256_set_m128i(a_high_times_b_high, a_low_times_b_low);
    // [AlBl || AhBh] ^ [0 || (AhBl ^ AlBh) || 0]
    c->value = _mm256_xor_si256(c->value, a_low_times_b_high_xor_a_high_times_b_low256);
}

void bf_unreduced_multiplication(bf_double_element_t *c, const bf_element_t *a,
                                 const bf_element_t *b) {

    __m128i a_low_times_b_low, a_low_times_b_high, a_high_times_b_low, a_high_times_b_high, a_low_times_b_high_xor_a_high_times_b_low;
    __m256i a_low_times_b_high_xor_a_high_times_b_low256;

    // 4 128x128 multiplications: AlBl, AlBh, AhBl, AhBh
    a_low_times_b_low = _mm_clmulepi64_si128((uint128_t) b->low,  (uint128_t) a->low, 0x00);
    a_low_times_b_high = _mm_clmulepi64_si128((uint128_t) b->high, (uint128_t) a->low, 0x00);
    a_high_times_b_low = _mm_clmulepi64_si128((uint128_t) b->low,  (uint128_t) a->high, 0x00);
    a_high_times_b_high = _mm_clmulepi64_si128((uint128_t) b->high, (uint128_t) a->high, 0x00);

    // 1 128x128 xor: AhBl ^ AlBh
    a_low_times_b_high_xor_a_high_times_b_low = _mm_xor_si128(a_low_times_b_high, a_high_times_b_low);

    // convert AlBl from 128 to 256
    __m128i zero = _mm_setzero_si128();

    // convert AhBl ^ AlBh from 128 to 256
    a_low_times_b_high_xor_a_high_times_b_low256 = _mm256_set_m128i(zero, a_low_times_b_high_xor_a_high_times_b_low);
    a_low_times_b_high_xor_a_high_times_b_low256 = _mm256_permute4x64_epi64(
            a_low_times_b_high_xor_a_high_times_b_low256, 0xD2); //0x93 works too

    // concatenate AlBl || AhBh
    c->value = _mm256_set_m128i(a_high_times_b_high, a_low_times_b_low);
    // [AlBl || AhBh] ^ [0 || (AhBl ^ AlBh) || 0]
    c->value = _mm256_xor_si256(c->value, a_low_times_b_high_xor_a_high_times_b_low256);
}

void bf_reduction(bf_element_t *reduced_a, bf_double_element_t *a) {

    // x^67 + x^5 + x^2 + x + 1

    uint64_t tmp = ((uint64_t) a->value[1] >> 62u) ^ ((uint64_t) a->value[2] << 2u);
    uint64_t  tmp1 = (uint64_t) a->value[1] ^ tmp ^ (tmp >> 3u) ^ (tmp >> 4u) ^ (tmp >> 5u);

    tmp = (tmp1 >> 3u) ^ ((uint64_t) a->value[2] << 61u);
    reduced_a->low = (uint64_t) a->value[0] ^ tmp ^ (tmp << 1u) ^ (tmp << 2u) ^ (tmp << 5u);

    reduced_a->high = tmp1 & ROLLO_I_BF_MASK_HIGH;
}

void bf_reduction_uint64(uint64_t *o, const uint64_t *e) {
    uint64_t tmp = (e[1] >> 62u) ^ (e[2] << 2u);
    o[1] = e[1] ^ tmp ^ (tmp >> 3u) ^ (tmp >> 4u) ^ (tmp >> 5u);

    tmp = (o[1] >> 3u) ^ (e[2] << 61u);
    o[0] = e[0] ^ tmp ^ (tmp << 1u) ^ (tmp << 2u) ^ (tmp << 5u);

    o[1] &= 0x0000000000000007u;
}

void bf_multiplication(bf_element_t *c, bf_element_t *a, bf_element_t *b) {

    bf_double_element_t axb;
    bf_unreduced_multiplication(&axb, a, b);
    bf_reduction(c, &axb);
}


void bf_unreduced_square(bf_double_element_t *c, bf_element_t *a) {
    static const __m128i squaring_left_mask = {0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F};
    static const __m128i squaring_lookup_table = {0x1514111005040100, 0x5554515045444140};
    __m128i a_left, a_right, tmp0, tmp1;
//    a_left = _mm_and_si128((__m128i) a->value, squaring_left_mask);
//    a_right = _mm_srli_epi64((__m128i) a->value, 4);
    a_left = _mm_and_si128((__m128i) {a->low, a->high}, squaring_left_mask);
    a_right = _mm_srli_epi64((__m128i) {a->low, a->high}, 4);
    a_right = _mm_and_si128(a_right, squaring_left_mask);

    a_left = _mm_shuffle_epi8(squaring_lookup_table, a_left);
    a_right = _mm_shuffle_epi8(squaring_lookup_table, a_right);

    tmp0 = _mm_unpacklo_epi8(a_left, a_right);
    tmp1 = _mm_unpackhi_epi8(a_left, a_right);
    c->value = _mm256_loadu2_m128i(&tmp0, &tmp1);
}

void bf_square(bf_element_t *c, bf_element_t *a) {

    bf_double_element_t c_tmp;
    bf_unreduced_square(&c_tmp, a);
    bf_reduction(c, &c_tmp);
}

void bf_inversion(bf_element_t *x_inverse, bf_element_t *x) {

    bf_element_t r0 = *x;
    bf_element_t r1, r3;

    bf_square(&r1, &r0); //r1 = r0^2 = a^2
    bf_multiplication(&r0, &r1, &r0); //r0 = r1 * r0 = a^3
    bf_square(&r1, &r0); //r1 = r0^2 = a^6
    bf_multiplication(&r3, &r1, x); //r3 = r0 * r1 = a^7

    //r1 = (r3) ^ (2 ^ 3) = a^(2^6 - 2^3)
    bf_square(&r1, &r3);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);

    bf_multiplication(&r0, &r1, &r3); //r0 = r1 * r3 = a^(2^6 - 1)

    //r1 = (r0) ^ (2 ^ 6) = a^((2^6 - 1) * 2^6)
    bf_square(&r1, &r0);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);

    bf_multiplication(&r0, &r1, &r0); //r0 = r0 * r1 = a^(2^12 - 1)

    //r1 = (r0) ^ (2 ^ 3) = a^(2^15 - 2^3)
    bf_square(&r1, &r0);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);

    bf_multiplication(&r0, &r1, &r3); //r0 = r0 * r1 = a^(2^15 - 1)

    //r1 = (r0) ^ (2 ^ 15) = a^((2^15 - 1) * 2^15)
    bf_square(&r1, &r0);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);

    bf_multiplication(&r0, &r1, &r0); //r0 = r0 * r1 = a^(2^30 - 1)

    //r1 = (r0) ^ (2 ^ 3) = a^(2^33 - 2^3)
    bf_square(&r1, &r0);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);

    bf_multiplication(&r0, &r1, &r3); //r0 = r0 * r1 = a^(2^33 - 1)

    //r1 = (r0) ^ (2 ^ 33) = a^((2^33 - 1) * 2^33)
    bf_square(&r1, &r0);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);
    bf_square(&r1, &r1);

    bf_multiplication(&r0, &r1, &r0); //r0 = r0 * r1 = a^(2^66 - 1)

    bf_square(x_inverse, &r0); //r0 = a^(2^67 - 2)
}

uint8_t uint8_t_ith_coefficient(bf_element_t *a, uint16_t i) {
    uint8_t pos = i / 64u;
    uint8_t bit =  ((uint64_t) (pos * (a->high >> (i - 64u))) ^ (1u-pos)*(a->low >> i)) & 0x1u;

    return bit;
}


void bf_and(bf_element_t *c, bf_element_t *b, bf_element_t *a) {
    c->low = b->low & a->low;
    c->high = b->high & a->high;
}

void bf_serialize_list(uint8_t * serialized_list, uint16_t serialized_list_len, bf_element_t *list, uint16_t number_of_elements) {
    uint16_t byte_counter = 0;
    for (uint16_t i = 0; i < number_of_elements; ++i) {
        for (uint16_t j = 0; j < ROLLO_I_BF_ELEMENT_BYTE_SIZE - 1; ++j) {
            serialized_list[byte_counter] = list[i].low >> (j*8u);
            byte_counter++;
        }
        serialized_list[byte_counter] = list[i].high;
        byte_counter++;
    }
}

#endif