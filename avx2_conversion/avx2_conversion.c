#include <stdio.h>
#include <string.h>
#include "avx2_conversion.h"
#include "lattices/kyber/kyber_constants.h"

#if defined(NEON)

__m256i _mm256_set1_epi16(const int16_t x) {

    __m256i vector;
    vector.high = vld1q_dup_s16(&x);
    vector.low = vld1q_dup_s16(&x);

    return vector;
}

__m128i _mm_set1_epi16(short a) {

    __m128i vector;
    vector = vld1q_dup_s16(&a);

    return vector;
}

__m256i _mm256_srli_epi16(__m256i a, const int imm8) {

    __m256i vector;
    vector.low = a.low;
    vector.high = a.high;
    switch (imm8) {
        case 0:
            break;
        case 1:
            vector.low = vshrq_n_u16(vector.low, 1);
            vector.high = vshrq_n_u16(vector.high, 1);
            break;
        case 2:
            vector.low = vshrq_n_u16(vector.low, 2);
            vector.high = vshrq_n_u16(vector.high, 2);
            break;
        case 3:
            vector.low = vshrq_n_u16(vector.low, 3);
            vector.high = vshrq_n_u16(vector.high, 3);
            break;
        case 4:
            vector.low = vshrq_n_u16(vector.low, 4);
            vector.high = vshrq_n_u16(vector.high, 4);
            break;
        case 5:
            vector.low = vshrq_n_u16(vector.low, 5);
            vector.high = vshrq_n_u16(vector.high, 5);
            break;
        case 6:
            vector.low = vshrq_n_u16(vector.low, 6);
            vector.high = vshrq_n_u16(vector.high, 6);
            break;
        case 7:
            vector.low = vshrq_n_u16(vector.low, 7);
            vector.high = vshrq_n_u16(vector.high, 7);
            break;
        case 8:
            vector.low = vshrq_n_u16(vector.low, 8);
            vector.high = vshrq_n_u16(vector.high, 8);
            break;
        case 9:
            vector.low = vshrq_n_u16(vector.low, 9);
            vector.high = vshrq_n_u16(vector.high, 9);
            break;
        case 10:
            vector.low = vshrq_n_u16(vector.low, 10);
            vector.high = vshrq_n_u16(vector.high, 10);
            break;
        case 11:
            vector.low = vshrq_n_u16(vector.low, 11);
            vector.high = vshrq_n_u16(vector.high, 11);
            break;
        case 12:
            vector.low = vshrq_n_u16(vector.low, 12);
            vector.high = vshrq_n_u16(vector.high, 12);
            break;
        case 13:
            vector.low = vshrq_n_u16(vector.low, 13);
            vector.high = vshrq_n_u16(vector.high, 13);
            break;
        case 14:
            vector.low = vshrq_n_u16(vector.low, 14);
            vector.high = vshrq_n_u16(vector.high, 14);
            break;
        case 15:
            vector.low = vshrq_n_u16(vector.low, 15);
            vector.high = vshrq_n_u16(vector.high, 15);
            break;
        default:
            break;
    }

    return vector;
}

__m256i _mm256_srai_epi16(__m256i a, const int imm8) {

    __m256i vector;
    vector.low = a.low;
    vector.high = a.high;

    switch (imm8) {
        case 0:
            break;
        case 1:
            vector.low = vshrq_n_s16(vector.low, 1);
            vector.high = vshrq_n_s16(vector.high, 1);
            break;
        case 2:
            vector.low = vshrq_n_s16(vector.low, 2);
            vector.high = vshrq_n_s16(vector.high, 2);
            break;
        case 3:
            vector.low = vshrq_n_s16(vector.low, 3);
            vector.high = vshrq_n_s16(vector.high, 3);
            break;
        case 4:
            vector.low = vshrq_n_s16(vector.low, 4);
            vector.high = vshrq_n_s16(vector.high, 4);
            break;
        case 5:
            vector.low = vshrq_n_s16(vector.low, 5);
            vector.high = vshrq_n_s16(vector.high, 5);
            break;
        case 6:
            vector.low = vshrq_n_s16(vector.low, 6);
            vector.high = vshrq_n_s16(vector.high, 6);
            break;
        case 7:
            vector.low = vshrq_n_s16(vector.low, 7);
            vector.high = vshrq_n_s16(vector.high, 7);
            break;
        case 8:
            vector.low = vshrq_n_s16(vector.low, 8);
            vector.high = vshrq_n_s16(vector.high, 8);
            break;
        case 9:
            vector.low = vshrq_n_s16(vector.low, 9);
            vector.high = vshrq_n_s16(vector.high, 9);
            break;
        case 10:
            vector.low = vshrq_n_s16(vector.low, 10);
            vector.high = vshrq_n_s16(vector.high, 10);
            break;
        case 11:
            vector.low = vshrq_n_s16(vector.low, 11);
            vector.high = vshrq_n_s16(vector.high, 11);
            break;
        case 12:
            vector.low = vshrq_n_s16(vector.low, 12);
            vector.high = vshrq_n_s16(vector.high, 12);
            break;
        case 13:
            vector.low = vshrq_n_s16(vector.low, 13);
            vector.high = vshrq_n_s16(vector.high, 13);
            break;
        case 14:
            vector.low = vshrq_n_s16(vector.low, 14);
            vector.high = vshrq_n_s16(vector.high, 14);
            break;
        case 15:
            vector.low = vshrq_n_s16(vector.low, 15);
            vector.high = vshrq_n_s16(vector.high, 15);
            break;
        default:
            break;
    }

    return vector;
}

__m256i _mm256_srai_epi32(__m256i a, int imm8){
    __m256i vector;
    int32x4_t low = vreinterpretq_s32_s16(a.low);
    int32x4_t high = vreinterpretq_s32_s16(a.high);

    switch (imm8) {
        case 0:
            break;
        case 1:
            low = vshrq_n_s32(low, 1);
            high = vshrq_n_s32(high, 1);
            break;
        case 2:
            low = vshrq_n_s32(low, 2);
            high = vshrq_n_s32(high, 2);
            break;
        case 3:
            low = vshrq_n_s32(low, 3);
            high = vshrq_n_s32(high, 3);
            break;
        case 4:
            low = vshrq_n_s32(low, 4);
            high = vshrq_n_s32(high, 4);
            break;
        case 5:
            low = vshrq_n_s32(low, 5);
            high = vshrq_n_s32(high, 5);
            break;
        case 6:
            low = vshrq_n_s32(low, 6);
            high = vshrq_n_s32(high, 6);
            break;
        case 7:
            low = vshrq_n_s32(low, 7);
            high = vshrq_n_s32(high, 7);
            break;
        case 8:
            low = vshrq_n_s32(low, 8);
            high = vshrq_n_s32(high, 8);
            break;
        case 9:
            low = vshrq_n_s32(low, 9);
            high = vshrq_n_s32(high, 9);
            break;
        case 10:
            low = vshrq_n_s32(low, 10);
            high = vshrq_n_s32(high, 10);
            break;
        case 11:
            low = vshrq_n_s32(low, 11);
            high = vshrq_n_s32(high, 11);
            break;
        case 12:
            low = vshrq_n_s32(low, 12);
            high = vshrq_n_s32(high, 12);
            break;
        case 13:
            low = vshrq_n_s32(low, 13);
            high = vshrq_n_s32(high, 13);
            break;
        case 14:
            low = vshrq_n_s32(low, 14);
            high = vshrq_n_s32(high, 14);
            break;
        case 15:
            low = vshrq_n_s32(low, 15);
            high = vshrq_n_s32(high, 15);
            break;
        case 16:
            low = vshrq_n_s32(low, 16);
            high = vshrq_n_s32(high, 16);
            break;
        case 17:
            low = vshrq_n_s32(low, 17);
            high = vshrq_n_s32(high, 17);
            break;
        case 18:
            low = vshrq_n_s32(low, 18);
            high = vshrq_n_s32(high, 18);
            break;
        case 19:
            low = vshrq_n_s32(low, 19);
            high = vshrq_n_s32(high, 19);
            break;
        case 20:
            low = vshrq_n_s32(low, 20);
            high = vshrq_n_s32(high, 20);
            break;
        case 21:
            low = vshrq_n_s32(low, 21);
            high = vshrq_n_s32(high, 21);
            break;
        case 22:
            low = vshrq_n_s32(low, 22);
            high = vshrq_n_s32(high, 22);
            break;
        case 23:
            low = vshrq_n_s32(low, 23);
            high = vshrq_n_s32(high, 23);
            break;
        case 24:
            low = vshrq_n_s32(low, 24);
            high = vshrq_n_s32(high, 24);
            break;
        case 25:
            low = vshrq_n_s32(low, 25);
            high = vshrq_n_s32(high, 25);
            break;
        case 26:
            low = vshrq_n_s32(low, 26);
            high = vshrq_n_s32(high, 26);
            break;
        case 27:
            low = vshrq_n_s32(low, 27);
            high = vshrq_n_s32(high, 27);
            break;
        case 28:
            low = vshrq_n_s32(low, 28);
            high = vshrq_n_s32(high, 28);
            break;
        case 29:
            low = vshrq_n_s32(low, 29);
            high = vshrq_n_s32(high, 29);
            break;
        case 30:
            low = vshrq_n_s32(low, 30);
            high = vshrq_n_s32(high, 30);
            break;
        case 31:
            low = vshrq_n_s32(low, 31);
            high = vshrq_n_s32(high, 31);
            break;
        default:
            break;
    }

    vector.low = vreinterpretq_s16_s32(low);
    vector.high = vreinterpretq_s16_s32(high);
    return vector;
}

__m256i _mm256_mullo_epi16(__m256i a, __m256i b) {

    __m256i vector;

    vector.high = vmulq_s16(a.high, b.high);
    vector.low = vmulq_s16(a.low, b.low);

    return vector;
}


__m256i _mm256_mulhi_epu16(__m256i a, __m256i b) {

    __m256i vector;

    uint16x4_t tmp_low_1 = {a.low[0], a.low[1], a.low[2], a.low[3]};
    uint16x4_t tmp_low_2 = {b.low[0], b.low[1], b.low[2], b.low[3]};
    uint16x4_t tmp_low_3 = {a.low[4], a.low[5], a.low[6], a.low[7]};
    uint16x4_t tmp_low_4 = {b.low[4], b.low[5], b.low[6], b.low[7]};

    uint32x4_t tmp_low_1_32 = vmull_u16(tmp_low_1, tmp_low_2);
    uint32x4_t tmp_low_2_32 = vmull_u16(tmp_low_3, tmp_low_4);

    tmp_low_1 = vshrn_n_u32(tmp_low_1_32, 16);
    tmp_low_2 = vshrn_n_u32(tmp_low_2_32, 16);

    vector.low = vcombine_s16(tmp_low_1, tmp_low_2);

    uint16x4_t tmp_high_1 = {a.high[0], a.high[1], a.high[2], a.high[3]};
    uint16x4_t tmp_high_2 = {b.high[0], b.high[1], b.high[2], b.high[3]};
    uint16x4_t tmp_high_3 = {a.high[4], a.high[5], a.high[6], a.high[7]};
    uint16x4_t tmp_high_4 = {b.high[4], b.high[5], b.high[6], b.high[7]};

    uint32x4_t tmp_high_1_32 = vmull_u16(tmp_high_1, tmp_high_2);
    uint32x4_t tmp_high_2_32 = vmull_u16(tmp_high_3, tmp_high_4);

    tmp_high_1 = vshrn_n_u32(tmp_high_1_32, 16);
    tmp_high_2 = vshrn_n_u32(tmp_high_2_32, 16);
    vector.high = vcombine_s16(tmp_high_1, tmp_high_2);

    return vector;
}


__m256i _mm256_mulhi_epi16(__m256i a, __m256i b) {

    __m256i vector;

    int16x4_t tmp_low_1 = {a.low[0], a.low[1], a.low[2], a.low[3]};
    int16x4_t tmp_low_2 = {b.low[0], b.low[1], b.low[2], b.low[3]};
    int16x4_t tmp_low_3 = {a.low[4], a.low[5], a.low[6], a.low[7]};
    int16x4_t tmp_low_4 = {b.low[4], b.low[5], b.low[6], b.low[7]};

    int32x4_t tmp_low_1_32 = vmull_s16(tmp_low_1, tmp_low_2);
    int32x4_t tmp_low_2_32 = vmull_s16(tmp_low_3, tmp_low_4);

    tmp_low_1 = vshrn_n_s32(tmp_low_1_32, 16);
    tmp_low_2 = vshrn_n_s32(tmp_low_2_32, 16);

    vector.low = vcombine_s16(tmp_low_1, tmp_low_2);

    int16x4_t tmp_high_1 = {a.high[0], a.high[1], a.high[2], a.high[3]};
    int16x4_t tmp_high_2 = {b.high[0], b.high[1], b.high[2], b.high[3]};
    int16x4_t tmp_high_3 = {a.high[4], a.high[5], a.high[6], a.high[7]};
    int16x4_t tmp_high_4 = {b.high[4], b.high[5], b.high[6], b.high[7]};

    int32x4_t tmp_high_1_32 = vmull_s16(tmp_high_1, tmp_high_2);
    int32x4_t tmp_high_2_32 = vmull_s16(tmp_high_3, tmp_high_4);

    tmp_high_1 = vshrn_n_s32(tmp_high_1_32, 16);
    tmp_high_2 = vshrn_n_s32(tmp_high_2_32, 16);
    vector.high = vcombine_s16(tmp_high_1, tmp_high_2);

    return vector;
}


__m256i _mm256_sub_epi16(__m256i a, __m256i b) {
    __m256i vector;

    vector.high = vsubq_s16(a.high, b.high);
    vector.low = vsubq_s16(a.low, b.low);

    return vector;

}

__m256i _mm256_add_epi16(__m256i a, __m256i b) {
    __m256i vector;

    vector.high = vaddq_s16(a.high, b.high);
    vector.low = vaddq_s16(a.low, b.low);

    return vector;

}

__m256i _mm256_unpacklo_epi16(__m256i a, __m256i b) {

    __m256i vector;

    vector.high = vzip1q_s16(a.high, b.high);
    vector.low = vzip1q_s16(a.low, b.low);

    return vector;

}

__m256i _mm256_unpackhi_epi16(__m256i a, __m256i b) {

    __m256i vector;

    vector.high = vzip2q_s16(a.high, b.high);
    vector.low = vzip2q_s16(a.low, b.low);

    return vector;

}

__m256i _mm256_xor_si256(__m256i a, __m256i b) {

    __m256i vector;

    vector.high = veorq_s16(a.high, b.high);
    vector.low = veorq_s16(a.low, b.low);

    return vector;

}

__m256i _mm256_and_si256(__m256i a, __m256i b) {

    __m256i vector;

    vector.high = vandq_s16(a.high, b.high);
    vector.low = vandq_s16(a.low, b.low);

    return vector;

}

__m256i _mm256_set1_epi32(int x) {

    __m256i vector;
    vector.high = vreinterpretq_s16_s32(vld1q_dup_s32(&x));
    vector.low = vreinterpretq_s16_s32(vld1q_dup_s32(&x));

    return vector;

}

__m256i _mm256_add_epi32(__m256i a, __m256i b) {

    __m256i vector;
    int32x4_t low = vreinterpretq_s32_s16(a.low);
    int32x4_t high = vreinterpretq_s32_s16(a.high);
    low = vaddq_s32(low, vreinterpretq_s32_s16(b.low));
    high = vaddq_s32(high, vreinterpretq_s32_s16(b.high));


    vector.low = vreinterpretq_s16_s32(low);
    vector.high = vreinterpretq_s16_s32(high);

    return vector;
}

__m256i _mm256_srli_epi32(__m256i a, int imm8) {

    __m256i vector;
    uint32x4_t low = vreinterpretq_u32_s16(a.low);
    uint32x4_t high = vreinterpretq_u32_s16(a.high);

    switch (imm8) {
        case 0:
            break;
        case 1:
            low = vshrq_n_u32(low, 1);
            high = vshrq_n_u32(high, 1);
            break;
        case 2:
            low = vshrq_n_u32(low, 2);
            high = vshrq_n_u32(high, 2);
            break;
        case 3:
            low = vshrq_n_u32(low, 3);
            high = vshrq_n_u32(high, 3);
            break;
        case 4:
            low = vshrq_n_u32(low, 4);
            high = vshrq_n_u32(high, 4);
            break;
        case 5:
            low = vshrq_n_u32(low, 5);
            high = vshrq_n_u32(high, 5);
            break;
        case 6:
            low = vshrq_n_u32(low, 6);
            high = vshrq_n_u32(high, 6);
            break;
        case 7:
            low = vshrq_n_u32(low, 7);
            high = vshrq_n_u32(high, 7);
            break;
        case 8:
            low = vshrq_n_u32(low, 8);
            high = vshrq_n_u32(high, 8);
            break;
        case 9:
            low = vshrq_n_u32(low, 9);
            high = vshrq_n_u32(high, 9);
            break;
        case 10:
            low = vshrq_n_u32(low, 10);
            high = vshrq_n_u32(high, 10);
            break;
        case 11:
            low = vshrq_n_u32(low, 11);
            high = vshrq_n_u32(high, 11);
            break;
        case 12:
            low = vshrq_n_u32(low, 12);
            high = vshrq_n_u32(high, 12);
            break;
        case 13:
            low = vshrq_n_u32(low, 13);
            high = vshrq_n_u32(high, 13);
            break;
        case 14:
            low = vshrq_n_u32(low, 14);
            high = vshrq_n_u32(high, 14);
            break;
        case 15:
            low = vshrq_n_u32(low, 15);
            high = vshrq_n_u32(high, 15);
            break;
        case 16:
            low = vshrq_n_u32(low, 16);
            high = vshrq_n_u32(high, 16);
            break;
        case 17:
            low = vshrq_n_u32(low, 17);
            high = vshrq_n_u32(high, 17);
            break;
        case 18:
            low = vshrq_n_u32(low, 18);
            high = vshrq_n_u32(high, 18);
            break;
        case 19:
            low = vshrq_n_u32(low, 19);
            high = vshrq_n_u32(high, 19);
            break;
        case 20:
            low = vshrq_n_u32(low, 20);
            high = vshrq_n_u32(high, 20);
            break;
        case 21:
            low = vshrq_n_u32(low, 21);
            high = vshrq_n_u32(high, 21);
            break;
        case 22:
            low = vshrq_n_u32(low, 22);
            high = vshrq_n_u32(high, 22);
            break;
        case 23:
            low = vshrq_n_u32(low, 23);
            high = vshrq_n_u32(high, 23);
            break;
        case 24:
            low = vshrq_n_u32(low, 24);
            high = vshrq_n_u32(high, 24);
            break;
        case 25:
            low = vshrq_n_u32(low, 25);
            high = vshrq_n_u32(high, 25);
            break;
        case 26:
            low = vshrq_n_u32(low, 26);
            high = vshrq_n_u32(high, 26);
            break;
        case 27:
            low = vshrq_n_u32(low, 27);
            high = vshrq_n_u32(high, 27);
            break;
        case 28:
            low = vshrq_n_u32(low, 28);
            high = vshrq_n_u32(high, 28);
            break;
        case 29:
            low = vshrq_n_u32(low, 29);
            high = vshrq_n_u32(high, 29);
            break;
        case 30:
            low = vshrq_n_u32(low, 30);
            high = vshrq_n_u32(high, 30);
            break;
        case 31:
            low = vshrq_n_u32(low, 31);
            high = vshrq_n_u32(high, 31);
            break;
        default:
            break;
    }

    vector.low = vreinterpretq_s16_u32(low);
    vector.high = vreinterpretq_s16_u32(high);
    return vector;
}

__m256i _mm256_sub_epi32(__m256i a, __m256i b) {

    __m256i vector;
    int32x4_t low = vreinterpretq_s32_s16(a.low);
    int32x4_t high = vreinterpretq_s32_s16(a.high);
    low = vsubq_s32(low, vreinterpretq_s32_s16(b.low));
    high = vsubq_s32(high, vreinterpretq_s32_s16(b.high));

    vector.low = vreinterpretq_s16_s32(low);
    vector.high = vreinterpretq_s16_s32(high);

    return vector;
}

__m256i _mm256_mullo_epi32(__m256i a, __m256i b) {

    __m256i vector;
    int32x4_t low = vreinterpretq_s32_s16(a.low);
    int32x4_t high = vreinterpretq_s32_s16(a.high);
    low = vmulq_s32(low, vreinterpretq_s32_s16(b.low));
    high = vmulq_s32(high, vreinterpretq_s32_s16(b.high));


    vector.low = vreinterpretq_s16_s32(low);
    vector.high = vreinterpretq_s16_s32(high);

    return vector;
}

__m256i _mm256_packus_epi32(__m256i a, __m256i b) {

    __m256i vector;
    int16x4_t low_first_half, low_second_half, high_first_half, high_second_half;

    low_first_half = vqmovn_u32(vreinterpretq_u32_s16(a.low));
    low_second_half = vqmovn_u32(vreinterpretq_u32_s16(b.low));
    high_first_half = vqmovn_u32(vreinterpretq_u32_s16(a.high));
    high_second_half = vqmovn_u32(vreinterpretq_u32_s16(b.high));

    vector.low = vcombine_s16(vreinterpret_u16_s16(low_first_half), vreinterpret_u16_s16(low_second_half));
    vector.high = vcombine_s16(vreinterpret_u16_s16(high_first_half), vreinterpret_u16_s16(high_second_half));

    return vector;
}


__m256i _mm256_packs_epi32(__m256i a, __m256i b) {

    __m256i vector;
    int16x4_t low_first_half, low_second_half, high_first_half, high_second_half;

    low_first_half = vqmovn_s32(vreinterpretq_s32_s16(a.low));
    low_second_half = vqmovn_s32(vreinterpretq_s32_s16(b.low));
    high_first_half = vqmovn_s32(vreinterpretq_s32_s16(a.high));
    high_second_half = vqmovn_s32(vreinterpretq_s32_s16(b.high));

    vector.low = vcombine_s16(low_first_half, low_second_half);
    vector.high = vcombine_s16(high_first_half, high_second_half);

    return vector;
}


__m256i _mm256_permute2f128_si256(__m256i a, __m256i b, int imm8) {

    uint16_t zero[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    __m256i vector;

    uint8_t control_first_half = (uint8_t) imm8 & (uint8_t) 0xBu; //discards the third bit
    uint8_t control_second_half = (uint8_t) (imm8 >> 4) & (uint8_t) 0xBu; //discards the third bit
    vector.low = vld1q_dup_u16(zero);
    vector.high = vld1q_dup_u16(zero);

    switch (control_first_half) {
        case 0:
            vector.low = a.low;
            break;

        case 1:
            vector.low = a.high;
            break;

        case 2:
            vector.low = b.low;
            break;

        case 3:
            vector.low = b.high;
            break;

        default:
            break;
    }

    switch (control_second_half) {
        case 0:
            vector.high = a.low;
            break;

        case 1:
            vector.high = a.high;
            break;

        case 2:
            vector.high = b.low;
            break;

        case 3:
            vector.high = b.high;
            break;

        default:
            break;
    }

    return vector;

}

__m256i _mm256_permute2x128_si256(__m256i a, __m256i b, const int imm8){
    return _mm256_permute2f128_si256(a,b,imm8);
}

__m256i _mm256_permute4x64_epi64(__m256i a, const int imm8) {

    __m256i vector;
    uint64x2_t low, high, tmp;
    low = vreinterpretq_u64_s16(a.low);
    high = vreinterpretq_u64_s16(a.high);

    int control = (imm8 & 0x3u);

    switch (control) {
        case 0:
            tmp[0] = low[0];
            break;

        case 1:
            tmp[0] = low[1];
            break;

        case 2:
            tmp[0] = high[0];
            break;

        case 3:
            tmp[0] = high[1];
            break;

        default:
            break;
    }

    control = ((imm8 >> 2) & 0x3u);

    switch (control) {
        case 0:
            tmp[1] = low[0];
            break;

        case 1:
            tmp[1] = low[1];
            break;

        case 2:
            tmp[1] = high[0];
            break;

        case 3:
            tmp[1] = high[1];
            break;

        default:
            break;
    }

    vector.low = vreinterpretq_s16_u64(tmp);

    control = (imm8 >> 4) & 0x3u;

    switch (control) {
        case 0:
            tmp[0] = low[0];
            break;

        case 1:
            tmp[0] = low[1];
            break;

        case 2:
            tmp[0] = high[0];
            break;

        case 3:
            tmp[0] = high[1];
            break;

        default:
            break;
    }

    control = (imm8 >> 6) & 0x3u;

    switch (control) {
        case 0:
            tmp[1] = low[0];
            break;

        case 1:
            tmp[1] = low[1];
            break;

        case 2:
            tmp[1] = high[0];
            break;

        case 3:
            tmp[1] = high[1];
            break;

        default:
            break;
    }
    vector.high = vreinterpretq_s16_u64(tmp);


    return vector;

}

__m256i _mm256_or_si256(__m256i a, __m256i b) {

    __m256i vector;

    vector.low = vorrq_s16(a.low, b.low);
    vector.high = vorrq_s16(a.high, b.high);

    return vector;
}

__m256i _mm256_cmpeq_epi16(__m256i a, __m256i b) {

    __m256i vector;

    vector.low = vceqq_s16(a.low, b.low);
    vector.high = vceqq_s16(a.high, b.high);

    return vector;

}

__m256i _mm256_cmpgt_epi16(__m256i a, __m256i b){
    __m256i returned_value;

    returned_value.low = vcgtq_s16(a.low, b.low);
    returned_value.high = vcgtq_s16(a.high, b.high);

    return returned_value;
}

__m256i _mm256_cvtepu16_epi32(__m128i a) {

    __m256i vector;
    int16x8_t zero = {0, 0, 0, 0, 0, 0, 0, 0};
    vector.high = vzip2q_s16(a, zero);
    vector.low = vzip1q_s16(a, zero);

    return vector;

}

__m128i _mm256_extractf128_si256(__m256i a, const int imm8) {

    __m128i vector = {0, 0, 0, 0};
    switch (imm8) {
        case 0:
            vector = a.low;
            break;
        case 1:
            vector = a.high;
            break;
        default:
            break;
    }
    return vector;
}

__m256i _mm256_set_m128i(__m128i hi, __m128i lo) {

    __m256i vector;
    vector.high = hi;
    vector.low = lo;

    return vector;
}

__m128i _mm_set1_epi32(int a) {

    __m128i vector;
    vector = vreinterpretq_s16_s32(vld1q_dup_s32(&a));

    return vector;
}

__m256i _mm256_set_epi32(int e7, int e6, int e5, int e4, int e3, int e2, int e1, int e0) {

    __m256i vector;

    uint32x4_t low, high;

    low[0] = e0;
    low[1] = e1;
    low[2] = e2;
    low[3] = e3;
    high[0] = e4;
    high[1] = e5;
    high[2] = e6;
    high[3] = e7;

    vector.low = vreinterpretq_s16_s32(low);
    vector.high = vreinterpretq_s16_s32(high);

    return vector;
}

__m256i _mm256_slli_epi32(__m256i a, int imm8) {

    __m256i vector;
    uint32x4_t low = vreinterpretq_u32_s16(a.low);
    uint32x4_t high = vreinterpretq_u32_s16(a.high);

    switch (imm8) {
        case 0:
            break;
        case 1:
            low = vshlq_n_u32(low, 1);
            high = vshlq_n_u32(high, 1);
            break;
        case 2:
            low = vshlq_n_u32(low, 2);
            high = vshlq_n_u32(high, 2);
            break;
        case 3:
            low = vshlq_n_u32(low, 3);
            high = vshlq_n_u32(high, 3);
            break;
        case 4:
            low = vshlq_n_u32(low, 4);
            high = vshlq_n_u32(high, 4);
            break;
        case 5:
            low = vshlq_n_u32(low, 5);
            high = vshlq_n_u32(high, 5);
            break;
        case 6:
            low = vshlq_n_u32(low, 6);
            high = vshlq_n_u32(high, 6);
            break;
        case 7:
            low = vshlq_n_u32(low, 7);
            high = vshlq_n_u32(high, 7);
            break;
        case 8:
            low = vshlq_n_u32(low, 8);
            high = vshlq_n_u32(high, 8);
            break;
        case 9:
            low = vshlq_n_u32(low, 9);
            high = vshlq_n_u32(high, 9);
            break;
        case 10:
            low = vshlq_n_u32(low, 10);
            high = vshlq_n_u32(high, 10);
            break;
        case 11:
            low = vshlq_n_u32(low, 11);
            high = vshlq_n_u32(high, 11);
            break;
        case 12:
            low = vshlq_n_u32(low, 12);
            high = vshlq_n_u32(high, 12);
            break;
        case 13:
            low = vshlq_n_u32(low, 13);
            high = vshlq_n_u32(high, 13);
            break;
        case 14:
            low = vshlq_n_u32(low, 14);
            high = vshlq_n_u32(high, 14);
            break;
        case 15:
            low = vshlq_n_u32(low, 15);
            high = vshlq_n_u32(high, 15);
            break;
        case 16:
            low = vshlq_n_u32(low, 16);
            high = vshlq_n_u32(high, 16);
            break;
        case 17:
            low = vshlq_n_u32(low, 17);
            high = vshlq_n_u32(high, 17);
            break;
        case 18:
            low = vshlq_n_u32(low, 18);
            high = vshlq_n_u32(high, 18);
            break;
        case 19:
            low = vshlq_n_u32(low, 19);
            high = vshlq_n_u32(high, 19);
            break;
        case 20:
            low = vshlq_n_u32(low, 20);
            high = vshlq_n_u32(high, 20);
            break;
        case 21:
            low = vshlq_n_u32(low, 21);
            high = vshlq_n_u32(high, 21);
            break;
        case 22:
            low = vshlq_n_u32(low, 22);
            high = vshlq_n_u32(high, 22);
            break;
        case 23:
            low = vshlq_n_u32(low, 23);
            high = vshlq_n_u32(high, 23);
            break;
        case 24:
            low = vshlq_n_u32(low, 24);
            high = vshlq_n_u32(high, 24);
            break;
        case 25:
            low = vshlq_n_u32(low, 25);
            high = vshlq_n_u32(high, 25);
            break;
        case 26:
            low = vshlq_n_u32(low, 26);
            high = vshlq_n_u32(high, 26);
            break;
        case 27:
            low = vshlq_n_u32(low, 27);
            high = vshlq_n_u32(high, 27);
            break;
        case 28:
            low = vshlq_n_u32(low, 28);
            high = vshlq_n_u32(high, 28);
            break;
        case 29:
            low = vshlq_n_u32(low, 29);
            high = vshlq_n_u32(high, 29);
            break;
        case 30:
            low = vshlq_n_u32(low, 30);
            high = vshlq_n_u32(high, 30);
            break;
        case 31:
            low = vshlq_n_u32(low, 31);
            high = vshlq_n_u32(high, 31);
            break;
        default:
            break;
    }

    vector.low = vreinterpretq_s16_u32(low);
    vector.high = vreinterpretq_s16_u32(high);
    return vector;
}

__m256i _mm256_permutevar8x32_epi32(__m256i a, __m256i idx) {

    __m256i vector;

    int32x4_t tmp_high = {0, 0, 0, 0};
    int32x4_t tmp_low = {0, 0, 0, 0};
    int32x4_t indexes, a_low, a_high;
    indexes = vreinterpretq_s32_s16(idx.low);
    a_low = vreinterpretq_s32_s16(a.low);
    a_high = vreinterpretq_s32_s16(a.high);

    if (indexes[0] >= 4 && indexes[0] < 8) {
        tmp_low[0] = a_high[indexes[0] - 4];
    } else if (indexes[0] >= 0) {
        tmp_low[0] = a_low[indexes[0]];
    }
    if (indexes[1] >= 4 && indexes[1] < 8) {
        tmp_low[1] = a_high[indexes[1] - 4];
    } else if (indexes[1] >= 0) {
        tmp_low[1] = a_low[indexes[1]];
    }
    if (indexes[2] >= 4 && indexes[2] < 8) {
        tmp_low[2] = a_high[indexes[2] - 4];
    } else if (indexes[2] >= 0) {
        tmp_low[2] = a_low[indexes[2]];
    }
    if (indexes[3] >= 4 && indexes[3] < 8) {
        tmp_low[3] = a_high[indexes[3] - 4];
    } else if (indexes[3] >= 0) {
        tmp_low[3] = a_low[indexes[3]];
    }

    indexes = vreinterpretq_s32_s16(idx.high);

    if (indexes[0] >= 4 && indexes[0] < 8) {
        tmp_high[0] = a_high[indexes[0] - 4];
    } else if (indexes[0] >= 0) {
        tmp_high[0] = a_low[indexes[0]];
    }
    if (indexes[1] >= 4 && indexes[1] < 8) {
        tmp_high[1] = a_high[indexes[1] - 4];
    } else if (indexes[1] >= 0) {
        tmp_high[1] = a_low[indexes[1]];
    }
    if (indexes[2] >= 4 && indexes[2] < 8) {
        tmp_high[2] = a_high[indexes[2] - 4];
    } else if (indexes[2] >= 0) {
        tmp_high[2] = a_low[indexes[2]];
    }
    if (indexes[3] >= 4 && indexes[3] < 8) {
        tmp_high[3] = a_high[indexes[3] - 4];
    } else if (indexes[3] >= 0) {
        tmp_high[3] = a_low[indexes[3]];
    }

    vector.high = vreinterpretq_s16_s32(tmp_high);
    vector.low = vreinterpretq_s16_s32(tmp_low);

    return vector;

}

__m256i _mm256_set_epi16(short e15, short e14, short e13, short e12, short e11, short e10, short e9, short e8, short e7,
                         short e6, short e5, short e4, short e3, short e2, short e1, short e0) {

    __m256i vector;

    vector.low[0] = e0;
    vector.low[1] = e1;
    vector.low[2] = e2;
    vector.low[3] = e3;
    vector.low[4] = e4;
    vector.low[5] = e5;
    vector.low[6] = e6;
    vector.low[7] = e7;
    vector.high[0] = e8;
    vector.high[1] = e9;
    vector.high[2] = e10;
    vector.high[3] = e11;
    vector.high[4] = e12;
    vector.high[5] = e13;
    vector.high[6] = e14;
    vector.high[7] = e15;

    return vector;
}

__m256i _mm256_set1_epi64x(long long a) {

    __m256i vector;

    int64x2_t tmp = {a, a};
    vector.high = vreinterpretq_s16_s64(tmp);
    vector.low = vreinterpretq_s16_s64(tmp);

    return vector;

}
//
//__m256i _mm256_set_epi64x(long long a, long long b, long long c, long long d ) {
//
//    __m256i vector;
//
//    int64x2_t tmp_high = {b, a};
//    int64x2_t tmp_low = {d, c};
//    vector.high = vreinterpretq_s16_s64(tmp_high);
//    vector.low = vreinterpretq_s16_s64(tmp_low);
//
//    return vector;
//
//}

__m256i _mm256_srli_epi64(__m256i a, int imm8) {

    __m256i vector;
    uint64x2_t low = vreinterpretq_u64_s16(a.low);
    uint64x2_t high = vreinterpretq_u64_s16(a.high);

    switch (imm8) {
        case 0:
            break;
        case 1 :
            low = vshrq_n_u64(low, 1);
            high = vshrq_n_u64(high, 1);
            break;
        case 2 :
            low = vshrq_n_u64(low, 2);
            high = vshrq_n_u64(high, 2);
            break;
        case 3 :
            low = vshrq_n_u64(low, 3);
            high = vshrq_n_u64(high, 3);
            break;
        case 4 :
            low = vshrq_n_u64(low, 4);
            high = vshrq_n_u64(high, 4);
            break;
        case 5 :
            low = vshrq_n_u64(low, 5);
            high = vshrq_n_u64(high, 5);
            break;
        case 6 :
            low = vshrq_n_u64(low, 6);
            high = vshrq_n_u64(high, 6);
            break;
        case 7 :
            low = vshrq_n_u64(low, 7);
            high = vshrq_n_u64(high, 7);
            break;
        case 8 :
            low = vshrq_n_u64(low, 8);
            high = vshrq_n_u64(high, 8);
            break;
        case 9 :
            low = vshrq_n_u64(low, 9);
            high = vshrq_n_u64(high, 9);
            break;
        case 10 :
            low = vshrq_n_u64(low, 10);
            high = vshrq_n_u64(high, 10);
            break;
        case 11 :
            low = vshrq_n_u64(low, 11);
            high = vshrq_n_u64(high, 11);
            break;
        case 12 :
            low = vshrq_n_u64(low, 12);
            high = vshrq_n_u64(high, 12);
            break;
        case 13 :
            low = vshrq_n_u64(low, 13);
            high = vshrq_n_u64(high, 13);
            break;
        case 14 :
            low = vshrq_n_u64(low, 14);
            high = vshrq_n_u64(high, 14);
            break;
        case 15 :
            low = vshrq_n_u64(low, 15);
            high = vshrq_n_u64(high, 15);
            break;
        case 16 :
            low = vshrq_n_u64(low, 16);
            high = vshrq_n_u64(high, 16);
            break;
        case 17 :
            low = vshrq_n_u64(low, 17);
            high = vshrq_n_u64(high, 17);
            break;
        case 18 :
            low = vshrq_n_u64(low, 18);
            high = vshrq_n_u64(high, 18);
            break;
        case 19 :
            low = vshrq_n_u64(low, 19);
            high = vshrq_n_u64(high, 19);
            break;
        case 20 :
            low = vshrq_n_u64(low, 20);
            high = vshrq_n_u64(high, 20);
            break;
        case 21 :
            low = vshrq_n_u64(low, 21);
            high = vshrq_n_u64(high, 21);
            break;
        case 22 :
            low = vshrq_n_u64(low, 22);
            high = vshrq_n_u64(high, 22);
            break;
        case 23 :
            low = vshrq_n_u64(low, 23);
            high = vshrq_n_u64(high, 23);
            break;
        case 24 :
            low = vshrq_n_u64(low, 24);
            high = vshrq_n_u64(high, 24);
            break;
        case 25 :
            low = vshrq_n_u64(low, 25);
            high = vshrq_n_u64(high, 25);
            break;
        case 26 :
            low = vshrq_n_u64(low, 26);
            high = vshrq_n_u64(high, 26);
            break;
        case 27 :
            low = vshrq_n_u64(low, 27);
            high = vshrq_n_u64(high, 27);
            break;
        case 28 :
            low = vshrq_n_u64(low, 28);
            high = vshrq_n_u64(high, 28);
            break;
        case 29 :
            low = vshrq_n_u64(low, 29);
            high = vshrq_n_u64(high, 29);
            break;
        case 30 :
            low = vshrq_n_u64(low, 30);
            high = vshrq_n_u64(high, 30);
            break;
        case 31 :
            low = vshrq_n_u64(low, 31);
            high = vshrq_n_u64(high, 31);
            break;
        case 32 :
            low = vshrq_n_u64(low, 32);
            high = vshrq_n_u64(high, 32);
            break;
        case 33 :
            low = vshrq_n_u64(low, 33);
            high = vshrq_n_u64(high, 33);
            break;
        case 34 :
            low = vshrq_n_u64(low, 34);
            high = vshrq_n_u64(high, 34);
            break;
        case 35 :
            low = vshrq_n_u64(low, 35);
            high = vshrq_n_u64(high, 35);
            break;
        case 36 :
            low = vshrq_n_u64(low, 36);
            high = vshrq_n_u64(high, 36);
            break;
        case 37 :
            low = vshrq_n_u64(low, 37);
            high = vshrq_n_u64(high, 37);
            break;
        case 38 :
            low = vshrq_n_u64(low, 38);
            high = vshrq_n_u64(high, 38);
            break;
        case 39 :
            low = vshrq_n_u64(low, 39);
            high = vshrq_n_u64(high, 39);
            break;
        case 40 :
            low = vshrq_n_u64(low, 40);
            high = vshrq_n_u64(high, 40);
            break;
        case 41 :
            low = vshrq_n_u64(low, 41);
            high = vshrq_n_u64(high, 41);
            break;
        case 42 :
            low = vshrq_n_u64(low, 42);
            high = vshrq_n_u64(high, 42);
            break;
        case 43 :
            low = vshrq_n_u64(low, 43);
            high = vshrq_n_u64(high, 43);
            break;
        case 44 :
            low = vshrq_n_u64(low, 44);
            high = vshrq_n_u64(high, 44);
            break;
        case 45 :
            low = vshrq_n_u64(low, 45);
            high = vshrq_n_u64(high, 45);
            break;
        case 46 :
            low = vshrq_n_u64(low, 46);
            high = vshrq_n_u64(high, 46);
            break;
        case 47 :
            low = vshrq_n_u64(low, 47);
            high = vshrq_n_u64(high, 47);
            break;
        case 48 :
            low = vshrq_n_u64(low, 48);
            high = vshrq_n_u64(high, 48);
            break;
        case 49 :
            low = vshrq_n_u64(low, 49);
            high = vshrq_n_u64(high, 49);
            break;
        case 50 :
            low = vshrq_n_u64(low, 50);
            high = vshrq_n_u64(high, 50);
            break;
        case 51 :
            low = vshrq_n_u64(low, 51);
            high = vshrq_n_u64(high, 51);
            break;
        case 52 :
            low = vshrq_n_u64(low, 52);
            high = vshrq_n_u64(high, 52);
            break;
        case 53 :
            low = vshrq_n_u64(low, 53);
            high = vshrq_n_u64(high, 53);
            break;
        case 54 :
            low = vshrq_n_u64(low, 54);
            high = vshrq_n_u64(high, 54);
            break;
        case 55 :
            low = vshrq_n_u64(low, 55);
            high = vshrq_n_u64(high, 55);
            break;
        case 56 :
            low = vshrq_n_u64(low, 56);
            high = vshrq_n_u64(high, 56);
            break;
        case 57 :
            low = vshrq_n_u64(low, 57);
            high = vshrq_n_u64(high, 57);
            break;
        case 58 :
            low = vshrq_n_u64(low, 58);
            high = vshrq_n_u64(high, 58);
            break;
        case 59 :
            low = vshrq_n_u64(low, 59);
            high = vshrq_n_u64(high, 59);
            break;
        case 60 :
            low = vshrq_n_u64(low, 60);
            high = vshrq_n_u64(high, 60);
            break;
        case 61 :
            low = vshrq_n_u64(low, 61);
            high = vshrq_n_u64(high, 61);
            break;
        case 62 :
            low = vshrq_n_u64(low, 62);
            high = vshrq_n_u64(high, 62);
            break;
        case 63 :
            low = vshrq_n_u64(low, 63);
            high = vshrq_n_u64(high, 63);
            break;
        default:
            break;
    }

    vector.low = vreinterpretq_s16_u64(low);
    vector.high = vreinterpretq_s16_u64(high);
    return vector;
}


__m256i _mm256_slli_epi64(__m256i a, int imm8) {

    __m256i vector;
    uint64x2_t low = vreinterpretq_u64_s16(a.low);
    uint64x2_t high = vreinterpretq_u64_s16(a.high);

    switch (imm8) {
        case 0:
            break;
        case 1 :
            low = vshlq_n_u64(low, 1);
            high = vshlq_n_u64(high, 1);
            break;
        case 2 :
            low = vshlq_n_u64(low, 2);
            high = vshlq_n_u64(high, 2);
            break;
        case 3 :
            low = vshlq_n_u64(low, 3);
            high = vshlq_n_u64(high, 3);
            break;
        case 4 :
            low = vshlq_n_u64(low, 4);
            high = vshlq_n_u64(high, 4);
            break;
        case 5 :
            low = vshlq_n_u64(low, 5);
            high = vshlq_n_u64(high, 5);
            break;
        case 6 :
            low = vshlq_n_u64(low, 6);
            high = vshlq_n_u64(high, 6);
            break;
        case 7 :
            low = vshlq_n_u64(low, 7);
            high = vshlq_n_u64(high, 7);
            break;
        case 8 :
            low = vshlq_n_u64(low, 8);
            high = vshlq_n_u64(high, 8);
            break;
        case 9 :
            low = vshlq_n_u64(low, 9);
            high = vshlq_n_u64(high, 9);
            break;
        case 10 :
            low = vshlq_n_u64(low, 10);
            high = vshlq_n_u64(high, 10);
            break;
        case 11 :
            low = vshlq_n_u64(low, 11);
            high = vshlq_n_u64(high, 11);
            break;
        case 12 :
            low = vshlq_n_u64(low, 12);
            high = vshlq_n_u64(high, 12);
            break;
        case 13 :
            low = vshlq_n_u64(low, 13);
            high = vshlq_n_u64(high, 13);
            break;
        case 14 :
            low = vshlq_n_u64(low, 14);
            high = vshlq_n_u64(high, 14);
            break;
        case 15 :
            low = vshlq_n_u64(low, 15);
            high = vshlq_n_u64(high, 15);
            break;
        case 16 :
            low = vshlq_n_u64(low, 16);
            high = vshlq_n_u64(high, 16);
            break;
        case 17 :
            low = vshlq_n_u64(low, 17);
            high = vshlq_n_u64(high, 17);
            break;
        case 18 :
            low = vshlq_n_u64(low, 18);
            high = vshlq_n_u64(high, 18);
            break;
        case 19 :
            low = vshlq_n_u64(low, 19);
            high = vshlq_n_u64(high, 19);
            break;
        case 20 :
            low = vshlq_n_u64(low, 20);
            high = vshlq_n_u64(high, 20);
            break;
        case 21 :
            low = vshlq_n_u64(low, 21);
            high = vshlq_n_u64(high, 21);
            break;
        case 22 :
            low = vshlq_n_u64(low, 22);
            high = vshlq_n_u64(high, 22);
            break;
        case 23 :
            low = vshlq_n_u64(low, 23);
            high = vshlq_n_u64(high, 23);
            break;
        case 24 :
            low = vshlq_n_u64(low, 24);
            high = vshlq_n_u64(high, 24);
            break;
        case 25 :
            low = vshlq_n_u64(low, 25);
            high = vshlq_n_u64(high, 25);
            break;
        case 26 :
            low = vshlq_n_u64(low, 26);
            high = vshlq_n_u64(high, 26);
            break;
        case 27 :
            low = vshlq_n_u64(low, 27);
            high = vshlq_n_u64(high, 27);
            break;
        case 28 :
            low = vshlq_n_u64(low, 28);
            high = vshlq_n_u64(high, 28);
            break;
        case 29 :
            low = vshlq_n_u64(low, 29);
            high = vshlq_n_u64(high, 29);
            break;
        case 30 :
            low = vshlq_n_u64(low, 30);
            high = vshlq_n_u64(high, 30);
            break;
        case 31 :
            low = vshlq_n_u64(low, 31);
            high = vshlq_n_u64(high, 31);
            break;
        case 32 :
            low = vshlq_n_u64(low, 32);
            high = vshlq_n_u64(high, 32);
            break;
        case 33 :
            low = vshlq_n_u64(low, 33);
            high = vshlq_n_u64(high, 33);
            break;
        case 34 :
            low = vshlq_n_u64(low, 34);
            high = vshlq_n_u64(high, 34);
            break;
        case 35 :
            low = vshlq_n_u64(low, 35);
            high = vshlq_n_u64(high, 35);
            break;
        case 36 :
            low = vshlq_n_u64(low, 36);
            high = vshlq_n_u64(high, 36);
            break;
        case 37 :
            low = vshlq_n_u64(low, 37);
            high = vshlq_n_u64(high, 37);
            break;
        case 38 :
            low = vshlq_n_u64(low, 38);
            high = vshlq_n_u64(high, 38);
            break;
        case 39 :
            low = vshlq_n_u64(low, 39);
            high = vshlq_n_u64(high, 39);
            break;
        case 40 :
            low = vshlq_n_u64(low, 40);
            high = vshlq_n_u64(high, 40);
            break;
        case 41 :
            low = vshlq_n_u64(low, 41);
            high = vshlq_n_u64(high, 41);
            break;
        case 42 :
            low = vshlq_n_u64(low, 42);
            high = vshlq_n_u64(high, 42);
            break;
        case 43 :
            low = vshlq_n_u64(low, 43);
            high = vshlq_n_u64(high, 43);
            break;
        case 44 :
            low = vshlq_n_u64(low, 44);
            high = vshlq_n_u64(high, 44);
            break;
        case 45 :
            low = vshlq_n_u64(low, 45);
            high = vshlq_n_u64(high, 45);
            break;
        case 46 :
            low = vshlq_n_u64(low, 46);
            high = vshlq_n_u64(high, 46);
            break;
        case 47 :
            low = vshlq_n_u64(low, 47);
            high = vshlq_n_u64(high, 47);
            break;
        case 48 :
            low = vshlq_n_u64(low, 48);
            high = vshlq_n_u64(high, 48);
            break;
        case 49 :
            low = vshlq_n_u64(low, 49);
            high = vshlq_n_u64(high, 49);
            break;
        case 50 :
            low = vshlq_n_u64(low, 50);
            high = vshlq_n_u64(high, 50);
            break;
        case 51 :
            low = vshlq_n_u64(low, 51);
            high = vshlq_n_u64(high, 51);
            break;
        case 52 :
            low = vshlq_n_u64(low, 52);
            high = vshlq_n_u64(high, 52);
            break;
        case 53 :
            low = vshlq_n_u64(low, 53);
            high = vshlq_n_u64(high, 53);
            break;
        case 54 :
            low = vshlq_n_u64(low, 54);
            high = vshlq_n_u64(high, 54);
            break;
        case 55 :
            low = vshlq_n_u64(low, 55);
            high = vshlq_n_u64(high, 55);
            break;
        case 56 :
            low = vshlq_n_u64(low, 56);
            high = vshlq_n_u64(high, 56);
            break;
        case 57 :
            low = vshlq_n_u64(low, 57);
            high = vshlq_n_u64(high, 57);
            break;
        case 58 :
            low = vshlq_n_u64(low, 58);
            high = vshlq_n_u64(high, 58);
            break;
        case 59 :
            low = vshlq_n_u64(low, 59);
            high = vshlq_n_u64(high, 59);
            break;
        case 60 :
            low = vshlq_n_u64(low, 60);
            high = vshlq_n_u64(high, 60);
            break;
        case 61 :
            low = vshlq_n_u64(low, 61);
            high = vshlq_n_u64(high, 61);
            break;
        case 62 :
            low = vshlq_n_u64(low, 62);
            high = vshlq_n_u64(high, 62);
            break;
        case 63 :
            low = vshlq_n_u64(low, 63);
            high = vshlq_n_u64(high, 63);
            break;
        default:
            break;
    }

    vector.low = vreinterpretq_s16_u64(low);
    vector.high = vreinterpretq_s16_u64(high);
    return vector;
}

__m256i _mm256_hadd_epi16(__m256i a, __m256i b) {


    __m256i vector;

    int16x8_t tmp1 = {a.low[0], a.low[2], a.low[4], a.low[6], b.low[0], b.low[2], b.low[4], b.low[6]};
    int16x8_t tmp2 = {a.low[1], a.low[3], a.low[5], a.low[7], b.low[1], b.low[3], b.low[5], b.low[7]};
    int16x8_t tmp3 = {a.high[0], a.high[2], a.high[4], a.high[6], b.high[0], b.high[2], b.high[4], b.high[6]};
    int16x8_t tmp4 = {a.high[1], a.high[3], a.high[5], a.high[7], b.high[1], b.high[3], b.high[5], b.high[7]};

    vector.low = vaddq_s16(tmp1, tmp2);
    vector.high = vaddq_s16(tmp3, tmp4);
    return vector;

}

__m256i _mm256_hadd_epi32(__m256i a, __m256i b) {

    __m256i vector;

    int32x4_t tmp_a = vreinterpretq_s32_s16(a.low);
    int32x4_t tmp_b = vreinterpretq_s32_s16(b.low);

    int32x4_t tmp1 = {tmp_a[0], tmp_a[2], tmp_b[0], tmp_b[2]};
    int32x4_t tmp2 = {tmp_a[1], tmp_a[3], tmp_b[1], tmp_b[3]};
    tmp_a = vreinterpretq_s32_s16(a.high);
    tmp_b = vreinterpretq_s32_s16(b.high);
    int32x4_t tmp3 = {tmp_a[0], tmp_a[2], tmp_b[0], tmp_b[2]};
    int32x4_t tmp4 = {tmp_a[1], tmp_a[3], tmp_b[1], tmp_b[3]};


    vector.low = vreinterpretq_s16_s32(vaddq_s32(tmp1, tmp2));
    vector.high = vreinterpretq_s16_s32(vaddq_s32(tmp3, tmp4));
    return vector;
}

__m256i _mm256_unpacklo_epi32(__m256i a, __m256i b) {

    __m256i vector;

    vector.low = vreinterpretq_s16_s32(vzip1q_s32(a.low, b.low));
    vector.high = vreinterpretq_s16_s32(vzip1q_s32(a.high, b.high));

    return vector;

}

__m256i _mm256_unpackhi_epi32(__m256i a, __m256i b) {

    __m256i vector;

    vector.low = vreinterpretq_s16_s32(vzip2q_s32(a.low, b.low));
    vector.high = vreinterpretq_s16_s32(vzip2q_s32(a.high, b.high));

    return vector;

}

__m256i _mm256_unpacklo_epi64(__m256i a, __m256i b) {

    __m256i vector;

    vector.low = vreinterpretq_s16_s64(vzip1q_s64(a.low, b.low));
    vector.high = vreinterpretq_s16_s64(vzip1q_s64(a.high, b.high));

    return vector;
}


__m256i _mm256_unpackhi_epi64(__m256i a, __m256i b) {

    __m256i vector;

    vector.low = vreinterpretq_s16_s64(vzip2q_s64(a.low, b.low));
    vector.high = vreinterpretq_s16_s64(vzip2q_s64(a.high, b.high));

    return vector;

}

__m256i _mm256_set_epi64x(__int64 e3, __int64 e2, __int64 e1, __int64 e0) {

    __m256i vector;

    int64x2_t low = {e0, e1};
    int64x2_t high = {e2, e3};

    vector.low = vreinterpretq_s16_s64(low);
    vector.high = vreinterpretq_s16_s64(high);

    return vector;
}


__m256i _mm256_set_epi8(char e31, char e30, char e29, char e28, char e27, char e26, char e25, char e24, char e23,
                        char e22, char e21, char e20, char e19, char e18, char e17, char e16, char e15, char e14,
                        char e13, char e12, char e11, char e10, char e9, char e8, char e7, char e6, char e5, char e4,
                        char e3, char e2, char e1, char e0) {

    __m256i vector;

    int8x16_t low = {e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15};
    int8x16_t high = {e16, e17, e18, e19, e20, e21, e22, e23, e24, e25, e26, e27, e28, e29, e30, e31};

    vector.low = vreinterpretq_s16_s8(low);
    vector.high = vreinterpretq_s16_s8(high);

    return vector;
}

__m256i _mm256_add_epi64(__m256i a, __m256i b) {

    __m256i vector;

    int64x2_t tmp1 = vreinterpretq_s64_s16(a.low);
    int64x2_t tmp2 = vreinterpretq_s64_s16(b.low);

    vector.low = vreinterpretq_s16_s64(vaddq_s64(tmp1, tmp2));

    tmp1 = vreinterpretq_s64_s16(a.high);
    tmp2 = vreinterpretq_s64_s16(b.high);

    vector.high = vreinterpretq_s16_s64(vaddq_s64(tmp1, tmp2));

    return vector;


}


__m256i _mm256_loadu2_m128i(__m128i const *hiaddr, __m128i const *loaddr) {

    __m256i vector;

    vector.high = *hiaddr;
    vector.low = *loaddr;

    return vector;

}

__m256i _mm256_loadu_si256(__m256i const *addr) {
    __m256i vector = *addr;
    return vector;
}

__m256i _mm256_load_si256(__m256i const *addr) {
    __m256i vector = *addr;
    return vector;
}

void _mm256_storeu_si256(__m256i *mem_addr, __m256i a) {
    *mem_addr = a;
}

void _mm256_store_si256(__m256i *mem_addr, __m256i a) {
    *mem_addr = a;
}



__m256i _mm256_setzero_si256(){
    __m256i returned_value;
    memset(&returned_value,0,KYBER_SYMMETRIC_KEY_LENGTH);
    return returned_value;
}


__m128i _mm_srli_epi64(__m128i a, int imm8) {

    uint64x2_t vector = vreinterpretq_u64_s16(a);

    switch (imm8) {
        case 0 :
            break;
        case 1 :
            vector = vshrq_n_u64(vector, 1);
            break;
        case 2 :
            vector = vshrq_n_u64(vector, 2);
            break;
        case 3 :
            vector = vshrq_n_u64(vector, 3);
            break;
        case 4 :
            vector = vshrq_n_u64(vector, 4);
            break;
        case 5 :
            vector = vshrq_n_u64(vector, 5);
            break;
        case 6 :
            vector = vshrq_n_u64(vector, 6);
            break;
        case 7 :
            vector = vshrq_n_u64(vector, 7);
            break;
        case 8 :
            vector = vshrq_n_u64(vector, 8);
            break;
        case 9 :
            vector = vshrq_n_u64(vector, 9);
            break;
        case 10 :
            vector = vshrq_n_u64(vector, 10);
            break;
        case 11 :
            vector = vshrq_n_u64(vector, 11);
            break;
        case 12 :
            vector = vshrq_n_u64(vector, 12);
            break;
        case 13 :
            vector = vshrq_n_u64(vector, 13);
            break;
        case 14 :
            vector = vshrq_n_u64(vector, 14);
            break;
        case 15 :
            vector = vshrq_n_u64(vector, 15);
            break;
        case 16 :
            vector = vshrq_n_u64(vector, 16);
            break;
        case 17 :
            vector = vshrq_n_u64(vector, 17);
            break;
        case 18 :
            vector = vshrq_n_u64(vector, 18);
            break;
        case 19 :
            vector = vshrq_n_u64(vector, 19);
            break;
        case 20 :
            vector = vshrq_n_u64(vector, 20);
            break;
        case 21 :
            vector = vshrq_n_u64(vector, 21);
            break;
        case 22 :
            vector = vshrq_n_u64(vector, 22);
            break;
        case 23 :
            vector = vshrq_n_u64(vector, 23);
            break;
        case 24 :
            vector = vshrq_n_u64(vector, 24);
            break;
        case 25 :
            vector = vshrq_n_u64(vector, 25);
            break;
        case 26 :
            vector = vshrq_n_u64(vector, 26);
            break;
        case 27 :
            vector = vshrq_n_u64(vector, 27);
            break;
        case 28 :
            vector = vshrq_n_u64(vector, 28);
            break;
        case 29 :
            vector = vshrq_n_u64(vector, 29);
            break;
        case 30 :
            vector = vshrq_n_u64(vector, 30);
            break;
        case 31 :
            vector = vshrq_n_u64(vector, 31);
            break;
        case 32 :
            vector = vshrq_n_u64(vector, 32);
            break;
        case 33 :
            vector = vshrq_n_u64(vector, 33);
            break;
        case 34 :
            vector = vshrq_n_u64(vector, 34);
            break;
        case 35 :
            vector = vshrq_n_u64(vector, 35);
            break;
        case 36 :
            vector = vshrq_n_u64(vector, 36);
            break;
        case 37 :
            vector = vshrq_n_u64(vector, 37);
            break;
        case 38 :
            vector = vshrq_n_u64(vector, 38);
            break;
        case 39 :
            vector = vshrq_n_u64(vector, 39);
            break;
        case 40 :
            vector = vshrq_n_u64(vector, 40);
            break;
        case 41 :
            vector = vshrq_n_u64(vector, 41);
            break;
        case 42 :
            vector = vshrq_n_u64(vector, 42);
            break;
        case 43 :
            vector = vshrq_n_u64(vector, 43);
            break;
        case 44 :
            vector = vshrq_n_u64(vector, 44);
            break;
        case 45 :
            vector = vshrq_n_u64(vector, 45);
            break;
        case 46 :
            vector = vshrq_n_u64(vector, 46);
            break;
        case 47 :
            vector = vshrq_n_u64(vector, 47);
            break;
        case 48 :
            vector = vshrq_n_u64(vector, 48);
            break;
        case 49 :
            vector = vshrq_n_u64(vector, 49);
            break;
        case 50 :
            vector = vshrq_n_u64(vector, 50);
            break;
        case 51 :
            vector = vshrq_n_u64(vector, 51);
            break;
        case 52 :
            vector = vshrq_n_u64(vector, 52);
            break;
        case 53 :
            vector = vshrq_n_u64(vector, 53);
            break;
        case 54 :
            vector = vshrq_n_u64(vector, 54);
            break;
        case 55 :
            vector = vshrq_n_u64(vector, 55);
            break;
        case 56 :
            vector = vshrq_n_u64(vector, 56);
            break;
        case 57 :
            vector = vshrq_n_u64(vector, 57);
            break;
        case 58 :
            vector = vshrq_n_u64(vector, 58);
            break;
        case 59 :
            vector = vshrq_n_u64(vector, 59);
            break;
        case 60 :
            vector = vshrq_n_u64(vector, 60);
            break;
        case 61 :
            vector = vshrq_n_u64(vector, 61);
            break;
        case 62 :
            vector = vshrq_n_u64(vector, 62);
            break;
        case 63 :
            vector = vshrq_n_u64(vector, 63);
            break;
        default:
            break;
    }

    return vreinterpretq_s16_u64(vector);
}

int _mm256_movemask_epi8(__m256i a){

    int mask = 0;

    int i;
    for (i = 0; i<16; i+=2){
        mask |= (((a.low[i/2] >> 7) & 0x1) << i);
        mask |= (((a.low[i/2] >> 15) & 0x1) << (i + 1));
    }

    for (i = 0; i<16; i+=2){
        mask |= (((a.high[i/2] >> 7) & 0x1) << (i + 16));
        mask |= (((a.high[i/2] >> 15) & 0x1) << (i + 1 + 16));
    }

    return mask;
}

__m256i _mm256_shuffle_epi8(__m256i a, __m256i b) {

    __m256i vector;

    int8x16_t tmp1, tmp2, tmp3, tmp4, tmp_index;
    int8x16_t shift_index = {16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16};

    tmp1 = vqtbl1q_s8(vreinterpretq_s8_s16(a.low), vreinterpretq_u8_s16(b.low));
    tmp2 = vqtbl1q_s8(vreinterpretq_s8_s16(a.low), vreinterpretq_u8_s16(b.high));
    tmp_index = vsubq_s8(vreinterpretq_s8_s16(b.low), shift_index);
    tmp3 = vqtbl1q_s8(vreinterpretq_s8_s16(a.high), tmp_index);
    tmp_index = vsubq_s8(vreinterpretq_s8_s16(b.high), shift_index);
    tmp4 = vqtbl1q_s8(vreinterpretq_s8_s16(a.high), tmp_index);

    vector.low = vreinterpretq_s16_s8(veorq_s8(tmp1, tmp3));
    vector.high = vreinterpretq_s16_s8(veorq_s8(tmp2, tmp4));

    return vector;

}

#endif

