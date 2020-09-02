#ifndef AVX2_CONVERSION_H
#define AVX2_CONVERSION_H

#if defined(AVX2)
#include <immintrin.h>


#ifndef _mm256_set_m128i
#define  _mm256_set_m128i(v0, v1) _mm256_insertf128_si256(_mm256_castsi128_si256(v1), (v0), 1)
#endif


#ifndef _mm256_loadu2_m128i
#define  _mm256_loadu2_m128i(low, high) _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((__m128i*)low)), _mm_loadu_si128((__m128i*)high), 1)
#endif

#endif
#if defined(NEON)

#include <arm_neon.h>

typedef struct __m256i {
    int16x8_t low;
    int16x8_t high;
} __m256i;

typedef int16x8_t __m128i;
typedef int64_t __int64;

//tested

__m256i _mm256_set1_epi16(int16_t x);

__m256i _mm256_set1_epi32(int x);

__m256i _mm256_srli_epi16(__m256i a, int imm8);

__m256i _mm256_srai_epi16(__m256i a, int imm8);

__m256i _mm256_mullo_epi16(__m256i a, __m256i b);

__m256i _mm256_mulhi_epu16(__m256i a, __m256i b);

__m256i _mm256_sub_epi16(__m256i a, __m256i b);

__m256i _mm256_add_epi16(__m256i a, __m256i b);

__m256i _mm256_unpacklo_epi16(__m256i a, __m256i b);

__m256i _mm256_unpackhi_epi16(__m256i a, __m256i b);

__m256i _mm256_and_si256(__m256i a, __m256i b);

__m256i _mm256_xor_si256(__m256i a, __m256i b);

__m256i _mm256_add_epi32(__m256i a, __m256i b);

__m256i _mm256_srli_epi32(__m256i a, int imm8);

__m256i _mm256_sub_epi32(__m256i a, __m256i b);

__m256i _mm256_mullo_epi32(__m256i a, __m256i b);

__m256i _mm256_packus_epi32(__m256i a, __m256i b);

__m256i _mm256_permute2f128_si256(__m256i a, __m256i b, int imm8);

__m256i _mm256_permute4x64_epi64(__m256i a, int imm8);

__m256i _mm256_shuffle_epi8(__m256i a, __m256i b);

__m256i _mm256_or_si256(__m256i a, __m256i b);

__m256i _mm256_cmpeq_epi16(__m256i a, __m256i b);

__m256i _mm256_cvtepu16_epi32(__m128i a);

__m128i _mm256_extractf128_si256(__m256i a, const int imm8);

__m128i _mm_set1_epi32(int a);

__m256i _mm256_set_epi32(int e7, int e6, int e5, int e4, int e3, int e2, int e1, int e0);

__m256i _mm256_slli_epi32(__m256i a, int imm8);

__m256i _mm256_permutevar8x32_epi32(__m256i a, __m256i idx);

__m256i _mm256_set_epi16(short e15, short e14, short e13, short e12, short e11, short e10, short e9, short e8, short e7,
                         short e6, short e5, short e4, short e3, short e2, short e1, short e0);

__m256i _mm256_set1_epi64x(long long a);

__m256i _mm256_srli_epi64(__m256i a, int imm8);

__m256i _mm256_slli_epi64(__m256i a, int imm8);

__m256i _mm256_hadd_epi16(__m256i a, __m256i b);

__m256i _mm256_hadd_epi32(__m256i a, __m256i b);

__m256i _mm256_unpacklo_epi32(__m256i a, __m256i b);

__m256i _mm256_unpackhi_epi32(__m256i a, __m256i b);

__m256i _mm256_unpacklo_epi64(__m256i a, __m256i b);

__m256i _mm256_unpackhi_epi64(__m256i a, __m256i b);

__m256i _mm256_set_epi64x(__int64 e3, __int64 e2, __int64 e1, __int64 e0);

__m256i _mm256_loadu2_m128i(__m128i const *hiaddr, __m128i const *loaddr);

__m128i _mm_srli_epi64(__m128i a, int imm8);



//not tested yet

__m256i _mm256_loadu_si256(__m256i const *addr);

void _mm256_storeu_si256(__m256i * mem_addr, __m256i a);

__m256i _mm256_setzero_si256();

__m256i _mm256_cmpgt_epi16(__m256i a, __m256i b);

__m256i _mm256_srai_epi32(__m256i a, int imm8);

__m256i _mm256_mulhi_epi16(__m256i a, __m256i b);

__m256i _mm256_packs_epi32(__m256i a, __m256i b);

__m256i _mm256_permute2x128_si256(__m256i a, __m256i b, const int imm8);

__m256i _mm256_load_si256(__m256i const *addr);

void _mm256_store_si256(__m256i *mem_addr, __m256i a);

int _mm256_movemask_epi8(__m256i a);

__m256i _mm256_set_m128i(__m128i hi, __m128i lo);

__m128i _mm_set1_epi16(short a);


__m256i _mm256_set_epi8(char e31, char e30, char e29, char e28, char e27, char e26, char e25, char e24, char e23,
                        char e22, char e21, char e20, char e19, char e18, char e17, char e16, char e15, char e14,
                        char e13, char e12, char e11, char e10, char e9, char e8, char e7, char e6, char e5, char e4,
                        char e3, char e2, char e1, char e0);

__m256i _mm256_add_epi64(__m256i a, __m256i b);


#endif


#endif //AVX2_CONVERSION_H
