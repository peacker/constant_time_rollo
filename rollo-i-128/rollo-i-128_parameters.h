#ifndef ROLLO_I_128_PARAMETERS_H
#define ROLLO_I_128_PARAMETERS_H


// =============================================
// ==               TYPEDEF                   ==
// =============================================

typedef __uint128_t uint128_t;

// =============================================
// ==               CODE PARAMETERS           ==
// =============================================

#define ROLLO_I_FINITE_FIELD_DEGREE 67u // m
#define ROLLO_I_BF_ELEMENT_BYTE_SIZE 9u //  (67 + 7) / 8
#define ROLLO_I_BF_MASK (uint128_t) ((((uint128_t) 0x7u) << 64u) ^(uint128_t) 0xffffffffffffffffu)
#define ROLLO_I_BF_MASK_LOW 0xffffffffffffffffu
#define ROLLO_I_BF_MASK_HIGH 0x7u
#define ROLLO_I_CODE_LENGTH 83u // n, P = X^83 + X^7 + X^4 + X^2 + 1
#define ROLLO_I_KEY_VECTORS_RANK_WEIGHT 8u // d
#define ROLLO_I_ERROR_VECTORS_RANK_WEIGHT 7u // r
#define ROLLO_I_SHARED_SECRET_BYTE_SIZE 32u
//#define ROLLO_I_SHARED_SECRET_BYTE_SIZE 64u

// =============================================
// ==               API PARAMETERS            ==
// =============================================


#endif //ROLLO_I_128_PARAMETERS_H
