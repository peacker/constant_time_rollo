#include <stdio.h>
#include <stdint.h>
#include <ctype.h>

#include "hexhelper.h"

static int x_hexdigit(int ch) {
    if (ch >= '0' && ch <= '9')
        return ch - '0';
    if (ch >= 'A' && ch <= 'F')
        return ch - 'A' + 10;
    if (ch >= 'a' && ch <= 'f')
        return ch - 'a' + 10;
    return -1;
}

size_t x_readhex(void *buf, size_t n, const char *str) {
    size_t i, j;
    int ch1, ch2;

    j = 0;
    i = 0;
    while (i < n) {
        ch1 = str[j++];
        if (ch1 == 0)
            break;
        if (isspace(ch1) || ch1 == ',' || ch1 == '\\')  // skip between bytes
            continue;

        if ((ch1 = x_hexdigit(ch1)) < 0)
            break;
        if ((ch2 = x_hexdigit(str[j++])) < 0)
            break;
        ((uint8_t *) buf)[i++] = (ch1 << 4) ^ ch2;
    }

    return i;
}

// dump hex

void x_hexdump(const void *dat, size_t n, const char *lab) {
    size_t i;

    printf("%6s = ", lab);

    for (i = 0; i < n; i++) {
        if ((i & 0x1F) == 0) {
            if (i > 0)
                printf("\n[%04X] = ", (unsigned) i);
        }
        printf("%02X", ((const uint8_t *) dat)[i]);
    }
    printf("\n");
}

