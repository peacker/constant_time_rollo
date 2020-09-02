#ifndef HEXHELPER_H
#define HEXHELPER_H

#include <stddef.h>

size_t x_readhex(void *buf, size_t n, const char *str);

void x_hexdump(const void *dat, size_t n, const char *lab);


#endif
