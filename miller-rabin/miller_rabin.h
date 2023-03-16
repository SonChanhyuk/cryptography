/*
 * SonChanhyuk
 * thscksgur0903@gmail.com
 */
#ifndef _MILLER_RABIN_H_
#define _MILLER_RABIN_H_

#include <stdint.h>

#define BASELEN 12
#define PRIME 1
#define COMPOSITE 0

uint64_t mod_add(uint64_t a, uint64_t b, uint64_t m);
uint64_t mod_sub(uint64_t a, uint64_t b, uint64_t m);
uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m);
uint64_t mod_pow(uint64_t a, uint64_t b, uint64_t m);
int miller_rabin(uint64_t n);

#endif
