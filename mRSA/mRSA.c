/*
 * SonChanhyuk
 * thscksgur0903@gmail.com
 */
#ifdef __linux__
#include <bsd/stdlib.h>
#elif __APPLE__
#include <stdlib.h>
#else
#include <stdlib.h>
#endif
#include "mRSA.h"

/*
 * mod_add() - computes a + b mod m
 */
static uint64_t mod_add(uint64_t a, uint64_t b, uint64_t m)
{
    a = a % m;
    b = b % m;
    return (a >= m-b ? a-(m-b) : a+b);
}

/*
 * mod_mul() - computes a * b mod m
 */
static uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m)
{
    a = a % m;
    b = b % m;
    uint64_t r = 0;
    while (b > 0) {
        if (b & 1) r = mod_add(r,a,m);
        b = b >> 1;
        a = mod_add(a,a,m);   // a = a << 1
    }
    return r;
}

/*
 * mod_pow() - computes a^b mod m
 */
static uint64_t mod_pow(uint64_t a, uint64_t b, uint64_t m)
{
    a = a % m;
    uint64_t r = 1;
    while (b > 0) {
        if (b & 1) r = mod_mul(r,a,m);
        b = b >> 1;
        a = mod_mul(a,a,m);
    }
    return r;
}

/*
 * gcd() - Euclidean algorithm
 */
static uint64_t gcd(uint64_t a, uint64_t b)
{
    uint64_t tmp;
    while (b != 0)
    {
        /* a = b, b = a % b */
        tmp = a % b;
        a = b;
        b = tmp;
    }
    return a;
}

/*
 * mul_inv() - computes multiplicative inverse a^-1 mod m
 * It returns 0 if no inverse exist.
 */
static uint64_t mul_inv(uint64_t a, uint64_t m)
{
    uint64_t d0 = a % m, d1 = m;
    uint64_t x0 = 1, x1 = 0, q, tmp, mul;

    while(d1 > 1)
    {
        q = d0 / d1;

        tmp = d0 % d1;
        d0 = d1;
        d1 = tmp;

        // q * x1
        mul = 0;
        tmp = x1;
        while (q > 0) 
        {
            if(q & 1) 
            {
                mul = mul + x1;
                if(mul < x1) mul += ((UINT64_MAX % m) + 1);
                mul = (mul % m);
            } 
            q = q >> 1;    
            x1 = ((x1 << 1) + ((x1 >> 63) & 1 ? (UINT64_MAX % m) + 1 : 0));
            x1 = (x1 % m);
        }

        x1 = x0 - mul;
        if(x0 < mul) x1 += m;
        x0 = tmp;
    }

    if(d1 == 1)
        return x1;

    return 0;
}

/*
 * Miller-Rabin Primality Testing against small sets of bases
 *
 * if n < 2^64,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, and 37.
 *
 * if n < 3317044064679887385961981,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, and 41.
 */
static const uint64_t a[BASELEN] = {2,3,5,7,11,13,17,19,23,29,31,37};

/*
 * miller_rabin() - Miller-Rabin Primality Test (deterministic version)
 *
 * n > 3, an odd integer to be tested for primality
 * It returns 1 if n is prime, 0 otherwise.
 */
static int miller_rabin(uint64_t n)
{
    if (n == 1 || (((n & 1) == 0) && (n != 2))) return COMPOSITE;

    // find k, q, k > 0, q = odd, n-1 = q * 2^k
    uint64_t q=n-1,k=0,aq; 
    while ((q & 1) == 0) {
        q = q >> 1;
        k++;
    }

    for (int i=0;i<BASELEN;i++) {
        // n이 배열 a에 속하면 소수
        if (a[i] == n) return PRIME;

        // aq = a ^ q mod n
        aq = mod_pow(a[i],q,n);
        if (aq == 1) continue;    

        for (int j=0;j<k;j++) {
            // aq^2^j == n-1 
            if (mod_pow(aq,1<<j,n) == n-1) break;
            if (j == k-1) return COMPOSITE;
        }
    }   
    // 위의 조건을 전부 통과하면 소수
    return PRIME;
}

/*
 * mRSA_generate_key() - generates mini RSA keys e, d and n
 *
 * Carmichael's totient function Lambda(n) is used.
 */
void mRSA_generate_key(uint64_t *e, uint64_t *d, uint64_t *n)
{
    uint64_t p, q, pi;
    *n = 0;

    // n = p*q의 길이가 64비트이고 소수 p와 q의 길이가 30이상이 되는 랜덤 p와q를 생성
    while(*n < MINIMUM_N){
        p = q = 0;
        while(p < 0x20000000 || miller_rabin(p) == COMPOSITE){
            p = arc4random();
            p = p | 1;
        }
        while(q < 0x20000000 || miller_rabin(q) == COMPOSITE){
            q = arc4random();
            q = q | 1;
        }
        *n = p * q;
    }

    // pi대신 카마이클함수 lambda를 사용, lambda(n) = lcm(p-1,q-1)
    pi = ((p-1)*(q-1))/gcd(p-1,q-1);

    // 랜덤 e, d를 설정한다.
    p = arc4random_uniform(pi);
    while(1){
        p = arc4random_uniform(pi);
        q = mul_inv(p,pi);
        if (gcd(p,pi) == 1 && p > 1 && q != 0) break;
    }
    
    *e = p;
    *d = q;
}

/*
 * mRSA_cipher() - compute m^k mod n
 *
 * If data >= n then returns 1 (error), otherwise 0 (success).
 */
int mRSA_cipher(uint64_t *m, uint64_t k, uint64_t n)
{
    // m mod n을 통해 m이 매우 작아지는 것을 방지
    if (*m > n) return 1;
    *m = mod_pow(*m,k,n);
    return 0;
}