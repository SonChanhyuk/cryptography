/*
 * SonChanhyuk
 * thscksgur0903@gmail.com
 */
#include "miller_rabin.h"

/*
 * mod_add() - computes a+b mod m
 * a와 b가 m보다 작다는 가정하에서 a+b >= m이면 결과에서 m을 빼줘야 하므로
 * 오버플로가 발생하지 않도록 a-(m-b)를 계산하고, 그렇지 않으면 그냥 a+b를 계산하면 된다.
 * a+b >= m을 검사하는 과정에서 오버플로가 발생할 수 있으므로 a >= m-b를 검사한다.
 */
uint64_t mod_add(uint64_t a, uint64_t b, uint64_t m)
{
    a = a % m;
    b = b % m;
    return (a >= m-b ? a-(m-b) : a+b);
}

/*
 * mod_sub() - computes a-b mod m
 * 만일 a < b이면 결과가 음수가 되므로 m을 더해서 양수로 만든다.
 */
uint64_t mod_sub(uint64_t a, uint64_t b, uint64_t m)
{
    a = a % m;
    b = b % m;
    return (a < b ? (m-b)+a : a-b);
}

/*
 * mod_mul() - computes a*b mod m
 * a*b에서 오버플로가 발생할 수 있기 때문에 덧셈을 사용하여 빠르게 계산할 수 있는
 * "double addition" 알고리즘을 사용한다. 그 알고리즘은 다음과 같다.
 *     r = 0;
 *     while (b > 0) {
 *         if (b & 1)
 *             r = mod_add(r, a, m);
 *         b = b >> 1;
 *         a = mod_add(a, a, m);
 *     }
 */
uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m)
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
 * a^b에서 오버플로가 발생할 수 있기 때문에 곱셈을 사용하여 빠르게 계산할 수 있는
 * "square multiplication" 알고리즘을 사용한다. 그 알고리즘은 다음과 같다.
 *     r = 1;
 *     while (b > 0) {
 *         if (b & 1)
 *             r = mod_mul(r, a, m);
 *         b = b >> 1;
 *         a = mod_mul(a, a, m);
 *     }
 */
uint64_t mod_pow(uint64_t a, uint64_t b, uint64_t m)
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
 * Miller-Rabin Primality Testing against small sets of bases
 *
 * if n < 2^64,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, and 37.
 *
 * if n < 3,317,044,064,679,887,385,961,981,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, and 41.
 */
const uint64_t a[BASELEN] = {2,3,5,7,11,13,17,19,23,29,31,37};

/*
 * miller_rabin() - Miller-Rabin Primality Test (deterministic version)
 *
 * n > 3, an odd integer to be tested for primality
 * It returns PRIME if n is prime, COMPOSITE otherwise.
 */
int miller_rabin(uint64_t n)
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