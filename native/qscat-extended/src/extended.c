#include <float.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

size_t qscat_ld_size(void) { return sizeof(long double); }
unsigned qscat_ld_bits(void) { return LDBL_MANT_DIG; }
unsigned qscat_ld_exponent(void) { return LDBL_MAX_EXP; }
int qscat_ld_supported(void) {
#if defined(__x86_64__) && defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return sizeof(long double) == 16 && LDBL_MANT_DIG == 64 && LDBL_MAX_EXP == 16384;
#else
    return 0;
#endif
}

/* Only called after the format handshake. Avoid unaligned loads and padding. */
static long double load(const unsigned char *p) {
    long double value;
    memcpy(&value, p, sizeof(value));
    return value;
}
static void store(unsigned char *p, long double value) {
    memset(p, 0, sizeof(value));
#if LDBL_MANT_DIG == 64 && defined(__x86_64__) && defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    memcpy(p, &value, 10);
#else
    memcpy(p, &value, sizeof(value));
#endif
}

/* Rust checks every product and byte length before entering this function. */
int qscat_ld_residual(const unsigned char *a_bytes, const unsigned char *b,
                     const unsigned char *x, unsigned char *out,
                     size_t m, size_t k, size_t n) {
    if (!qscat_ld_supported()) return -2;
    const size_t s = sizeof(long double);
    long double *a = malloc(m * k * s);
    long double *xt = malloc(n * k * s);
    if (!a || !xt) { free(a); free(xt); return -1; }
    memcpy(a, a_bytes, m * k * s);
    for (size_t j = 0; j < n; ++j)
        for (size_t l = 0; l < k; ++l)
            xt[j * k + l] = load(x + (l * n + j) * s);
    for (size_t i = 0; i < m; ++i) {
        size_t j = 0;
        for (; j + 3 < n; j += 4) {
            long double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
            for (size_t l = 0; l < k; ++l) {
                long double v = a[i * k + l];
                s0 += v * xt[j * k + l];
                s1 += v * xt[(j + 1) * k + l];
                s2 += v * xt[(j + 2) * k + l];
                s3 += v * xt[(j + 3) * k + l];
            }
            store(out + (i * n + j) * s, load(b + (i * n + j) * s) - s0);
            store(out + (i * n + j + 1) * s, load(b + (i * n + j + 1) * s) - s1);
            store(out + (i * n + j + 2) * s, load(b + (i * n + j + 2) * s) - s2);
            store(out + (i * n + j + 3) * s, load(b + (i * n + j + 3) * s) - s3);
        }
        for (; j < n; ++j) {
            long double sum = 0;
            for (size_t l = 0; l < k; ++l)
                sum += a[i * k + l] * xt[j * k + l];
            store(out + (i * n + j) * s, load(b + (i * n + j) * s) - sum);
        }
    }
    free(a); free(xt);
    return 0;
}

void qscat_ld_fixture(unsigned char *out, size_t entries) {
    for (size_t i = 0; i < entries; ++i)
        store(out + i * sizeof(long double), ((long double)(i % 17) - 8) / 8);
}
