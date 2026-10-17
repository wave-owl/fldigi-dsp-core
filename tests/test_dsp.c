/*
 * fldigi-dsp-core tests: every primitive against golden data produced by fldigi
 * 4.2.13's own headers (tests/golden/make_golden.cxx). Comparisons are bit-exact.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "fldigi_dsp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures;
static const char *golden_dir = "tests/golden";

#define CHECK(cond, ...)                                                                 \
    do {                                                                                 \
        if (!(cond)) {                                                                   \
            failures++;                                                                  \
            fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__);                         \
            fprintf(stderr, __VA_ARGS__);                                                \
            fputc('\n', stderr);                                                         \
        }                                                                                \
    } while (0)

/* Same generator as make_golden.cxx. */
static uint32_t lcg_state;
static double lcg_next(void)
{
    lcg_state = lcg_state * 1664525u + 1013904223u;
    return ((double)(int32_t)lcg_state) / 2147483648.0;
}

static void *read_file(const char *name, size_t *len)
{
    char path[512];
    FILE *f;
    long n;
    void *buf;

    snprintf(path, sizeof path, "%s/%s", golden_dir, name);
    f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "cannot open %s\n", path);
        exit(2);
    }
    fseek(f, 0, SEEK_END);
    n = ftell(f);
    fseek(f, 0, SEEK_SET);
    buf = malloc((size_t)n + 1);
    if (!buf || fread(buf, 1, (size_t)n, f) != (size_t)n) {
        fprintf(stderr, "cannot read %s\n", path);
        exit(2);
    }
    fclose(f);
    *len = (size_t)n;
    return buf;
}

static int same_doubles(const double *a, const double *b, size_t n)
{
    return memcmp(a, b, n * sizeof(double)) == 0;
}

static void test_version(void)
{
    CHECK(strcmp(fdsp_version(), "fldigi-dsp-core " FDSP_VERSION_STRING
                                 " contract=dsp-c1 abi=1") == 0, "version %s", fdsp_version());
}

static void test_fft(void)
{
    static const size_t sizes[] = {4, 8, 16, 64, 256, 512, 1024, 2048, 4096, 8192};
    fdsp_r2fft bad;
    size_t s, i, len;
    char name[64];

    CHECK(fdsp_r2fft_init(&bad, 2) < 0, "size 2 accepted");
    CHECK(fdsp_r2fft_init(&bad, 96) < 0, "size 96 accepted");
    CHECK(bad.bitrev == NULL && bad.twiddle == NULL && bad.size == 0, "failed init not zeroed");
    fdsp_r2fft_free(&bad);

    for (s = 0; s < sizeof sizes / sizeof sizes[0]; s++) {
        size_t n = sizes[s];
        fdsp_r2fft fft;
        fdsp_cmpx *buf, *s0, *s1;
        double *tw, *out, *sep, *join;
        uint32_t *br;

        CHECK(fdsp_r2fft_init(&fft, n) == 0, "init %zu", n);
        snprintf(name, sizeof name, "fft%zu_twiddle.f64", n);
        tw = (double *)read_file(name, &len);
        /* The twiddles come from the platform's cos()/sin(), which may differ from the
         * golden machine's libm (or its sincos) in the last bit: compare within 2^-52,
         * then run the transforms on the golden twiddles, which must match bit for bit. */
        CHECK(len == 2 * n * sizeof(double), "twiddle size %zu", n);
        for (i = 0; i < 2 * n; i++) {
            double d = tw[i] - ((const double *)fft.twiddle)[i];
            CHECK(d <= 2.3e-16 && d >= -2.3e-16, "twiddle %zu[%zu] off by %g", n, i, d);
        }
        memcpy(fft.twiddle, tw, 2 * n * sizeof(double));
        snprintf(name, sizeof name, "fft%zu_bitrev.u32", n);
        br = (uint32_t *)read_file(name, &len);
        for (i = 0; i < n; i++)
            CHECK(br[i] == fft.bitrev[i], "bitrev %zu[%zu]", n, i);

        buf = (fdsp_cmpx *)malloc(n * sizeof *buf);
        s0 = (fdsp_cmpx *)malloc(n / 2 * sizeof *s0);
        s1 = (fdsp_cmpx *)malloc(n / 2 * sizeof *s1);
        lcg_state = (uint32_t)n;
        for (i = 0; i < n; i++) {
            buf[i].re = lcg_next();
            buf[i].im = lcg_next();
        }
        fdsp_r2fft_process(&fft, buf);
        snprintf(name, sizeof name, "fft%zu_out.f64", n);
        out = (double *)read_file(name, &len);
        CHECK(len == 2 * n * sizeof(double) && same_doubles(out, (const double *)buf, 2 * n),
              "fft output %zu", n);

        fdsp_r2fft_separ_two_reals(&fft, buf, s0, s1);
        snprintf(name, sizeof name, "fft%zu_separ.f64", n);
        sep = (double *)read_file(name, &len);
        CHECK(same_doubles(sep, (const double *)s0, n) &&
              same_doubles(sep + n, (const double *)s1, n), "separ %zu", n);

        fdsp_r2fft_join_two_reals(&fft, s0, s1, buf);
        fdsp_r2fft_process(&fft, buf);
        snprintf(name, sizeof name, "fft%zu_join.f64", n);
        join = (double *)read_file(name, &len);
        CHECK(same_doubles(join, (const double *)buf, 2 * n), "join %zu", n);

        free(tw);
        free(br);
        free(out);
        free(sep);
        free(join);
        free(buf);
        free(s0);
        free(s1);
        fdsp_r2fft_free(&fft);
        fdsp_r2fft_free(&fft);
    }
}

static void test_fft_math(void)
{
    /* Independent of fldigi: a single tone lands in its bin with the kernel e^{-j...}. */
    fdsp_r2fft fft;
    fdsp_cmpx x[64];
    size_t i;
    fdsp_r2fft_init(&fft, 64);
    for (i = 0; i < 64; i++) {
        x[i].re = fft.twiddle[(5 * i) % 64].re;  /* e^{+j 2 pi 5 i / 64} */
        x[i].im = fft.twiddle[(5 * i) % 64].im;
    }
    fdsp_r2fft_process(&fft, x);
    for (i = 0; i < 64; i++) {
        double e = x[i].re * x[i].re + x[i].im * x[i].im;
        if (i == 5)
            CHECK(e > 4095.0 && e < 4097.0, "tone bin energy %g", e);
        else
            CHECK(e < 1e-18, "leakage bin %zu %g", i, e);
    }
    fdsp_r2fft_free(&fft);
}

static void test_fht(void)
{
    size_t n, i, len;
    char name[64];
    for (n = 2; n <= 128; n *= 2) {
        double a[128], b[128], *g;
        lcg_state = (uint32_t)(1000 + n);
        for (i = 0; i < n; i++)
            a[i] = b[i] = lcg_next();
        fdsp_fht(a, n);
        fdsp_ifht(b, n);
        snprintf(name, sizeof name, "fht%zu.f64", n);
        g = (double *)read_file(name, &len);
        CHECK(len == 2 * n * sizeof(double) && same_doubles(g, a, n) &&
              same_doubles(g + n, b, n), "fht %zu", n);
        /* FHT(IFHT(x)) = n x for these power-of-two lengths (independent check). */
        fdsp_fht(b, n);
        lcg_state = (uint32_t)(1000 + n);
        for (i = 0; i < n; i++) {
            double ref = lcg_next() * (double)n, d = b[i] - ref;
            CHECK(d < 1e-9 && d > -1e-9, "fht roundtrip %zu[%zu]", n, i);
        }
        free(g);
    }
    for (n = 32; n <= 64; n *= 2) {
        int8_t *g, v[64];
        size_t c;
        snprintf(name, sizeof name, "ifht_i8_%zu.i8", n);
        g = (int8_t *)read_file(name, &len);
        CHECK(len == 2 * n * n, "ifht_i8 size %zu", n);
        for (c = 0; c < 2 * n; c++) {
            memset(v, 0, sizeof v);
            v[c % n] = (int8_t)(c < n ? 1 : -1);
            fdsp_ifht_i8(v, n);
            CHECK(memcmp(v, g + c * n, n) == 0, "ifht_i8 %zu unit %zu", n, c);
        }
        free(g);
    }
}

static void test_gray(void)
{
    size_t len;
    uint8_t *g = (uint8_t *)read_file("gray.u8", &len);
    int i;
    CHECK(len == 512, "gray size");
    for (i = 0; i < 256; i++) {
        CHECK(fdsp_gray_code((uint8_t)i) == g[i], "gray %d", i);
        CHECK(fdsp_binary_code((uint8_t)i) == g[256 + i], "binary %d", i);
        CHECK(fdsp_binary_code(fdsp_gray_code((uint8_t)i)) == i, "gray roundtrip %d", i);
    }
    free(g);
}

static void test_lowpass3(void)
{
    static const double weights[] = {1.0 / 4, 1.0 / 32, 1.0 / 100};
    static const double feedbacks[] = {0.1, 0.5};
    size_t len, w, fb, k = 0;
    double *g = (double *)read_file("lowpass3.f64", &len);
    CHECK(len == 6 * 2002 * sizeof(double), "lowpass3 size");
    for (w = 0; w < 3; w++) {
        for (fb = 0; fb < 2; fb++) {
            fdsp_lowpass3 f;
            int i, ok = 1;
            fdsp_lowpass3_set(&f, 0.0);
            lcg_state = (uint32_t)(77 + 10 * w + fb);
            for (i = 0; i < 2000; i++) {
                fdsp_lowpass3_process(&f, lcg_next(), weights[w], feedbacks[fb]);
                ok &= memcmp(&f.output, &g[k++], sizeof(double)) == 0;
            }
            ok &= memcmp(&f.out1, &g[k++], sizeof(double)) == 0;
            ok &= memcmp(&f.out2, &g[k++], sizeof(double)) == 0;
            CHECK(ok, "lowpass3 weight %zu feedback %zu", w, fb);
        }
    }
    /* A constant input settles at the input (independent check). */
    {
        fdsp_lowpass3 f;
        int i;
        fdsp_lowpass3_set(&f, 0.0);
        for (i = 0; i < 400; i++)
            fdsp_lowpass3_process(&f, 3.0, 0.25, 0.1);
        CHECK(f.output > 2.999999 && f.output < 3.000001, "lowpass3 settle %g", f.output);
    }
    free(g);
}

int main(int argc, char **argv)
{
    if (argc > 1)
        golden_dir = argv[1];
    test_version();
    test_fft();
    test_fft_math();
    test_fht();
    test_gray();
    test_lowpass3();
    if (failures) {
        fprintf(stderr, "%d failure(s)\n", failures);
        return 1;
    }
    printf("fldigi-dsp-core: all tests passed (%s)\n", fdsp_version());
    return 0;
}
