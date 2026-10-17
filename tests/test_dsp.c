/*
 * fldigi-dsp-core tests: every primitive against golden data produced by fldigi
 * 4.2.13's own headers (tests/golden/make_golden.cxx). Comparisons are bit-exact.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "fldigi_dsp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

/* ------------------------------------------------------------------------------------
 * 0.2.0 primitives against make_golden_fldigi.cxx (fldigi's own sources).
 * Tables computed with libm (cosine table, filter shapes, FIR taps, sliding-DFT
 * rotations) are compared within a few ulp, then the golden tables are injected so the
 * processing itself is compared bit for bit on every platform.
 * ------------------------------------------------------------------------------------ */

static int close_doubles(const double *a, const double *b, size_t n, double tol)
{
    size_t i;
    for (i = 0; i < n; i++) {
        double d = a[i] - b[i];
        if (d > tol || d < -tol)
            return 0;
    }
    return 1;
}

static void test_gfft(void)
{
    size_t n, i, len;
    char name[64];
    fdsp_gfft bad;
    CHECK(fdsp_gfft_init(&bad, 8) < 0, "gfft size 8 accepted");
    CHECK(fdsp_gfft_init(&bad, 100) < 0, "gfft size 100 accepted");
    for (n = 16; n <= 8192; n *= 2) {
        fdsp_gfft fft;
        fdsp_cmpx *buf = (fdsp_cmpx *)malloc(n * sizeof *buf);
        double *gold, *utbl;
        CHECK(fdsp_gfft_init(&fft, n) == 0, "gfft init %zu", n);
        snprintf(name, sizeof name, "gfft%zu_utbl.f64", n);
        utbl = (double *)read_file(name, &len);
        CHECK(len == (n / 4 + 1) * sizeof(double) &&
              close_doubles(utbl, fft.utbl, n / 4 + 1, 2.3e-16), "gfft utbl %zu", n);
        memcpy(fft.utbl, utbl, len);
        lcg_state = (uint32_t)(5000 + n);
        for (i = 0; i < n; i++) {
            buf[i].re = lcg_next();
            buf[i].im = lcg_next();
        }
        snprintf(name, sizeof name, "gfft%zu.f64", n);
        gold = (double *)read_file(name, &len);
        fdsp_gfft_forward(&fft, buf);
        CHECK(same_doubles(gold, (const double *)buf, 2 * n), "gfft forward %zu", n);
        fdsp_gfft_inverse(&fft, buf);
        CHECK(same_doubles(gold + 2 * n, (const double *)buf, 2 * n), "gfft inverse %zu", n);
        free(gold);
        free(utbl);
        free(buf);
        fdsp_gfft_free(&fft);
        fdsp_gfft_free(&fft);
    }
}

static void test_fftfilt(void)
{
    static const double cases[3][2] = {
        {(2000.0 - 0.5 * 2.0 * 281.25) / 8000.0, (2000.0 + 0.5 * 2.0 * 281.25) / 8000.0},
        {(1500.0 - 0.5 * 2.0 * 193.798828125) / 11025.0,
         (1500.0 + 0.5 * 2.0 * 193.798828125) / 11025.0},
        {0.0, 0.1},
    };
    int c, i, k;
    size_t len, pos;
    char name[64];
    double *utbl = (double *)read_file("gfft1024_utbl.f64", &len);
    for (c = 0; c < 3; c++) {
        fdsp_fftfilt f;
        double *shape, *gold;
        size_t glen;
        CHECK(fdsp_fftfilt_init(&f, cases[c][0], cases[c][1], 1024) == 0, "fftfilt init");
        snprintf(name, sizeof name, "fftfilt%d_shape.f64", c);
        shape = (double *)read_file(name, &len);
        CHECK(len == 2048 * sizeof(double) &&
              close_doubles(shape, (const double *)f.filter, 2048, 1e-14), "fftfilt shape %d", c);
        memcpy(f.filter, shape, len);
        memcpy(f.fft.utbl, utbl, 257 * sizeof(double));
        snprintf(name, sizeof name, "fftfilt%d_out.f64", c);
        gold = (double *)read_file(name, &glen);
        lcg_state = (uint32_t)(77 + c);
        pos = 0;
        for (i = 0; i < 6000; i++) {
            const fdsp_cmpx *zp = NULL;
            fdsp_cmpx in;
            int n;
            in.re = lcg_next();
            in.im = lcg_next();
            n = fdsp_fftfilt_run(&f, in, &zp);
            for (k = 0; k < n; k++, pos += 2) {
                if (pos + 1 >= glen / sizeof(double) ||
                    gold[pos] != zp[k].re || gold[pos + 1] != zp[k].im) {
                    CHECK(0, "fftfilt %d sample %zu", c, pos / 2);
                    i = 6000;
                    break;
                }
            }
        }
        CHECK(pos == glen / sizeof(double), "fftfilt %d output count %zu", c, pos);
        free(shape);
        free(gold);
        fdsp_fftfilt_free(&f);
    }
    free(utbl);
}

static void test_fir(void)
{
    int c, i;
    size_t len, glen, pos;
    char name[64];
    for (c = 0; c < 3; c++) {
        fdsp_fir f;
        int flen = c == 0 ? 37 : c == 1 ? 127 : 65, rc;
        double *taps, *gold;
        if (c == 0)
            rc = fdsp_fir_init_hilbert(&f, flen, 1);
        else if (c == 1)
            rc = fdsp_fir_init_lowpass(&f, flen, 1, 0.0625);
        else
            rc = fdsp_fir_init_bandpass(&f, flen, 2, 0.1, 0.2);
        CHECK(rc == 0, "fir init %d", c);
        snprintf(name, sizeof name, "fir%d_taps.f64", c);
        taps = (double *)read_file(name, &len);
        CHECK(close_doubles(taps, f.ifilter, (size_t)flen, 1e-16) &&
              close_doubles(taps + flen, f.qfilter, (size_t)flen, 1e-16), "fir taps %d", c);
        memcpy(f.ifilter, taps, (size_t)flen * sizeof(double));
        memcpy(f.qfilter, taps + flen, (size_t)flen * sizeof(double));
        snprintf(name, sizeof name, "fir%d_out.f64", c);
        gold = (double *)read_file(name, &glen);
        lcg_state = (uint32_t)(300 + c);
        pos = 0;
        for (i = 0; i < 9000; i++) {
            fdsp_cmpx in, out;
            in.re = lcg_next();
            in.im = lcg_next();
            if (fdsp_fir_run(&f, in, &out)) {
                if (gold[pos] != out.re || gold[pos + 1] != out.im) {
                    CHECK(0, "fir %d output %zu", c, pos / 2);
                    break;
                }
                pos += 2;
            }
        }
        CHECK(pos == glen / sizeof(double), "fir %d output count", c);
        free(taps);
        free(gold);
        fdsp_fir_free(&f);
    }
}

static void test_fir_irun(void)
{
    int c, i;
    size_t len, glen, pos;
    char name[64];
    for (c = 0; c < 2; c++) {
        fdsp_fir f;
        int flen = c == 0 ? 256 : 64, dec = c == 0 ? 1 : 2;
        double *taps, *gold;
        snprintf(name, sizeof name, "firirun%d_taps.f64", c);
        taps = (double *)read_file(name, &len);
        CHECK(len == (size_t)flen * sizeof(double), "fir irun taps %d", c);
        CHECK(fdsp_fir_init(&f, flen, dec, taps, NULL) == 0, "fir irun init %d", c);
        CHECK(f.qfilter[0] == 0.0 && f.qfilter[flen - 1] == 0.0, "fir irun zero Q %d", c);
        snprintf(name, sizeof name, "firirun%d_out.f64", c);
        gold = (double *)read_file(name, &glen);
        lcg_state = (uint32_t)(320 + c);
        pos = 0;
        for (i = 0; i < 9000; i++) {
            double out = 0.0, x = fabs(lcg_next());
            if (fdsp_fir_irun(&f, x, &out)) {
                if (pos >= glen / sizeof(double) || gold[pos] != out) {
                    CHECK(0, "fir irun %d output %zu", c, pos);
                    break;
                }
                pos++;
            }
        }
        CHECK(pos == glen / sizeof(double), "fir irun %d output count", c);
        free(taps);
        free(gold);
        fdsp_fir_free(&f);
    }
    {
        fdsp_fir f;
        double t[3] = {1, 2, 3};
        CHECK(fdsp_fir_init(&f, 3, 1, NULL, NULL) == -1, "fir init without taps");
        CHECK(fdsp_fir_init(&f, 0, 1, t, NULL) == -1, "fir init length 0");
    }
}

static void test_cmac(void)
{
    const int lens[2] = {256, 11};
    size_t glen, pos = 0;
    double *gold = (double *)read_file("cmac.f64", &glen);
    int c, i, p;
    for (c = 0; c < 2; c++) {
        int len = lens[c];
        fdsp_cmpx *a = (fdsp_cmpx *)malloc((size_t)len * sizeof *a);
        fdsp_cmpx *b = (fdsp_cmpx *)malloc((size_t)len * sizeof *b);
        lcg_state = (uint32_t)(340 + c);
        for (i = 0; i < len; i++) {
            a[i].re = lcg_next();
            a[i].im = lcg_next();
        }
        for (i = 0; i < len; i++) {
            b[i].re = lcg_next();
            b[i].im = lcg_next();
        }
        for (p = 0; p < len + 3; p++, pos += 2) {
            fdsp_cmpx z = fdsp_cmac(a, b, p, len);
            if (pos + 1 >= glen / sizeof(double) || gold[pos] != z.re || gold[pos + 1] != z.im) {
                CHECK(0, "cmac len %d ptr %d", len, p);
                break;
            }
        }
        free(a);
        free(b);
    }
    CHECK(pos == glen / sizeof(double), "cmac count");
    free(gold);
}

static void test_movavg(void)
{
    size_t len, pos = 0;
    int l, i;
    double *gold = (double *)read_file("movavg.f64", &len);
    for (l = 8; l <= 16; l += 8) {
        fdsp_movavg m;
        CHECK(fdsp_movavg_init(&m, l) == 0, "movavg init");
        lcg_state = (uint32_t)(400 + l);
        for (i = 0; i < 500; i++)
            CHECK(fdsp_movavg_run(&m, 1000.0 * lcg_next()) == gold[pos++], "movavg %d/%d", l, i);
        fdsp_movavg_reset(&m);
        CHECK(fdsp_movavg_run(&m, 3.25) == gold[pos++], "movavg reset");
        CHECK(fdsp_movavg_run(&m, -1.5) == gold[pos++], "movavg after reset");
    }
    CHECK(pos == len / sizeof(double), "movavg count");
    free(gold);
}

static void test_sfft(void)
{
    static const int cases[3][3] = {{1024, 184, 248}, {512, 91, 121}, {64, 3, 20}};
    int c, i, k;
    size_t len, glen, pos;
    char name[64];
    for (c = 0; c < 3; c++) {
        int n = cases[c][0], first = cases[c][1], last = cases[c][2];
        fdsp_sfft s;
        fdsp_cmpx *res = (fdsp_cmpx *)calloc((size_t)(3 * (last - first)), sizeof *res);
        double *vr, *gold;
        CHECK(fdsp_sfft_init(&s, n, first, last) == 0, "sfft init");
        snprintf(name, sizeof name, "sfft%d_vrot.f64", c);
        vr = (double *)read_file(name, &len);
        CHECK(len == (size_t)(2 * n + 1) * sizeof(double) &&
              close_doubles(vr, (const double *)s.vrot, (size_t)(2 * n), 2.3e-16) &&
              close_doubles(vr + 2 * n, &s.k2, 1, 1e-15), "sfft rotations %d", c);
        memcpy(s.vrot, vr, (size_t)(2 * n) * sizeof(double));
        s.k2 = vr[2 * n];
        snprintf(name, sizeof name, "sfft%d.f64", c);
        gold = (double *)read_file(name, &glen);
        lcg_state = (uint32_t)(500 + c);
        pos = 0;
        for (i = 0; i < 3000; i++) {
            fdsp_cmpx in;
            in.re = lcg_next();
            in.im = lcg_next();
            fdsp_sfft_run(&s, in, res + 1, 3);
            if (i % 7 == 0 || i > 2990) {
                for (k = 0; k < last - first; k++, pos += 2) {
                    if (gold[pos] != res[1 + 3 * k].re || gold[pos + 1] != res[1 + 3 * k].im) {
                        CHECK(0, "sfft %d step %d bin %d", c, i, k);
                        i = 3000;
                        break;
                    }
                }
            }
        }
        CHECK(pos == glen / sizeof(double), "sfft %d count", c);
        fdsp_sfft_reset(&s);
        CHECK(s.count == 0 && s.bins[first].re == 0.0, "sfft reset");
        free(vr);
        free(gold);
        free(res);
        fdsp_sfft_free(&s);
    }
}

static void test_conv(void)
{
    static const int ks[2][3] = {{7, 0x6d, 0x4f}, {15, 044735, 063057}};
    int c, i;
    size_t len;
    char name[64];
    fdsp_viterbi bad;
    CHECK(fdsp_viterbi_init(&bad, 16, 1, 1) < 0, "viterbi k 16 accepted");
    for (c = 0; c < 2; c++) {
        fdsp_encoder enc;
        fdsp_viterbi v1, v8;
        int32_t *codes, *dec1, *met1, *dec8;
        unsigned char *sym;
        int k = ks[c][0];
        CHECK(fdsp_encoder_init(&enc, k, ks[c][1], ks[c][2]) == 0, "encoder init");
        snprintf(name, sizeof name, "conv%d_enc.i32", k);
        codes = (int32_t *)read_file(name, &len);
        lcg_state = (uint32_t)(600 + c);
        for (i = 0; i < 4000; i++) {
            int bit = lcg_next() > 0.0;
            CHECK(fdsp_encoder_encode(&enc, bit) == codes[i], "encode k=%d %d", k, i);
            lcg_next();
            lcg_next();
        }
        snprintf(name, sizeof name, "conv%d_sym.u8", k);
        sym = (unsigned char *)read_file(name, &len);
        snprintf(name, sizeof name, "conv%d_dec1.i32", k);
        dec1 = (int32_t *)read_file(name, &len);
        snprintf(name, sizeof name, "conv%d_met1.i32", k);
        met1 = (int32_t *)read_file(name, &len);
        snprintf(name, sizeof name, "conv%d_dec8.i32", k);
        dec8 = (int32_t *)read_file(name, &len);
        CHECK(fdsp_viterbi_init(&v1, k, ks[c][1], ks[c][2]) == 0, "viterbi init");
        CHECK(fdsp_viterbi_set_traceback(&v1, c == 0 ? 45 : 180) == 0, "traceback");
        CHECK(fdsp_viterbi_set_chunksize(&v1, 1) == 0, "chunksize");
        CHECK(fdsp_viterbi_init(&v8, k, ks[c][1], ks[c][2]) == 0, "viterbi init 8");
        for (i = 0; i < 4000; i++) {
            int met = 0, r = fdsp_viterbi_decode(&v1, sym + 2 * i, &met);
            CHECK(r == dec1[i] && met == met1[i], "viterbi k=%d chunk 1 pair %d", k, i);
            met = 0;
            r = fdsp_viterbi_decode(&v8, sym + 2 * i, &met);
            CHECK(r == dec8[2 * i] && met == dec8[2 * i + 1], "viterbi k=%d chunk 8 pair %d", k, i);
        }
        free(codes);
        free(sym);
        free(dec1);
        free(met1);
        free(dec8);
        fdsp_encoder_free(&enc);
        fdsp_viterbi_free(&v1);
        fdsp_viterbi_free(&v8);
    }
}

static void test_viterbi_long_run(void)
{
    /* Independent of fldigi: 9 million clean pairs (past the int overflow point of
     * fldigi's metrics with chunk size 1) still decode, and UBSan stays quiet. */
    fdsp_encoder enc;
    fdsp_viterbi v;
    long i, errors = 0;
    int delay = 45, hist[64] = {0};
    fdsp_encoder_init(&enc, 7, 0x6d, 0x4f);
    fdsp_viterbi_init(&v, 7, 0x6d, 0x4f);
    fdsp_viterbi_set_traceback(&v, delay);
    fdsp_viterbi_set_chunksize(&v, 1);
    lcg_state = 4242;
    for (i = 0; i < 9000000L; i++) {
        int bit = (lcg_state >> 13) & 1, d, r;
        unsigned char s[2];
        lcg_next();
        d = fdsp_encoder_encode(&enc, bit);
        s[0] = (d & 1) ? 255 : 0;
        s[1] = (d & 2) ? 255 : 0;
        r = fdsp_viterbi_decode(&v, s, NULL);
        if (i > delay + 1 && r != hist[(i - delay) & 63])
            errors++;
        hist[i & 63] = bit;
    }
    CHECK(errors == 0, "long run: %ld bit errors", errors);
    fdsp_encoder_free(&enc);
    fdsp_viterbi_free(&v);
}

static void test_interleave(void)
{
    static const int cases[4][2] = {{4, 10}, {4, 50}, {4, 4}, {5, 6}};
    size_t len, pos = 0;
    int c, i, k;
    unsigned char *gold = (unsigned char *)read_file("interleave.u8", &len);
    for (c = 0; c < 4; c++) {
        int size = cases[c][0], depth = cases[c][1];
        fdsp_interleave fwd, rev, bfwd;
        CHECK(fdsp_interleave_init(&fwd, size, depth, FDSP_INTERLEAVE_FWD) == 0, "il init");
        fdsp_interleave_init(&rev, size, depth, FDSP_INTERLEAVE_REV);
        fdsp_interleave_init(&bfwd, size, depth, FDSP_INTERLEAVE_FWD);
        lcg_state = (uint32_t)(700 + c);
        for (i = 0; i < 400; i++) {
            unsigned char s[8], r[8];
            unsigned int bits;
            for (k = 0; k < size; k++)
                s[k] = r[k] = (unsigned char)(256.0 * (0.5 + 0.5 * lcg_next()));
            fdsp_interleave_symbols(&fwd, s);
            fdsp_interleave_symbols(&rev, r);
            bits = (unsigned int)(1024.0 * (0.5 + 0.5 * lcg_next())) & ((1u << size) - 1);
            fdsp_interleave_bits(&bfwd, &bits);
            CHECK(memcmp(gold + pos, s, (size_t)size) == 0, "il fwd %d/%d", c, i);
            pos += (size_t)size;
            CHECK(memcmp(gold + pos, r, (size_t)size) == 0, "il rev %d/%d", c, i);
            pos += (size_t)size;
            CHECK(gold[pos] == (unsigned char)bits, "il bits %d/%d", c, i);
            pos++;
        }
        fdsp_interleave_free(&fwd);
        fdsp_interleave_free(&rev);
        fdsp_interleave_free(&bfwd);
    }
    CHECK(pos == len, "interleave count");
    free(gold);
}

static void test_varicode(void)
{
    size_t len, pos = 0;
    int c;
    unsigned int s;
    char *enc = (char *)read_file("mfsk_varienc.txt", &len);
    int32_t *dec;
    for (c = -1; c <= 256; c++) {
        CHECK(strcmp(enc + pos, fdsp_mfsk_varienc(c)) == 0, "varienc %d", c);
        pos += strlen(enc + pos) + 1;
    }
    CHECK(pos == len, "varienc count");
    dec = (int32_t *)read_file("mfsk_varidec.i32", &len);
    for (s = 0; s < 4096; s++)
        CHECK(fdsp_mfsk_varidec(s) == dec[s], "varidec %u", s);
    free(dec);
    dec = (int32_t *)read_file("parity.i32", &len);
    for (s = 0; s < 10000; s++)
        CHECK(fdsp_parity((unsigned long)s * 7) == dec[s], "parity %u", s * 7);
    CHECK(fdsp_decayavg(1.0, 5.0, 4) == 2.0 && fdsp_decayavg(1.0, 5.0, 1) == 5.0, "decayavg");
    free(enc);
    free(dec);
}

/* Rational resampler (not an fldigi port): unit gain, a tone keeps its frequency and
 * amplitude, an out-of-band tone is suppressed, the output does not depend on the
 * chunking, and max_out bounds every call. */
static double tone_power(const double *x, size_t n, double f, double rate)
{
    double re = 0, im = 0;
    size_t i;
    for (i = 0; i < n; i++) {
        re += x[i] * cos(2 * M_PI * f * (double)i / rate);
        im += x[i] * sin(2 * M_PI * f * (double)i / rate);
    }
    return (re * re + im * im) / ((double)n * (double)n) * 4.0;
}

static void test_resampler(void)
{
    static const int rates[][2] = {{8000, 11025}, {11025, 8000}, {48000, 11025},
                                   {8000, 8000}, {12000, 8000}};
    size_t c;
    fdsp_resampler bad;
    CHECK(fdsp_resampler_init(&bad, 0, 8000, 32) == -1, "resampler rate 0");
    CHECK(fdsp_resampler_init(&bad, 8000, 11025, 2) == -1, "resampler taps 2");
    CHECK(fdsp_resampler_init(&bad, 8000, 8191 * 7, 32) == -1, "resampler factor");
    for (c = 0; c < sizeof rates / sizeof rates[0]; c++) {
        int in = rates[c][0], out = rates[c][1];
        size_t n = (size_t)in * 2, i, k1 = 0, k2 = 0, off;
        double *x = (double *)malloc(n * sizeof(double));
        double *y1, *y2, p, pa;
        fdsp_resampler r1, r2;
        CHECK(fdsp_resampler_init(&r1, in, out, 32) == 0, "resampler init %d/%d", in, out);
        CHECK(fdsp_resampler_init(&r2, in, out, 32) == 0, "resampler init %d/%d", in, out);
        y1 = (double *)malloc(fdsp_resampler_max_out(&r1, n) * sizeof(double));
        y2 = (double *)malloc(fdsp_resampler_max_out(&r1, n) * sizeof(double));
        /* 1500 Hz in band + 3.9 kHz below the lower Nyquist only when it is not aliased */
        for (i = 0; i < n; i++)
            x[i] = 0.5 * sin(2 * M_PI * 1500.0 * (double)i / in);
        k1 = fdsp_resampler_run(&r1, x, n, y1);
        CHECK(k1 <= fdsp_resampler_max_out(&r1, n), "resampler bound %d/%d", in, out);
        lcg_state = 7;
        for (off = 0; off < n;) {
            size_t m = 1 + (size_t)(fabs(lcg_next()) * 700), got;
            if (m > n - off)
                m = n - off;
            got = fdsp_resampler_run(&r2, x + off, m, y2 + k2);
            CHECK(got <= fdsp_resampler_max_out(&r2, m), "resampler chunk bound");
            k2 += got;
            off += m;
        }
        CHECK(k1 == k2 && same_doubles(y1, y2, k1), "resampler chunking %d/%d", in, out);
        CHECK(labs((long)k1 - (long)(n * (size_t)out / (size_t)in)) <= 1, "resampler count");
        /* skip the filter's settling time, then measure the tone */
        p = tone_power(y1 + out / 4, k1 - (size_t)out / 4, 1500.0, out);
        CHECK(fabs(10 * log10(p / 0.25)) < 0.05, "resampler gain %d/%d: %g", in, out, p);
        pa = tone_power(y1 + out / 4, k1 - (size_t)out / 4, 1700.0, out);
        CHECK(pa < p * 1e-4, "resampler leakage %d/%d", in, out);
        if (in > out) {
            /* a tone above the output Nyquist must not alias into the band */
            double fa = 0.5 * out + 900.0;
            for (i = 0; i < n; i++)
                x[i] = 0.5 * sin(2 * M_PI * fa * (double)i / in);
            fdsp_resampler_reset(&r1);
            k1 = fdsp_resampler_run(&r1, x, n, y1);
            p = tone_power(y1 + out / 4, k1 - (size_t)out / 4, out - fa, out);
            CHECK(p < 0.25 * 1e-6, "resampler alias %d/%d: %g", in, out, p);
        }
        fdsp_resampler_free(&r1);
        fdsp_resampler_free(&r2);
        free(x);
        free(y1);
        free(y2);
    }
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
    test_gfft();
    test_fftfilt();
    test_fir();
    test_fir_irun();
    test_cmac();
    test_movavg();
    test_sfft();
    test_conv();
    test_viterbi_long_run();
    test_interleave();
    test_varicode();
    test_resampler();
    if (failures) {
        fprintf(stderr, "%d failure(s)\n", failures);
        return 1;
    }
    printf("fldigi-dsp-core: all tests passed (%s)\n", fdsp_version());
    return 0;
}
