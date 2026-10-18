/*
 * fldigi-dsp-core: fldigi 4.2.13's FIR filter, moving average, sliding DFT and
 * overlap-add FFT filter (src/filters/filters.cxx, src/filters/fftfilt.cxx,
 * src/include/filters.h, src/include/fftfilt.h, src/include/misc.h).
 *
 * Copyright (C) 2026 Gerolf Ziegenhain
 * Original C++: (c) 2006-2008 Dave Freese W1HKJ and the fldigi authors.
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * std::complex<double> arithmetic is written out as fldigi's compiler evaluates it
 * ((a + jb)(c + jd) = (ac - bd) + j(ad + bc) with the products in separate statements,
 * |z| = hypot). Expressions that are one statement in fldigi stay one statement here, so
 * a build with floating-point contraction (-ffp-contract=on, as fldigi's own) fuses the
 * same operations as fldigi's binary. The sliding DFT keeps
 * fldigi's long double damping constant K1, so its results follow the platform's
 * long double like fldigi's do.
 */
#include "fldigi_dsp.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static fdsp_cmpx cmul(fdsp_cmpx a, fdsp_cmpx b)
{
    /* std::complex<double> operator* as libc++/libstdc++ evaluate it: four products in
     * separate statements, so a compiler that contracts a * b + c within one expression
     * (clang and gcc by default) does not fuse them into FMAs here either */
    fdsp_cmpx r;
    double ac = a.re * b.re, bd = a.im * b.im, ad = a.re * b.im, bc = a.im * b.re;
    r.re = ac - bd;
    r.im = ad + bc;
    return r;
}

double fdsp_decayavg(double average, double input, int weight)
{
    /* include/misc.h:59-63 */
    if (weight <= 1)
        return input;
    return ((input - average) / (double)weight) + average;
}

int fdsp_parity(unsigned long w)
{
    /* hweight32 + parity, misc/misc.cxx:34-42, 67-70 */
    w = (w & 0x55555555UL) + ((w >> 1) & 0x55555555UL);
    w = (w & 0x33333333UL) + ((w >> 2) & 0x33333333UL);
    w = (w & 0x0F0F0F0FUL) + ((w >> 4) & 0x0F0F0F0FUL);
    w = (w & 0x00FF00FFUL) + ((w >> 8) & 0x00FF00FFUL);
    w = (w & 0x0000FFFFUL) + ((w >> 16) & 0x0000FFFFUL);
    return (int)(w & 1);
}

double fdsp_blackman(double x)
{
    /* include/misc.h:66-69: 0.42 - 0.50 * cos(2 * M_PI * x) + 0.08 * cos(4 * M_PI * x),
     * with the two multiply-adds fldigi's reference build (clang -O2, FP contraction on)
     * fuses, written out so the result does not depend on this build's contraction. */
    double c2 = cos(2 * M_PI * x);
    double c4 = cos(4 * M_PI * x);
    return fma(0.08, c4, fma(-0.50, c2, 0.42));
}

void fdsp_sincos(double x, double *s, double *c)
{
#if defined(__APPLE__)
    __sincos(x, s, c);
#else
    *s = sin(x);
    *c = cos(x);
#endif
}

/* ---------------------------------------------------------------------------------- */
/* fftfilt, filters/fftfilt.cxx:55-237, include/fftfilt.h:47-55                                       */
/* ---------------------------------------------------------------------------------- */

static double ff_fsinc(double fc, int i, int len)
{
    /* include/fftfilt.h:47-50 */
    return (i == len / 2) ? 2.0 * fc :
           sin(2 * M_PI * fc * (i - len / 2)) / (M_PI * (i - len / 2));
}

static double ff_blackman(int i, int len)
{
    /* include/fftfilt.h:51-55 */
    return (0.42 -
            0.50 * cos(2.0 * M_PI * i / len) +
            0.08 * cos(4.0 * M_PI * i / len));
}

void fdsp_fftfilt_free(fdsp_fftfilt *f)
{
    /* ~fftfilt, filters/fftfilt.cxx:112-122 */
    if (!f)
        return;
    fdsp_gfft_free(&f->fft);
    free(f->filter);
    free(f->timedata);
    free(f->freqdata);
    free(f->output);
    free(f->ovlbuf);
    memset(f, 0, sizeof *f);
}

static void ff_clear(fdsp_fftfilt *f)
{
    /* clear_filter, filters/fftfilt.cxx:55-67 */
    memset(f->filter, 0, (size_t)f->flen * sizeof(fdsp_cmpx));
    memset(f->timedata, 0, (size_t)f->flen * sizeof(fdsp_cmpx));
    memset(f->freqdata, 0, (size_t)f->flen * sizeof(fdsp_cmpx));
    memset(f->output, 0, (size_t)f->flen * sizeof(fdsp_cmpx));
    memset(f->ovlbuf, 0, (size_t)f->flen2 * sizeof(fdsp_cmpx));
    f->inptr = 0;
}

int fdsp_fftfilt_init(fdsp_fftfilt *f, double f1, double f2, int len)
{
    if (!f)
        return -1;
    memset(f, 0, sizeof *f);
    if (len < 32)
        return -1;
    /* init_filter, filters/fftfilt.cxx:69-80 */
    f->flen = len;
    f->flen2 = len >> 1;
    if (fdsp_gfft_init(&f->fft, (size_t)len) < 0)
        return -1;
    f->filter = (fdsp_cmpx *)calloc((size_t)len, sizeof(fdsp_cmpx));
    f->timedata = (fdsp_cmpx *)calloc((size_t)len, sizeof(fdsp_cmpx));
    f->freqdata = (fdsp_cmpx *)calloc((size_t)len, sizeof(fdsp_cmpx));
    f->output = (fdsp_cmpx *)calloc((size_t)len, sizeof(fdsp_cmpx));
    f->ovlbuf = (fdsp_cmpx *)calloc((size_t)f->flen2, sizeof(fdsp_cmpx));
    if (!f->filter || !f->timedata || !f->freqdata || !f->output || !f->ovlbuf) {
        fdsp_fftfilt_free(f);
        return -1;
    }
    fdsp_fftfilt_create(f, f1, f2);
    return 0;
}

void fdsp_fftfilt_create(fdsp_fftfilt *f, double f1, double f2)
{
    /* create_filter, filters/fftfilt.cxx:124-196. ht lives in filter[] before the transform
     * (fldigi copies ht to filter with memcpy, :151). */
    int i, b_lowpass, b_highpass;
    double scale = 0, mag;
    fdsp_cmpx *ht = f->filter;

    ff_clear(f);
    b_lowpass = (f2 != 0);
    b_highpass = (f1 != 0);
    for (i = 0; i < f->flen2; i++) {
        ht[i].re = 0;
        ht[i].im = 0;
        if (b_lowpass)
            ht[i].re += ff_fsinc(f2, i, f->flen2);
        if (b_highpass)
            ht[i].re -= ff_fsinc(f1, i, f->flen2);
    }
    if (b_highpass && f2 < f1)
        ht[f->flen2 / 2].re += 1;
    for (i = 0; i < f->flen2; i++)
        ht[i].re *= ff_blackman(i, f->flen2);

    fdsp_gfft_forward(&f->fft, f->filter);

    for (i = 0; i < f->flen2; i++) {
        mag = hypot(f->filter[i].re, f->filter[i].im);
        if (mag > scale)
            scale = mag;
    }
    if (scale != 0) {
        for (i = 0; i < f->flen; i++) {
            f->filter[i].re /= scale;
            f->filter[i].im /= scale;
        }
    }
    f->pass = 1;
}

int fdsp_fftfilt_run(fdsp_fftfilt *f, fdsp_cmpx in, const fdsp_cmpx **out)
{
    /* run, filters/fftfilt.cxx:202-237 */
    int i;
    f->timedata[f->inptr++] = in;
    if (f->inptr < f->flen2)
        return 0;
    if (f->pass)
        --f->pass;
    memcpy(f->freqdata, f->timedata, (size_t)f->flen * sizeof(fdsp_cmpx));
    fdsp_gfft_forward(&f->fft, f->freqdata);
    for (i = 0; i < f->flen; i++)
        f->freqdata[i] = cmul(f->freqdata[i], f->filter[i]);
    fdsp_gfft_inverse(&f->fft, f->freqdata);
    for (i = 0; i < f->flen2; i++) {
        f->output[i].re = f->ovlbuf[i].re + f->freqdata[i].re;
        f->output[i].im = f->ovlbuf[i].im + f->freqdata[i].im;
        f->ovlbuf[i] = f->freqdata[i + f->flen2];
    }
    f->inptr = 0;
    if (f->pass)
        return 0;
    *out = f->output;
    return f->flen2;
}

static double ff_sinc(double x)
{
    /* sinc, include/misc.h:43-46 */
    return (fabs(x) < 1e-10) ? 1.0 : (sin(M_PI * x) / (M_PI * x));
}

void fdsp_fftfilt_rtty_filter(fdsp_fftfilt *f, double fc)
{
    /* rtty_filter, filters/fftfilt.cxx:243-290: raised-cosine (cos^2) response, amplitude
     * equalised by 1 / sinc, with a linear phase of -i * pi / 2. Overwrites filter[0 ..
     * flen2 - 1] and filter[flen2 + 1 .. flen - 1]; filter[flen2] keeps the value of the
     * previous shape (fldigi's quirk: the loop never reaches it). The cos/sin pairs of one
     * argument are one sincos call each, as fldigi's reference build evaluates them. */
    int i;
    double dht, x, s, c;
    fc *= 1.4;
    for (i = 0; i < f->flen2; ++i) {
        x = (double)i / (double)(f->flen2);
        dht = x <= 0 ? 1.0 :
              x > 2.0 * fc ? 0.0 :
              cos((M_PI * x) / (fc * 4.0));
        dht *= dht;
        dht /= ff_sinc(2.0 * i * fc);
        fdsp_sincos((double)i * -0.5 * M_PI, &s, &c);
        f->filter[i].re = dht * c;
        f->filter[i].im = dht * s;
        fdsp_sincos((double)i * +0.5 * M_PI, &s, &c);
        f->filter[(f->flen - i) % f->flen].re = dht * c;
        f->filter[(f->flen - i) % f->flen].im = dht * s;
    }
    f->pass = 1;
}

/* ---------------------------------------------------------------------------------- */
/* C_FIR_filter, filters/filters.cxx:46-187, include/filters.h:40-110                                 */
/* ---------------------------------------------------------------------------------- */

static double fir_sinc(double x)
{
    /* include/filters.h:59-64 */
    if (fabs(x) < 1e-10)
        return 1.0;
    return sin(M_PI * x) / (M_PI * x);
}

static double fir_cosc(double x)
{
    /* include/filters.h:65-70 */
    if (fabs(x) < 1e-10)
        return 0.0;
    return (1.0 - cos(M_PI * x)) / (M_PI * x);
}

static double fir_hamming(double x)
{
    /* include/filters.h:71-73 */
    return 0.54 - 0.46 * cos(2 * M_PI * x);
}

static double fir_mac(const double *a, const double *b, unsigned int size)
{
    /* mac, include/filters.h:74-91: four partial sums. fldigi's reference build (Apple
     * clang -O2 on arm64, FP contraction on) vectorises the loop over groups of four
     * iterations without fusing product and sum (the in-order reduction is done on
     * extracted lanes), and runs the remaining iterations and the tail as fused
     * multiply-adds. That evaluation is written out here, so every compiler and
     * optimisation level gives the reference result. */
    double sum = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, p;
    unsigned int iters = size > 3 ? (size - 4) / 4 + 1 : 0;
    unsigned int vec = size >= 16 ? iters & ~3u : 0;
    unsigned int it;
    for (it = 0; it < vec; it++, a += 4, b += 4) {
        p = a[0] * b[0];
        sum += p;
        p = a[1] * b[1];
        sum2 += p;
        p = a[2] * b[2];
        sum3 += p;
        p = a[3] * b[3];
        sum4 += p;
    }
    for (; it < iters; it++, a += 4, b += 4) {
        sum = fma(a[0], b[0], sum);
        sum2 = fma(a[1], b[1], sum2);
        sum3 = fma(a[2], b[2], sum3);
        sum4 = fma(a[3], b[3], sum4);
    }
    for (size -= 4 * iters; size; --size)
        sum = fma(*a++, *b++, sum);
    return sum + sum2 + sum3 + sum4;
}

void fdsp_fir_bp_taps(double *fir, int len, int hilbert, double f1, double f2)
{
    /* bp_FIR, filters/filters.cxx:92-119 */
    int i;
    double t, h, x;
    for (i = 0; i < len; i++) {
        t = i - (len - 1.0) / 2.0;
        h = i * (1.0 / (len - 1.0));
        if (!hilbert) {
            x = (2 * f2 * fir_sinc(2 * f2 * t) -
                 2 * f1 * fir_sinc(2 * f1 * t)) * fir_hamming(h);
        } else {
            x = (2 * f2 * fir_cosc(2 * f2 * t) -
                 2 * f1 * fir_cosc(2 * f1 * t)) * fir_hamming(h);
            x = -x;
        }
        fir[i] = x;
    }
}

void fdsp_fir_free(fdsp_fir *f)
{
    if (!f)
        return;
    free(f->ifilter);
    free(f->qfilter);
    f->ifilter = f->qfilter = NULL;
    f->length = 0;
}

static int fir_init(fdsp_fir *f, int len, int dec, const double *itaps, const double *qtaps)
{
    /* init, filters/filters.cxx:58-84 */
    if (!f || len < 1 || len >= FDSP_FIR_BUFLEN || dec < 1)
        return -1;
    memset(f, 0, sizeof *f);
    f->length = len;
    f->decimateratio = dec;
    f->ifilter = (double *)malloc((size_t)len * sizeof(double));
    f->qfilter = (double *)malloc((size_t)len * sizeof(double));
    if (!f->ifilter || !f->qfilter) {
        fdsp_fir_free(f);
        return -1;
    }
    memcpy(f->ifilter, itaps, (size_t)len * sizeof(double));
    memcpy(f->qfilter, qtaps, (size_t)len * sizeof(double));
    f->pointer = len;
    f->counter = 0;
    return 0;
}

int fdsp_fir_init(fdsp_fir *f, int len, int dec, const double *itaps, const double *qtaps)
{
    /* init, filters/filters.cxx:58-84 with the caller's taps; a missing Q branch is zero */
    int rc;
    double *zero;
    if (!itaps || len < 1 || len >= FDSP_FIR_BUFLEN)
        return -1;
    if (qtaps)
        return fir_init(f, len, dec, itaps, qtaps);
    zero = (double *)calloc((size_t)len, sizeof(double));
    if (!zero)
        return -1;
    rc = fir_init(f, len, dec, itaps, zero);
    free(zero);
    return rc;
}

int fdsp_fir_init_hilbert(fdsp_fir *f, int len, int dec)
{
    /* init_hilbert, filters/filters.cxx:152-158 */
    double *fi, *fq;
    int rc;
    if (len < 2 || len >= FDSP_FIR_BUFLEN)
        return -1;
    fi = (double *)malloc((size_t)len * sizeof(double));
    fq = (double *)malloc((size_t)len * sizeof(double));
    if (!fi || !fq) {
        free(fi);
        free(fq);
        return -1;
    }
    fdsp_fir_bp_taps(fi, len, 0, 0.05, 0.45);
    fdsp_fir_bp_taps(fq, len, 1, 0.05, 0.45);
    rc = fir_init(f, len, dec, fi, fq);
    free(fi);
    free(fq);
    return rc;
}

int fdsp_fir_init_lowpass(fdsp_fir *f, int len, int dec, double freq)
{
    /* init_lowpass, filters/filters.cxx:128-133 */
    return fdsp_fir_init_bandpass(f, len, dec, 0.0, freq);
}

int fdsp_fir_init_bandpass(fdsp_fir *f, int len, int dec, double f1, double f2)
{
    /* init_bandpass, filters/filters.cxx:142-146 */
    double *fi;
    int rc;
    if (len < 2 || len >= FDSP_FIR_BUFLEN)
        return -1;
    fi = (double *)malloc((size_t)len * sizeof(double));
    if (!fi)
        return -1;
    fdsp_fir_bp_taps(fi, len, 0, f1, f2);
    rc = fir_init(f, len, dec, fi, fi);
    free(fi);
    return rc;
}

int fdsp_fir_run(fdsp_fir *f, fdsp_cmpx in, fdsp_cmpx *out)
{
    /* run, filters/filters.cxx:168-187 */
    f->ibuffer[f->pointer] = in.re;
    f->qbuffer[f->pointer] = in.im;
    f->counter++;
    if (f->counter == f->decimateratio) {
        out->re = fir_mac(&f->ibuffer[f->pointer - f->length], f->ifilter,
                          (unsigned int)f->length);
        out->im = fir_mac(&f->qbuffer[f->pointer - f->length], f->qfilter,
                          (unsigned int)f->length);
    }
    f->pointer++;
    if (f->pointer == FDSP_FIR_BUFLEN) {
        memmove(f->ibuffer, f->ibuffer + FDSP_FIR_BUFLEN - f->length,
                (size_t)f->length * sizeof(double));
        memmove(f->qbuffer, f->qbuffer + FDSP_FIR_BUFLEN - f->length,
                (size_t)f->length * sizeof(double));
        f->pointer = f->length;
    }
    if (f->counter == f->decimateratio) {
        f->counter = 0;
        return 1;
    }
    return 0;
}

int fdsp_fir_irun(fdsp_fir *f, double in, double *out)
{
    /* Irun, filters/filters.cxx:193-217: the dot product covers the length samples before
     * the new one, as run() does */
    double *iptr = f->ibuffer + f->pointer;
    f->pointer++;
    f->counter++;
    *iptr = in;
    if (f->counter == f->decimateratio)
        *out = fir_mac(iptr - f->length, f->ifilter, (unsigned int)f->length);
    if (f->pointer == FDSP_FIR_BUFLEN) {
        iptr = f->ibuffer + FDSP_FIR_BUFLEN - f->length;
        memcpy(f->ibuffer, iptr, (size_t)f->length * sizeof(double));
        f->pointer = f->length;
    }
    if (f->counter == f->decimateratio) {
        f->counter = 0;
        return 1;
    }
    return 0;
}

fdsp_cmpx fdsp_cmac(const fdsp_cmpx *a, const fdsp_cmpx *b, int ptr, int len)
{
    /* cmac, include/complex.h:33-41 */
    fdsp_cmpx z, p;
    int i;
    z.re = z.im = 0.0;
    ptr %= len;
    for (i = 0; i < len; i++) {
        p = cmul(a[i], b[ptr]);
        z.re += p.re;
        z.im += p.im;
        ptr = (ptr + 1) % len;
    }
    return z;
}

/* ---------------------------------------------------------------------------------- */
/* Cmovavg, filters/filters.cxx:261-307                                                       */
/* ---------------------------------------------------------------------------------- */

int fdsp_movavg_init(fdsp_movavg *m, int len)
{
    if (!m || len < 1 || len > FDSP_MOVAVG_MAX)
        return -1;
    memset(m, 0, sizeof *m);
    m->len = len;
    m->empty = 1;
    return 0;
}

void fdsp_movavg_reset(fdsp_movavg *m)
{
    m->empty = 1;
}

double fdsp_movavg_run(fdsp_movavg *m, double a)
{
    int i;
    if (m->empty) {
        m->empty = 0;
        m->out = 0;
        for (i = 0; i < m->len; i++) {
            m->in[i] = a;
            m->out += a;
        }
        m->pint = 0;
        return a;
    }
    m->out = m->out - m->in[m->pint] + a;
    m->in[m->pint] = a;
    if (++m->pint >= m->len)
        m->pint = 0;
    return m->out / m->len;
}

/* ---------------------------------------------------------------------------------- */
/* Cmovavg with a heap window, filters/filters.cxx:261-307 (setLength 294-302)          */
/* ---------------------------------------------------------------------------------- */

int fdsp_movavgh_init(fdsp_movavgh *m, int len)
{
    if (!m)
        return -1;
    memset(m, 0, sizeof *m);
    if (len < 1)
        return -1;
    m->in = (double *)malloc((size_t)len * sizeof(double));
    if (!m->in)
        return -1;
    m->len = len;
    m->empty = 1;
    return 0;
}

void fdsp_movavgh_free(fdsp_movavgh *m)
{
    if (!m)
        return;
    free(m->in);
    m->in = NULL;
    m->len = 0;
}

void fdsp_movavgh_reset(fdsp_movavgh *m)
{
    m->empty = 1;
}

double fdsp_movavgh_run(fdsp_movavgh *m, double a)
{
    int i;
    if (!m->in)
        return a;
    if (m->empty) {
        m->empty = 0;
        m->out = 0;
        for (i = 0; i < m->len; i++) {
            m->in[i] = a;
            m->out += a;
        }
        m->pint = 0;
        return a;
    }
    m->out = m->out - m->in[m->pint] + a;
    m->in[m->pint] = a;
    if (++m->pint >= m->len)
        m->pint = 0;
    return m->out / m->len;
}

/* ---------------------------------------------------------------------------------- */
/* sfft, filters/filters.cxx:402-465, include/filters.h:136-151                                       */
/* ---------------------------------------------------------------------------------- */

#define SFFT_K1 0.99999999999L   /* include/filters.h:137 */

void fdsp_sfft_free(fdsp_sfft *s)
{
    if (!s)
        return;
    free(s->vrot);
    free(s->bins);
    free(s->delay);
    memset(s, 0, sizeof *s);
}

int fdsp_sfft_init(fdsp_sfft *s, int len, int first, int last)
{
    /* sfft::sfft, filters/filters.cxx:407-424 */
    int i;
    double phi = 0.0, tau;
    if (!s)
        return -1;
    memset(s, 0, sizeof *s);
    if (len < 1 || first < 0 || last < first || last > len)
        return -1;
    s->vrot = (fdsp_cmpx *)calloc((size_t)len, sizeof(fdsp_cmpx));
    s->bins = (fdsp_cmpx *)calloc((size_t)len, sizeof(fdsp_cmpx));
    s->delay = (fdsp_cmpx *)calloc((size_t)len, sizeof(fdsp_cmpx));
    if (!s->vrot || !s->bins || !s->delay) {
        fdsp_sfft_free(s);
        return -1;
    }
    s->fftlen = len;
    s->first = first;
    s->last = last;
    s->ptr = 0;
    tau = 2.0 * M_PI / len;
    s->k2 = 1.0;
    for (i = 0; i < len; i++) {
        double sn, cs;
        fdsp_sincos(phi, &sn, &cs);
        s->vrot[i].re = (double)(SFFT_K1 * cs);
        s->vrot[i].im = (double)(SFFT_K1 * sn);
        phi += tau;
        s->k2 = (double)(s->k2 * SFFT_K1);
    }
    s->count = 0;
    return 0;
}

void fdsp_sfft_reset(fdsp_sfft *s)
{
    memset(s->delay, 0, (size_t)s->fftlen * sizeof(fdsp_cmpx));
    memset(s->bins, 0, (size_t)s->fftlen * sizeof(fdsp_cmpx));
    s->count = 0;
}

void fdsp_sfft_run(fdsp_sfft *s, fdsp_cmpx input, fdsp_cmpx *result, int stride)
{
    /* sfft::run, filters/filters.cxx:447-465 */
    int i;
    fdsp_cmpx *de = &s->delay[s->ptr];
    fdsp_cmpx z, a, b;
    /* fused in fldigi's reference build (fmsub), explicit here for every build */
    z.re = fma(-s->k2, de->re, input.re);
    z.im = fma(-s->k2, de->im, input.im);
    *de = input;
    ++s->ptr;
    if (s->ptr >= s->fftlen)
        s->ptr = 0;
    for (i = s->first; i < s->last; i++, result += stride) {
        a = cmul(s->bins[i], s->vrot[i]);
        b = cmul(z, s->vrot[i]);
        s->bins[i].re = a.re + b.re;
        s->bins[i].im = a.im + b.im;
        *result = s->bins[i];
    }
    if (s->count < s->fftlen)
        s->count++;
}
