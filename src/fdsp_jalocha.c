/*
 * fldigi-dsp-core: Pawel Jalocha's generic DSP helpers as used by fldigi 4.2.13
 * (src/include/jalocha/pj_fft.h, pj_fht.h, pj_gray.h, pj_lowpass3.h).
 *
 * Copyright (C) 2026 Gerolf Ziegenhain
 * Original C++: (c) 1999-2004 Pawel Jalocha SP9VRC, (c) 2019 David Freese W1HKJ.
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Every function keeps the operation order of the C++ templates instantiated with
 * double, so the results are bit-identical to fldigi's code compiled without
 * floating-point contraction.
 */
#include "fldigi_dsp.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const char *fdsp_version(void)
{
    return "fldigi-dsp-core " FDSP_VERSION_STRING " contract=dsp-c1 abi=1";
}

/* ---------------------------------------------------------------------------------- */
/* r2FFT, pj_fft.h:63-317                                                             */
/* ---------------------------------------------------------------------------------- */

void fdsp_r2fft_free(fdsp_r2fft *fft)
{
    /* r2FFT::Free, pj_fft.h:84-89 */
    if (!fft)
        return;
    free(fft->bitrev);
    free(fft->twiddle);
    fft->bitrev = NULL;
    fft->twiddle = NULL;
    fft->size = 0;
}

int fdsp_r2fft_init(fdsp_r2fft *fft, size_t size)
{
    size_t idx, ridx, mask, rmask, size4, rest;
    double phase;

    if (!fft)
        return -1;
    fft->size = 0;
    fft->bitrev = NULL;
    fft->twiddle = NULL;
    /* r2FFT::Preset, pj_fft.h:92-100: a power of two, at least 4 */
    if (size < 4)
        return -1;
    for (rest = size; (rest & 1) == 0; rest >>= 1)
        ;
    if (rest != 1)
        return -1;
    fft->bitrev = (size_t *)malloc(size * sizeof(size_t));
    fft->twiddle = (fdsp_cmpx *)malloc(size * sizeof(fdsp_cmpx));
    if (!fft->bitrev || !fft->twiddle) {
        fdsp_r2fft_free(fft);
        return -1;
    }
    fft->size = size;
    /* pj_fft.h:105-113: first quarter by cos/sin (Cmpx::SetPhase, pj_cmpx.h:49-52),
     * the other three quarters by rotation by +90 degrees */
    size4 = size / 4;
    for (idx = 0; idx < size4; idx++) {
        phase = (2 * M_PI * idx) / size;
        fft->twiddle[idx].re = 1.0 * cos(phase);
        fft->twiddle[idx].im = 1.0 * sin(phase);
    }
    for (; idx < size; idx++) {
        fft->twiddle[idx].re = (-fft->twiddle[idx - size4].im);
        fft->twiddle[idx].im = fft->twiddle[idx - size4].re;
    }
    /* pj_fft.h:114-120: bit-reverse table */
    for (idx = 0; idx < size; idx++) {
        for (ridx = 0, mask = size / 2, rmask = 1; mask; mask >>= 1, rmask <<= 1) {
            if (idx & mask)
                ridx |= rmask;
        }
        fft->bitrev[idx] = ridx;
    }
    return 0;
}

/* Butterfly with the conjugate twiddle, pj_fft.h:277-285 */
static void fft_bf(fdsp_cmpx *x0, fdsp_cmpx *x1, const fdsp_cmpx *w)
{
    fdsp_cmpx x1w;
    x1w.re = x1->re * w->re + x1->im * w->im;
    x1w.im = (-x1->re * w->im) + x1->im * w->re;
    x1->re = x0->re - x1w.re;
    x1->im = x0->im - x1w.im;
    x0->re = x0->re + x1w.re;
    x0->im = x0->im + x1w.im;
}

/* 2-point FFT of the first pass, pj_fft.h:289-297 */
static void fft_2(fdsp_cmpx *x0, fdsp_cmpx *x1)
{
    fdsp_cmpx x1w;
    x1w.re = x1->re;
    x1w.im = x1->im;
    x1->re = x0->re - x1->re;
    x1->im = x0->im - x1->im;
    x0->re += x1w.re;
    x0->im += x1w.im;
}

void fdsp_r2fft_process(const fdsp_r2fft *fft, fdsp_cmpx *x)
{
    size_t idx, ridx, groups, groupsize2, group, bf, twid;
    size_t size = fft->size, size2 = fft->size / 2;
    fdsp_cmpx tmp;

    /* Scramble, pj_fft.h:129-140 */
    for (idx = 0; idx < size; idx++) {
        if ((ridx = fft->bitrev[idx]) > idx) {
            tmp = x[idx];
            x[idx] = x[ridx];
            x[ridx] = tmp;
        }
    }
    /* CoreProc, pj_fft.h:210-224 */
    for (bf = 0; bf < size; bf += 2)
        fft_2(&x[bf], &x[bf + 1]);
    for (groups = size2 / 2, groupsize2 = 2; groups; groups >>= 1, groupsize2 <<= 1) {
        for (group = 0, bf = 0; group < groups; group++, bf += groupsize2) {
            for (twid = 0; twid < size2; twid += groups, bf++)
                fft_bf(&x[bf], &x[bf + groupsize2], &fft->twiddle[twid]);
        }
    }
}

void fdsp_r2fft_separ_two_reals(const fdsp_r2fft *fft, const fdsp_cmpx *buf,
                                fdsp_cmpx *out0, fdsp_cmpx *out1)
{
    /* SeparTwoReals, pj_fft.h:161-178 */
    int idx, size = (int)fft->size, half = (int)(fft->size / 2);
    out0[0].re = buf[0].re;
    out1[0].re = buf[0].im;
    for (idx = 1; idx < half; idx++) {
        out0[idx].re = buf[idx].re + buf[size - idx].re;
        out0[idx].im = buf[idx].im - buf[size - idx].im;
        out1[idx].re = buf[idx].im + buf[size - idx].im;
        out1[idx].im = (-buf[idx].re) + buf[size - idx].re;
    }
    out0[0].im = buf[half].re;
    out1[0].im = buf[half].im;
}

void fdsp_r2fft_join_two_reals(const fdsp_r2fft *fft, const fdsp_cmpx *in0,
                               const fdsp_cmpx *in1, fdsp_cmpx *buf)
{
    /* JoinTwoReals, pj_fft.h:186-203 */
    int idx, size = (int)fft->size, half = (int)(fft->size / 2);
    buf[0].re = 2 * in0[0].re;
    buf[0].im = (2 * in1[0].re);
    for (idx = 1; idx < half; idx++) {
        buf[idx].re = in0[idx].re + in1[idx].im;
        buf[idx].im = (-in0[idx].im) + in1[idx].re;
        buf[size - idx].re = in0[idx].re - in1[idx].im;
        buf[size - idx].im = in0[idx].im + in1[idx].re;
    }
    buf[half].re = 2 * in0[0].im;
    buf[half].im = 2 * in1[0].im;
}

/* ---------------------------------------------------------------------------------- */
/* Fast Hadamard transform, pj_fht.h:11-48                                            */
/* ---------------------------------------------------------------------------------- */

void fdsp_fht(double *data, size_t len)
{
    /* FHT, pj_fht.h:11-27 */
    size_t step, ptr, ptr2;
    double bit1, bit2, nb1, nb2;
    for (step = 1; step < len; step *= 2) {
        for (ptr = 0; ptr < len; ptr += 2 * step) {
            for (ptr2 = ptr; (ptr2 - ptr) < step; ptr2 += 1) {
                bit1 = data[ptr2];
                bit2 = data[ptr2 + step];
                nb1 = bit2;
                nb1 += bit1;
                nb2 = bit2;
                nb2 -= bit1;
                data[ptr2] = nb1;
                data[ptr2 + step] = nb2;
            }
        }
    }
}

void fdsp_ifht(double *data, size_t len)
{
    /* IFHT, pj_fht.h:31-48 */
    size_t step, ptr, ptr2;
    double bit1, bit2, nb1, nb2;
    for (step = len / 2; step; step /= 2) {
        for (ptr = 0; ptr < len; ptr += 2 * step) {
            for (ptr2 = ptr; (ptr2 - ptr) < step; ptr2 += 1) {
                bit1 = data[ptr2];
                bit2 = data[ptr2 + step];
                nb1 = bit1;
                nb1 -= bit2;
                nb2 = bit1;
                nb2 += bit2;
                data[ptr2] = nb1;
                data[ptr2 + step] = nb2;
            }
        }
    }
}

void fdsp_ifht_i8(int8_t *data, size_t len)
{
    /* IFHT instantiated with int8_t (MFSK_Encoder::EncodeCharacter, pj_mfsk.h:1165) */
    size_t step, ptr, ptr2;
    int8_t bit1, bit2, nb1, nb2;
    for (step = len / 2; step; step /= 2) {
        for (ptr = 0; ptr < len; ptr += 2 * step) {
            for (ptr2 = ptr; (ptr2 - ptr) < step; ptr2 += 1) {
                bit1 = data[ptr2];
                bit2 = data[ptr2 + step];
                nb1 = bit1;
                nb1 = (int8_t)(nb1 - bit2);
                nb2 = bit1;
                nb2 = (int8_t)(nb2 + bit2);
                data[ptr2] = nb1;
                data[ptr2 + step] = nb2;
            }
        }
    }
}

/* ---------------------------------------------------------------------------------- */
/* Gray code, pj_gray.h:11-21                                                         */
/* ---------------------------------------------------------------------------------- */

uint8_t fdsp_gray_code(uint8_t binary)
{
    return (uint8_t)(binary ^ (binary >> 1));
}

uint8_t fdsp_binary_code(uint8_t gray)
{
    gray ^= (uint8_t)(gray >> 4);
    gray ^= (uint8_t)(gray >> 2);
    gray ^= (uint8_t)(gray >> 1);
    return gray;
}

/* ---------------------------------------------------------------------------------- */
/* LowPass3_Filter, pj_lowpass3.h:13-41                                               */
/* ---------------------------------------------------------------------------------- */

void fdsp_lowpass3_set(fdsp_lowpass3 *f, double level)
{
    f->out1 = level;
    f->out2 = level;
    f->output = level;
}

void fdsp_lowpass3_process(fdsp_lowpass3 *f, double input, double weight, double feedback)
{
    double d_i1, d_12, d_23;
    weight *= 2.0;
    d_i1 = input;
    d_i1 -= f->out1;
    d_12 = f->out1;
    d_12 -= f->out2;
    d_23 = f->out2;
    d_23 -= f->output;
    d_i1 *= weight;
    f->out1 += d_i1;
    d_12 *= weight;
    f->out2 += d_12;
    d_23 *= weight;
    f->output += d_23;
    d_23 *= feedback;
    f->out2 += d_23;
}
