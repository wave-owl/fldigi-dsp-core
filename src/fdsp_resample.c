/*
 * fldigi-dsp-core: rational polyphase resampler between modem and station rates.
 *
 * Copyright (C) 2026 Gerolf Ziegenhain
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Not a port: fldigi converts between the sound card rate and the modem rate with
 * libsamplerate (src/soundcard/sound.cxx, SRC_SINC_*), which the cores do not link. This
 * is the shared replacement for the cores whose modems run at 11025 Hz (MFSK11/22,
 * THOR/DominoEX 5/11/22/44/88) while the station delivers 8 kHz (wave-owl R-832.1,
 * R-834.3; umbrella add-fldigi-native-modes, open question 2).
 *
 * out_rate / in_rate = up / down after reducing by the greatest common divisor. The
 * prototype low-pass g(t) is a Kaiser-windowed sinc with cutoff 0.45 x min(in, out) and
 * taps_per_phase x ceil(down / up) input samples of support, sampled at up x in_rate; output k is
 *   y[k] = sum_j g(j + p / up) x[n - j],  n = floor(k down / up),  p = k down mod up,
 * evaluated in the same order for every chunking, so the output does not depend on how
 * the input is split.
 */
#include "fldigi_dsp.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define KAISER_BETA 8.0
#define CUTOFF 0.45

static long gcd_l(long a, long b)
{
    while (b) {
        long t = a % b;
        a = b;
        b = t;
    }
    return a;
}

/* Modified Bessel function of the first kind, order 0 (power series). */
static double bessel_i0(double x)
{
    double sum = 1.0, term = 1.0, q = x * x / 4.0;
    int k;
    for (k = 1; k < 64; k++) {
        term *= q / ((double)k * (double)k);
        sum += term;
        if (term < sum * 1e-17)
            break;
    }
    return sum;
}

void fdsp_resampler_free(fdsp_resampler *r)
{
    if (!r)
        return;
    free(r->h);
    free(r->hist);
    memset(r, 0, sizeof *r);
}

int fdsp_resampler_init(fdsp_resampler *r, int in_rate, int out_rate, int taps_per_phase)
{
    long g;
    int up, down, taps, n, total;
    double fc, half, norm, i0b;
    memset(r, 0, sizeof *r);
    if (in_rate <= 0 || out_rate <= 0 || taps_per_phase < 4 || taps_per_phase > 256)
        return -1;
    g = gcd_l(in_rate, out_rate);
    up = (int)(out_rate / g);
    down = (int)(in_rate / g);
    if (up > FDSP_RESAMPLE_MAX_UP || down > FDSP_RESAMPLE_MAX_UP)
        return -1;
    /* a decimating resampler needs the support in output samples: scale the taps */
    taps = taps_per_phase * ((down + up - 1) / up);
    total = up * taps;
    r->h = (double *)calloc((size_t)total, sizeof(double));
    r->hist = (double *)calloc((size_t)taps, sizeof(double));
    if (!r->h || !r->hist) {
        fdsp_resampler_free(r);
        return -1;
    }
    r->in_rate = in_rate;
    r->out_rate = out_rate;
    r->up = up;
    r->down = down;
    r->taps = taps;
    /* cutoff in cycles per input sample */
    fc = CUTOFF * (in_rate < out_rate ? in_rate : out_rate) / (double)in_rate;
    half = (double)taps / 2.0;
    i0b = bessel_i0(KAISER_BETA);
    for (n = 0; n < total; n++) {
        double t = (double)n / (double)up - half;   /* input samples from the centre */
        double x = 2.0 * fc * t, s, w, u = t / half;
        s = fabs(x) < 1e-12 ? 1.0 : sin(M_PI * x) / (M_PI * x);
        w = fabs(u) >= 1.0 ? 0.0 : bessel_i0(KAISER_BETA * sqrt(1.0 - u * u)) / i0b;
        r->h[n] = 2.0 * fc * s * w;
    }
    /* unit gain at DC for every phase */
    for (n = 0; n < up; n++) {
        int j;
        norm = 0.0;
        for (j = 0; j < taps; j++)
            norm += r->h[n + j * up];
        if (norm != 0.0)
            for (j = 0; j < taps; j++)
                r->h[n + j * up] /= norm;
    }
    fdsp_resampler_reset(r);
    return 0;
}

void fdsp_resampler_reset(fdsp_resampler *r)
{
    if (r->hist)
        memset(r->hist, 0, (size_t)r->taps * sizeof(double));
    r->pos = 0;
    r->acc = 0;
}

size_t fdsp_resampler_max_out(const fdsp_resampler *r, size_t n_in)
{
    /* each input emits ceil((up - acc) / down) outputs with 0 <= acc < down */
    return (n_in * (size_t)r->up + (size_t)r->down - 1) / (size_t)r->down + 1;
}

size_t fdsp_resampler_run(fdsp_resampler *r, const double *in, size_t n, double *out)
{
    size_t i, k = 0;
    const int taps = r->taps, up = r->up;
    for (i = 0; i < n; i++) {
        /* newest sample at hist[pos], older ones at pos + 1, pos + 2, ... (mod taps) */
        r->pos = (r->pos == 0 ? taps : r->pos) - 1;
        r->hist[r->pos] = in[i];
        while (r->acc < up) {
            const double *h = r->h + r->acc;
            double y = 0.0;
            int j, idx = r->pos;
            for (j = 0; j < taps; j++) {
                y += h[j * up] * r->hist[idx];
                if (++idx == taps)
                    idx = 0;
            }
            out[k++] = y;
            r->acc += r->down;
        }
        r->acc -= up;
    }
    return k;
}
