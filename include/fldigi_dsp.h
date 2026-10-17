/*
 * fldigi-dsp-core -- shared DSP primitives of fldigi 4.2.13, ported to C99.
 *
 * Copyright (C) 2026 Gerolf Ziegenhain
 * Ported from fldigi 4.2.13 (https://github.com/wave-owl/fldigi, tag v4.2.13):
 *   src/include/jalocha/pj_fft.h, pj_fht.h, pj_gray.h, pj_lowpass3.h
 *   (c) 1999-2004 Pawel Jalocha SP9VRC, (c) 2019 David Freese W1HKJ.
 *
 * This program is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this
 * program. If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef FLDIGI_DSP_H
#define FLDIGI_DSP_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Version of this library. Family cores pin it (dsp.lock) and refuse to build
 * against a different version or file checksum (wave-owl R-920.3, INV-920). */
#define FDSP_VERSION_MAJOR 0
#define FDSP_VERSION_MINOR 1
#define FDSP_VERSION_PATCH 0
#define FDSP_VERSION_STRING "0.1.0"

/* "fldigi-dsp-core 0.1.0 contract=dsp-c1 abi=1". Static string. */
const char *fdsp_version(void);

/* ------------------------------------------------------------------------------------
 * Complex number (pj_cmpx.h: Cmpx<double>, fields Re, Im)
 * ------------------------------------------------------------------------------------ */
typedef struct fdsp_cmpx {
    double re;
    double im;
} fdsp_cmpx;

/* ------------------------------------------------------------------------------------
 * Radix-2 FFT of Pawel Jalocha (pj_fft.h:63-317, class r2FFT).
 *
 * Forward transform with kernel e^{-j 2 pi k n / N}, in place, unscrambling included.
 * No normalisation. The two-real-sequences helpers follow pj_fft.h:159-203.
 * ------------------------------------------------------------------------------------ */
typedef struct fdsp_r2fft {
    size_t size;            /* FFT length, a power of two >= 4 */
    size_t *bitrev;         /* bit-reverse index table (BitRevIdx) */
    fdsp_cmpx *twiddle;     /* full circle of twiddles: twiddle[k] = e^{j 2 pi k / size} */
} fdsp_r2fft;

/* Allocate and fill the tables (r2FFT::Preset, pj_fft.h:92-125). Returns 0, or -1 for a
 * size that is not a power of two >= 4 or when allocation fails (then *fft is zeroed). */
int fdsp_r2fft_init(fdsp_r2fft *fft, size_t size);
/* Free the tables (r2FFT::Free). Safe on a zeroed or already freed object. */
void fdsp_r2fft_free(fdsp_r2fft *fft);
/* In-place complex FFT of fft->size points (r2FFT::Process, pj_fft.h:249-250). */
void fdsp_r2fft_process(const fdsp_r2fft *fft, fdsp_cmpx *x);
/* Spectra of two real sequences packed as re/im of buf (SeparTwoReals, pj_fft.h:161-178).
 * out0/out1 hold size/2 points each; out[0].im carries the Nyquist bin. */
void fdsp_r2fft_separ_two_reals(const fdsp_r2fft *fft, const fdsp_cmpx *buf,
                                fdsp_cmpx *out0, fdsp_cmpx *out1);
/* Inverse of separ_two_reals prepared for a forward process() (JoinTwoReals,
 * pj_fft.h:186-203): after fdsp_r2fft_process(buf), buf[].re is sequence 0 and buf[].im
 * sequence 1, both scaled by 2 * size. */
void fdsp_r2fft_join_two_reals(const fdsp_r2fft *fft, const fdsp_cmpx *in0,
                               const fdsp_cmpx *in1, fdsp_cmpx *buf);

/* ------------------------------------------------------------------------------------
 * Fast (Walsh-)Hadamard transform (pj_fht.h:11-48). len is a power of two.
 * ------------------------------------------------------------------------------------ */
void fdsp_fht(double *data, size_t len);           /* FHT  (pj_fht.h:11-27) */
void fdsp_ifht(double *data, size_t len);          /* IFHT (pj_fht.h:31-48) */
void fdsp_ifht_i8(int8_t *data, size_t len);       /* IFHT on int8_t (MFSK_Encoder) */

/* ------------------------------------------------------------------------------------
 * Gray code (pj_gray.h:11-21).
 * ------------------------------------------------------------------------------------ */
uint8_t fdsp_gray_code(uint8_t binary);            /* Binary ^ (Binary >> 1) */
uint8_t fdsp_binary_code(uint8_t gray);            /* inverse, BinaryCode(uint8_t) */

/* ------------------------------------------------------------------------------------
 * Third-order IIR integrator (pj_lowpass3.h:13-41, class LowPass3_Filter).
 * ------------------------------------------------------------------------------------ */
typedef struct fdsp_lowpass3 {
    double out1;
    double out2;
    double output;
} fdsp_lowpass3;

/* Set all three stages to level (LowPass3_Filter::Set / operator=). */
void fdsp_lowpass3_set(fdsp_lowpass3 *f, double level);
/* One step with weight (1 / peaking time) and feedback (fldigi default 0.1). */
void fdsp_lowpass3_process(fdsp_lowpass3 *f, double input, double weight, double feedback);

#ifdef __cplusplus
}
#endif

#endif /* FLDIGI_DSP_H */
