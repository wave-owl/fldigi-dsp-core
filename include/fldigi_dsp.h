/*
 * fldigi-dsp-core -- shared DSP primitives of fldigi 4.2.13, ported to C99.
 *
 * Copyright (C) 2026 Gerolf Ziegenhain
 * Ported from fldigi 4.2.13 (https://github.com/wave-owl/fldigi, tag v4.2.13):
 *   src/include/jalocha/pj_fft.h, pj_fht.h, pj_gray.h, pj_lowpass3.h
 *   (c) 1999-2004 Pawel Jalocha SP9VRC, (c) 2019 David Freese W1HKJ;
 *   src/include/gfft.h (John Green, public domain; C++ (c) 2013 Dave Freese W1HKJ);
 *   src/filters/filters.cxx, fftfilt.cxx, viterbi.cxx, src/mfsk/interleave.cxx,
 *   src/mfsk/mfskvaricode.cxx (c) 2006-2010 Dave Freese W1HKJ and the fldigi authors,
 *   based on gmfsk (c) Tomi Manninen OH2BNS.
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
#define FDSP_VERSION_MINOR 2
#define FDSP_VERSION_PATCH 0
#define FDSP_VERSION_STRING "0.2.0"

/* "fldigi-dsp-core 0.2.0 contract=dsp-c1 abi=1". Static string. */
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

/* ------------------------------------------------------------------------------------
 * Small helpers of src/include/misc.h.
 * ------------------------------------------------------------------------------------ */
/* decayavg, include/misc.h:59-63 */
double fdsp_decayavg(double average, double input, int weight);
/* parity, misc/misc.cxx:34-42, 67-70 (hweight32 & 1) */
int fdsp_parity(unsigned long w);
/* sin(x) and cos(x) as an optimising compiler evaluates a cos/sin pair of one argument:
 * one sincos call (macOS __sincos, whose last bit can differ from separate sin and cos
 * calls), so the result does not depend on the optimisation level of the build. */
void fdsp_sincos(double x, double *s, double *c);

/* ------------------------------------------------------------------------------------
 * Complex FFT of fldigi (src/include/gfft.h, g_fft<double>, John Green's algorithm).
 * Forward kernel e^{-j 2 pi k n / N}, unnormalised; the inverse is scaled by 1/N.
 * ------------------------------------------------------------------------------------ */
typedef struct fdsp_gfft {
    size_t size;            /* a power of two, 16 .. 2^28 */
    int log2n;
    double *utbl;           /* cosine table, size/4 + 1 entries (fftCosInit) */
    short *brlow;           /* bit-reverse table (fftBRInit) */
} fdsp_gfft;

/* g_fft(size) / fftInit, include/gfft.h:120-124, 3242-3265. Returns 0, or -1 (size not a power
 * of two in 16 .. 2^28, out of memory; *fft is then zeroed). */
int fdsp_gfft_init(fdsp_gfft *fft, size_t size);
void fdsp_gfft_free(fdsp_gfft *fft);
/* ComplexFFT / InverseComplexFFT, include/gfft.h:3315-3337: in place on fft->size points. */
void fdsp_gfft_forward(const fdsp_gfft *fft, fdsp_cmpx *buf);
void fdsp_gfft_inverse(const fdsp_gfft *fft, fdsp_cmpx *buf);

/* ------------------------------------------------------------------------------------
 * fftfilt: overlap-add FFT convolution filter (src/filters/fftfilt.cxx:55-237).
 * ------------------------------------------------------------------------------------ */
typedef struct fdsp_fftfilt {
    int flen;               /* FFT length */
    int flen2;              /* flen / 2: block size and impulse response length */
    fdsp_gfft fft;
    fdsp_cmpx *filter;      /* H(w), flen points */
    fdsp_cmpx *timedata;    /* flen points, second half stays zero */
    fdsp_cmpx *freqdata;
    fdsp_cmpx *output;      /* flen2 output samples of the last block */
    fdsp_cmpx *ovlbuf;      /* flen2 overlap samples */
    int inptr;
    int pass;
} fdsp_fftfilt;

/* fftfilt(f1, f2, len), filters/fftfilt.cxx:69-100, 124-196. f1 < f2: band pass; f1 == 0: low
 * pass; f2 == 0: high pass; f1 > f2: band reject (fractions of the sample rate). len is
 * a power of two >= 32. Returns 0, or -1 (bad length, out of memory; zeroed). */
int fdsp_fftfilt_init(fdsp_fftfilt *f, double f1, double f2, int len);
void fdsp_fftfilt_free(fdsp_fftfilt *f);
/* create_filter, filters/fftfilt.cxx:124-196: new shape, state cleared. */
void fdsp_fftfilt_create(fdsp_fftfilt *f, double f1, double f2);
/* run, filters/fftfilt.cxx:202-237: one input sample; returns flen2 when a block of output is
 * ready in *out (valid until the next call), else 0. */
int fdsp_fftfilt_run(fdsp_fftfilt *f, fdsp_cmpx in, const fdsp_cmpx **out);

/* ------------------------------------------------------------------------------------
 * C_FIR_filter (src/filters/filters.cxx:46-187, src/include/filters.h:40-110).
 * ------------------------------------------------------------------------------------ */
#define FDSP_FIR_BUFLEN 4096    /* FIRBufferLen, include/filters.h:41 */
typedef struct fdsp_fir {
    int length;
    int decimateratio;
    double *ifilter;
    double *qfilter;
    double ibuffer[FDSP_FIR_BUFLEN];
    double qbuffer[FDSP_FIR_BUFLEN];
    int pointer;
    int counter;
} fdsp_fir;

/* bp_FIR, filters/filters.cxx:92-119: len Hamming-windowed taps into taps[], f1 < f2 fractions
 * of the sample rate; hilbert != 0 gives the (negated) quadrature branch. */
void fdsp_fir_bp_taps(double *taps, int len, int hilbert, double f1, double f2);
/* init_hilbert (filters/filters.cxx:152-158), init_lowpass (:128-133), init_bandpass (:142-146).
 * Return 0, or -1 (len < 1 or >= FDSP_FIR_BUFLEN, out of memory). */
int fdsp_fir_init_hilbert(fdsp_fir *f, int len, int dec);
int fdsp_fir_init_lowpass(fdsp_fir *f, int len, int dec, double freq);
int fdsp_fir_init_bandpass(fdsp_fir *f, int len, int dec, double f1, double f2);
void fdsp_fir_free(fdsp_fir *f);
/* run, filters/filters.cxx:168-187: returns 1 when *out holds a (decimated) output. */
int fdsp_fir_run(fdsp_fir *f, fdsp_cmpx in, fdsp_cmpx *out);

/* ------------------------------------------------------------------------------------
 * Cmovavg: moving average (filters/filters.cxx:261-307). The first value fills the window.
 * ------------------------------------------------------------------------------------ */
#define FDSP_MOVAVG_MAX 2048    /* MAXMOVAVG, include/filters.h:116 */
typedef struct fdsp_movavg {
    double in[FDSP_MOVAVG_MAX];
    double out;
    int len;
    int pint;
    int empty;
} fdsp_movavg;

/* Cmovavg(len), filters/filters.cxx:261-266; len 1 .. FDSP_MOVAVG_MAX. Returns 0 or -1. */
int fdsp_movavg_init(fdsp_movavg *m, int len);
void fdsp_movavg_reset(fdsp_movavg *m);        /* filters/filters.cxx:304-307 */
double fdsp_movavg_run(fdsp_movavg *m, double a);   /* filters/filters.cxx:273-292 */

/* ------------------------------------------------------------------------------------
 * sfft: sliding DFT over len samples for bins first .. last-1 (filters/filters.cxx:402-465),
 * with fldigi's damping K1 = 0.99999999999L per sample.
 * ------------------------------------------------------------------------------------ */
typedef struct fdsp_sfft {
    int fftlen;
    int first;
    int last;
    int ptr;
    int count;
    double k2;
    fdsp_cmpx *vrot;        /* len rotation factors */
    fdsp_cmpx *bins;        /* len bins (only first .. last-1 used) */
    fdsp_cmpx *delay;       /* len past input samples */
} fdsp_sfft;

/* sfft(len, first, last), filters/filters.cxx:407-424. Returns 0, or -1 (bad range, memory). */
int fdsp_sfft_init(fdsp_sfft *s, int len, int first, int last);
void fdsp_sfft_free(fdsp_sfft *s);
void fdsp_sfft_reset(fdsp_sfft *s);            /* filters/filters.cxx:432-436 */
/* run, filters/filters.cxx:447-465: result[k * stride] = bin first + k for k < last - first. */
void fdsp_sfft_run(fdsp_sfft *s, fdsp_cmpx input, fdsp_cmpx *result, int stride);

/* ------------------------------------------------------------------------------------
 * Convolutional code (src/filters/viterbi.cxx, src/include/viterbi.h): rate 1/2,
 * constraint length k <= 15, soft symbols 0 .. 255 (128 = erasure).
 * ------------------------------------------------------------------------------------ */
#define FDSP_VITERBI_PATHMEM 256    /* PATHMEM, include/viterbi.h:28 */
typedef struct fdsp_viterbi {
    int traceback_len;
    int chunksize;
    int nstates;
    int outsize;
    int k;
    int poly1;
    int poly2;
    int *output;            /* outsize entries */
    int *metrics;           /* PATHMEM x nstates */
    int *history;           /* PATHMEM x nstates */
    int sequence[FDSP_VITERBI_PATHMEM];
    int mettab[2][256];
    unsigned int ptr;
} fdsp_viterbi;

/* viterbi(k, poly1, poly2) + init(), filters/viterbi.cxx:36-82. Returns 0, or -1 (k outside
 * 3 .. 15, out of memory; zeroed). Default traceback k * 12, chunk size 8. */
int fdsp_viterbi_init(fdsp_viterbi *v, int k, int poly1, int poly2);
void fdsp_viterbi_free(fdsp_viterbi *v);
void fdsp_viterbi_reset(fdsp_viterbi *v);              /* filters/viterbi.cxx:85-92 */
int fdsp_viterbi_set_traceback(fdsp_viterbi *v, int t);  /* filters/viterbi.cxx:94-100 */
int fdsp_viterbi_set_chunksize(fdsp_viterbi *v, int c);  /* filters/viterbi.cxx:102-108 */
/* decode, filters/viterbi.cxx:154-214: one symbol pair; returns the decoded bits of a chunk (and
 * the path metric difference in *metric when not NULL) or -1 while none is ready. */
int fdsp_viterbi_decode(fdsp_viterbi *v, const unsigned char sym[2], int *metric);

typedef struct fdsp_encoder {
    int *output;
    unsigned int shreg;
    unsigned int shregmask;
    int k;
    int poly1;
    int poly2;
} fdsp_encoder;

/* encoder(k, poly1, poly2), filters/viterbi.cxx:218-257. Returns 0 or -1. */
int fdsp_encoder_init(fdsp_encoder *e, int k, int poly1, int poly2);
void fdsp_encoder_free(fdsp_encoder *e);
/* encode, filters/viterbi.cxx:234-239: two code bits (poly1 in bit 0, poly2 in bit 1). */
int fdsp_encoder_encode(fdsp_encoder *e, int bit);

/* ------------------------------------------------------------------------------------
 * Block interleaver (src/mfsk/interleave.cxx, src/include/interleave.h):
 * size x size x depth, forward (TX) or reverse (RX, filled with erasures).
 * ------------------------------------------------------------------------------------ */
#define FDSP_INTERLEAVE_FWD 0
#define FDSP_INTERLEAVE_REV 1
#define FDSP_PUNCTURE 128
typedef struct fdsp_interleave {
    int size;
    int depth;
    int len;
    int direction;
    unsigned char *table;
} fdsp_interleave;

/* interleave(size, depth, dir), mfsk/interleave.cxx:34-42; size 1 .. 32. Returns 0 or -1. */
int fdsp_interleave_init(fdsp_interleave *il, int size, int depth, int direction);
void fdsp_interleave_free(fdsp_interleave *il);
void fdsp_interleave_flush(fdsp_interleave *il);       /* mfsk/interleave.cxx:92-100 */
/* symbols, mfsk/interleave.cxx:57-76: size soft symbols in place. */
void fdsp_interleave_symbols(fdsp_interleave *il, unsigned char *psyms);
/* bits, mfsk/interleave.cxx:78-90: the low size bits of *pbits, MSB first, in place. */
void fdsp_interleave_bits(fdsp_interleave *il, unsigned int *pbits);

/* ------------------------------------------------------------------------------------
 * IZ8BLY MFSK varicode (src/mfsk/mfskvaricode.cxx), also THOR's primary set and the
 * MultiPSK DominoEX FEC code.
 * ------------------------------------------------------------------------------------ */
/* varienc, mfsk/mfskvaricode.cxx:332-337: the code of byte c as a '0'/'1' string (c outside
 * 0 .. 255 gives the code of NUL). */
const char *fdsp_mfsk_varienc(int c);
/* varidec, mfsk/mfskvaricode.cxx:339-349: the byte of a code (bits MSB first, without the
 * separating "001"), or -1. */
int fdsp_mfsk_varidec(unsigned int symbol);

#ifdef __cplusplus
}
#endif

#endif /* FLDIGI_DSP_H */
