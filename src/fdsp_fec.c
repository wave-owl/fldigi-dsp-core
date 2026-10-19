/*
 * fldigi-dsp-core: rate 1/2 convolutional encoder and soft-decision Viterbi decoder
 * (src/filters/viterbi.cxx, src/include/viterbi.h), the block interleaver
 * (src/mfsk/interleave.cxx, src/include/interleave.h) and the IZ8BLY MFSK varicode
 * (src/mfsk/mfskvaricode.cxx) of fldigi 4.2.13.
 *
 * Copyright (C) 2026 Gerolf Ziegenhain
 * Original C++: (c) 2006-2010 Dave Freese W1HKJ and the fldigi authors, after gmfsk
 * (c) Tomi Manninen OH2BNS; varicode by Nino Porcino IZ8BLY.
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * One deliberate difference: fldigi renormalises the path metrics only on calls that do
 * not end in a traceback (filters/viterbi.cxx:198-211). With a chunk size of 1 (THOR, MFSK)
 * every call ends in one, so the int metrics grow by up to 256 per symbol pair and
 * overflow after about 8.4 million pairs (signed overflow, undefined in C). This port
 * applies the same renormalisation on every call. It subtracts one constant from all
 * metrics, so decisions, tracebacks and the reported metric difference are unchanged.
 */
#include "fldigi_dsp.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------------------- */
/* viterbi, filters/viterbi.cxx:36-214                                                         */
/* ---------------------------------------------------------------------------------- */
/*
 * Storage (0.9.0). fldigi keeps PATHMEM rows of path metrics and of predecessor states
 * (2 x 256 x 2^(k-1) ints: 64 MB for K=16). Its results only need:
 *
 * - the path metrics of the previous step (the add-compare-select) and of the current one
 *   (the best state of a traceback): two rows;
 * - the decision of every state and step: the predecessor is s >> 1 or (s + nstates) >> 1,
 *   one bit (fldigi's history[p][s] is p0 or p1);
 * - the metric difference of a traceback, metrics[end][seq[end]] - metrics[start][seq[start]]
 *   (viterbi.cxx:138-149). Along the traced path every step adds exactly one branch metric
 *   met[output[s]] (viterbi.cxx:184-193), so the difference is the sum of the branch metrics
 *   of the chunk's steps; the four met[] values of every step are kept for that. The
 *   renormalisation shifts all of fldigi's rows by one constant and leaves the difference
 *   unchanged.
 *
 * Rows not yet written since the last reset hold metric 0 and predecessor 0 in fldigi
 * (memset in reset, viterbi.cxx:85-92); a traceback into them continues at state 0 and adds
 * nothing. `filled` counts the written rows to reproduce that. Decisions, decoded bits and
 * metric differences are identical to fldigi's (goldens K=5, 7, 13, 15, 16).
 */

#define VROW(v, p) ((v)->metrics + (size_t)((p) & 1u) * (size_t)(v)->nstates)

static int vhis_bit(const fdsp_viterbi *v, unsigned int p, int s)
{
    size_t i = (size_t)p * (size_t)v->nstates + (size_t)s;
    return (int)((v->history[i >> 5] >> (i & 31u)) & 1u);
}

/* fldigi's history[p][s]: the predecessor of state s at step p (0 in unwritten rows). */
static int vpred(const fdsp_viterbi *v, unsigned int p, int s)
{
    if (p >= v->filled)
        return 0;
    return vhis_bit(v, p, s) ? (s + v->nstates) >> 1 : s >> 1;
}

/* The branch metric of step p into state s (0 in unwritten rows). */
static int vbranch(const fdsp_viterbi *v, unsigned int p, int s)
{
    int full;
    if (p >= v->filled)
        return 0;
    full = vhis_bit(v, p, s) ? s + v->nstates : s;
    return v->branch[4 * (size_t)p + (size_t)v->output[full]];
}

void fdsp_viterbi_free(fdsp_viterbi *v)
{
    /* viterbi::~viterbi, filters/viterbi.cxx:54-61 */
    if (!v)
        return;
    free(v->output);
    free(v->metrics);
    free(v->history);
    free(v->branch);
    memset(v, 0, sizeof *v);
}

int fdsp_viterbi_init(fdsp_viterbi *v, int k, int poly1, int poly2)
{
    /* viterbi::viterbi + init, filters/viterbi.cxx:36-82 */
    int i;
    size_t bits;
    if (!v)
        return -1;
    memset(v, 0, sizeof *v);
    if (k < 3 || k > 16)
        return -1;
    v->outsize = 1 << k;
    v->nstates = 1 << (k - 1);
    v->k = k;
    v->poly1 = poly1;
    v->poly2 = poly2;
    bits = (size_t)FDSP_VITERBI_PATHMEM * (size_t)v->nstates;
    v->output = (int *)malloc((size_t)v->outsize * sizeof(int));
    v->metrics = (int *)calloc(2 * (size_t)v->nstates, sizeof(int));
    v->history = (uint32_t *)calloc((bits + 31) / 32, sizeof(uint32_t));
    v->branch = (int *)calloc(4 * (size_t)FDSP_VITERBI_PATHMEM, sizeof(int));
    if (!v->output || !v->metrics || !v->history || !v->branch) {
        fdsp_viterbi_free(v);
        return -1;
    }
    v->traceback_len = k * 12;
    v->chunksize = 8;
    for (i = 0; i < v->outsize; i++)
        v->output[i] = fdsp_parity((unsigned long)(poly1 & i)) |
                       (fdsp_parity((unsigned long)(poly2 & i)) << 1);
    for (i = 0; i < 256; i++) {
        v->mettab[0][i] = 128 - i;
        v->mettab[1][i] = i - 128;
    }
    memset(v->sequence, 0, sizeof v->sequence);
    fdsp_viterbi_reset(v);
    return 0;
}

void fdsp_viterbi_reset(fdsp_viterbi *v)
{
    /* viterbi::reset, filters/viterbi.cxx:85-92 */
    size_t bits = (size_t)FDSP_VITERBI_PATHMEM * (size_t)v->nstates;
    memset(v->metrics, 0, 2 * (size_t)v->nstates * sizeof(int));
    memset(v->history, 0, (bits + 31) / 32 * sizeof(uint32_t));
    memset(v->branch, 0, 4 * (size_t)FDSP_VITERBI_PATHMEM * sizeof(int));
    v->filled = 0;
    v->ptr = 0;
}

size_t fdsp_viterbi_footprint(const fdsp_viterbi *v)
{
    size_t bits = (size_t)FDSP_VITERBI_PATHMEM * (size_t)v->nstates;
    return (size_t)v->outsize * sizeof(int) + 2 * (size_t)v->nstates * sizeof(int)
           + (bits + 31) / 32 * sizeof(uint32_t) + 4 * (size_t)FDSP_VITERBI_PATHMEM * sizeof(int);
}

int fdsp_viterbi_set_traceback(fdsp_viterbi *v, int trace)
{
    /* filters/viterbi.cxx:94-100 */
    if (trace < 0 || trace > FDSP_VITERBI_PATHMEM - 1)
        return -1;
    v->traceback_len = trace;
    return 0;
}

int fdsp_viterbi_set_chunksize(fdsp_viterbi *v, int chunk)
{
    /* filters/viterbi.cxx:102-108 */
    if (chunk < 1 || chunk > v->traceback_len)
        return -1;
    v->chunksize = chunk;
    return 0;
}

static int viterbi_traceback(fdsp_viterbi *v, int *metric)
{
    /* viterbi::traceback, filters/viterbi.cxx:110-152 */
    int bestmetric, beststate, i, sum = 0;
    unsigned int p, c = 0;
    const int *row;

    p = (v->ptr - 1) % FDSP_VITERBI_PATHMEM;
    row = VROW(v, p);
    bestmetric = INT_MIN;
    beststate = 0;
    for (i = 0; i < v->nstates; i++) {
        if (row[i] > bestmetric) {
            bestmetric = row[i];
            beststate = i;
        }
    }
    v->sequence[p] = beststate;
    for (i = 0; i < v->traceback_len; i++) {
        unsigned int prev = (p - 1) % FDSP_VITERBI_PATHMEM;
        v->sequence[prev] = vpred(v, p, v->sequence[p]);
        p = prev;
    }
    for (i = 0; i < v->chunksize; i++) {
        c = (c << 1) | (unsigned int)(v->sequence[p] & 1);
        p = (p + 1) % FDSP_VITERBI_PATHMEM;
        /* metrics[p][seq[p]] - metrics[p - 1][seq[p - 1]] along the traced path */
        sum += vbranch(v, p, v->sequence[p]);
    }
    if (metric)
        *metric = sum;
    return (int)c;
}

static void viterbi_renormalise(fdsp_viterbi *v, unsigned int currptr)
{
    /* filters/viterbi.cxx:201-211, on the current row (the previous one is not read again) */
    int *row = VROW(v, currptr), i;
    if (row[0] > INT_MAX / 2) {
        for (i = 0; i < v->nstates; i++)
            row[i] -= INT_MAX / 2;
    }
    if (row[0] < INT_MIN / 2) {
        for (i = 0; i < v->nstates; i++)
            row[i] += INT_MIN / 2;
    }
}

int fdsp_viterbi_decode(fdsp_viterbi *v, const unsigned char sym[2], int *metric)
{
    /* viterbi::decode, filters/viterbi.cxx:154-214 */
    unsigned int currptr, prevptr;
    int met[4], n, *cur, *branch;
    const int *prev;
    uint32_t *his;
    size_t base;

    currptr = v->ptr;
    prevptr = (currptr - 1) % FDSP_VITERBI_PATHMEM;
    met[0] = v->mettab[0][sym[1]] + v->mettab[0][sym[0]];
    met[1] = v->mettab[0][sym[1]] + v->mettab[1][sym[0]];
    met[2] = v->mettab[1][sym[1]] + v->mettab[0][sym[0]];
    met[3] = v->mettab[1][sym[1]] + v->mettab[1][sym[0]];
    branch = v->branch + 4 * (size_t)currptr;
    for (n = 0; n < 4; n++)
        branch[n] = met[n];

    cur = VROW(v, currptr);
    prev = VROW(v, prevptr);
    his = v->history;
    base = (size_t)currptr * (size_t)v->nstates;
    for (n = 0; n < v->nstates; n++) {
        int p0, p1, s0, s1, m0, m1;
        size_t bit = base + (size_t)n;
        s0 = n;
        s1 = n + v->nstates;
        p0 = s0 >> 1;
        p1 = s1 >> 1;
        m0 = prev[p0] + met[v->output[s0]];
        m1 = prev[p1] + met[v->output[s1]];
        if (m0 > m1) {
            cur[n] = m0;
            his[bit >> 5] &= ~((uint32_t)1 << (bit & 31u));
        } else {
            cur[n] = m1;
            his[bit >> 5] |= (uint32_t)1 << (bit & 31u);
        }
    }
    if (v->filled < FDSP_VITERBI_PATHMEM)
        v->filled++;
    v->ptr = (v->ptr + 1) % FDSP_VITERBI_PATHMEM;
    /* see the file comment: renormalised before a traceback as well */
    viterbi_renormalise(v, currptr);
    if ((v->ptr % (unsigned int)v->chunksize) == 0)
        return viterbi_traceback(v, metric);
    return -1;
}

/* ---------------------------------------------------------------------------------- */
/* encoder, filters/viterbi.cxx:218-257                                                        */
/* ---------------------------------------------------------------------------------- */

void fdsp_encoder_free(fdsp_encoder *e)
{
    if (!e)
        return;
    free(e->output);
    memset(e, 0, sizeof *e);
}

int fdsp_encoder_init(fdsp_encoder *e, int k, int poly1, int poly2)
{
    int i, size;
    if (!e)
        return -1;
    memset(e, 0, sizeof *e);
    if (k < 3 || k > 16)
        return -1;
    size = 1 << k;
    e->output = (int *)malloc((size_t)size * sizeof(int));
    if (!e->output)
        return -1;
    e->k = k;
    e->poly1 = poly1;
    e->poly2 = poly2;
    for (i = 0; i < size; i++)
        e->output[i] = fdsp_parity((unsigned long)(poly1 & i)) |
                       (fdsp_parity((unsigned long)(poly2 & i)) << 1);
    e->shreg = 0;
    e->shregmask = (unsigned int)size - 1;
    return 0;
}

int fdsp_encoder_encode(fdsp_encoder *e, int bit)
{
    e->shreg = (e->shreg << 1) | (unsigned int)!!bit;
    return e->output[e->shreg & e->shregmask];
}

/* ---------------------------------------------------------------------------------- */
/* interleave, mfsk/interleave.cxx:34-100, include/interleave.h:31-50                               */
/* ---------------------------------------------------------------------------------- */

#define ITAB(il, i, j, k) (&(il)->table[((il)->size * (il)->size * (i)) + ((il)->size * (j)) + (k)])

int fdsp_interleave_init(fdsp_interleave *il, int size, int depth, int direction)
{
    if (!il)
        return -1;
    memset(il, 0, sizeof *il);
    if (size < 1 || size > 32 || depth < 1)
        return -1;
    il->size = size;
    il->depth = depth;
    il->direction = direction;
    il->len = size * size * depth;
    il->table = (unsigned char *)malloc((size_t)il->len);
    if (!il->table)
        return -1;
    fdsp_interleave_flush(il);
    return 0;
}

void fdsp_interleave_free(fdsp_interleave *il)
{
    if (!il)
        return;
    free(il->table);
    memset(il, 0, sizeof *il);
}

void fdsp_interleave_flush(fdsp_interleave *il)
{
    if (il->direction == FDSP_INTERLEAVE_REV)
        memset(il->table, FDSP_PUNCTURE, (size_t)il->len);
    else
        memset(il->table, 0, (size_t)il->len);
}

void fdsp_interleave_symbols(fdsp_interleave *il, unsigned char *psyms)
{
    int i, j, k, size = il->size;
    for (k = 0; k < il->depth; k++) {
        for (i = 0; i < size; i++)
            for (j = 0; j < size - 1; j++)
                *ITAB(il, k, i, j) = *ITAB(il, k, i, j + 1);
        for (i = 0; i < size; i++)
            *ITAB(il, k, i, size - 1) = psyms[i];
        for (i = 0; i < size; i++) {
            if (il->direction == FDSP_INTERLEAVE_FWD)
                psyms[i] = *ITAB(il, k, i, size - i - 1);
            else
                psyms[i] = *ITAB(il, k, i, i);
        }
    }
}

void fdsp_interleave_bits(fdsp_interleave *il, unsigned int *pbits)
{
    unsigned char syms[32];
    int i, size = il->size;
    for (i = 0; i < size; i++)
        syms[i] = (unsigned char)((*pbits >> (size - i - 1)) & 1);
    fdsp_interleave_symbols(il, syms);
    for (*pbits = 0, i = 0; i < size; i++)
        *pbits = (*pbits << 1) | syms[i];
}

/* ---------------------------------------------------------------------------------- */
/* IZ8BLY MFSK varicode, mfskvaricode.cxx                                              */
/* ---------------------------------------------------------------------------------- */

/* varicode[], mfsk/mfskvaricode.cxx:35-292 (IZ8BLY MFSK varicode, 256 codes) */
static const char *const mfsk_varicode[256] = {
    "11101011100", "11101100000", "11101101000", "11101101100",
    "11101110000", "11101110100", "11101111000", "11101111100",
    "10101000", "11110000000", "11110100000", "11110101000",
    "11110101100", "10101100", "11110110000", "11110110100",
    "11110111000", "11110111100", "11111000000", "11111010000",
    "11111010100", "11111011000", "11111011100", "11111100000",
    "11111101000", "11111101100", "11111110000", "11111110100",
    "11111111000", "11111111100", "100000000000", "101000000000",
    "100", "111000000", "111111100", "1011011000",
    "1010101000", "1010100000", "1000000000", "110111100",
    "111110100", "111110000", "1010110100", "111100000",
    "10100000", "111011000", "111010100", "111101000",
    "11100000", "11110000", "101000000", "101010100",
    "101110100", "101100000", "101101100", "110100000",
    "110000000", "110101100", "111101100", "111111000",
    "1011000000", "111011100", "1010111100", "111010000",
    "1010000000", "10111100", "100000000", "11010100",
    "11011100", "10111000", "11111000", "101010000",
    "101011000", "11000000", "110110100", "101111100",
    "11110100", "11101000", "11111100", "11010000",
    "11101100", "110110000", "11011000", "10110100",
    "10110000", "101011100", "110101000", "101101000",
    "101110000", "101111000", "110111000", "1011101000",
    "1011010000", "1011101100", "1011010100", "1010110000",
    "1010101100", "10100", "1100000", "111000",
    "110100", "1000", "1010000", "1011000",
    "110000", "11000", "10000000", "1110000",
    "101100", "1000000", "11100", "10000",
    "1010100", "1111000", "100000", "101000",
    "1100", "111100", "1101100", "1101000",
    "1110100", "1011100", "1111100", "1011011100",
    "1010111000", "1011100000", "1011110000", "101010000000",
    "101010100000", "101010101000", "101010101100", "101010110000",
    "101010110100", "101010111000", "101010111100", "101011000000",
    "101011010000", "101011010100", "101011011000", "101011011100",
    "101011100000", "101011101000", "101011101100", "101011110000",
    "101011110100", "101011111000", "101011111100", "101100000000",
    "101101000000", "101101010000", "101101010100", "101101011000",
    "101101011100", "101101100000", "101101101000", "101101101100",
    "101101110000", "101101110100", "101101111000", "101101111100",
    "1011110100", "1011111000", "1011111100", "1100000000",
    "1101000000", "1101010000", "1101010100", "1101011000",
    "1101011100", "1101100000", "1101101000", "1101101100",
    "1101110000", "1101110100", "1101111000", "1101111100",
    "1110000000", "1110100000", "1110101000", "1110101100",
    "1110110000", "1110110100", "1110111000", "1110111100",
    "1111000000", "1111010000", "1111010100", "1111011000",
    "1111011100", "1111100000", "1111101000", "1111101100",
    "1111110000", "1111110100", "1111111000", "1111111100",
    "10000000000", "10100000000", "10101000000", "10101010000",
    "10101010100", "10101011000", "10101011100", "10101100000",
    "10101101000", "10101101100", "10101110000", "10101110100",
    "10101111000", "10101111100", "10110000000", "10110100000",
    "10110101000", "10110101100", "10110110000", "10110110100",
    "10110111000", "10110111100", "10111000000", "10111010000",
    "10111010100", "10111011000", "10111011100", "10111100000",
    "10111101000", "10111101100", "10111110000", "10111110100",
    "10111111000", "10111111100", "11000000000", "11010000000",
    "11010100000", "11010101000", "11010101100", "11010110000",
    "11010110100", "11010111000", "11010111100", "11011000000",
    "11011010000", "11011010100", "11011011000", "11011011100",
    "11011100000", "11011101000", "11011101100", "11011110000",
    "11011110100", "11011111000", "11011111100", "11100000000",
    "11101000000", "11101010000", "11101010100", "11101011000",
};

/* varidecode[], mfsk/mfskvaricode.cxx:297-330 */
static const unsigned int mfsk_varidecode[256] = {
    0x75C, 0x760, 0x768, 0x76C, 0x770, 0x774, 0x778, 0x77C,
    0x0A8, 0x780, 0x7A0, 0x7A8, 0x7AC, 0x0AC, 0x7B0, 0x7B4,
    0x7B8, 0x7BC, 0x7C0, 0x7D0, 0x7D4, 0x7D8, 0x7DC, 0x7E0,
    0x7E8, 0x7EC, 0x7F0, 0x7F4, 0x7F8, 0x7FC, 0x800, 0xA00,
    0x004, 0x1C0, 0x1FC, 0x2D8, 0x2A8, 0x2A0, 0x200, 0x1BC,
    0x1F4, 0x1F0, 0x2B4, 0x1E0, 0x0A0, 0x1D8, 0x1D4, 0x1E8,
    0x0E0, 0x0F0, 0x140, 0x154, 0x174, 0x160, 0x16C, 0x1A0,
    0x180, 0x1AC, 0x1EC, 0x1F8, 0x2C0, 0x1DC, 0x2BC, 0x1D0,
    0x280, 0x0BC, 0x100, 0x0D4, 0x0DC, 0x0B8, 0x0F8, 0x150,
    0x158, 0x0C0, 0x1B4, 0x17C, 0x0F4, 0x0E8, 0x0FC, 0x0D0,
    0x0EC, 0x1B0, 0x0D8, 0x0B4, 0x0B0, 0x15C, 0x1A8, 0x168,
    0x170, 0x178, 0x1B8, 0x2E8, 0x2D0, 0x2EC, 0x2D4, 0x2B0,
    0x2AC, 0x014, 0x060, 0x038, 0x034, 0x008, 0x050, 0x058,
    0x030, 0x018, 0x080, 0x070, 0x02C, 0x040, 0x01C, 0x010,
    0x054, 0x078, 0x020, 0x028, 0x00C, 0x03C, 0x06C, 0x068,
    0x074, 0x05C, 0x07C, 0x2DC, 0x2B8, 0x2E0, 0x2F0, 0xA80,
    0xAA0, 0xAA8, 0xAAC, 0xAB0, 0xAB4, 0xAB8, 0xABC, 0xAC0,
    0xAD0, 0xAD4, 0xAD8, 0xADC, 0xAE0, 0xAE8, 0xAEC, 0xAF0,
    0xAF4, 0xAF8, 0xAFC, 0xB00, 0xB40, 0xB50, 0xB54, 0xB58,
    0xB5C, 0xB60, 0xB68, 0xB6C, 0xB70, 0xB74, 0xB78, 0xB7C,
    0x2F4, 0x2F8, 0x2FC, 0x300, 0x340, 0x350, 0x354, 0x358,
    0x35C, 0x360, 0x368, 0x36C, 0x370, 0x374, 0x378, 0x37C,
    0x380, 0x3A0, 0x3A8, 0x3AC, 0x3B0, 0x3B4, 0x3B8, 0x3BC,
    0x3C0, 0x3D0, 0x3D4, 0x3D8, 0x3DC, 0x3E0, 0x3E8, 0x3EC,
    0x3F0, 0x3F4, 0x3F8, 0x3FC, 0x400, 0x500, 0x540, 0x550,
    0x554, 0x558, 0x55C, 0x560, 0x568, 0x56C, 0x570, 0x574,
    0x578, 0x57C, 0x580, 0x5A0, 0x5A8, 0x5AC, 0x5B0, 0x5B4,
    0x5B8, 0x5BC, 0x5C0, 0x5D0, 0x5D4, 0x5D8, 0x5DC, 0x5E0,
    0x5E8, 0x5EC, 0x5F0, 0x5F4, 0x5F8, 0x5FC, 0x600, 0x680,
    0x6A0, 0x6A8, 0x6AC, 0x6B0, 0x6B4, 0x6B8, 0x6BC, 0x6C0,
    0x6D0, 0x6D4, 0x6D8, 0x6DC, 0x6E0, 0x6E8, 0x6EC, 0x6F0,
    0x6F4, 0x6F8, 0x6FC, 0x700, 0x740, 0x750, 0x754, 0x758,
};

const char *fdsp_mfsk_varienc(int c)
{
    /* varienc, mfsk/mfskvaricode.cxx:332-337 */
    if (c >= 0 && c < 256)
        return mfsk_varicode[c];
    return mfsk_varicode[0];
}

int fdsp_mfsk_varidec(unsigned int symbol)
{
    /* varidec, mfsk/mfskvaricode.cxx:339-349 */
    int i;
    for (i = 0; i < 256; i++)
        if (symbol == mfsk_varidecode[i])
            return i;
    return -1;
}
