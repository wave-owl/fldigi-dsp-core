// Golden data generator for the fldigi-dsp-core 0.2 primitives.
//
// Compiles fldigi 4.2.13's own sources (g_fft, fftfilt, C_FIR_filter, Cmovavg, sfft,
// viterbi/encoder, interleave, MFSK varicode, parity) and writes their results for
// fixed inputs to tests/golden/. The C tests compare the port against these files.
//
//   make golden FLDIGI=/path/to/fldigi
//
// $FLDIGI is a checkout of https://github.com/wave-owl/fldigi at tag v4.2.13
// (commit e7c3ab3709ad62a7d91829364eb48a05fa325c0d).
//
// SPDX-License-Identifier: GPL-3.0-or-later
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <string>
#include <vector>

#include "complex.h"
#include "gfft.h"
#include "fftfilt.h"
#include "filters.h"
#include "viterbi.h"
#include "interleave.h"
#include "mfskvaricode.h"
#include "misc.h"

// fldigi's .cxx files are compiled into this program (see the Makefile target).

static uint32_t lcg_state;
static double lcg_next(void)
{
    lcg_state = lcg_state * 1664525u + 1013904223u;
    return ((double)(int32_t)lcg_state) / 2147483648.0;
}

static std::string gdir;
static void write_file(const char *name, const void *data, size_t n)
{
    std::string path = gdir + "/" + name;
    FILE *f = fopen(path.c_str(), "wb");
    if (!f || fwrite(data, 1, n, f) != n) {
        fprintf(stderr, "cannot write %s\n", path.c_str());
        exit(1);
    }
    fclose(f);
}

static void put_c(std::vector<double> &v, const cmplx &z)
{
    v.push_back(z.real());
    v.push_back(z.imag());
}

// Access to the protected/private state the C tests inject (tables computed with libm).
struct fftfilt_peek : public fftfilt {
    fftfilt_peek(double f1, double f2, int len) : fftfilt(f1, f2, len) {}
    fftfilt_peek(double f, int len) : fftfilt(f, len) {}
    const cmplx *shape() const { return filter; }
};

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "usage: %s <golden dir>\n", argv[0]);
        return 2;
    }
    gdir = argv[1];
    char name[64];

    // g_fft: forward and inverse of LCG data, sizes 16 .. 8192.
    for (int n = 16; n <= 8192; n *= 2) {
        g_fft<double> fft(n);
        std::vector<cmplx> buf(n);
        std::vector<double> out;
        lcg_state = (uint32_t)(5000 + n);
        for (int i = 0; i < n; i++)
            buf[i] = cmplx(lcg_next(), lcg_next());
        fft.ComplexFFT(buf.data());
        for (int i = 0; i < n; i++)
            put_c(out, buf[i]);
        fft.InverseComplexFFT(buf.data());
        for (int i = 0; i < n; i++)
            put_c(out, buf[i]);
        // cosine table: recomputed exactly as fftCosInit does
        std::vector<double> utbl(n / 4 + 1);
        utbl[0] = 1.0;
        for (int i = 1; i < n / 4; i++)
            utbl[i] = cos((2.0 * FFT_PI * (double)i) / (double)n);
        utbl[n / 4] = 0.0;
        snprintf(name, sizeof name, "gfft%d.f64", n);
        write_file(name, out.data(), out.size() * sizeof(double));
        snprintf(name, sizeof name, "gfft%d_utbl.f64", n);
        write_file(name, utbl.data(), utbl.size() * sizeof(double));
    }

    // fftfilt: the THOR16 / DominoEX11 band-pass shapes and a low pass, 1024 points,
    // 6000 LCG input samples; outputs of every completed block.
    {
        const double cases[][2] = {
            {(2000.0 - 0.5 * 2.0 * 281.25) / 8000.0, (2000.0 + 0.5 * 2.0 * 281.25) / 8000.0},
            {(1500.0 - 0.5 * 2.0 * 193.798828125) / 11025.0,
             (1500.0 + 0.5 * 2.0 * 193.798828125) / 11025.0},
            {0.0, 0.1},
        };
        for (int c = 0; c < 3; c++) {
            fftfilt_peek f(cases[c][0], cases[c][1], 1024);
            std::vector<double> shape, out;
            for (int i = 0; i < 1024; i++)
                put_c(shape, f.shape()[i]);
            lcg_state = (uint32_t)(77 + c);
            for (int i = 0; i < 6000; i++) {
                cmplx *zp;
                double re = lcg_next(), im = lcg_next();
                int n = f.run(cmplx(re, im), &zp);
                for (int k = 0; k < n; k++)
                    put_c(out, zp[k]);
            }
            snprintf(name, sizeof name, "fftfilt%d_shape.f64", c);
            write_file(name, shape.data(), shape.size() * sizeof(double));
            snprintf(name, sizeof name, "fftfilt%d_out.f64", c);
            write_file(name, out.data(), out.size() * sizeof(double));
        }
    }

    // fftfilt::rtty_filter (0.7.0): NAVTEX's mark/space low passes (navtex.cxx:998-1005,
    // fftfilt(100/11025, 512) then rtty_filter(100/11025)) and an RTTY 45.45 Bd shape at
    // 8000 Hz, 1024 points; 6000 LCG input samples, outputs of every completed block.
    {
        const double cases[][2] = {{100.0 / 11025.0, 512}, {45.45 / 8000.0, 1024}};
        for (int c = 0; c < 2; c++) {
            int len = (int)cases[c][1];
            fftfilt_peek f(cases[c][0], len);
            f.rtty_filter(cases[c][0]);
            std::vector<double> shape, out;
            for (int i = 0; i < len; i++)
                put_c(shape, f.shape()[i]);
            lcg_state = (uint32_t)(177 + c);
            for (int i = 0; i < 6000; i++) {
                cmplx *zp;
                double re = lcg_next(), im = lcg_next();
                int n = f.run(cmplx(re, im), &zp);
                for (int k = 0; k < n; k++)
                    put_c(out, zp[k]);
            }
            snprintf(name, sizeof name, "fftfilt_rtty%d_shape.f64", c);
            write_file(name, shape.data(), shape.size() * sizeof(double));
            snprintf(name, sizeof name, "fftfilt_rtty%d_out.f64", c);
            write_file(name, out.data(), out.size() * sizeof(double));
        }
    }

    // C_FIR_filter: Hilbert 37 (THOR/DominoEX), low pass 127, band pass 65 with
    // decimation 2; 9000 samples (two wraps of the 4096-sample buffer).
    {
        for (int c = 0; c < 3; c++) {
            C_FIR_filter fir;
            int len = c == 0 ? 37 : c == 1 ? 127 : 65;
            int dec = c == 2 ? 2 : 1;
            double *ti, *tq;
            if (c == 0) {
                fir.init_hilbert(len, dec);
                ti = fir.bp_FIR(len, 0, 0.05, 0.45);
                tq = fir.bp_FIR(len, 1, 0.05, 0.45);
            } else if (c == 1) {
                fir.init_lowpass(len, dec, 0.0625);
                ti = fir.bp_FIR(len, 0, 0.0, 0.0625);
                tq = fir.bp_FIR(len, 0, 0.0, 0.0625);
            } else {
                fir.init_bandpass(len, dec, 0.1, 0.2);
                ti = fir.bp_FIR(len, 0, 0.1, 0.2);
                tq = fir.bp_FIR(len, 0, 0.1, 0.2);
            }
            std::vector<double> taps(ti, ti + len), out;
            taps.insert(taps.end(), tq, tq + len);
            delete [] ti;
            delete [] tq;
            lcg_state = (uint32_t)(300 + c);
            for (int i = 0; i < 9000; i++) {
                cmplx o;
                double re = lcg_next(), im = lcg_next();
                if (fir.run(cmplx(re, im), o))
                    put_c(out, o);
            }
            snprintf(name, sizeof name, "fir%d_taps.f64", c);
            write_file(name, taps.data(), taps.size() * sizeof(double));
            snprintf(name, sizeof name, "fir%d_out.f64", c);
            write_file(name, out.data(), out.size() * sizeof(double));
        }
    }

    // C_FIR_filter with own taps, I branch only (Irun, 0.4.0): Throb's sync filter, the
    // semi pulse of 256 taps (THROB1) with decimation 1, and a full pulse of 64 taps
    // (THROB4) with decimation 2; 9000 non-negative LCG samples (magnitudes).
    {
        for (int c = 0; c < 2; c++) {
            int len = c == 0 ? 256 : 64, dec = c == 0 ? 1 : 2;
            std::vector<double> taps(len), out;
            for (int i = 0; i < len; i++)
                taps[i] = 0.5 * (1 - cos(2 * M_PI * i / len)) + 0.001 * (c + 1) * i;
            C_FIR_filter fir;
            fir.init(len, dec, taps.data(), NULL);
            lcg_state = (uint32_t)(320 + c);
            for (int i = 0; i < 9000; i++) {
                double o = 0.0, x = fabs(lcg_next());
                if (fir.Irun(x, o))
                    out.push_back(o);
            }
            snprintf(name, sizeof name, "firirun%d_taps.f64", c);
            write_file(name, taps.data(), taps.size() * sizeof(double));
            snprintf(name, sizeof name, "firirun%d_out.f64", c);
            write_file(name, out.data(), out.size() * sizeof(double));
        }
    }

    // cmac (complex.h:33-41, 0.4.0): circular correlation of two LCG vectors of length 256
    // and 11 at every start pointer 0 .. len + 2.
    {
        std::vector<double> out;
        const int lens[2] = {256, 11};
        for (int c = 0; c < 2; c++) {
            int len = lens[c];
            std::vector<cmplx> a(len), b(len);
            lcg_state = (uint32_t)(340 + c);
            for (int i = 0; i < len; i++) {
                double re = lcg_next(), im = lcg_next();
                a[i] = cmplx(re, im);
            }
            for (int i = 0; i < len; i++) {
                double re = lcg_next(), im = lcg_next();
                b[i] = cmplx(re, im);
            }
            for (int p = 0; p < len + 3; p++)
                put_c(out, cmac(a.data(), b.data(), p, len));
        }
        write_file("cmac.f64", out.data(), out.size() * sizeof(double));
    }

    // Cmovavg: lengths 8 and 16, 500 values.
    {
        std::vector<double> out;
        for (int len = 8; len <= 16; len += 8) {
            Cmovavg m(len);
            lcg_state = (uint32_t)(400 + len);
            for (int i = 0; i < 500; i++)
                out.push_back(m.run(1000.0 * lcg_next()));
            m.reset();
            out.push_back(m.run(3.25));
            out.push_back(m.run(-1.5));
        }
        write_file("movavg.f64", out.data(), out.size() * sizeof(double));
    }

    // Cmovavg grown with setLength (feld.cxx:236): lengths 2049 and 11428, 12000 values.
    {
        std::vector<double> out;
        const int lens[2] = {2049, 11428};
        for (int k = 0; k < 2; k++) {
            Cmovavg m(200);
            m.setLength(lens[k]);
            lcg_state = (uint32_t)(900 + lens[k]);
            for (int i = 0; i < 12000; i++)
                out.push_back(m.run(1000.0 * lcg_next()));
            m.reset();
            out.push_back(m.run(3.25));
            out.push_back(m.run(-1.5));
        }
        write_file("movavgh.f64", out.data(), out.size() * sizeof(double));
    }

    // sfft: (len, first, last) of THOR8 (1024, 184, 248), THOR16 slow paths (512, 91, 121)
    // and a small one; every step of 3000 LCG samples, stride 3 as THOR's paths use it.
    {
        const int cases[][3] = {{1024, 184, 248}, {512, 91, 121}, {64, 3, 20}};
        for (int c = 0; c < 3; c++) {
            int len = cases[c][0], first = cases[c][1], last = cases[c][2];
            sfft s(len, first, last);
            std::vector<cmplx> res(3 * (last - first));
            std::vector<double> out;
            lcg_state = (uint32_t)(500 + c);
            for (int i = 0; i < 3000; i++) {
                double re = lcg_next(), im = lcg_next();
                s.run(cmplx(re, im), res.data() + 1, 3);
                if (i % 7 == 0 || i > 2990)
                    for (int k = 0; k < last - first; k++)
                        put_c(out, res[1 + 3 * k]);
            }
            snprintf(name, sizeof name, "sfft%d.f64", c);
            write_file(name, out.data(), out.size() * sizeof(double));
            // rotation table and k2 as sfft::sfft computes them
            std::vector<double> vr;
            double phi = 0.0, tau = 2.0 * M_PI / len, k2 = 1.0;
            for (int i = 0; i < len; i++) {
                vr.push_back((double)(K1 * cos(phi)));
                vr.push_back((double)(K1 * sin(phi)));
                phi += tau;
                k2 *= K1;
            }
            vr.push_back(k2);
            snprintf(name, sizeof name, "sfft%d_vrot.f64", c);
            write_file(name, vr.data(), vr.size() * sizeof(double));
        }
    }

    // Convolutional code: K=7 (THOR, MFSK, DominoEX FEC, PSK-R), K=15 (THOR 25/50/100,
    // OFDM-500F), and since 0.8.0 K=5 (QPSK), K=13 (8PSK FL/500F) and K=16 (8PSK125F/250F).
    // Encoder output for 4000 LCG bits; decoder over noisy soft symbols with erasures,
    // traceback ks[c][3], chunk size 1 (outputs and metrics), then chunk size 8.
    {
        const int ks[][4] = {{7, 0x6d, 0x4f, 45}, {15, 044735, 063057, 15 * 12},
                             {5, 0x17, 0x19, 5 * 12}, {13, 016461, 012767, 13 * 16},
                             {16, 0152711, 0126723, 16 * 12}};
        for (int c = 0; c < 5; c++) {
            encoder enc(ks[c][0], ks[c][1], ks[c][2]);
            std::vector<int32_t> codes, dec1, met1, dec8;
            std::vector<unsigned char> sym;
            lcg_state = (uint32_t)(600 + c);
            for (int i = 0; i < 4000; i++) {
                int bit = lcg_next() > 0.0;
                int d = enc.encode(bit);
                codes.push_back(d);
                for (int b = 0; b < 2; b++) {
                    double soft = ((d >> b) & 1) ? 200.0 : 55.0;
                    soft += 90.0 * lcg_next();
                    if (soft < 0) soft = 0;
                    if (soft > 255) soft = 255;
                    sym.push_back(i % 37 == 5 ? 128 : (unsigned char)soft);
                }
            }
            viterbi v1(ks[c][0], ks[c][1], ks[c][2]);
            v1.settraceback(ks[c][3]);
            v1.setchunksize(1);
            for (int i = 0; i < 4000; i++) {
                int met = 0;
                int r = v1.decode(&sym[2 * i], &met);
                dec1.push_back(r);
                met1.push_back(met);
            }
            viterbi v8(ks[c][0], ks[c][1], ks[c][2]);
            for (int i = 0; i < 4000; i++) {
                int met = 0;
                dec8.push_back(v8.decode(&sym[2 * i], &met));
                dec8.push_back(met);
            }
            snprintf(name, sizeof name, "conv%d_enc.i32", ks[c][0]);
            write_file(name, codes.data(), codes.size() * sizeof(int32_t));
            snprintf(name, sizeof name, "conv%d_sym.u8", ks[c][0]);
            write_file(name, sym.data(), sym.size());
            snprintf(name, sizeof name, "conv%d_dec1.i32", ks[c][0]);
            write_file(name, dec1.data(), dec1.size() * sizeof(int32_t));
            snprintf(name, sizeof name, "conv%d_met1.i32", ks[c][0]);
            write_file(name, met1.data(), met1.size() * sizeof(int32_t));
            snprintf(name, sizeof name, "conv%d_dec8.i32", ks[c][0]);
            write_file(name, dec8.data(), dec8.size() * sizeof(int32_t));
        }
    }

    // Interleavers: 4x4x10 (THOR), 4x4x50 (THOR 50/100), 4x4x4 (DominoEX FEC), 5x5x6
    // (MFSK32-like), forward and reverse on symbols, forward on bits.
    {
        const int cases[][2] = {{4, 10}, {4, 50}, {4, 4}, {5, 6}};
        std::vector<unsigned char> out;
        for (int c = 0; c < 4; c++) {
            int size = cases[c][0], depth = cases[c][1];
            interleave fwd(size, depth, INTERLEAVE_FWD), rev(size, depth, INTERLEAVE_REV);
            interleave bfwd(size, depth, INTERLEAVE_FWD);
            lcg_state = (uint32_t)(700 + c);
            for (int i = 0; i < 400; i++) {
                unsigned char s[8], r[8];
                for (int k = 0; k < size; k++)
                    s[k] = r[k] = (unsigned char)(256.0 * (0.5 + 0.5 * lcg_next()));
                fwd.symbols(s);
                rev.symbols(r);
                unsigned int bits = (unsigned int)(1024.0 * (0.5 + 0.5 * lcg_next())) &
                                    ((1u << size) - 1);
                bfwd.bits(&bits);
                out.insert(out.end(), s, s + size);
                out.insert(out.end(), r, r + size);
                out.push_back((unsigned char)bits);
            }
        }
        write_file("interleave.u8", out.data(), out.size());
    }

    // MFSK varicode: every code (NUL-separated) and the decode of every 12-bit value.
    {
        std::string enc;
        for (int c = -1; c <= 256; c++) {
            enc += varienc(c);
            enc += '\0';
        }
        write_file("mfsk_varienc.txt", enc.data(), enc.size());
        std::vector<int32_t> dec;
        for (unsigned int s = 0; s < 4096; s++)
            dec.push_back(varidec(s));
        write_file("mfsk_varidec.i32", dec.data(), dec.size() * sizeof(int32_t));
        std::vector<int32_t> par;
        for (unsigned long w = 0; w < 70000; w += 7)
            par.push_back(parity(w));
        write_file("parity.i32", par.data(), par.size() * sizeof(int32_t));
    }
    return 0;
}
