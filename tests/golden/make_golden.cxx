// Golden data generator for fldigi-dsp-core.
//
// Compiles fldigi 4.2.13's own Jalocha headers (src/include/jalocha) and writes their
// results for fixed inputs to tests/golden/*.f64|*.u8|*.i8. The C tests compare the
// port against these files bit for bit.
//
//   c++ -std=c++11 -O2 -ffp-contract=off -I$FLDIGI/src/include/jalocha \
//       tests/golden/make_golden.cxx -o /tmp/make_golden
//   /tmp/make_golden tests/golden
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

#include "pj_struc.h"
#include "pj_cmpx.h"
#include "pj_fft.h"
#include "pj_fht.h"
#include "pj_gray.h"
#include "pj_lowpass3.h"

// Same generator as tests/test_dsp.c: 32-bit LCG (Numerical Recipes), mapped to [-1, 1).
static uint32_t lcg_state;
static double lcg_next(void)
{
    lcg_state = lcg_state * 1664525u + 1013904223u;
    return ((double)(int32_t)lcg_state) / 2147483648.0;
}

static void write_file(const std::string &dir, const char *name, const void *data, size_t n)
{
    std::string path = dir + "/" + name;
    FILE *f = fopen(path.c_str(), "wb");
    if (!f || fwrite(data, 1, n, f) != n) {
        fprintf(stderr, "cannot write %s\n", path.c_str());
        exit(1);
    }
    fclose(f);
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "usage: %s <golden dir>\n", argv[0]);
        return 2;
    }
    std::string dir = argv[1];
    char name[64];

    // r2FFT: tables and forward transforms of LCG data, sizes 4 .. 8192.
    const size_t sizes[] = {4, 8, 16, 64, 256, 512, 1024, 2048, 4096, 8192};
    for (size_t s = 0; s < sizeof sizes / sizeof sizes[0]; s++) {
        size_t n = sizes[s];
        r2FFT< Cmpx<double> > fft;
        if (fft.Preset((int)n) < 0)
            return 1;
        std::vector<double> tw(2 * n), x(2 * n);
        for (size_t i = 0; i < n; i++) {
            tw[2 * i] = fft.Twiddle[i].Re;
            tw[2 * i + 1] = fft.Twiddle[i].Im;
        }
        snprintf(name, sizeof name, "fft%zu_twiddle.f64", n);
        write_file(dir, name, tw.data(), tw.size() * sizeof(double));
        std::vector<uint32_t> br(n);
        for (size_t i = 0; i < n; i++)
            br[i] = (uint32_t)fft.BitRevIdx[i];
        snprintf(name, sizeof name, "fft%zu_bitrev.u32", n);
        write_file(dir, name, br.data(), br.size() * sizeof(uint32_t));

        std::vector< Cmpx<double> > buf(n);
        lcg_state = (uint32_t)n;
        for (size_t i = 0; i < n; i++) {
            buf[i].Re = lcg_next();
            buf[i].Im = lcg_next();
        }
        fft.Process(buf.data());
        for (size_t i = 0; i < n; i++) {
            x[2 * i] = buf[i].Re;
            x[2 * i + 1] = buf[i].Im;
        }
        snprintf(name, sizeof name, "fft%zu_out.f64", n);
        write_file(dir, name, x.data(), x.size() * sizeof(double));

        // SeparTwoReals of the transform, then JoinTwoReals + Process.
        std::vector< Cmpx<double> > s0(n / 2), s1(n / 2);
        fft.SeparTwoReals(buf.data(), s0.data(), s1.data());
        std::vector<double> sep(2 * n);
        for (size_t i = 0; i < n / 2; i++) {
            sep[2 * i] = s0[i].Re;
            sep[2 * i + 1] = s0[i].Im;
            sep[n + 2 * i] = s1[i].Re;
            sep[n + 2 * i + 1] = s1[i].Im;
        }
        snprintf(name, sizeof name, "fft%zu_separ.f64", n);
        write_file(dir, name, sep.data(), sep.size() * sizeof(double));
        fft.JoinTwoReals(s0.data(), s1.data(), buf.data());
        fft.Process(buf.data());
        for (size_t i = 0; i < n; i++) {
            x[2 * i] = buf[i].Re;
            x[2 * i + 1] = buf[i].Im;
        }
        snprintf(name, sizeof name, "fft%zu_join.f64", n);
        write_file(dir, name, x.data(), x.size() * sizeof(double));
    }

    // FHT / IFHT on LCG doubles, lengths 2 .. 128.
    for (size_t n = 2; n <= 128; n *= 2) {
        std::vector<double> a(n), b(n), out(2 * n);
        lcg_state = (uint32_t)(1000 + n);
        for (size_t i = 0; i < n; i++)
            a[i] = b[i] = lcg_next();
        FHT(a.data(), n);
        IFHT(b.data(), n);
        for (size_t i = 0; i < n; i++) {
            out[i] = a[i];
            out[n + i] = b[i];
        }
        snprintf(name, sizeof name, "fht%zu.f64", n);
        write_file(dir, name, out.data(), out.size() * sizeof(double));
    }

    // IFHT on int8_t: every signed unit vector of length 32 and 64 (MFSK_Encoder input).
    for (size_t n = 32; n <= 64; n *= 2) {
        std::vector<int8_t> all;
        for (size_t c = 0; c < 2 * n; c++) {
            std::vector<int8_t> v(n, 0);
            if (c < n)
                v[c] = 1;
            else
                v[c - n] = -1;
            IFHT(v.data(), n);
            all.insert(all.end(), v.begin(), v.end());
        }
        snprintf(name, sizeof name, "ifht_i8_%zu.i8", n);
        write_file(dir, name, all.data(), all.size());
    }

    // Gray code and its inverse for every byte.
    {
        uint8_t g[512];
        for (int i = 0; i < 256; i++) {
            g[i] = GrayCode((uint8_t)i);
            g[256 + i] = BinaryCode((uint8_t)i);
        }
        write_file(dir, "gray.u8", g, sizeof g);
    }

    // LowPass3_Filter: 2000 LCG steps for three weights, feedback 0.1 (default) and 0.5.
    {
        const double weights[] = {1.0 / 4, 1.0 / 32, 1.0 / 100};
        const double feedbacks[] = {0.1, 0.5};
        std::vector<double> out;
        for (size_t w = 0; w < 3; w++) {
            for (size_t fb = 0; fb < 2; fb++) {
                LowPass3_Filter<double> f;
                f.Set(0.0);
                lcg_state = (uint32_t)(77 + 10 * w + fb);
                for (int i = 0; i < 2000; i++) {
                    f.Process(lcg_next(), weights[w], feedbacks[fb]);
                    out.push_back(f.Output);
                }
                out.push_back(f.Out1);
                out.push_back(f.Out2);
            }
        }
        write_file(dir, "lowpass3.f64", out.data(), out.size() * sizeof(double));
    }
    return 0;
}
