# fldigi-dsp-core

[![CI](https://github.com/wave-owl/fldigi-dsp-core/actions/workflows/ci.yml/badge.svg)](https://github.com/wave-owl/fldigi-dsp-core/actions/workflows/ci.yml)

Shared DSP primitives of [fldigi](https://github.com/wave-owl/fldigi) 4.2.13, ported to C99,
for the native mode cores of [wave-owl](https://github.com/wave-owl/wave-owl). Each family
core (`fldigi-olivia-core`, `fldigi-ifk-core` for THOR/DominoEX, later MFSK, Throb, FSQ/IFKP, NAVTEX,
Hellschreiber) consumes this repository as the git submodule `dsp/`, compiles the files it
needs and refuses to build when the submodule is not the pinned version (see
[Pinning](#pinning-and-the-version-gate)).

| | |
|---|---|
| **Is** | A set of small, allocation-explicit C99 functions that reproduce fldigi's own DSP building blocks bit for bit (same operation order, golden tests generated from fldigi's own headers). |
| **Is not** | A modem. Mode logic (framing, synchronisation, character sets) lives in the family cores. Not a general DSP library either: only what fldigi's modems use, in fldigi's arithmetic. |
| | Threaded or stateful behind your back: no globals, no threads, no I/O. Objects that need tables (`fdsp_r2fft`) allocate in their `*_init` and nowhere else. |

## Contents

- [Scope](#scope)
- [API reference](#api-reference)
- [Usage examples](#usage-examples)
- [Pinning and the version gate](#pinning-and-the-version-gate)
- [Verification](#verification)
- [Performance and memory](#performance-and-memory)
- [Building](#building)
- [Maintenance and roadmap](#maintenance-and-roadmap)
- [License and credits](#license-and-credits)

## Scope

The scope follows the shared-DSP table of wave-owl's umbrella change
`add-fldigi-native-modes` (design.md, "fldigi-dsp-core: scope sketch", requirement R-920).
Version 0.1.0 implemented what the first consumer, `fldigi-olivia-core`, needs; 0.2.0 adds
what `fldigi-ifk-core` (THOR and DominoEX) needs, 0.3.0 the rate converter of
`fldigi-mfsk-core` (MFSK11/22 at the station's 8 kHz). The remaining rows are reserved and get their
own golden tests when the first family that needs them lands (wave-owl D-926: the first family
change that needs a primitive adds it here).

| Primitive | fldigi 4.2.13 source | Status | Consumers |
|---|---|---|---|
| Radix-2 FFT, two-real-sequences split/join (`fdsp_r2fft_*`) | `src/include/jalocha/pj_fft.h:63-317`, `pj_cmpx.h:49-52` | **0.1.0** | Olivia/Contestia |
| Fast Hadamard transform, forward and inverse (`fdsp_fht`, `fdsp_ifht`, `fdsp_ifht_i8`) | `src/include/jalocha/pj_fht.h:11-48` | **0.1.0** | Olivia/Contestia |
| Gray code (`fdsp_gray_code`, `fdsp_binary_code`) | `src/include/jalocha/pj_gray.h:11-21` | **0.1.0** | Olivia/Contestia |
| Third-order IIR integrator (`fdsp_lowpass3_*`) | `src/include/jalocha/pj_lowpass3.h:13-41` | **0.1.0** | Olivia/Contestia |
| Generic complex FFT (`fdsp_gfft_*`, John Green's algorithm) | `src/include/gfft.h:142-2290, 3241-3345` | **0.2.0** (complex only) | THOR/DominoEX via fftfilt; later MT63 (FFT only), FSQ/IFKP |
| Overlap-add FFT filter (`fdsp_fftfilt_*`: low/high/band-pass, band reject) | `src/filters/fftfilt.cxx:55-237` | **0.2.0** (raised-cosine and RTTY shapes planned) | THOR/DominoEX; later NAVTEX, Throb |
| FIR filter and Hilbert (`fdsp_fir_*`, `C_FIR_filter`) | `src/filters/filters.cxx:46-187` | **0.2.0** | THOR/DominoEX; later MFSK, Throb, FSQ/IFKP |
| Moving average (`fdsp_movavg_*`, `Cmovavg`) | `src/filters/filters.cxx:261-307` | **0.2.0** | THOR/DominoEX; later MFSK, Hellschreiber, Throb, FSQ/IFKP |
| Sliding DFT (`fdsp_sfft_*`) | `src/filters/filters.cxx:402-465` | **0.2.0** | THOR/DominoEX; later MFSK |
| Convolutional encoder and soft Viterbi decoder, K <= 15 (`fdsp_encoder_*`, `fdsp_viterbi_*`) | `src/filters/viterbi.cxx:36-257` | **0.2.0** | THOR (K=7, K=15), DominoEX MultiPSK FEC; later MFSK |
| Block interleaver (`fdsp_interleave_*`) | `src/mfsk/interleave.cxx:34-100` | **0.2.0** | THOR, DominoEX FEC; later MFSK |
| IZ8BLY MFSK varicode (`fdsp_mfsk_varienc`, `fdsp_mfsk_varidec`) | `src/mfsk/mfskvaricode.cxx:35-349` | **0.2.0** | THOR (primary set), DominoEX FEC; later MFSK |
| `decayavg`, `parity` (`fdsp_decayavg`, `fdsp_parity`) | `src/include/misc.h:59-63`, `src/misc/misc.cxx:34-70` | **0.2.0** | THOR/DominoEX |
| Rational polyphase resampler (`fdsp_resampler_*`) | – (fldigi uses libsamplerate between the sound card and the modem rate) | **0.3.0** | MFSK11/22 (8000 <-> 11025 Hz); later THOR/DominoEX 11025 Hz submodes |
| Varicode tables of FSQ and IFKP | `src/fsq/fsq_varicode.cxx`, `src/ifkp/ifkp_varicode.cxx` | planned | FSQ, IFKP |

THOR's secondary character set (`thorvaricode.cxx`) and the DominoEX nibble varicode
(`dominovar.cxx`) are family-specific and live in `fldigi-ifk-core` (wave-owl D-820).

Not in scope (R-920.4): MT63's own DSP (vendored by `fldigi-mt63-core`) and anything SSTV
(no fldigi lineage). RTTY and PSK keep their already shipped filters; migrating them here is an
optional follow-up.

**API rules for every primitive** (so that one gate and one test style fit all consumers):
the prefix `fdsp_`; plain structs owned by the caller; tables allocated only in
`fdsp_<x>_init()` and released by `fdsp_<x>_free()`; processing functions never allocate; no
global state; every function cites the fldigi lines it ports; every primitive has a golden
test generated from fldigi's own code.

## API reference

All declarations are in `include/fldigi_dsp.h`.

| Function | fldigi | Description |
|---|---|---|
| `const char *fdsp_version(void)` | – | `"fldigi-dsp-core 0.3.0 contract=dsp-c1 abi=1"`. `FDSP_VERSION_STRING` holds `"0.3.0"`. |
| `int fdsp_r2fft_init(fdsp_r2fft *f, size_t n)` | `r2FFT::Preset`, `pj_fft.h:92-125` | Tables for an `n`-point FFT, `n` a power of two `>= 4`. Twiddles cover the full circle, `twiddle[k] = e^{+j 2 pi k / n}` (first quarter from `cos`/`sin`, the rest by rotation, as fldigi). Returns 0, or -1 (bad size, out of memory; the struct is then zeroed). |
| `void fdsp_r2fft_free(fdsp_r2fft *f)` | `r2FFT::Free` | Releases the tables; safe to call twice. |
| `void fdsp_r2fft_process(const fdsp_r2fft *f, fdsp_cmpx *x)` | `r2FFT::Process`, `pj_fft.h:249-250` | In-place forward FFT `X[k] = sum x[n] e^{-j 2 pi k n / N}`, bit reversal first, radix-2 butterflies with the conjugated twiddle (`pj_fft.h:210-224, 277-297`). Not normalised. |
| `void fdsp_r2fft_separ_two_reals(f, buf, out0, out1)` | `SeparTwoReals`, `pj_fft.h:161-178` | Spectra (`n/2` points each) of the two real sequences packed as `re`/`im` of `buf`; `out[0].im` holds the Nyquist bin. Scaled like fldigi (the energy of each sequence times `n`, times 2 for the cross terms). |
| `void fdsp_r2fft_join_two_reals(f, in0, in1, buf)` | `JoinTwoReals`, `pj_fft.h:186-203` | The inverse arrangement: after `fdsp_r2fft_process(buf)` the real parts are sequence 0 and the imaginary parts sequence 1, both scaled by `2 n`. |
| `void fdsp_fht(double *x, size_t n)` / `fdsp_ifht(double *x, size_t n)` | `FHT`/`IFHT`, `pj_fht.h:11-48` | In-place fast Walsh-Hadamard transform and its inverse (fldigi's butterfly order); `fdsp_fht(fdsp_ifht(x)) = n x`. |
| `void fdsp_ifht_i8(int8_t *x, size_t n)` | `IFHT<int8_t>`, `pj_mfsk.h:1165` | The inverse transform on `int8_t`, as the Olivia encoder uses it. |
| `uint8_t fdsp_gray_code(uint8_t b)` / `fdsp_binary_code(uint8_t g)` | `pj_gray.h:11-21` | Gray code and its inverse. |
| `double fdsp_decayavg(double avg, double in, int w)` / `int fdsp_parity(unsigned long w)` | `misc.h:59-63`, `misc.cxx:34-70` | fldigi's exponential average and bit parity. |
| `void fdsp_sincos(double x, double *s, double *c)` | – | sine and cosine with one `sincos` call (macOS `__sincos`), as fldigi's optimised build evaluates a cos/sin pair. |
| `fdsp_gfft_init/free/forward/inverse` | `g_fft<double>`, `gfft.h` | Complex FFT of 16 .. 2^28 points, in place; the inverse is scaled by 1/N (as fldigi's). Tables in `init`. |
| `fdsp_fftfilt_init(f, f1, f2, len)`, `_create`, `_run(f, in, &out)`, `_free` | `fftfilt.cxx:55-237` | Overlap-add FFT filter (`len` FFT points, `len/2` taps of Blackman-windowed sinc, unit peak gain). `run` takes one sample and returns `len/2` when a block of output is ready. |
| `fdsp_fir_init_hilbert/_lowpass/_bandpass`, `fdsp_fir_run`, `fdsp_fir_bp_taps`, `fdsp_fir_free` | `filters.cxx:46-187` | `C_FIR_filter` with optional decimation; the output lags one sample as in fldigi (dot product over the previous `len` samples). |
| `fdsp_movavg_init/reset/run` | `filters.cxx:261-307` | `Cmovavg`; the first value fills the window. |
| `fdsp_sfft_init/reset/run/free` | `filters.cxx:402-465` | Sliding DFT over `len` samples for bins `first .. last-1`, damping `K1 = 0.99999999999L`; `run` writes bin `k` to `result[k * stride]`. |
| `fdsp_encoder_init/encode/free`, `fdsp_viterbi_init/reset/set_traceback/set_chunksize/decode/free` | `viterbi.cxx:36-257` | Rate-1/2 code, K 3..15, soft symbols 0..255 (128 = erasure). `decode` returns the bits of a finished chunk or -1. Path metrics are renormalised on every call (fldigi skips it on traceback calls and overflows after ~8.4 M pairs with chunk size 1); decisions are identical. |
| `fdsp_interleave_init/flush/symbols/bits/free` | `interleave.cxx:34-100` | `size x size x depth` interleaver, forward (TX, zero-filled) or reverse (RX, erasure-filled). |
| `fdsp_mfsk_varienc(c)` / `fdsp_mfsk_varidec(code)` | `mfskvaricode.cxx:332-349` | IZ8BLY varicode: the code string of a byte, the byte of a code (or -1). |
| `fdsp_resampler_init(r, in_rate, out_rate, taps)`, `_run(r, in, n, out)`, `_max_out(r, n)`, `_reset`, `_free` | – (not a port) | Streaming rational resampler: `out/in` reduced to `up/down` (each at most 4096), Kaiser-windowed sinc (beta 8, cutoff 0.45 x the lower rate, `taps x ceil(down/up)` input samples of support), unit DC gain per phase. The output does not depend on the chunking; `_max_out` bounds one call. 32 taps give a flat pass band (< 0.05 dB at 1500 Hz) and > 60 dB alias rejection. |
| `void fdsp_lowpass3_set(fdsp_lowpass3 *f, double level)` | `LowPass3_Filter::Set`, `pj_lowpass3.h:29-41` | All three stages to `level`. |
| `void fdsp_lowpass3_process(fdsp_lowpass3 *f, double in, double w, double fb)` | `LowPass3_Filter::Process`, `pj_lowpass3.h:17-27` | One step of the three-stage integrator; `w` is 1 / peaking time, `fb` the feedback (fldigi's default 0.1: overshoot about 1e-6). The result is `f->output`. |

## Usage examples

### C

```c
/* Spectrum of two real blocks with one FFT, as fldigi's MFSK demodulator does. */
#include <stdio.h>
#include "fldigi_dsp.h"

int main(void)
{
    fdsp_r2fft fft;
    fdsp_cmpx buf[256], s0[128], s1[128];
    int i;
    if (fdsp_r2fft_init(&fft, 256) < 0)
        return 1;
    for (i = 0; i < 256; i++) {
        buf[i].re = fft.twiddle[(10 * i) % 256].re;  /* cos, bin 10 */
        buf[i].im = fft.twiddle[(20 * i) % 256].im;  /* sin, bin 20 */
    }
    fdsp_r2fft_process(&fft, buf);
    fdsp_r2fft_separ_two_reals(&fft, buf, s0, s1);
    printf("|S0[10]|^2 = %g, |S1[20]|^2 = %g\n",
           s0[10].re * s0[10].re + s0[10].im * s0[10].im,
           s1[20].re * s1[20].re + s1[20].im * s1[20].im);
    fdsp_r2fft_free(&fft);
    return 0;
}
```

```
cc -std=c99 -Iinclude example.c src/fdsp_jalocha.c -lm && ./a.out
```

### Python (ctypes)

The library is linked statically into the family cores; for experiments build a shared
object yourself (`cc -shared -fPIC -Iinclude src/fdsp_jalocha.c -o libfdsp.so -lm`):

```python
import ctypes
lib = ctypes.CDLL("./libfdsp.so")
lib.fdsp_version.restype = ctypes.c_char_p
print(lib.fdsp_version().decode())
x = (ctypes.c_double * 8)(1, 0, 0, 0, 0, 0, 0, 0)
lib.fdsp_fht(x, ctypes.c_size_t(8))
print(list(x))          # [1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0] (fldigi's FHT)
```

## Pinning and the version gate

wave-owl requirement R-920.3 (invariant INV-920): a family core checks the `dsp/` submodule
before it compiles anything. Each family core carries a `dsp.lock`:

```
version=0.3.0
commit=<git commit of this repository>
checksum=<content checksum of the linked files>
files=include/fldigi_dsp.h src/fdsp_jalocha.c src/fdsp_gfft.c src/fdsp_filters.c src/fdsp_fec.c src/fdsp_resample.c
```

and its `tools/check_dsp.sh` refuses the build when `dsp/` reports another
`FDSP_VERSION_STRING`, is checked out at another commit, or when the content checksum
differs. The checksum is the sha256 of the sorted list of `"<sha256>  <path>"` lines of the
linked files; `make checksum` prints it for the files in `MANIFEST`
(`tools/dsp_checksum.sh`). Bumping the pin is a deliberate change in the family core:
update the submodule, `dsp.lock` and re-run that core's goldens.

The contract tag `dsp-c1` changes only with an incompatible API change; additions keep it and
raise the minor version.

## Verification

`make test` compares every primitive with golden data produced by fldigi's own code:
`tests/golden/make_golden.cxx` includes `src/include/jalocha/pj_*.h` of fldigi 4.2.13
(mirror `wave-owl/fldigi`, commit `e7c3ab3709ad62a7d91829364eb48a05fa325c0d`), and
`tests/golden/make_golden_fldigi.cxx` compiles fldigi's `fftfilt.cxx`, `filters.cxx`,
`viterbi.cxx`, `interleave.cxx`, `mfskvaricode.cxx` and `misc.cxx` with `gfft.h`. Both run on
fixed pseudo-random inputs (a 32-bit LCG shared with the C test) and write the results;
`make golden FLDIGI=/path/to/fldigi` regenerates them.

| Test | Golden | Comparison |
|---|---|---|
| FFT tables, sizes 4 .. 8192 | `fft<n>_twiddle.f64`, `fft<n>_bitrev.u32` | bit-reverse table exact; twiddles within 2^-52 (they come from the platform's `cos`/`sin`, and an optimising compiler may turn a `cos`/`sin` pair into `sincos`) |
| FFT, separate, join, sizes 4 .. 8192 | `fft<n>_out.f64`, `_separ`, `_join` | bit for bit, computed on the golden twiddles |
| FHT / IFHT, lengths 2 .. 128 | `fht<n>.f64` | bit for bit, plus `FHT(IFHT(x)) = n x` |
| IFHT on `int8_t`, every signed unit vector of length 32 and 64 | `ifht_i8_<n>.i8` | exact |
| Gray code | `gray.u8` | exact for all 256 bytes, plus the round trip |
| LowPass3, 3 weights x 2 feedbacks x 2000 steps | `lowpass3.f64` | bit for bit, plus settling to a constant |
| g_fft forward and inverse, sizes 16 .. 8192 | `gfft<n>.f64`, `gfft<n>_utbl.f64` | cosine table within 2^-52, transforms bit for bit on the golden table |
| fftfilt: THOR16 and DominoEX11 band passes, a low pass, 6000 samples | `fftfilt<c>_shape.f64`, `_out.f64` | shape within 1e-14, every output block bit for bit on the golden shape |
| C_FIR_filter: Hilbert 37, low pass 127, band pass 65 with decimation 2, 9000 samples (buffer wraps) | `fir<c>_taps.f64`, `_out.f64` | taps within 1e-16, outputs bit for bit (golden from `filters.cxx` built like fldigi's reference build, see below) |
| Cmovavg, lengths 8 and 16, with reset | `movavg.f64` | bit for bit |
| sfft (1024/184..248, 512/91..121, 64/3..20), stride 3, 3000 steps | `sfft<c>.f64`, `_vrot.f64` | rotations within 2^-52 (`long double` K1), bins bit for bit on the golden rotations |
| Encoder and Viterbi, K=7 (0x6d, 0x4f) and K=15 (044735, 063057), noisy soft symbols with erasures, chunk 1 (traceback 45/180) and 8 | `conv<k>_*.i32`, `conv<k>_sym.u8` | encoder output, decoded bits and metrics exact; plus 9 million clean pairs without error (the overflow case) |
| Interleaver 4x4x10, 4x4x50, 4x4x4, 5x5x6, symbols both directions and bits | `interleave.u8` | exact |
| MFSK varicode, every byte and every 12-bit code; parity | `mfsk_varienc.txt`, `mfsk_varidec.i32`, `parity.i32` | exact |
| Resampler 8000/11025, 11025/8000, 48000/11025, 8000/8000, 12000/8000 (no fldigi golden: not a port) | – | gain of a 1500 Hz tone within 0.05 dB, leakage and aliasing below -40/-60 dB, output identical for random chunkings, `_max_out` holds |

Both sides are compiled without floating-point contraction (`-ffp-contract=off`), so the
comparison is between identical operation sequences, with one exception: fldigi's
`filters.cxx` (FIR, `Cmovavg`, `sfft`) is built as fldigi's reference build compiles it
(Apple clang `-O2`, contraction on, arm64; generate the goldens there). In that build the
sliding DFT's input step is one fused multiply-add and the FIR dot product is vectorised over
groups of four iterations without fusion, the rest fused. `fdsp_sfft_run` and the FIR's `mac`
write exactly that evaluation out with `fma()`, so they give the reference result with any
compiler, optimisation level and contraction setting.

**Floating-point notes for family cores.** fldigi multiplies `std::complex<double>` with the
four products in separate statements, so they are never fused; this library does the same.
Family cores compile with `-ffp-contract=on` to fuse the same scalar expressions as fldigi's
reference build, and `fdsp_sincos` computes a sine/cosine pair with one `sincos` call as an
optimising compiler does (on macOS its last bit can differ from separate `sin` and `cos`). CI runs the tests on Ubuntu and macOS,
and again under AddressSanitizer and UndefinedBehaviorSanitizer (`make sanitize`).

## Performance and memory

The primitives are fldigi's code paths; there is nothing to tune. An `n`-point
`fdsp_r2fft` holds `n` complex twiddles and `n` indices (24 `n` bytes); 8192 points (the
largest Olivia window) take 192 KiB and one transform runs in about 45 µs on an Apple M4 (`-O3`).

## Building

```
make            # build/libfldigi_dsp.a
make test       # golden tests
make sanitize   # the same under ASan + UBSan
make checksum   # content checksum for a family core's dsp.lock
```

Plain make, C99, `-Wall -Wextra -Werror -pedantic`, no dependencies besides libm. Family
cores do not use this Makefile: they compile the files listed in their `dsp.lock`.

## Maintenance and roadmap

- **fldigi updates**: wave-owl's `scripts/fldigi_upstream_check.py` (skill
  `fldigi-upstream-sync`) checks every `file:line` citation of this repository against a new
  fldigi tag, because the family cores carry it as `third_party/fldigi-*-core/dsp`.
- **New primitives** follow the scope table above: the family change that first needs one
  ports it here with a golden test from fldigi's code, raises the minor version, and pins
  the new version in its own `dsp.lock`. Existing pins are not touched.
- **API status**: the 0.1 API was reviewed against the scope table of
  `add-fldigi-native-modes` (umbrella task 3.2) when `fldigi-olivia-core` adopted it.

## License and credits

GPL-3.0-or-later, as fldigi. The ported code is by Pawel Jalocha SP9VRC (FFT, FHT, Gray code,
integrator, 1999-2004) with changes by Dave Freese W1HKJ (2019) in fldigi; fldigi's filters,
convolutional code, interleaver and MFSK varicode by Dave Freese W1HKJ and the fldigi authors
(after gmfsk by Tomi Manninen OH2BNS; varicode by Nino Porcino IZ8BLY). The complex FFT
(`src/fdsp_gfft.c`) ports fldigi's `gfft.h`, John Green's public-domain FFT as a C++ template by
Dave Freese, LGPL-3.0-or-later in fldigi; LGPL-3.0 permits distributing it as part of this
GPL-3.0-or-later library. The C port and the
tests: Gerolf Ziegenhain DG6FL, 2026. See `CHANGELOG.md` for the version history.
