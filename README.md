# fldigi-dsp-core

[![CI](https://github.com/wave-owl/fldigi-dsp-core/actions/workflows/ci.yml/badge.svg)](https://github.com/wave-owl/fldigi-dsp-core/actions/workflows/ci.yml)

Shared DSP primitives of [fldigi](https://github.com/wave-owl/fldigi) 4.2.13, ported to C99,
for the native mode cores of [wave-owl](https://github.com/wave-owl/wave-owl). Each family
core (`fldigi-olivia-core`, later THOR/DominoEX, MFSK, Throb, FSQ/IFKP, NAVTEX,
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
Version 0.1.0 implements what the first consumer, `fldigi-olivia-core`, needs; the other rows
are reserved and get their own golden tests when the first family that needs them lands
(wave-owl D-926: the first family change that needs a primitive adds it here).

| Primitive | fldigi 4.2.13 source | Status | Consumers |
|---|---|---|---|
| Radix-2 FFT, two-real-sequences split/join (`fdsp_r2fft_*`) | `src/include/jalocha/pj_fft.h:63-317`, `pj_cmpx.h:49-52` | **0.1.0** | Olivia/Contestia |
| Fast Hadamard transform, forward and inverse (`fdsp_fht`, `fdsp_ifht`, `fdsp_ifht_i8`) | `src/include/jalocha/pj_fht.h:11-48` | **0.1.0** | Olivia/Contestia |
| Gray code (`fdsp_gray_code`, `fdsp_binary_code`) | `src/include/jalocha/pj_gray.h:11-21` | **0.1.0** | Olivia/Contestia |
| Third-order IIR integrator (`fdsp_lowpass3_*`) | `src/include/jalocha/pj_lowpass3.h:13-41` | **0.1.0** | Olivia/Contestia |
| Generic FFT (`g_fft`) | `src/include/gfft.h` | planned | MT63 (FFT only), FSQ/IFKP, fftfilt |
| Overlap-add FFT filter (`fftfilt`: low/band-pass, raised cosine, RTTY shape) | `src/filters/fftfilt.cxx` | planned | NAVTEX, Throb, RTTY (migration, optional) |
| FIR filter and Hilbert (`C_FIR_filter`) | `src/filters/filters.cxx:46-250` | planned | MFSK, Throb, FSQ/IFKP |
| Moving average (`Cmovavg`) | `src/filters/filters.cxx:261-310` | planned | MFSK, Hellschreiber, Throb, FSQ/IFKP |
| Sliding DFT | `src/filters/filters.cxx:402-465` | planned | MFSK, THOR/DominoEX |
| K=7 convolutional encoder and soft Viterbi decoder | `src/filters/viterbi.cxx:36-257` | planned | THOR, MFSK |
| Block interleaver | `src/mfsk/interleave.cxx` | planned | THOR, MFSK |
| Varicode tables (MFSK, THOR, DominoEX, FSQ, IFKP) | `src/mfsk/mfskvaricode.cxx`, `src/thor/thorvaricode.cxx`, `src/dominoex/dominovar.cxx`, `src/fsq/fsq_varicode.cxx`, `src/ifkp/ifkp_varicode.cxx` | planned | THOR, DominoEX, MFSK, FSQ, IFKP |

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
| `const char *fdsp_version(void)` | – | `"fldigi-dsp-core 0.1.0 contract=dsp-c1 abi=1"`. `FDSP_VERSION_STRING` holds `"0.1.0"`. |
| `int fdsp_r2fft_init(fdsp_r2fft *f, size_t n)` | `r2FFT::Preset`, `pj_fft.h:92-125` | Tables for an `n`-point FFT, `n` a power of two `>= 4`. Twiddles cover the full circle, `twiddle[k] = e^{+j 2 pi k / n}` (first quarter from `cos`/`sin`, the rest by rotation, as fldigi). Returns 0, or -1 (bad size, out of memory; the struct is then zeroed). |
| `void fdsp_r2fft_free(fdsp_r2fft *f)` | `r2FFT::Free` | Releases the tables; safe to call twice. |
| `void fdsp_r2fft_process(const fdsp_r2fft *f, fdsp_cmpx *x)` | `r2FFT::Process`, `pj_fft.h:249-250` | In-place forward FFT `X[k] = sum x[n] e^{-j 2 pi k n / N}`, bit reversal first, radix-2 butterflies with the conjugated twiddle (`pj_fft.h:210-224, 277-297`). Not normalised. |
| `void fdsp_r2fft_separ_two_reals(f, buf, out0, out1)` | `SeparTwoReals`, `pj_fft.h:161-178` | Spectra (`n/2` points each) of the two real sequences packed as `re`/`im` of `buf`; `out[0].im` holds the Nyquist bin. Scaled like fldigi (the energy of each sequence times `n`, times 2 for the cross terms). |
| `void fdsp_r2fft_join_two_reals(f, in0, in1, buf)` | `JoinTwoReals`, `pj_fft.h:186-203` | The inverse arrangement: after `fdsp_r2fft_process(buf)` the real parts are sequence 0 and the imaginary parts sequence 1, both scaled by `2 n`. |
| `void fdsp_fht(double *x, size_t n)` / `fdsp_ifht(double *x, size_t n)` | `FHT`/`IFHT`, `pj_fht.h:11-48` | In-place fast Walsh-Hadamard transform and its inverse (fldigi's butterfly order); `fdsp_fht(fdsp_ifht(x)) = n x`. |
| `void fdsp_ifht_i8(int8_t *x, size_t n)` | `IFHT<int8_t>`, `pj_mfsk.h:1165` | The inverse transform on `int8_t`, as the Olivia encoder uses it. |
| `uint8_t fdsp_gray_code(uint8_t b)` / `fdsp_binary_code(uint8_t g)` | `pj_gray.h:11-21` | Gray code and its inverse. |
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
version=0.1.0
commit=<git commit of this repository>
checksum=<content checksum of the linked files>
files=include/fldigi_dsp.h src/fdsp_jalocha.c
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

`make test` compares every primitive with golden data produced by fldigi's own headers:
`tests/golden/make_golden.cxx` includes `src/include/jalocha/pj_*.h` of fldigi 4.2.13
(mirror `wave-owl/fldigi`, commit `e7c3ab3709ad62a7d91829364eb48a05fa325c0d`), runs them on
fixed pseudo-random inputs (a 32-bit LCG shared with the C test) and writes the results;
`make golden FLDIGI=/path/to/fldigi` regenerates them.

| Test | Golden | Comparison |
|---|---|---|
| FFT tables, sizes 4 .. 8192 | `fft<n>_twiddle.f64`, `fft<n>_bitrev.u32` | bit-reverse table exact; twiddles within 2^-52 (they come from the platform's `cos`/`sin`, and an optimising compiler may turn a `cos`/`sin` pair into `sincos`) |
| FFT, separate, join, sizes 4 .. 8192 | `fft<n>_out.f64`, `_separ`, `_join` | bit for bit, computed on the golden twiddles |
| FHT / IFHT, lengths 2 .. 128 | `fht<n>.f64` | bit for bit, plus `FHT(IFHT(x)) = n x` |
| IFHT on `int8_t`, every signed unit vector of length 32 and 64 | `ifht_i8_<n>.i8` | exact |
| Gray code | `gray.u8` | exact for all 256 bytes, plus the round trip |
| LowPass3, 3 weights x 2 feedbacks x 2000 steps | `lowpass3.f64` | bit for bit, plus settling to a constant |

Both sides are compiled without floating-point contraction (`-ffp-contract=off`), so the
comparison is between identical operation sequences. CI runs the tests on Ubuntu and macOS,
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
integrator, 1999-2004) with changes by Dave Freese W1HKJ (2019) in fldigi. The C port and the
tests: Gerolf Ziegenhain DG6FL, 2026. See `CHANGELOG.md` for the version history.
