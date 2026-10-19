# Changelog

All notable changes to fldigi-dsp-core. Versions follow semantic versioning; the contract
tag `dsp-c1` changes only with an incompatible API change.

## 0.9.1 - 2026-10-19

- The K=16 footprint check of the tests expects the real figure, 1.5 MB (output table
  256 kB, two metric rows 256 kB, decision bits 1 MB); 0.9.0 checked against 1.4 MB and
  failed. No change in the library.

## 0.9.0 - 2026-10-19

Memory for the PSK variants core (wave-owl change `add-fldigi-psk-variants`, R-931.6). The
API is extended; every result is unchanged.

- `fdsp_viterbi` stores what fldigi's results depend on instead of fldigi's two
  PATHMEM x 2^(k-1) int tables: two path-metric rows, one decision bit per state and step and
  the four branch metrics of every step. The metric difference of a traceback is the sum of
  the branch metrics along the traced path (what fldigi's two table reads amount to); rows
  not written since the reset behave as fldigi's zeroed rows. K=16 needs about 1.5 MB instead
  of 64 MB, K=13 160 kB instead of 8 MB. All `conv*` goldens (K=5, 7, 13, 15, 16) and the
  9-million-pair run pass unchanged. `fdsp_viterbi_footprint` reports the allocation. The
  struct fields `metrics`/`history` change meaning (family cores do not read them).
- `fdsp_firh_*`: `C_FIR_filter` with heap buffers of a chosen length (fldigi's
  `FIRBufferLen` only decides when the history is moved, `filters.cxx:180-185`), and
  `fdsp_firh_clear` to restart a filter without allocating. Same outputs as `fdsp_fir` for
  every buffer length (tested on fldigi's FIR goldens with length + 1, length + 7 and 4096).

## 0.8.0 - 2026-10-19

The FEC range `fldigi-psk-core` needs for fldigi's PSK variants (wave-owl change
`add-fldigi-psk-variants`). Additive: the 0.7 API and its results are unchanged.

- `fdsp_viterbi_init` and `fdsp_encoder_init` accept constraint length 16 (fldigi's `K16`
  code of 8PSK125F/8PSK250F, `src/psk/psk.cxx:91-94`, 996-1001); fldigi's `viterbi` class
  has no upper limit, the port stopped at 15.
- Golden data from fldigi's own `encoder`/`viterbi` for K=5 (QPSK, `psk.cxx:66-68`), K=13
  (8PSK FL/500F and xPSK, `psk.cxx:83-85`, traceback 13 x 16) and K=16 (traceback 16 x 12),
  chunk size 1 with metrics and chunk size 8.

## 0.7.0 - 2026-10-18

The primitive `fldigi-navtex-core` needs (wave-owl change `add-fldigi-navtex-core`).
Additive: the 0.6 API and its results are unchanged.

- `fdsp_fftfilt_rtty_filter`: fldigi's raised-cosine RTTY filter shape of
  `src/filters/fftfilt.cxx:243-290` (`f *= 1.4`, cos^2 response equalised by 1 / sinc, linear
  phase -i pi / 2), the mark and space low passes of the NAVTEX/SITOR-B receiver
  (`navtex.cxx:998-1005`). Like fldigi it leaves `filter[flen2]` at the value of the previous
  shape, and evaluates each cos/sin pair with one `sincos` call as the reference build does.
- Golden data from fldigi's own `fftfilt::rtty_filter` (`fftfilt_rtty0/1_shape.f64`,
  `_out.f64`).

## 0.6.1 - 2026-10-18

- The `blackman.f64` golden test compares within 1e-15 off macOS (glibc's `cos()` differs in
  the last bit of about 5 % of the values); bit for bit with the macOS libm as before. No
  functional change.

## 0.6.0 - 2026-10-18

The primitive `fldigi-fsq-core` needs (wave-owl change `add-fldigi-fsq-ifkp-core`).
Additive: the 0.5 API and its results are unchanged.

- `fdsp_blackman`: fldigi's `blackman()` window of `src/include/misc.h:66-69`, with the two
  multiply-adds fldigi's reference build fuses written out as `fma`, so the window is the
  same whatever contraction the consumer's build uses (FSQ and IFKP fill their 4096-point
  FFT window with it, `fsq.cxx:206-207`, `ifkp.cxx:162-163`).
- Golden data from fldigi's own `blackman()` compiled like the reference build
  (`blackman.f64`, `tests/golden/make_golden_window.cxx`).

## 0.5.0 - 2026-10-17

The primitive `fldigi-feldhell-core` needs (wave-owl change `add-fldigi-hellschreiber-core`).
Additive: the 0.4 API and its results are unchanged.

- `fdsp_movavgh_*`: fldigi's `Cmovavg` with a heap window of any length, as Hellschreiber's
  receiver grows it with `setLength(500 * samplerate / rxpixrate)` (`feld.cxx` restart, up to
  91428 samples, beyond the fixed `FDSP_MOVAVG_MAX` window of `fdsp_movavg`). Allocation only
  in `fdsp_movavgh_init`.
- Golden data from fldigi's own `Cmovavg` grown with `setLength` (`movavgh.f64`).

## 0.4.2 - 2026-10-17

- The 0.4.1 entry below no longer quotes the old comment. No functional change.

## 0.4.1 - 2026-10-17

- The `cmac` comment no longer puts the version right after the fldigi line range, which
  wave-owl's citation check parsed as a second range starting at line 0. No functional
  change.

## 0.4.0 - 2026-10-17

The primitives `fldigi-throb-core` needs (wave-owl change `add-fldigi-throb-core`). Additive:
the 0.3 API and its results are unchanged.

- `fdsp_fir_init`: `C_FIR_filter::init` with the caller's taps (`filters.cxx:58-84`); a missing
  Q branch is zero (Throb's sync filter, `throb.cxx:231-232`).
- `fdsp_fir_irun`: `C_FIR_filter::Irun` (`filters.cxx:193-217`), the I branch alone.
- `fdsp_cmac`: the circular complex multiply-accumulate of `complex.h:33-41` (Throb's tone
  correlation), products evaluated as `std::complex` does.
- Golden data from fldigi's own `filters.cxx` and `complex.h` (`firirun*.f64`, `cmac.f64`).

## 0.3.0 - 2026-10-17

The rate converter for `fldigi-mfsk-core` (wave-owl change `add-fldigi-mfsk-core`). Additive:
the 0.2 API and its results are unchanged.

- `fdsp_resampler_*`: streaming rational polyphase resampler (Kaiser-windowed sinc, unit DC
  gain per phase, chunking-independent output). Not an fldigi port: fldigi converts between
  the sound card and the modem rate with libsamplerate. It lets the 11025 Hz modems (MFSK11/22,
  later THOR/DominoEX 5/11/22/44/88) run behind the station's 8 kHz audio (umbrella
  `add-fldigi-native-modes`, open question 2: one shared copy instead of one per core).
- New file `src/fdsp_resample.c` (in `MANIFEST`, so family cores pin it in `dsp.lock`).

## 0.2.0 - 2026-10-17

Primitives for `fldigi-ifk-core` (THOR and DominoEX; wave-owl change
`add-fldigi-thor-dominoex-core`). Additive: the 0.1 API and its results are unchanged.

- fldigi's complex FFT `g_fft<double>` (John Green's algorithm, `gfft.h`), converted from the
  template mechanically with the same operation order.
- `fftfilt` overlap-add filter (`fftfilt.cxx:55-237`), `C_FIR_filter` incl. Hilbert pair
  (`filters.cxx:46-187`), `Cmovavg` (`:261-307`), sliding DFT `sfft` (`:402-465`).
- Rate-1/2 convolutional encoder and soft-decision Viterbi decoder for K up to 15
  (`viterbi.cxx`), with the path-metric renormalisation applied on every call (fldigi's
  metrics overflow after ~8.4 M symbol pairs with chunk size 1; decisions unchanged).
- Block interleaver (`interleave.cxx`) and the IZ8BLY MFSK varicode (`mfskvaricode.cxx`);
  `decayavg` and `parity`.
- Golden data from fldigi's own sources (`tests/golden/make_golden_fldigi.cxx`); tables that
  come from libm are compared within a few ulp and injected, the processing bit for bit.
- `std::complex` products are never fused (as in fldigi); the FIR dot product and the sliding
  DFT input step follow fldigi's reference build (Apple clang -O2, arm64) explicitly with
  `fma()`; `fdsp_sincos` evaluates sine/cosine pairs with one `sincos` call. Family cores built
  with contraction on reproduce fldigi's binary bit for bit at every optimisation level.

## 0.1.0 - 2026-10-17

First release, for `fldigi-olivia-core` (wave-owl change `add-fldigi-olivia-core`, umbrella
`add-fldigi-native-modes`, R-920).

- Pawel Jalocha's radix-2 FFT with the two-real-sequences split and join
  (`pj_fft.h:63-317`), fast Hadamard transform and inverse including the `int8_t` variant
  (`pj_fht.h:11-48`), Gray code (`pj_gray.h:11-21`) and the three-stage integrator
  `LowPass3_Filter` (`pj_lowpass3.h:13-41`), all from fldigi 4.2.13.
- Golden tests generated from fldigi's own headers (`tests/golden/make_golden.cxx`), bit for
  bit; CI on Ubuntu and macOS with AddressSanitizer and UndefinedBehaviorSanitizer.
- `MANIFEST` and `tools/dsp_checksum.sh`: the content checksum that family cores pin in
  their `dsp.lock`.
- Scope table for the primitives the other fldigi families will add.
