# Changelog

All notable changes to fldigi-dsp-core. Versions follow semantic versioning; the contract
tag `dsp-c1` changes only with an incompatible API change.

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
