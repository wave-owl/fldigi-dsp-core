# Changelog

All notable changes to fldigi-dsp-core. Versions follow semantic versioning; the contract
tag `dsp-c1` changes only with an incompatible API change.

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
