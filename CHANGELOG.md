# Changelog

All notable changes to fldigi-dsp-core. Versions follow semantic versioning; the contract
tag `dsp-c1` changes only with an incompatible API change.

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
