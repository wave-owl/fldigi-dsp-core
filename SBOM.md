# Source SBOM

Source repository inventory only; completeness is UNKNOWN. Syft inventories tracked source/package manifests. Gitlinks and declared vendored/model evidence are inspected at immutable pins. CMake declarations are conditional requirements, not resolved installations. This is NOT proof of operating-system packages, static-link closure, build-time resolution, or the contents of release binaries.

Generation: **success**

## Source identity

wave-owl/fldigi-dsp-core `b2af2a5b55c6ce73c63506226dbf5d179090c1e8`

- sbom:role: source-root
- sbom:git-commit: b2af2a5b55c6ce73c63506226dbf5d179090c1e8
- sbom:git-tree: 858206f8ff19db5d7e4bb8f11ef763e38b1d67da
- sbom:source-repository: https://github.com/wave-owl/fldigi-dsp-core
- sbom:release-tags: ["v0.9.2"]
- sbom:license-evidence: [{"path": "LICENSE", "sha256": "3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986"}]
- sbom:license-resolution: unresolved; manual review required

## Generation errors

None.

## Components and unresolved licenses

| Component | Version | Role | License |
| --- | --- | --- | --- |
| wave-owl/fldigi-dsp-core | b2af2a5b55c6ce73c63506226dbf5d179090c1e8 | source-root | unresolved; manual review required |
| src/fdsp_jalocha.c | e7c3ab3709ad62a7d91829364eb48a05fa325c0d | declared-vendored-source | unresolved; manual review required |
| src/fdsp_gfft.c | e7c3ab3709ad62a7d91829364eb48a05fa325c0d | declared-vendored-source | unresolved; manual review required |
| src/fdsp_filters.c | e7c3ab3709ad62a7d91829364eb48a05fa325c0d | declared-vendored-source | unresolved; manual review required |
| src/fdsp_fec.c | e7c3ab3709ad62a7d91829364eb48a05fa325c0d | declared-vendored-source | unresolved; manual review required |
| actions/checkout | v4 | Syft source/package discovery | unresolved |
| /Users/gerolfziegenhain/src/wave-owl/local/sbom-core-preflight-v2/fldigi-dsp-core/.source-input/.github/workflows/ci.yml | unresolved | Syft source/package discovery | unresolved |

License evidence paths/hashes, source manifests, model hashes and conditional
CMake declarations are recorded in SBOM.cdx.json. Declared upstream pins do not
prove vendored equivalence. No dependency installation or build was performed.
