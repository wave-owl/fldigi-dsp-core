#!/bin/sh
# Content checksum of fldigi-dsp-core as pinned by family cores (wave-owl R-920.3).
#
#   sh tools/dsp_checksum.sh [DIR]
#
# Lists every file named in DIR/MANIFEST (sorted) as "<sha256>  <path>", then prints the
# sha256 of that listing. Family cores record this value in their dsp.lock and refuse to
# build when it differs. SPDX-License-Identifier: GPL-3.0-or-later
set -eu
dir=${1:-.}
if command -v sha256sum >/dev/null 2>&1; then
    sha() { sha256sum "$@"; }
else
    sha() { shasum -a 256 "$@"; }
fi
cd "$dir"
LC_ALL=C sort MANIFEST | while read -r f; do
    [ -n "$f" ] || continue
    sha "$f"
done | sha | cut -d' ' -f1
