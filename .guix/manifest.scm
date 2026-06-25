;; Guix manifest for the tools needed to build LikeGT and run the `geno`
;; external-alignment smoke test.
;;
;; `gfainject` is also required by `likegt geno`, but the currently fetchable
;; guix-bioinformatics channel state does not package it.  Keep checking for it
;; with scripts/check-geno-tools.sh until we add a local package definition or
;; the upstream channel exposes it again.

(specifications->manifest
 (list
  "gcc-toolchain"
  "cmake"
  "pkg-config"
  "clang"
  "jemalloc"
  "zlib"
  "htslib"
  "samtools"
  "minimap2"
  "bwa"
  "odgi"
  "gafpack"))
