;; Guix manifest for the tools needed to build LikeGT and run the `geno`
;; external-alignment smoke test.
;;
;; `gfainject` is also required by `likegt geno`, but the currently fetchable
;; guix-bioinformatics channel state does not package it.  Install the pinned
;; Cargo-built binary with scripts/install-geno-cargo-tools.sh.

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
