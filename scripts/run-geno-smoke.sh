#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_dir"

if [[ "${LIKEGT_SKIP_ENV_SH:-0}" != 1 && -f "$repo_dir/env.sh" ]]; then
    # Octopus needs this wrapper environment for the Rust/C++ build stack.
    # Set LIKEGT_SKIP_ENV_SH=1 when running in a different prepared shell.
    # shellcheck source=/dev/null
    source "$repo_dir/env.sh"
fi

"$repo_dir/scripts/check-geno-tools.sh" --required-only

cargo test --locked test_geno_pipeline_with_external_bam_if_tools_available -- --nocapture
