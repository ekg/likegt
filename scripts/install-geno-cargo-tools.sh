#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tool_root="${LIKEGT_TOOL_ROOT:-$repo_dir/target/tools}"
force=0

usage() {
    cat <<'USAGE'
Usage: scripts/install-geno-cargo-tools.sh [--force]

Installs Rust helper binaries used by `likegt geno` into:
  ${LIKEGT_TOOL_ROOT:-target/tools}/bin

Installed tools:
  gfainject
  gafpack
USAGE
}

for arg in "$@"; do
    case "$arg" in
        --force)
            force=1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "unknown argument: $arg" >&2
            usage >&2
            exit 2
            ;;
    esac
done

mkdir -p "$tool_root"

install_tool() {
    local name=$1
    local url=$2
    local rev=$3
    local bin_path="$tool_root/bin/$name"
    local install_force=$force

    if [[ "$force" != 1 && -x "$bin_path" ]]; then
        if compatible_tool "$name" "$bin_path"; then
            echo "found $name at $bin_path"
            return
        fi
        echo "replacing incompatible $name at $bin_path"
        install_force=1
    fi

    echo "installing $name into $tool_root"
    local force_args=()
    if [[ "$install_force" == 1 ]]; then
        force_args=(--force)
    fi
    CARGO_NET_GIT_FETCH_WITH_CLI="${CARGO_NET_GIT_FETCH_WITH_CLI:-true}" \
        cargo install "${force_args[@]}" --locked --git "$url" --rev "$rev" --root "$tool_root" "$name"
}

compatible_tool() {
    local name=$1
    local bin_path=$2

    case "$name" in
        gfainject)
            "$bin_path" --help 2>&1 | grep -q -- 'gfainject --gfa'
            ;;
        gafpack)
            "$bin_path" --help 2>&1 | grep -q -- '--gfa <GFA>'
            ;;
        *)
            return 0
            ;;
    esac
}

install_tool \
    gfainject \
    https://github.com/ekg/gfainject.git \
    f5feb7b6218a7885aba11fa62de111d941cd61f3

install_tool \
    gafpack \
    https://github.com/ekg/gafpack.git \
    b539b0f507698779c89827c82d5979c1b3146333

echo
echo "Add these tools to PATH with:"
echo "  export PATH=\"$tool_root/bin:\$PATH\""
