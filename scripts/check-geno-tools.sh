#!/usr/bin/env bash
set -u

show_optional=1
strict_optional=0

usage() {
    cat <<'USAGE'
Usage: scripts/check-geno-tools.sh [--required-only] [--strict-optionals]

Checks the external tools used by `likegt geno`.

Required:
  samtools minimap2 gfainject gafpack odgi

Optional:
  bwa        needed for --aligner bwa-mem
  panplexity needed for build-time graph annotation
USAGE
}

for arg in "$@"; do
    case "$arg" in
        --required-only)
            show_optional=0
            ;;
        --strict-optionals)
            strict_optional=1
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

required_tools=(samtools minimap2 gfainject gafpack odgi)
optional_tools=(bwa panplexity)
status=0

version_line() {
    local tool=$1
    case "$tool" in
        samtools)
            "$tool" --version 2>&1 | head -n 1
            ;;
        minimap2)
            printf 'minimap2 '
            "$tool" --version 2>&1 | head -n 1
            ;;
        gfainject)
            "$tool" --help 2>&1 | head -n 1
            ;;
        gafpack)
            "$tool" --version 2>&1 | head -n 1
            ;;
        odgi)
            "$tool" 2>&1 | head -n 1
            ;;
        bwa)
            "$tool" 2>&1 | head -n 1
            ;;
        panplexity)
            "$tool" --version 2>&1 | head -n 1
            ;;
        *)
            "$tool" --version 2>&1 | head -n 1
            ;;
    esac
}

check_tool() {
    local tool=$1
    local required=$2
    local path

    path=$(command -v "$tool" 2>/dev/null || true)
    if [[ -z "$path" ]]; then
        if [[ "$required" == 1 ]]; then
            printf 'missing required: %s\n' "$tool"
            status=1
        else
            printf 'missing optional: %s\n' "$tool"
            if [[ "$strict_optional" == 1 ]]; then
                status=1
            fi
        fi
        return
    fi

    printf 'found %-10s %s\n' "$tool" "$path"
    version_line "$tool" | sed 's/^/  /'
}

echo "Checking required likegt geno tools"
for tool in "${required_tools[@]}"; do
    check_tool "$tool" 1
done

if [[ "$show_optional" == 1 ]]; then
    echo
    echo "Checking optional tools"
    for tool in "${optional_tools[@]}"; do
        check_tool "$tool" 0
    done
fi

if command -v guix >/dev/null 2>&1; then
    echo
    echo "Guix is available at $(command -v guix)"
    if [[ -f .guix/manifest.scm ]]; then
        echo "Guix manifest: .guix/manifest.scm"
        echo "Pinned shell:"
        echo "  guix time-machine -C .guix/channels.scm -- shell -m .guix/manifest.scm"
    fi
fi

exit "$status"
