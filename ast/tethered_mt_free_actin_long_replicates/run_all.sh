#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
SYSTEM_FILTER="${SYSTEM_FILTER:-small_system}"

case "$SYSTEM_FILTER" in
    all)
        SEARCH_DIRS=(small_system large_system)
        ;;
    small|small_system)
        SEARCH_DIRS=(small_system)
        ;;
    large|large_system)
        SEARCH_DIRS=(large_system)
        ;;
    *)
        echo "Unknown SYSTEM_FILTER='$SYSTEM_FILTER'. Use all, small, or large." >&2
        exit 1
        ;;
esac

cd "$HERE"
for search_dir in "${SEARCH_DIRS[@]}"; do
    while IFS= read -r config; do
        run_dir="$(dirname "$config")"
        echo "Running ${run_dir}"
        (cd "$run_dir" && "$ROOT/sim" config.cym)
    done < <(find "$search_dir" -path '*/r[0-9][0-9][0-9][0-9]/config.cym' | sort)
done
