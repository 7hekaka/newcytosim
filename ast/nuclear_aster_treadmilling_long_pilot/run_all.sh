#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
MODE="${1:-production}"

cd "$HERE"

if [[ "$MODE" == "smoke" ]]; then
    CONFIGS=(00_smoke_parse/config.cym)
else
    mapfile -t CONFIGS < <(find production -path '*/r[0-9][0-9][0-9][0-9]/config.cym' | sort)
fi

if [[ "${#CONFIGS[@]}" -eq 0 ]]; then
    echo "No configs found. Generate first:" >&2
    echo "  python3 ast/nuclear_aster_treadmilling_long_pilot/generate_treadmilling_long_pilot.py" >&2
    exit 1
fi

for config in "${CONFIGS[@]}"; do
    run_dir="$(dirname "$config")"
    if [[ -s "$run_dir/objects.cmo" ]]; then
        echo "Skipping $run_dir (objects.cmo exists)"
        continue
    fi
    echo "Running $run_dir"
    (cd "$run_dir" && "$ROOT/sim" config.cym)
done
