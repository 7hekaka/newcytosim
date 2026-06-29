#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

cd "$HERE"
for config in */r*/config.cym; do
    run_dir="$(dirname "$config")"
    if [[ -s "$run_dir/objects.cmo" ]]; then
        echo "Skipping $run_dir (objects.cmo exists)"
        continue
    fi
    echo "Running $run_dir"
    (cd "$run_dir" && "$ROOT/sim" config.cym)
done
