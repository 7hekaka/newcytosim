#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../" && pwd)"
here="$root/ast/nuclear_propulsion_mechanism_mockups"
sim="$root/build_mini_turnover/bin/sim"

if [[ ! -x "$sim" ]]; then
    echo "Missing executable: $sim" >&2
    exit 1
fi

mode="${1:-smoke}"

case "$mode" in
    smoke)
        find "$here/00_smoke_all_cases" -mindepth 1 -maxdepth 1 -type d | sort | while read -r run; do
            echo "SMOKE $(basename "$run")"
            (cd "$run" && "$sim" config.cym > sim.out 2> sim.err)
        done
        ;;
    quick)
        find "$here/01_quick_20s" -mindepth 2 -maxdepth 2 -type d -name 'r0001' | sort | while read -r run; do
            echo "QUICK $(realpath --relative-to="$here" "$run")"
            (cd "$run" && "$sim" config.cym > sim.out 2> sim.err)
        done
        ;;
    first-reps)
        find "$here/production" -mindepth 2 -maxdepth 2 -type d -name 'r0001' | sort | while read -r run; do
            echo "RUN $(realpath --relative-to="$here" "$run")"
            (cd "$run" && "$sim" config.cym > sim.out 2> sim.err)
        done
        ;;
    all)
        find "$here/production" -mindepth 2 -maxdepth 2 -type d -name 'r*' | sort | while read -r run; do
            echo "RUN $(realpath --relative-to="$here" "$run")"
            (cd "$run" && "$sim" config.cym > sim.out 2> sim.err)
        done
        ;;
    *)
        echo "Usage: $0 [smoke|quick|first-reps|all]" >&2
        exit 2
        ;;
esac
