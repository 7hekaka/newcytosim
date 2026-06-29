#!/usr/bin/env bash
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
sim="${SIM_EXE:-$root/build_cluster/bin/sim}"
submit="${SUBMIT:-python3 submit_slurm.py}"

if [[ ! -x "$sim" ]]; then
    echo "Missing cluster sim executable: $sim" >&2
    echo "Set SIM_EXE=/path/to/sim if needed." >&2
    exit 1
fi

mapfile -t configs < <(find "$here/production" -mindepth 3 -maxdepth 3 -name config.cym | sort)

if [[ "${#configs[@]}" -eq 0 ]]; then
    echo "No production configs found. Run generate_mechanism_mockups.py first." >&2
    exit 1
fi

echo "Submitting ${#configs[@]} configs"
if [[ "${DRY_RUN:-0}" == "1" ]]; then
    printf '%s\n' "${configs[@]}"
    exit 0
fi

cd "$here"
$submit "$sim" mem="${MEM_MB:-8192}" cpu="${CPU:-4}" time="${TIME_LIMIT:-1-00:00:00}" \
    queue="${QUEUE:-condo-sabel1}" account="${ACCOUNT:-ACF-UTK0049}" qos="${QOS:-condo}" \
    "${configs[@]}"
