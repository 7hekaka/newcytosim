#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$HERE"

SIM="${SIM:-$ROOT/build_cluster/bin/sim}"
QUEUE="${QUEUE:-condo-sabel1}"
ACCOUNT="${ACCOUNT:-ACF-UTK0049}"
QOS="${QOS:-condo}"
NODELIST="${NODELIST:-ber1528,ber1529}"
HOURS="${HOURS:-72}"
MEM="${MEM:-8192}"
CPU="${CPU:-4}"
DRY_RUN="${DRY_RUN:-0}"
WRAPPER="$HERE/run_sim_with_omp.sh"

mapfile -t CONFIGS < <(find production -path '*/r[0-9][0-9][0-9][0-9]/config.cym' | sort)
if [[ "${#CONFIGS[@]}" -eq 0 ]]; then
    echo "No production configs found under $HERE/production" >&2
    echo "Generate them first:" >&2
    echo "  cd $ROOT && python3 ast/nuclear_aster_treadmilling_long_pilot/generate_treadmilling_long_pilot.py" >&2
    exit 1
fi

if [[ "$DRY_RUN" == "1" ]]; then
    echo "Dry run: would submit ${#CONFIGS[@]} treadmilling-only nuclear-aster jobs"
    printf '%s\n' "${CONFIGS[@]}"
    exit 0
fi

if [[ ! -x "$SIM" ]]; then
    echo "Missing executable: $SIM" >&2
    echo "Compile first, usually from an allocated Bergamo node:" >&2
    echo "  cd $ROOT" >&2
    echo "  ./python/run/compile_cytosim_cluster.sh" >&2
    exit 1
fi

cat > "$WRAPPER" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS="\${OMP_NUM_THREADS:-\${SLURM_CPUS_PER_TASK:-$CPU}}"
export OPENBLAS_NUM_THREADS=1
export OMP_PROC_BIND="\${OMP_PROC_BIND:-spread}"
export OMP_PLACES="\${OMP_PLACES:-cores}"
exec "$SIM" "\$@"
EOF
chmod +x "$WRAPPER"

ARGS=(
    "$WRAPPER"
    "queue=$QUEUE"
    "account=$ACCOUNT"
    "qos=$QOS"
    "hours=$HOURS"
    "mem=$MEM"
    "cpu=$CPU"
    "nodelist=$NODELIST"
)

echo "Submitting ${#CONFIGS[@]} treadmilling-only nuclear-aster jobs"
echo "sim: $SIM"
echo "queue=$QUEUE account=$ACCOUNT qos=$QOS nodelist=$NODELIST hours=$HOURS mem=$MEM cpu=$CPU"

python3 "$ROOT/python/run/submit.py" "${ARGS[@]}" "${CONFIGS[@]}"
