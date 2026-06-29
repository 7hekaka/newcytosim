#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$HERE"

SIM="${SIM:-$ROOT/build_cluster/bin/sim}"
QUEUE="${QUEUE:-condo-sabel1}"
ACCOUNT="${ACCOUNT:-ACF-UTK0049}"
QOS="${QOS:-condo}"
NODELIST="${NODELIST:-ber1528,ber1529}"
HOURS="${HOURS:-48}"
MEM="${MEM:-8192}"
CPU="${CPU:-4}"
WRAPPER="$HERE/run_sim_with_omp.sh"

if [[ ! -x "$SIM" ]]; then
    echo "Missing executable: $SIM" >&2
    echo "Compile first on a Bergamo node, then rerun this submit script." >&2
    exit 1
fi

mapfile -t CONFIGS < <(find cluster_runs -path '*/config.cym' | sort)
if [[ "${#CONFIGS[@]}" -eq 0 ]]; then
    echo "No production configs found under $HERE/cluster_runs" >&2
    echo "Generate them first:" >&2
    echo "  cd $HERE && python3 python3 analysis/generate_turnover_full_comparable.py" >&2
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

echo "Submitting ${#CONFIGS[@]} full-comparable aggressive turnover jobs"
echo "sim: $SIM"
echo "queue=$QUEUE account=$ACCOUNT qos=$QOS nodelist=$NODELIST hours=$HOURS mem=$MEM cpu=$CPU"

python3 "$ROOT/python/run/submit.py" "${ARGS[@]}" "${CONFIGS[@]}"
