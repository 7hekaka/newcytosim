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
SYSTEM_FILTER="${SYSTEM_FILTER:-all}"
DRY_RUN="${DRY_RUN:-0}"
WRAPPER="$HERE/run_sim_with_omp.sh"

if [[ ! -x "$SIM" ]]; then
    echo "Missing executable: $SIM" >&2
    echo "Compile first, usually from an allocated Bergamo node:" >&2
    echo "  cd $ROOT" >&2
    echo "  ./python/run/compile_cytosim_cluster.sh" >&2
    exit 1
fi

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

CONFIGS=()
for search_dir in "${SEARCH_DIRS[@]}"; do
    if [[ -d "$search_dir" ]]; then
        while IFS= read -r config; do
            CONFIGS+=("$config")
        done < <(find "$search_dir" -path '*/r[0-9][0-9][0-9][0-9]/config.cym' | sort)
    fi
done

if [[ "${#CONFIGS[@]}" -eq 0 ]]; then
    echo "No replicated configs found under $HERE" >&2
    echo "Generate them first:" >&2
    echo "  cd $ROOT && python3 ast/tethered_mt_free_actin_long_replicates/generate_tethered_mt_free_actin_long_replicates.py" >&2
    exit 1
fi

if [[ "$DRY_RUN" == "1" ]]; then
    echo "Dry run: would submit ${#CONFIGS[@]} replicated tethered-MT/free-actin jobs"
    echo "system_filter=$SYSTEM_FILTER"
    printf '%s\n' "${CONFIGS[@]}"
    exit 0
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

echo "Submitting ${#CONFIGS[@]} replicated tethered-MT/free-actin jobs"
echo "system_filter=$SYSTEM_FILTER"
echo "sim: $SIM"
echo "queue=$QUEUE account=$ACCOUNT qos=$QOS nodelist=$NODELIST hours=$HOURS mem=$MEM cpu=$CPU"

python3 "$ROOT/python/run/submit.py" "${ARGS[@]}" "${CONFIGS[@]}"
