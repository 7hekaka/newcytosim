#!/usr/bin/env bash
set -Eeuo pipefail

# Run this inside an allocated compute node, for example:
#   salloc -p bergamo -w ber1528 --cpus-per-task=8 --mem=8G --time=02:00:00
#   ./python/run/compile_cytosim_cluster.sh

ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$ROOT"

STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="${LOG:-$ROOT/compile_diagnostics_${STAMP}.log}"
BUILD_DIR="${BUILD_DIR:-$ROOT/build_cluster}"
JOBS="${JOBS:-${SLURM_CPUS_PER_TASK:-4}}"
TARGETS="${TARGETS:-sim}"
MAKE_PLAY="${MAKE_PLAY:-OFF}"
MAKE_TOOLS="${MAKE_TOOLS:-OFF}"
MAKE_TESTS="${MAKE_TESTS:-OFF}"
CYTOSIM_NATIVE_ARCH="${CYTOSIM_NATIVE_ARCH:-OFF}"
CYTOSIM_ARCH_FLAGS="${CYTOSIM_ARCH_FLAGS:-}"

exec > >(tee "$LOG") 2>&1

fail()
{
    local status=$?
    echo
    echo "FAILED with status $status"
    echo "Log: $LOG"
    exit "$status"
}
trap fail ERR

echo "=== Cytosim cluster compile diagnostic ==="
echo "date: $(date -Is)"
echo "repo: $ROOT"
echo "git_head: $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "hostname: $(hostname)"
echo "build_dir: $BUILD_DIR"
echo "targets: $TARGETS"
echo "jobs: $JOBS"
echo "native_arch: $CYTOSIM_NATIVE_ARCH"
echo "arch_flags: ${CYTOSIM_ARCH_FLAGS:-<none>}"
echo

echo "=== Node ==="
uname -a || true
lscpu || true
grep -m1 '^flags' /proc/cpuinfo || true
echo

echo "=== Modules ==="
if command -v module >/dev/null 2>&1; then
    if [[ -n "${MODULES:-}" ]]; then
        for mod in $MODULES; do
            echo "module load $mod"
            module load "$mod"
        done
    fi
    module list || true
else
    echo "module command not available"
fi
echo

echo "=== Toolchain ==="
command -v cmake || true
cmake --version || true
command -v g++ || true
g++ --version || true
command -v make || true
make --version | head -3 || true
echo

echo "=== Configure ==="
cmake_args=(
    -S "$ROOT"
    -B "$BUILD_DIR"
    -DCMAKE_BUILD_TYPE=Release
    -DMAKE_SIM=ON
    -DMAKE_PLAY="$MAKE_PLAY"
    -DMAKE_TOOLS="$MAKE_TOOLS"
    -DMAKE_TESTS="$MAKE_TESTS"
    -DCYTOSIM_NATIVE_ARCH="$CYTOSIM_NATIVE_ARCH"
)
if [[ -n "$CYTOSIM_ARCH_FLAGS" ]]; then
    cmake_args+=("-DCYTOSIM_ARCH_FLAGS=$CYTOSIM_ARCH_FLAGS")
fi
printf ' %q' cmake "${cmake_args[@]}"
echo
cmake "${cmake_args[@]}"
echo

echo "=== Build ==="
for target in $TARGETS; do
    echo "--- target: $target ---"
    cmake --build "$BUILD_DIR" --target "$target" --parallel "$JOBS" --verbose
done
echo

echo "=== Binary sanity check ==="
if [[ -x "$BUILD_DIR/bin/sim" ]]; then
    file "$BUILD_DIR/bin/sim" || true
    ldd "$BUILD_DIR/bin/sim" || true
    timeout 20s "$BUILD_DIR/bin/sim" help >/tmp/cytosim_help_${STAMP}.txt
    head -20 /tmp/cytosim_help_${STAMP}.txt
    rm -f /tmp/cytosim_help_${STAMP}.txt
else
    echo "No sim binary found at $BUILD_DIR/bin/sim"
fi

echo
echo "SUCCESS"
echo "Log: $LOG"
