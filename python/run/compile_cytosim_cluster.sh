#!/usr/bin/env bash
set -Eeuo pipefail

# Run this inside an allocated Bergamo compute node, for example:
#   salloc --account=ACF-UTK0049 --partition=condo-sabel1 --qos=condo \
#       --nodes=1 --ntasks=1 --cpus-per-task=8 --time=03:00:00 \
#       --nodelist=ber1528
#   ./python/run/compile_cytosim_cluster.sh

ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$ROOT"

STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="${LOG:-$ROOT/compile_diagnostics_${STAMP}.log}"
BUILD_DIR="${BUILD_DIR:-$ROOT/build_cluster}"
JOBS="${JOBS:-${SLURM_CPUS_PER_TASK:-4}}"
TARGETS="${TARGETS:-sim report}"
FRESH_BUILD="${FRESH_BUILD:-1}"
MAKE_PLAY="${MAKE_PLAY:-OFF}"
MAKE_TOOLS="${MAKE_TOOLS:-ON}"
MAKE_TESTS="${MAKE_TESTS:-OFF}"
CYTOSIM_NATIVE_ARCH="${CYTOSIM_NATIVE_ARCH:-OFF}"
CYTOSIM_ARCH_FLAGS="${CYTOSIM_ARCH_FLAGS:-}"
CYTOSIM_ENABLE_OPENMP="${CYTOSIM_ENABLE_OPENMP:-ON}"
CHEW_MODE="${CHEW_MODE:-2}"
OPENBLAS_ROOT="${OPENBLAS_ROOT:-$HOME/OpenBLAS}"
USER_CMAKE="${USER_CMAKE:-$HOME/.local/opt/cmake-3.30.5-linux-x86_64/bin/cmake}"

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
echo "fresh_build: $FRESH_BUILD"
echo "native_arch: $CYTOSIM_NATIVE_ARCH"
echo "arch_flags: ${CYTOSIM_ARCH_FLAGS:-<none>}"
echo "openmp: $CYTOSIM_ENABLE_OPENMP"
echo "chew_mode: $CHEW_MODE"
echo

if [[ "$(hostname)" == login* && "${ALLOW_LOGIN_COMPILE:-0}" != "1" ]]; then
    cat <<EOF
This is a login node. Compile from an allocated compute shell instead:

  salloc --account=ACF-UTK0049 --partition=condo-sabel1 --qos=condo \\
    --nodes=1 --ntasks=1 --cpus-per-task=8 --time=03:00:00 \\
    --nodelist=ber1528
  srun --pty bash
  ./python/run/compile_cytosim_cluster.sh

Set ALLOW_LOGIN_COMPILE=1 only if you deliberately want to override this guard.
EOF
    exit 2
fi

echo "=== Node ==="
uname -a || true
lscpu || true
grep -m1 '^flags' /proc/cpuinfo || true
echo

echo "=== Modules ==="
if command -v module >/dev/null 2>&1; then
    MODULES="${MODULES:-gcc netlib-lapack openblas}"
    if [[ -n "$MODULES" ]]; then
        for mod in $MODULES; do
            echo "module load $mod"
            module load "$mod" || echo "warning: could not load module $mod"
        done
    fi
    module list || true
else
    echo "module command not available"
fi
echo

echo "=== Enable requested chewer mode ==="
python3 - <<PY
from pathlib import Path
import re
path = Path("src/sim/fiber_prop.h")
text = path.read_text()
new = re.sub(r"#define\\s+NEW_FIBER_END_CHEW\\s+\\d+", "#define NEW_FIBER_END_CHEW   ${CHEW_MODE}", text, count=1)
if new == text:
    raise SystemExit("could not find NEW_FIBER_END_CHEW in src/sim/fiber_prop.h")
path.write_text(new)
PY
grep -n "NEW_FIBER_END_CHEW" src/sim/fiber_prop.h
echo

echo "=== BLAS/LAPACK paths ==="
mkdir -p "$HOME/lib"
BLAS_CMAKE_ARGS=()
if [[ -f "$OPENBLAS_ROOT/libopenblas.a" ]]; then
    OPENBLAS_LIB="$OPENBLAS_ROOT/libopenblas.a"
elif [[ -f "$OPENBLAS_ROOT/lib/libopenblas.a" ]]; then
    OPENBLAS_LIB="$OPENBLAS_ROOT/lib/libopenblas.a"
elif [[ -f "$OPENBLAS_ROOT/libopenblas.so" ]]; then
    OPENBLAS_LIB="$OPENBLAS_ROOT/libopenblas.so"
elif [[ -f "$OPENBLAS_ROOT/lib/libopenblas.so" ]]; then
    OPENBLAS_LIB="$OPENBLAS_ROOT/lib/libopenblas.so"
else
    OPENBLAS_LIB=""
fi

if [[ -n "$OPENBLAS_LIB" ]]; then
    echo "OpenBLAS: $OPENBLAS_LIB"
    ln -sfn "$OPENBLAS_LIB" "$HOME/lib/libblas.${OPENBLAS_LIB##*.}"
    ln -sfn "$OPENBLAS_LIB" "$HOME/lib/liblapack.${OPENBLAS_LIB##*.}"
    BLAS_CMAKE_ARGS+=("-DBLAS_LIB:FILEPATH=$OPENBLAS_LIB" "-DLAPACK_LIB:FILEPATH=$OPENBLAS_LIB")
else
    echo "OpenBLAS not found under $OPENBLAS_ROOT; CMake will search system/module libraries"
fi
export LD_LIBRARY_PATH="$OPENBLAS_ROOT:$OPENBLAS_ROOT/lib:$HOME/lib:${LD_LIBRARY_PATH:-}"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo

echo "=== Toolchain ==="
if [[ -z "${CC:-}" ]]; then
    if CC_BIN="$(command -v gcc 2>/dev/null)"; then
        export CC="$CC_BIN"
    fi
fi
if [[ -z "${CXX:-}" ]]; then
    if CXX_BIN="$(command -v g++ 2>/dev/null)"; then
        export CXX="$CXX_BIN"
    fi
fi
if [[ -n "${CMAKE:-}" && -x "${CMAKE:-}" ]]; then
    CMAKE_BIN="$CMAKE"
elif [[ -x "$USER_CMAKE" ]]; then
    CMAKE_BIN="$USER_CMAKE"
else
    CMAKE_BIN="$(command -v cmake)"
fi
echo "CMAKE_BIN=$CMAKE_BIN"
"$CMAKE_BIN" --version || true
echo "CC=${CC:-<unset>}"
if [[ -n "${CC:-}" ]]; then "$CC" --version || true; fi
echo "CXX=${CXX:-<unset>}"
if [[ -n "${CXX:-}" ]]; then "$CXX" --version || true; fi
command -v g++ || true
g++ --version || true
command -v make || true
make --version | head -3 || true
echo

echo "=== Configure ==="
if [[ "$FRESH_BUILD" == "1" ]]; then
    echo "Removing old build directory: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
fi

cmake_args=(
    -S "$ROOT"
    -B "$BUILD_DIR"
    -DCMAKE_BUILD_TYPE=Release
    "-DCMAKE_C_COMPILER=${CC:-cc}"
    "-DCMAKE_CXX_COMPILER=${CXX:-c++}"
    -DMAKE_SIM=ON
    -DMAKE_PLAY="$MAKE_PLAY"
    -DMAKE_TOOLS="$MAKE_TOOLS"
    -DMAKE_TESTS="$MAKE_TESTS"
    -DCYTOSIM_NATIVE_ARCH="$CYTOSIM_NATIVE_ARCH"
    -DCYTOSIM_ENABLE_OPENMP="$CYTOSIM_ENABLE_OPENMP"
    "${BLAS_CMAKE_ARGS[@]}"
)
if [[ -n "$CYTOSIM_ARCH_FLAGS" ]]; then
    cmake_args+=("-DCYTOSIM_ARCH_FLAGS=$CYTOSIM_ARCH_FLAGS")
fi
printf ' %q' "$CMAKE_BIN" "${cmake_args[@]}"
echo
"$CMAKE_BIN" "${cmake_args[@]}"
echo

echo "=== Build ==="
for target in $TARGETS; do
    echo "--- target: $target ---"
    "$CMAKE_BIN" --build "$BUILD_DIR" --target "$target" --parallel "$JOBS" --verbose
done
echo

echo "=== Binary sanity check ==="
if [[ -x "$BUILD_DIR/bin/sim" ]]; then
    file "$BUILD_DIR/bin/sim" || true
    ldd "$BUILD_DIR/bin/sim" || true
    ldd "$BUILD_DIR/bin/sim" | grep -Ei 'gomp|omp|blas|lapack|openblas|gfortran|stdc\\+\\+' || true
    timeout 20s "$BUILD_DIR/bin/sim" help >/tmp/cytosim_help_${STAMP}.txt
    head -20 /tmp/cytosim_help_${STAMP}.txt
    rm -f /tmp/cytosim_help_${STAMP}.txt
else
    echo "No sim binary found at $BUILD_DIR/bin/sim"
fi

echo
echo "=== Optional minus-chew sanity run ==="
SANITY_CONF="$ROOT/turnover_mini_minus_chew/validation/minus_end_chew/config.cym"
if [[ -x "$BUILD_DIR/bin/sim" && -x "$BUILD_DIR/bin/report" && -f "$SANITY_CONF" ]]; then
    SANITY_DIR="$BUILD_DIR/sanity_minus_chew_${STAMP}"
    mkdir -p "$SANITY_DIR"
    cp "$SANITY_CONF" "$SANITY_DIR/config.cym"
    (
        cd "$SANITY_DIR"
        export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
        export OPENBLAS_NUM_THREADS=1
        "$BUILD_DIR/bin/sim" config.cym
        "$BUILD_DIR/bin/report" fiber:length frame=0,30 > fiber_length_check.txt
        cat fiber_length_check.txt
    )
else
    echo "skipping sanity run; missing sim/report or $SANITY_CONF"
fi

echo
echo "SUCCESS"
echo "Log: $LOG"
