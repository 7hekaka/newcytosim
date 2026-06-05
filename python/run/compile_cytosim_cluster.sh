#!/usr/bin/env bash
set -Eeuo pipefail

# Run this inside an allocated compute node, for example:
#   salloc --account=ACF-UTK0049 --partition=condo-sabel1 --qos=condo \
#       --nodes=1 --ntasks=1 --cpus-per-task=8 --time=02:00:00
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
MODULES="${MODULES:-netlib-lapack openblas}"
CMAKE_MODULES="${CMAKE_MODULES:-cmake/3.30.5-gcc}"
CMAKE_EXE="${CMAKE_EXE:-${CMAKE:-}}"
PRIVATE_CMAKE="${PRIVATE_CMAKE:-$HOME/.local/opt/cmake-3.30.5-linux-x86_64/bin/cmake}"
LEGACY_PRIVATE_CMAKE="${LEGACY_PRIVATE_CMAKE:-/nfs/home/kacheamp/.local/opt/cmake-3.30.5-linux-x86_64/bin/cmake}"
OPENBLAS_ROOT="${OPENBLAS_ROOT:-$HOME/OpenBLAS}"
LOCAL_LIB_DIR="${LOCAL_LIB_DIR:-$HOME/lib}"
SETUP_OPENBLAS_SHIMS="${SETUP_OPENBLAS_SHIMS:-ON}"

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
echo "requested_modules: ${MODULES:-<none>}"
echo "private_cmake: $PRIVATE_CMAKE"
echo "legacy_private_cmake: $LEGACY_PRIVATE_CMAKE"
echo "openblas_root: $OPENBLAS_ROOT"
echo "local_lib_dir: $LOCAL_LIB_DIR"
echo

echo "=== Node ==="
uname -a || true
lscpu || true
grep -m1 '^flags' /proc/cpuinfo || true
echo

echo "=== Modules ==="
if ! command -v module >/dev/null 2>&1 && [[ -r /etc/profile.d/modules.sh ]]; then
    set +u
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh || true
    set -u
fi
if command -v module >/dev/null 2>&1; then
    if [[ -n "$MODULES" ]]; then
        for mod in $MODULES; do
            echo "module load $mod"
            if ! module load "$mod"; then
                echo "WARNING: failed to load module '$mod'; continuing so local libraries can still be used"
            fi
        done
    fi
    module list || true
else
    echo "module command not available"
fi
echo

echo "=== BLAS/LAPACK setup ==="
prepend_path()
{
    local var_name="$1"
    local dir="$2"
    [[ -d "$dir" ]] || return 0
    local current="${!var_name:-}"
    case ":$current:" in
        *":$dir:"*) ;;
        *) export "$var_name=$dir${current:+:$current}" ;;
    esac
}

find_openblas()
{
    local candidate
    for candidate in \
        "$OPENBLAS_ROOT/libopenblas.a" \
        "$OPENBLAS_ROOT/libopenblas.so" \
        "$OPENBLAS_ROOT/lib/libopenblas.a" \
        "$OPENBLAS_ROOT/lib/libopenblas.so"
    do
        if [[ -e "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

openblas_lib="$(find_openblas || true)"
if [[ "$SETUP_OPENBLAS_SHIMS" != "OFF" && -n "$openblas_lib" ]]; then
    mkdir -p "$LOCAL_LIB_DIR"
    case "$openblas_lib" in
        *.a)
            ln -sfn "$openblas_lib" "$LOCAL_LIB_DIR/libblas.a"
            ln -sfn "$openblas_lib" "$LOCAL_LIB_DIR/liblapack.a"
            ;;
        *.so)
            ln -sfn "$openblas_lib" "$LOCAL_LIB_DIR/libblas.so"
            ln -sfn "$openblas_lib" "$LOCAL_LIB_DIR/liblapack.so"
            ;;
        *)
            ln -sfn "$openblas_lib" "$LOCAL_LIB_DIR/libblas"
            ln -sfn "$openblas_lib" "$LOCAL_LIB_DIR/liblapack"
            ;;
    esac
    echo "openblas_lib: $openblas_lib"
else
    echo "openblas_lib: <not found under $OPENBLAS_ROOT>"
fi

prepend_path LD_LIBRARY_PATH "$LOCAL_LIB_DIR"
prepend_path LD_LIBRARY_PATH "$OPENBLAS_ROOT"
prepend_path LD_LIBRARY_PATH "$OPENBLAS_ROOT/lib"
prepend_path LIBRARY_PATH "$LOCAL_LIB_DIR"
prepend_path LIBRARY_PATH "$OPENBLAS_ROOT"
prepend_path LIBRARY_PATH "$OPENBLAS_ROOT/lib"
prepend_path CMAKE_LIBRARY_PATH "$LOCAL_LIB_DIR"
prepend_path CMAKE_LIBRARY_PATH "$OPENBLAS_ROOT"
prepend_path CMAKE_LIBRARY_PATH "$OPENBLAS_ROOT/lib"
ls -l "$LOCAL_LIB_DIR"/libblas* "$LOCAL_LIB_DIR"/liblapack* 2>/dev/null || true
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-<unset>}"
echo "LIBRARY_PATH=${LIBRARY_PATH:-<unset>}"
echo "CMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH:-<unset>}"
echo

echo "=== Toolchain ==="
if [[ -z "$CMAKE_EXE" ]]; then
    for candidate in "$PRIVATE_CMAKE" "$LEGACY_PRIVATE_CMAKE"; do
        if [[ -x "$candidate" ]]; then
            CMAKE_EXE="$candidate"
            break
        fi
    done
fi

if [[ -z "$CMAKE_EXE" && -n "$CMAKE_MODULES" && "$(command -v module || true)" ]]; then
    for mod in $CMAKE_MODULES; do
        echo "fallback module load $mod"
        if module load "$mod" && command -v cmake >/dev/null 2>&1; then
            CMAKE_EXE="$(command -v cmake)"
            break
        fi
    done
fi

if [[ -z "$CMAKE_EXE" && "$(command -v cmake || true)" ]]; then
    CMAKE_EXE="$(command -v cmake)"
fi

if [[ -z "$CMAKE_EXE" || ! -x "$CMAKE_EXE" ]]; then
    echo "ERROR: no usable CMake found."
    echo "Set CMAKE_EXE=/path/to/cmake or PRIVATE_CMAKE=/path/to/cmake and run again."
    exit 2
fi

echo "cmake_exe: $CMAKE_EXE"
"$CMAKE_EXE" --version || true
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
printf ' %q' "$CMAKE_EXE" "${cmake_args[@]}"
echo
"$CMAKE_EXE" "${cmake_args[@]}"
echo

echo "=== Build ==="
for target in $TARGETS; do
    echo "--- target: $target ---"
    "$CMAKE_EXE" --build "$BUILD_DIR" --target "$target" --parallel "$JOBS" --verbose
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
