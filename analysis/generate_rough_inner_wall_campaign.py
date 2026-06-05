#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUT = ROOT / "clu_rough_inner_wall_init_minusz"

SOURCE_ROOTS = {
    "rotatable": ROOT / "clu_init_minusz",
    "fixed_global": ROOT / "clu_fixed_global_init_minusz",
}
DEFAULT_SOURCE_SUBDIR = Path("total480") / "c40_m12"

RUN_RE = re.compile(r"^r\d{4}$")
SPACE_SHAPE_RE = re.compile(r"(set\s+space\s+stripbox\s*\{\s*shape\s*=\s*)\S+", re.S)
STRIPBOX_RE = re.compile(r"new\s+stripbox\s*\{.*?\n\}", re.S)
POSITION_RE = re.compile(
    r"new\s+1\s+myosin1\s*\{\s*position\s*=\s*"
    r"([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\}"
)

DEFAULT_INNER = 10.5
DEFAULT_OUTER = 11.0
DEFAULT_BOTTOM = -20.0
DEFAULT_TOP = 20.0
DEFAULT_ATTACH_OFFSET = 0.01
DEFAULT_AMPLITUDE = 0.25
DEFAULT_THETA_MODE = 11.0
DEFAULT_Z_MODE = 7.0
DEFAULT_ROUGH_SEED = 17
DEFAULT_ROUGH_COMPONENTS = 12
UINT32_MASK = 0xFFFFFFFF


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def parse_param(block: str, name: str, default: float) -> float:
    match = re.search(rf"\b{name}\s*=\s*([-+0-9.eE]+)", block)
    if not match:
        return default
    return float(match.group(1))


def hash_u32(value: int) -> int:
    value &= UINT32_MASK
    value ^= value >> 16
    value = (value * 0x7FEB352D) & UINT32_MASK
    value ^= value >> 15
    value = (value * 0x846CA68B) & UINT32_MASK
    value ^= value >> 16
    return value & UINT32_MASK


def hash_unit(seed: int, channel: int) -> float:
    mixed = (seed ^ ((0x9E3779B9 * (channel + 1)) & UINT32_MASK)) & UINT32_MASK
    return (hash_u32(mixed) & 0x00FFFFFF) / float(0x01000000)


def bounded_mode(mode: float, fallback: int) -> int:
    if mode < 0:
        return fallback
    return max(0, int(math.floor(mode + 0.5)))


def rough_inner_radius(
    theta: float,
    z: float,
    *,
    inner: float,
    bottom: float,
    top: float,
    amplitude: float,
    theta_mode: float,
    z_mode: float,
    rough_seed: int,
    rough_components: int,
    phase: float,
) -> float:
    if amplitude <= 0:
        return inner

    height = top - bottom
    z_phase = 0.0 if abs(height) < 1.0e-12 else 2.0 * math.pi * (z - bottom) / height
    theta_max = max(1, bounded_mode(theta_mode, 11))
    theta_min = min(theta_max, 2)
    z_max = bounded_mode(z_mode, 7)
    count = max(1, min(int(rough_components), 64))

    total = 0.0
    norm = 0.0
    for k in range(count):
        base = (int(rough_seed) + 0x85EBCA6B * (k + 1)) & UINT32_MASK
        mt = theta_min + int(hash_unit(base, 0) * (theta_max - theta_min + 1))
        mz = int(hash_unit(base, 1) * (z_max + 1)) if z_max else 0
        direction = -1.0 if hash_unit(base, 2) < 0.5 else 1.0
        component_phase = 2.0 * math.pi * hash_unit(base, 3) + phase
        weight = (0.65 + 0.70 * hash_unit(base, 4)) / math.sqrt(mt + mz + 1.0)
        total += weight * math.sin(mt * theta + direction * mz * z_phase + component_phase)
        norm += weight * weight

    wave = total / math.sqrt(max(norm, 1.0e-12))
    return inner + amplitude * math.tanh(0.9 * wave)


def replace_stripbox(
    text: str,
    *,
    amplitude: float,
    theta_mode: float,
    z_mode: float,
    rough_seed: int,
    rough_components: int,
    phase: float,
) -> tuple[str, dict[str, float]]:
    match = STRIPBOX_RE.search(text)
    if not match:
        raise RuntimeError("could not find `new stripbox` block")

    old_block = match.group(0)
    params = {
        "outer": parse_param(old_block, "outer", DEFAULT_OUTER),
        "inner": parse_param(old_block, "inner", DEFAULT_INNER),
        "top": parse_param(old_block, "top", DEFAULT_TOP),
        "bottom": parse_param(old_block, "bottom", DEFAULT_BOTTOM),
    }
    if params["outer"] <= params["inner"] + amplitude:
        raise RuntimeError(
            "rough wall amplitude is too large for the annulus thickness: "
            f"outer={params['outer']}, inner={params['inner']}, amplitude={amplitude}"
        )

    new_block = (
        "new stripbox {\n"
        f"    outer = {params['outer']:g}\n"
        f"    inner = {params['inner']:g}\n"
        f"    top = {params['top']:g}\n"
        f"    bottom = {params['bottom']:g}\n"
        f"    amplitude = {amplitude:g}\n"
        f"    theta_mode = {theta_mode:g}\n"
        f"    z_mode = {z_mode:g}\n"
        f"    rough_seed = {rough_seed:d}\n"
        f"    rough_components = {rough_components:d}\n"
        f"    phase = {phase:g}\n"
        "}"
    )

    text = text[: match.start()] + new_block + text[match.end() :]
    text, shape_count = SPACE_SHAPE_RE.subn(r"\1rough_annulus", text, count=1)
    if shape_count != 1:
        raise RuntimeError("could not replace stripbox space shape")
    text, display_count = re.subn(
        r"(set\s+space\s+stripbox\s*\{[^{}]*display\s*=\s*\(\s*color\s*=\s*)[^;)]+(\s*;?\s*\)[^{}]*\})",
        r"\g<1>0x9AA0A866\2",
        text,
        count=1,
        flags=re.S,
    )
    if display_count != 1:
        raise RuntimeError("could not replace stripbox display color")
    return text, params


def reproject_motors(
    text: str,
    *,
    inner: float,
    bottom: float,
    top: float,
    amplitude: float,
    theta_mode: float,
    z_mode: float,
    rough_seed: int,
    rough_components: int,
    phase: float,
    attach_offset: float,
) -> tuple[str, int]:
    def repl(match: re.Match[str]) -> str:
        x = float(match.group(1))
        y = float(match.group(2))
        z = float(match.group(3))
        theta = math.atan2(y, x)
        radius = rough_inner_radius(
            theta,
            z,
            inner=inner,
            bottom=bottom,
            top=top,
            amplitude=amplitude,
            theta_mode=theta_mode,
            z_mode=z_mode,
            rough_seed=rough_seed,
            rough_components=rough_components,
            phase=phase,
        ) + attach_offset
        return f"new 1 myosin1 {{ position = {radius * math.cos(theta):.3f} {radius * math.sin(theta):.3f} {z:.3f} }}"

    return POSITION_RE.subn(repl, text)


def update_cluster_comments(
    text: str,
    *,
    inner: float,
    outer: float,
    bottom: float,
    top: float,
    amplitude: float,
    theta_mode: float,
    z_mode: float,
    rough_seed: int,
    rough_components: int,
    phase: float,
    attach_offset: float,
) -> str:
    text = re.sub(
        r"^% inner_radius=.*$",
        (
            f"% rough_inner_wall mean_inner={inner:g} outer={outer:g} "
            f"bottom={bottom:g} top={top:g} amplitude={amplitude:g} "
            f"theta_mode={theta_mode:g} z_mode={z_mode:g} "
            f"rough_seed={rough_seed:d} rough_components={rough_components:d} phase={phase:g} "
            f"attach_offset={attach_offset:g}"
        ),
        text,
        count=1,
        flags=re.M,
    )
    text = text.replace(
        "% stationary surface clusters on INNER wall",
        "% stationary surface clusters on rough INNER wall",
        1,
    )
    text = re.sub(
        r"^% r_place=.*$",
        f"% motor heads projected to rough inner radius + {attach_offset:g} um",
        text,
        count=1,
        flags=re.M,
    )
    return text


def transform_config(
    src: Path,
    dst: Path,
    *,
    amplitude: float,
    theta_mode: float,
    z_mode: float,
    rough_seed: int,
    rough_components: int,
    phase: float,
    attach_offset: float,
) -> dict[str, object]:
    text = read_text(src)
    text, params = replace_stripbox(
        text,
        amplitude=amplitude,
        theta_mode=theta_mode,
        z_mode=z_mode,
        rough_seed=rough_seed,
        rough_components=rough_components,
        phase=phase,
    )
    text, motor_count = reproject_motors(
        text,
        inner=params["inner"],
        bottom=params["bottom"],
        top=params["top"],
        amplitude=amplitude,
        theta_mode=theta_mode,
        z_mode=z_mode,
        rough_seed=rough_seed,
        rough_components=rough_components,
        phase=phase,
        attach_offset=attach_offset,
    )
    if motor_count <= 0:
        raise RuntimeError(f"found no myosin1 positions in {src}")
    text = update_cluster_comments(
        text,
        inner=params["inner"],
        outer=params["outer"],
        bottom=params["bottom"],
        top=params["top"],
        amplitude=amplitude,
        theta_mode=theta_mode,
        z_mode=z_mode,
        rough_seed=rough_seed,
        rough_components=rough_components,
        phase=phase,
        attach_offset=attach_offset,
    )
    write_text(dst, text)
    return {
        "source": src,
        "config": dst,
        "motors": motor_count,
        **params,
    }


def discover_runs(case_dir: Path) -> list[Path]:
    if not case_dir.is_dir():
        raise RuntimeError(f"missing source case directory: {case_dir}")
    runs = sorted(child for child in case_dir.iterdir() if child.is_dir() and RUN_RE.match(child.name))
    if not runs:
        raise RuntimeError(f"no run directories found in {case_dir}")
    return runs


def write_submit_script(out_dir: Path) -> None:
    script = """#!/usr/bin/env bash
set -euo pipefail

cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

QUEUE="${QUEUE:-bergamo}"
HOURS="${HOURS:-24}"
MEM="${MEM:-8192}"
CPU="${CPU:-1}"
NODELIST="${NODELIST:-ber1528,ber1529}"
EXE="${EXE:-./sim}"

mapfile -t CONFIGS < <(find rotatable fixed_global -path '*/config.cym' | sort)
if [[ "${#CONFIGS[@]}" -eq 0 ]]; then
    echo "No configs found under $(pwd)" >&2
    exit 1
fi
if [[ ! -x "$EXE" ]]; then
    echo "Missing executable $EXE; copy or build the rough_annulus sim binary first." >&2
    exit 1
fi

ARGS=("$EXE" "queue=$QUEUE" "hours=$HOURS" "mem=$MEM" "cpu=$CPU" "nodelist=$NODELIST")
if [[ -n "${ACCOUNT:-}" ]]; then
    ARGS+=("account=$ACCOUNT")
fi

python3 ../python/run/submit.py "${ARGS[@]}" "${CONFIGS[@]}"
"""
    path = out_dir / "submit_cluster.sh"
    write_text(path, script)
    path.chmod(0o755)


def write_readme(
    out_dir: Path,
    *,
    amplitude: float,
    theta_mode: float,
    z_mode: float,
    rough_seed: int,
    rough_components: int,
    phase: float,
    source_subdir: Path,
) -> None:
    readme = f"""# Rough Inner Wall Initially Minus-Z Campaign

This campaign is generated from the canonical initially minus-z flat-wall runs:

- `rotatable`: `clu_init_minusz/{source_subdir.as_posix()}`
- `fixed_global`: `clu_fixed_global_init_minusz/{source_subdir.as_posix()}`

Only the domain geometry and motor-head radial placement are changed. The outer wall is smooth,
the inner wall is rigid and spatially imperfect, and motor heads are projected to the local rough
inner radius plus the original attachment offset.

Rough wall parameters:

- amplitude: `{amplitude:g}` um
- theta_mode: `{theta_mode:g}` maximum angular mode
- z_mode: `{z_mode:g}` maximum axial mode
- rough_seed: `{rough_seed:d}`
- rough_components: `{rough_components:d}`
- phase: `{phase:g}`

The wall is a seeded irregular Fourier field, not a single repeated sine wave. New configs also set
the space display color to a translucent grey (`0x9AA0A866`) so `play` shows the wall as a shaded
surface rather than blue contour rings.

The submit helper pins Slurm to the Bergamo nodes by default through `nodelist=ber1528,ber1529`.
Override with `NODELIST=...` if needed.

It expects a rough-annulus-aware executable at `./sim`; override with `EXE=...` if you want to use another path.
The local CMake build uses `-march=native`, so for cluster submission use a binary compiled on, or known to run on,
the target node family. Given the current cluster behavior, that means using a Bergamo-safe `sim`.
"""
    write_text(out_dir / "README.md", readme)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--source-subdir", type=Path, default=DEFAULT_SOURCE_SUBDIR)
    parser.add_argument("--out-subdir", type=Path, default=None)
    parser.add_argument("--amplitude", type=float, default=DEFAULT_AMPLITUDE)
    parser.add_argument("--theta-mode", type=float, default=DEFAULT_THETA_MODE)
    parser.add_argument("--z-mode", type=float, default=DEFAULT_Z_MODE)
    parser.add_argument("--rough-seed", type=int, default=DEFAULT_ROUGH_SEED)
    parser.add_argument("--rough-components", type=int, default=DEFAULT_ROUGH_COMPONENTS)
    parser.add_argument("--phase", type=float, default=0.0)
    parser.add_argument("--attach-offset", type=float, default=DEFAULT_ATTACH_OFFSET)
    args = parser.parse_args()
    out_subdir = args.out_subdir if args.out_subdir is not None else args.source_subdir

    records: list[dict[str, object]] = []
    for label, source_root in SOURCE_ROOTS.items():
        source_dir = source_root / args.source_subdir
        for run_dir in discover_runs(source_dir):
            src = run_dir / "config.cym"
            dst = args.out / label / out_subdir / run_dir.name / "config.cym"
            records.append(
                {
                    "condition": label,
                    "run": run_dir.name,
                    **transform_config(
                        src,
                        dst,
                        amplitude=args.amplitude,
                        theta_mode=args.theta_mode,
                        z_mode=args.z_mode,
                        rough_seed=args.rough_seed,
                        rough_components=args.rough_components,
                        phase=args.phase,
                        attach_offset=args.attach_offset,
                    ),
                }
            )

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "manifest.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "condition",
                "run",
                "motors",
                "inner",
                "outer",
                "bottom",
                "top",
                "source",
                "config",
            ],
        )
        writer.writeheader()
        for record in records:
            writer.writerow({key: record[key] for key in writer.fieldnames})

    write_submit_script(args.out)
    write_readme(
        args.out,
        amplitude=args.amplitude,
        theta_mode=args.theta_mode,
        z_mode=args.z_mode,
        rough_seed=args.rough_seed,
        rough_components=args.rough_components,
        phase=args.phase,
        source_subdir=args.source_subdir,
    )
    print(f"Wrote {len(records)} configs to {args.out}")


if __name__ == "__main__":
    main()
