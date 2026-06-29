#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from generate_rough_inner_wall_campaign import (
    DEFAULT_AMPLITUDE,
    DEFAULT_ATTACH_OFFSET,
    DEFAULT_ROUGH_COMPONENTS,
    DEFAULT_ROUGH_SEED,
    DEFAULT_THETA_MODE,
    DEFAULT_Z_MODE,
    ROOT,
    SOURCE_ROOTS,
    discover_runs,
    filter_runs,
    parse_run_names,
    transform_config,
)


DEFAULT_OUT = ROOT / "clu_mobile_inner_wall_c80_m12_init_minusz"
DEFAULT_SOURCE_SUBDIR = Path("mpc12") / "c80_m12"
MYOSIN_BLOCK_RE = re.compile(r"(set\s+single\s+myosin1\s*\{\n.*?\n\})", re.S)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def set_or_add(block: str, name: str, value: str) -> str:
    pattern = re.compile(rf"(^\s*{re.escape(name)}\s*=\s*).*$", re.M)
    block, count = pattern.subn(rf"\g<1>{value}", block, count=1)
    if count:
        return block

    activity_match = re.search(r"^\s*activity\s*=.*$", block, flags=re.M)
    if activity_match:
        insert_at = activity_match.end()
        return block[:insert_at] + f"\n    {name} = {value}" + block[insert_at:]

    return block[:-2].rstrip() + f"\n    {name} = {value}\n}}"


def make_mobile_anchor(text: str, *, anchor_d: float) -> str:
    def replace_block(match: re.Match[str]) -> str:
        block = match.group(1)
        block = set_or_add(block, "diffusion", "0")
        block = set_or_add(block, "activity", "fixed")
        block = set_or_add(block, "anchor_mode", "slide")
        block = set_or_add(block, "anchor_D", f"{anchor_d:g}")
        block = set_or_add(block, "confine", "surface, 0, stripbox")
        return block

    text, count = MYOSIN_BLOCK_RE.subn(replace_block, text, count=1)
    if count != 1:
        raise RuntimeError("could not find `set single myosin1` block")

    replacements = {
        "% stationary surface clusters on INNER wall": "% initially clustered sliding surface motors on INNER wall",
        "% stationary surface clusters on rough INNER wall": "% initially clustered sliding surface motors on rough INNER wall",
    }
    for old, new in replacements.items():
        text = text.replace(old, new, 1)
    return text


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

mapfile -t CONFIGS < <(find smooth_rotatable smooth_fixed_global rough_rotatable rough_fixed_global -path '*/config.cym' | sort)
if [[ "${#CONFIGS[@]}" -eq 0 ]]; then
    echo "No configs found under $(pwd)" >&2
    exit 1
fi
if [[ ! -x "$EXE" ]]; then
    echo "Missing executable $EXE; copy or build the sliding-anchor/rough-annulus sim binary first." >&2
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
    source_subdir: Path,
    run_start: int | None,
    run_stop: int | None,
    run_names: set[str] | None,
    anchor_d: float,
    amplitude: float,
    theta_mode: float,
    z_mode: float,
    rough_seed: int,
    rough_components: int,
    phase: float,
) -> None:
    if run_names:
        run_text = ", ".join(f"`{name}`" for name in sorted(run_names))
    elif run_start is None and run_stop is None:
        run_text = "all discovered runs"
    elif run_start is None:
        run_text = f"runs up to `r{run_stop:04d}`"
    elif run_stop is None:
        run_text = f"runs from `r{run_start:04d}` onward"
    else:
        run_text = f"`r{run_start:04d}` through `r{run_stop:04d}`"

    readme = f"""# Mobile Inner Wall Initially Minus-Z Campaign

Canonical c80/m12, initially minus-z runs with membrane-localized mobile motor anchors.

Included run range: {run_text}

Source configs:

- rotatable orientation: `clu_init_minusz/{source_subdir.as_posix()}`
- fixed-global orientation: `clu_fixed_global_init_minusz/{source_subdir.as_posix()}`

Generated cases:

- `smooth_rotatable`
- `smooth_fixed_global`
- `rough_rotatable`
- `rough_fixed_global`

The motor single remains `activity=fixed` so it still acts as a force-generating Picket, but it
uses `anchor_mode=slide`, `anchor_D={anchor_d:g}`, and `confine=surface, 0, stripbox` so the anchor
point diffuses on the membrane instead of staying at one permanent wall position.

Rough wall parameters:

- amplitude: `{amplitude:g}` um
- theta_mode: `{theta_mode:g}`
- z_mode: `{z_mode:g}`
- rough_seed: `{rough_seed:d}`
- rough_components: `{rough_components:d}`
- phase: `{phase:g}`

These configs require a Cytosim build that includes both `rough_annulus` and sliding Picket anchors.
The submit helper pins Slurm to Bergamo nodes by default through `nodelist=ber1528,ber1529`.
"""
    write_text(out_dir / "README.md", readme)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--source-subdir", type=Path, default=DEFAULT_SOURCE_SUBDIR)
    parser.add_argument("--out-subdir", type=Path, default=None)
    parser.add_argument("--run-start", type=int, default=None)
    parser.add_argument("--run-stop", type=int, default=None)
    parser.add_argument("--runs", default=None, help="comma-separated explicit run list, e.g. r0002,r0005,r0008")
    parser.add_argument("--anchor-d", type=float, default=0.05)
    parser.add_argument("--amplitude", type=float, default=DEFAULT_AMPLITUDE)
    parser.add_argument("--theta-mode", type=float, default=DEFAULT_THETA_MODE)
    parser.add_argument("--z-mode", type=float, default=DEFAULT_Z_MODE)
    parser.add_argument("--rough-seed", type=int, default=DEFAULT_ROUGH_SEED)
    parser.add_argument("--rough-components", type=int, default=DEFAULT_ROUGH_COMPONENTS)
    parser.add_argument("--phase", type=float, default=0.0)
    parser.add_argument("--attach-offset", type=float, default=DEFAULT_ATTACH_OFFSET)
    args = parser.parse_args()

    if args.anchor_d <= 0:
        raise RuntimeError("--anchor-d must be > 0 for mobile-anchor configs")

    out_subdir = args.out_subdir if args.out_subdir is not None else args.source_subdir
    run_names = parse_run_names(args.runs)
    records: list[dict[str, object]] = []

    for orientation, source_root in SOURCE_ROOTS.items():
        source_dir = source_root / args.source_subdir
        for run_dir in filter_runs(discover_runs(source_dir), args.run_start, args.run_stop, run_names):
            src = run_dir / "config.cym"

            smooth_label = f"smooth_{orientation}"
            smooth_dst = args.out / smooth_label / out_subdir / run_dir.name / "config.cym"
            write_text(smooth_dst, make_mobile_anchor(read_text(src), anchor_d=args.anchor_d))
            records.append(
                {
                    "geometry": "smooth",
                    "orientation": orientation,
                    "condition": smooth_label,
                    "run": run_dir.name,
                    "anchor_D": args.anchor_d,
                    "source": src,
                    "config": smooth_dst,
                }
            )

            rough_label = f"rough_{orientation}"
            rough_dst = args.out / rough_label / out_subdir / run_dir.name / "config.cym"
            rough_record = transform_config(
                src,
                rough_dst,
                amplitude=args.amplitude,
                theta_mode=args.theta_mode,
                z_mode=args.z_mode,
                rough_seed=args.rough_seed,
                rough_components=args.rough_components,
                phase=args.phase,
                attach_offset=args.attach_offset,
            )
            write_text(rough_dst, make_mobile_anchor(read_text(rough_dst), anchor_d=args.anchor_d))
            records.append(
                {
                    "geometry": "rough",
                    "orientation": orientation,
                    "condition": rough_label,
                    "run": run_dir.name,
                    "anchor_D": args.anchor_d,
                    "source": src,
                    "config": rough_dst,
                    "inner": rough_record["inner"],
                    "outer": rough_record["outer"],
                    "bottom": rough_record["bottom"],
                    "top": rough_record["top"],
                }
            )

    args.out.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "geometry",
        "orientation",
        "condition",
        "run",
        "anchor_D",
        "inner",
        "outer",
        "bottom",
        "top",
        "source",
        "config",
    ]
    with (args.out / "manifest.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow({key: record.get(key, "") for key in fieldnames})

    write_submit_script(args.out)
    write_readme(
        args.out,
        source_subdir=args.source_subdir,
        run_start=args.run_start,
        run_stop=args.run_stop,
        run_names=run_names,
        anchor_d=args.anchor_d,
        amplitude=args.amplitude,
        theta_mode=args.theta_mode,
        z_mode=args.z_mode,
        rough_seed=args.rough_seed,
        rough_components=args.rough_components,
        phase=args.phase,
    )
    print(f"Wrote {len(records)} configs to {args.out}")


if __name__ == "__main__":
    main()
