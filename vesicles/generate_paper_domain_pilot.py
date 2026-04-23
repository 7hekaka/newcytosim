#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "vesicles" / "paper_domain_pilot_v2"
RUNS = OUT / "cluster_runs"

N_FILAMENTS = 256
FILAMENT_LENGTH_UM = 15.0
VESICLE_RADIUS_UM = 0.4
PILOT_VESICLE_COUNT = 20
MOTORS_PER_VESICLE = 32
PILOT_XLINK_RATIO = 8
SWEEP_XLINK_RATIOS = (4, 8, 16)
SWEEP_VESICLE_COUNTS = (10, 20, 50, 100)
SWEEP_MOTORS_PER_VESICLE = (8, 16, 32, 64)

ANNULUS_INNER_UM = 10.0
ANNULUS_OUTER_UM = 11.0
ANNULUS_HEIGHT_UM = 40.0
ANNULUS_VOLUME_UM3 = math.pi * (ANNULUS_OUTER_UM**2 - ANNULUS_INNER_UM**2) * ANNULUS_HEIGHT_UM

DISC_RADIUS_UM = 0.5 * (ANNULUS_INNER_UM + ANNULUS_OUTER_UM)
DISC_HEIGHT_UM = 1.0
DISC_VOLUME_UM3 = math.pi * DISC_RADIUS_UM**2 * DISC_HEIGHT_UM
DISC_DENSITY_SCALE = DISC_VOLUME_UM3 / ANNULUS_VOLUME_UM3

SHELL_THICKNESS_UM = 1.0
SHELL_TARGET = 3.0 * ANNULUS_VOLUME_UM3 / (4.0 * math.pi)
SHELL_OUTER_UM = (3.0 + math.sqrt(9.0 - 12.0 * (1.0 - SHELL_TARGET))) / 6.0
SHELL_INNER_UM = SHELL_OUTER_UM - SHELL_THICKNESS_UM


DOMAINS = [
    {
        "domain": "annulus",
        "domain_label": "Annulus slab",
        "shape": "annulus",
        "inner": ANNULUS_INNER_UM,
        "outer": ANNULUS_OUTER_UM,
        "bottom": -0.5 * ANNULUS_HEIGHT_UM,
        "top": 0.5 * ANNULUS_HEIGHT_UM,
        "volume_um3": ANNULUS_VOLUME_UM3,
        "density_scale": 1.0,
    },
    {
        "domain": "disc",
        "domain_label": "Thin density-matched cylinder",
        "shape": "cylinderZ",
        "radius": DISC_RADIUS_UM,
        "bottom": -0.5 * DISC_HEIGHT_UM,
        "top": 0.5 * DISC_HEIGHT_UM,
        "volume_um3": DISC_VOLUME_UM3,
        "density_scale": DISC_DENSITY_SCALE,
    },
    {
        "domain": "spherical_shell",
        "domain_label": "Matched-volume spherical shell",
        "shape": "spherical_shell",
        "outer_radius": SHELL_OUTER_UM,
        "inner_radius": SHELL_INNER_UM,
        "volume_um3": ANNULUS_VOLUME_UM3,
        "density_scale": 1.0,
    },
]

CONDITIONS = [
    {
        "condition": "c00_actin_only",
        "condition_label": "Actin only",
        "include_xlinks": 0,
        "vesicle_count": 0,
        "motors_per_vesicle": 0,
    },
    {
        "condition": "c01_inert_vesicles",
        "condition_label": "Inert vesicles, no crosslinkers",
        "include_xlinks": 0,
        "vesicle_count": PILOT_VESICLE_COUNT,
        "motors_per_vesicle": 0,
    },
    {
        "condition": "c02_xlink_inert_vesicles",
        "condition_label": "Crosslinkers plus inert vesicles",
        "include_xlinks": 1,
        "xlink_ratio": PILOT_XLINK_RATIO,
        "vesicle_count": PILOT_VESICLE_COUNT,
        "motors_per_vesicle": 0,
    },
    {
        "condition": "c03_motor_no_xlink",
        "condition_label": "Motorized vesicles, no crosslinkers",
        "include_xlinks": 0,
        "vesicle_count": PILOT_VESICLE_COUNT,
        "motors_per_vesicle": MOTORS_PER_VESICLE,
    },
    {
        "condition": "c04_motor_xlink",
        "condition_label": "Motorized vesicles plus crosslinkers",
        "include_xlinks": 1,
        "xlink_ratio": PILOT_XLINK_RATIO,
        "vesicle_count": PILOT_VESICLE_COUNT,
        "motors_per_vesicle": MOTORS_PER_VESICLE,
    },
]


def block(title: str, body: str) -> str:
    indented = "\n".join(f"    {line}" if line.strip() else "" for line in body.splitlines())
    return f"{title}\n{{\n{indented}\n}}\n"


def render_header() -> str:
    return (
        block(
            "set simul system",
            """steric = 1, 500
time_step = 0.004
kT = 0.0042
viscosity = 0.1""",
        )
        + "set system display { back_color=white }\n"
    )


def render_space(domain: dict) -> str:
    if domain["domain"] == "annulus":
        return block("set space cell", "shape = annulus\ndisplay = ( color = blue; )") + block(
            "new cell",
            f"""outer = {domain['outer']:.6f}
inner = {domain['inner']:.6f}
top = {domain['top']:.6f}
bottom = {domain['bottom']:.6f}""",
        )
    if domain["domain"] == "disc":
        return block("set space cell", "shape = cylinderZ") + block(
            "new cell",
            f"""radius = {domain['radius']:.6f}
top = {domain['top']:.6f}
bottom = {domain['bottom']:.6f}""",
        )
    if domain["domain"] == "spherical_shell":
        return (
            block("set space cell", "shape = spherical_shell")
            + block(
                "new cell",
                f"""outer = {domain['outer_radius']:.6f}
inner = {domain['inner_radius']:.6f}""",
            )
        )
    raise ValueError(f"unknown domain {domain['domain']}")


def placement(domain: dict, radius_margin: float = 0.0) -> str:
    if domain["domain"] == "spherical_shell":
        inner = domain["inner_radius"] + radius_margin
        return f"placement = inside, cell, ( R > {inner:.6f} )"
    return "placement = inside, cell"


def xlink_ratio_for(condition: dict) -> int:
    if not condition["include_xlinks"]:
        return 0
    return int(condition.get("xlink_ratio", PILOT_XLINK_RATIO))


def scaled_count(domain: dict, nominal: int, minimum: int = 0) -> int:
    nominal = int(nominal)
    if nominal <= 0:
        return 0
    count = int(round(nominal * float(domain.get("density_scale", 1.0))))
    return max(minimum, count)


def filament_count_for(domain: dict) -> int:
    return scaled_count(domain, N_FILAMENTS, minimum=1)


def vesicle_count_for(domain: dict, condition: dict) -> int:
    return scaled_count(domain, int(condition["vesicle_count"]), minimum=1)


def crosslinker_count_for(domain: dict, condition: dict) -> int:
    return filament_count_for(domain) * xlink_ratio_for(condition)


def total_motor_count_for(domain: dict, condition: dict) -> int:
    return vesicle_count_for(domain, condition) * int(condition["motors_per_vesicle"])


def render_biology_sets(condition: dict) -> str:
    parts = [
        block(
            "set hand actin_binder",
            """binding_rate = 10
binding_range = 0.150
unbinding_rate = 0.08
unbinding_force = 5
display = ( color=cyan; )""",
        ),
        block(
            "set couple crosslinker",
            """hand1 = actin_binder
hand2 = actin_binder
stiffness = 10
diffusion = 10
specificity = parallel""",
        ),
        block(
            "set fiber actin1",
            """rigidity = 0.075
segmentation = 0.18
steric = 1, 0.025
confine = inside, 200, cell
display = ( color=black; )""",
        ),
    ]
    if condition["motors_per_vesicle"] > 0:
        parts.extend(
            [
                block(
                    "set hand motor",
                    """activity = move
binding = 10, 0.05
unbinding = 0.1, 3
unloaded_speed = 2.0
stall_force = 5
display = ( color=red; size=2 )""",
                ),
                block(
                    "set single motor_on_surface",
                    """activity = wrist
hand = motor
stiffness = 10""",
                ),
            ]
        )
    parts.append(
        block(
            "set solid vesicle",
            """viscosity = 1
steric = 1, 0.025
confine = all_inside, 300, cell
display = ( color=blue; )""",
        )
    )
    return "\n".join(parts)


def render_fibers(domain: dict) -> str:
    n_filaments = filament_count_for(domain)
    return (
        block(
            f"new {n_filaments} actin1",
            f"""length = {FILAMENT_LENGTH_UM:.6f}
position = inside
{placement(domain)}""",
        )
        + "\n"
        + block("run 10000 system", "nb_frames = 100")
    )


def render_xlinks(domain: dict, condition: dict) -> str:
    n_xlinks = crosslinker_count_for(domain, condition)
    if n_xlinks <= 0:
        return block("run 25000 system", "nb_frames = 250")
    return (
        block(
            f"new {n_xlinks} crosslinker",
            f"""position = inside
{placement(domain)}""",
        )
        + "\n"
        + block("run 25000 system", "nb_frames = 250")
    )


def render_vesicles(domain: dict, condition: dict) -> str:
    count = vesicle_count_for(domain, condition)
    if count <= 0:
        return block("run 200000 system", "nb_frames = 1000")
    motors = int(condition["motors_per_vesicle"])
    if motors > 0:
        sphere_line = f"sphere1 = center, {VESICLE_RADIUS_UM:.6f}, {motors} motor_on_surface"
    else:
        sphere_line = f"sphere1 = center, {VESICLE_RADIUS_UM:.6f}"
    return (
        block(
            f"new {count} solid vesicle",
            f"""{sphere_line}
position = inside
{placement(domain, VESICLE_RADIUS_UM)}""",
        )
        + "\n"
        + block("run 200000 system", "nb_frames = 1000")
    )


def render_config(domain: dict, condition: dict) -> str:
    lines = [
        "% Vesicle paper: motor clusters in matched-volume confinement",
        f"% domain={domain['domain']} condition={condition['condition']}",
        f"% volume_um3={domain['volume_um3']:.6f}",
        f"% density_scale={float(domain.get('density_scale', 1.0)):.6f}",
        f"% n_filaments={filament_count_for(domain)} vesicle_count={vesicle_count_for(domain, condition)}",
        "",
        render_header(),
        render_space(domain),
        render_biology_sets(condition),
        render_fibers(domain),
        render_xlinks(domain, condition),
        render_vesicles(domain, condition),
    ]
    return "\n".join(lines)


def write_submit_script() -> None:
    script = """#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/../.. && pwd)"
cd "$ROOT/vesicles/paper_domain_pilot_v2/cluster_runs"

QUEUE="${QUEUE:-icelake}"
HOURS="${HOURS:-24}"
MEM="${MEM:-8192}"
CPU="${CPU:-1}"

if [[ -z "${SIM_EXE:-}" ]]; then
    if [[ -x "$ROOT/build/bin/sim" ]]; then
        SIM_EXE="$ROOT/build/bin/sim"
    else
        SIM_EXE="sim"
    fi
fi

ARGS=("$SIM_EXE" "queue=$QUEUE" "hours=$HOURS" "mem=$MEM" "cpu=$CPU")
if [[ -n "${ACCOUNT:-}" ]]; then
    ARGS+=("account=$ACCOUNT")
fi

python3 "$ROOT/python/run/submit.py" "${ARGS[@]}" run????/config.cym
"""
    path = OUT / "submit_cluster.sh"
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def write_readme(rows: list[dict]) -> None:
    pilot_table = "\n".join(
        f"- `run{int(row['index']):04d}`: `{row['domain']}` / `{row['condition']}`"
        for row in rows
    )
    readme = f"""# Vesicle Paper Domain Pilot

Single-trajectory sanity campaign for the vesicle-bound motor-cluster paper.

## Standardized geometry
- Annulus standard: inner `{ANNULUS_INNER_UM:.3f} um`, outer `{ANNULUS_OUTER_UM:.3f} um`, height `{ANNULUS_HEIGHT_UM:.3f} um`, volume `{ANNULUS_VOLUME_UM3:.3f} um^3`.
- Disc: clear cylinderZ radius `{DISC_RADIUS_UM:.3f} um`, height `{DISC_HEIGHT_UM:.3f} um`, volume `{DISC_VOLUME_UM3:.3f} um^3`, density scale `{DISC_DENSITY_SCALE:.5f}` relative to annulus.
- Spherical shell pilot: native clear `spherical_shell`, outer radius `{SHELL_OUTER_UM:.3f} um`, inner radius `{SHELL_INNER_UM:.3f} um`, thickness `{SHELL_THICKNESS_UM:.3f} um`, matched volume `{ANNULUS_VOLUME_UM3:.3f} um^3`.

## Shared biology
- Filaments: nominal `{N_FILAMENTS}` actin filaments, `{FILAMENT_LENGTH_UM:.1f} um` length. Annulus/shell use `{N_FILAMENTS}`; the density-matched disc uses `{filament_count_for(DOMAINS[1])}`.
- Vesicles: radius `{VESICLE_RADIUS_UM:.1f} um`; nominal pilot count `{PILOT_VESICLE_COUNT}` and `{MOTORS_PER_VESICLE}` motors per vesicle. Annulus/shell use `{PILOT_VESICLE_COUNT}` vesicles; the disc uses `{vesicle_count_for(DOMAINS[1], CONDITIONS[1])}`.
- Crosslink canonical condition: `{PILOT_XLINK_RATIO}` crosslinkers per filament. Annulus/shell use `{N_FILAMENTS * PILOT_XLINK_RATIO}` crosslinkers; the disc uses `{crosslinker_count_for(DOMAINS[1], CONDITIONS[2])}`.

## Conditions
- `c00_actin_only`: actin-only baseline.
- `c01_inert_vesicles`: vesicle steric/crowding control without crosslinkers or motors.
- `c02_xlink_inert_vesicles`: crosslinker-only network control with inert vesicles.
- `c03_motor_no_xlink`: motorized vesicles without crosslinkers.
- `c04_motor_xlink`: full motorized and crosslinked condition.

## Pilot run map
{pilot_table}

## Visual checks after cluster completion
- First compare full-condition runs across domains: annulus `run0004`, disc `run0009`, spherical shell `run0014`.
- Then compare controls within each domain to isolate mechanism: actin-only, inert-vesicle-only, crosslinker-plus-inert-vesicles, motor-without-crosslinkers, and motor-plus-crosslinkers.
- The shell now uses native hard inner/outer shell confinement rather than a visible steric core.

## Fair comparison metrics
- Use domain-normalized compaction metrics instead of raw center-of-mass motion: occupied-volume fraction, radius of gyration normalized by domain scale, pairwise-distance contraction, and density heterogeneity.
- Report active transport relative to matched passive controls: motorized minus inert-vesicle and crosslinker-only baselines.
- Separate connectivity and motor effects using no-crosslinker motorized runs and crosslinker-plus-inert-vesicle runs.
- Track engagement of vesicle-bound motors using bound motor count, vesicle-filament contact count, and fraction of vesicles with at least one bound motor.

## Production plan after this pilot
- If all geometries look numerically stable, expand to the initial single-trajectory matrix: 3 domains, crosslink ratios `1:4`, `1:8`, `1:16`, vesicle counts `10`, `20`, `50`, `100`, and motors per vesicle `8`, `16`, `32`, `64`.
- If the motor-per-vesicle axis remains weaker than vesicle count, downselect the replicated paper matrix to fixed `32` motors per vesicle.
- Keep controls in every domain: actin-only, inert vesicles, crosslinkers plus inert vesicles, motorized vesicles without crosslinkers, and no-motor crosslinked networks.
- Use single trajectories only for triage; the paper-scale run should use replicates and medoid selection for representative movies.

## Submission
Run from repo root:

```bash
vesicles/paper_domain_pilot_v2/submit_cluster.sh
```

This submits `{len(rows)}` configs from `cluster_runs/run????/config.cym`.
Set `SIM_EXE=/path/to/rebuilt/sim` if the cluster should use a specific rebuilt executable.

## Notes
- Only the annulus space is colored blue because its custom drawer shows outlines. Disc and shell spaces are intentionally left uncolored/clear so filaments remain visible.
- The native `spherical_shell` shape requires a rebuilt Cytosim binary containing `SpaceShell`; older `sim`/`play` binaries will not parse shell configs.
- The current annulus example in repo uses `256` filaments of length `15 um`; this pilot preserves that density convention for annulus/shell and density-scales the smaller disc.
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")


def main() -> None:
    RUNS.mkdir(parents=True, exist_ok=True)
    rows = []
    idx = 0
    for domain in DOMAINS:
        for condition in CONDITIONS:
            run_dir = RUNS / f"run{idx:04d}"
            run_dir.mkdir(parents=True, exist_ok=True)
            (run_dir / "config.cym").write_text(render_config(domain, condition), encoding="utf-8")
            xlink_ratio = xlink_ratio_for(condition)
            n_filaments = filament_count_for(domain)
            n_xlinks = crosslinker_count_for(domain, condition)
            vesicle_count = vesicle_count_for(domain, condition)
            row = {
                "index": idx,
                "run_dir": f"run{idx:04d}",
                "domain": domain["domain"],
                "domain_label": domain["domain_label"],
                "condition": condition["condition"],
                "condition_label": condition["condition_label"],
                "volume_um3": f"{domain['volume_um3']:.6f}",
                "density_scale": f"{float(domain.get('density_scale', 1.0)):.6f}",
                "nominal_n_filaments": N_FILAMENTS,
                "n_filaments": n_filaments,
                "filament_length_um": FILAMENT_LENGTH_UM,
                "include_xlinks": condition["include_xlinks"],
                "xlink_ratio_per_filament": xlink_ratio,
                "n_crosslinkers": n_xlinks,
                "vesicle_radius_um": VESICLE_RADIUS_UM,
                "nominal_vesicle_count": condition["vesicle_count"],
                "vesicle_count": vesicle_count,
                "motors_per_vesicle": condition["motors_per_vesicle"],
                "total_motors": vesicle_count * condition["motors_per_vesicle"],
                "config_path": str(run_dir / "config.cym"),
            }
            if domain["domain"] == "annulus":
                row.update({"inner_um": domain["inner"], "outer_um": domain["outer"], "bottom_um": domain["bottom"], "top_um": domain["top"]})
            elif domain["domain"] == "disc":
                row.update({"radius_um": domain["radius"], "bottom_um": domain["bottom"], "top_um": domain["top"]})
            else:
                row.update({"outer_radius_um": domain["outer_radius"], "inner_radius_um": domain["inner_radius"]})
            rows.append(row)
            idx += 1

    fieldnames = list(rows[0].keys())
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with (OUT / "manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    write_submit_script()
    write_readme(rows)
    print(OUT / "README.md")
    print(OUT / "manifest.csv")
    print(f"wrote {len(rows)} configs under {RUNS}")


if __name__ == "__main__":
    main()
