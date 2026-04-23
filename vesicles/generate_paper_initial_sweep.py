#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import generate_paper_domain_pilot as base


OUT = base.ROOT / "vesicles" / "paper_initial_sweep_v2"
RUNS = OUT / "cluster_runs"


def make_condition(
    category: str,
    condition: str,
    condition_label: str,
    *,
    xlink_ratio: int = 0,
    vesicle_count: int = 0,
    motors_per_vesicle: int = 0,
) -> dict:
    data = {
        "category": category,
        "condition": condition,
        "condition_label": condition_label,
        "include_xlinks": 1 if xlink_ratio else 0,
        "xlink_ratio": xlink_ratio,
        "vesicle_count": vesicle_count,
        "motors_per_vesicle": motors_per_vesicle,
    }
    return data


def make_conditions() -> list[dict]:
    conditions = [
        make_condition(
            "actin_only",
            "c00_actin_only",
            "Actin only",
        )
    ]

    for xlink_ratio in base.SWEEP_XLINK_RATIOS:
        conditions.append(
            make_condition(
                "xlink_only",
                f"c01_xlink_only_x{xlink_ratio:02d}",
                f"Crosslinkers only, ratio 1:{xlink_ratio}",
                xlink_ratio=xlink_ratio,
            )
        )

    for vesicle_count in base.SWEEP_VESICLE_COUNTS:
        conditions.append(
            make_condition(
                "inert_vesicles",
                f"c02_inert_vesicles_v{vesicle_count:03d}",
                f"Inert vesicles only, vesicles {vesicle_count}",
                vesicle_count=vesicle_count,
            )
        )

    for xlink_ratio in base.SWEEP_XLINK_RATIOS:
        for vesicle_count in base.SWEEP_VESICLE_COUNTS:
            conditions.append(
                make_condition(
                    "xlink_inert_vesicles",
                    f"c03_xlink_inert_x{xlink_ratio:02d}_v{vesicle_count:03d}",
                    f"Crosslinkers plus inert vesicles, ratio 1:{xlink_ratio}, vesicles {vesicle_count}",
                    xlink_ratio=xlink_ratio,
                    vesicle_count=vesicle_count,
                )
            )

    for vesicle_count in base.SWEEP_VESICLE_COUNTS:
        for motors_per_vesicle in base.SWEEP_MOTORS_PER_VESICLE:
            conditions.append(
                make_condition(
                    "motor_no_xlink",
                    f"c04_motor_no_xlink_v{vesicle_count:03d}_m{motors_per_vesicle:02d}",
                    f"Motorized vesicles without crosslinkers, vesicles {vesicle_count}, motors/vesicle {motors_per_vesicle}",
                    vesicle_count=vesicle_count,
                    motors_per_vesicle=motors_per_vesicle,
                )
            )

    for xlink_ratio in base.SWEEP_XLINK_RATIOS:
        for vesicle_count in base.SWEEP_VESICLE_COUNTS:
            for motors_per_vesicle in base.SWEEP_MOTORS_PER_VESICLE:
                conditions.append(
                    make_condition(
                        "motor_xlink",
                        f"c05_motor_xlink_x{xlink_ratio:02d}_v{vesicle_count:03d}_m{motors_per_vesicle:02d}",
                        (
                            "Motorized vesicles plus crosslinkers, "
                            f"ratio 1:{xlink_ratio}, vesicles {vesicle_count}, motors/vesicle {motors_per_vesicle}"
                        ),
                        xlink_ratio=xlink_ratio,
                        vesicle_count=vesicle_count,
                        motors_per_vesicle=motors_per_vesicle,
                    )
                )

    return conditions


def geometry_fields(domain: dict) -> dict:
    if domain["domain"] == "annulus":
        return {
            "inner_um": domain["inner"],
            "outer_um": domain["outer"],
            "bottom_um": domain["bottom"],
            "top_um": domain["top"],
        }
    if domain["domain"] == "disc":
        return {
            "radius_um": domain["radius"],
            "bottom_um": domain["bottom"],
            "top_um": domain["top"],
        }
    return {
        "outer_radius_um": domain["outer_radius"],
        "inner_radius_um": domain["inner_radius"],
    }


def make_manifest_row(idx: int, run_dir: Path, domain: dict, condition: dict) -> dict:
    xlink_ratio = base.xlink_ratio_for(condition)
    n_filaments = base.filament_count_for(domain)
    n_crosslinkers = base.crosslinker_count_for(domain, condition)
    vesicle_count = base.vesicle_count_for(domain, condition)
    row = {
        "index": idx,
        "run_dir": f"run{idx:04d}",
        "domain": domain["domain"],
        "domain_label": domain["domain_label"],
        "category": condition["category"],
        "condition": condition["condition"],
        "condition_label": condition["condition_label"],
        "volume_um3": f"{domain['volume_um3']:.6f}",
        "density_scale": f"{float(domain.get('density_scale', 1.0)):.6f}",
        "nominal_n_filaments": base.N_FILAMENTS,
        "n_filaments": n_filaments,
        "filament_length_um": base.FILAMENT_LENGTH_UM,
        "include_xlinks": condition["include_xlinks"],
        "xlink_ratio_per_filament": xlink_ratio,
        "n_crosslinkers": n_crosslinkers,
        "vesicle_radius_um": base.VESICLE_RADIUS_UM,
        "nominal_vesicle_count": condition["vesicle_count"],
        "vesicle_count": vesicle_count,
        "motors_per_vesicle": condition["motors_per_vesicle"],
        "total_motors": vesicle_count * condition["motors_per_vesicle"],
        "config_path": str(run_dir / "config.cym"),
    }
    row.update(geometry_fields(domain))
    return row


def write_submit_script() -> None:
    script = """#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/../.. && pwd)"
cd "$ROOT/vesicles/paper_initial_sweep_v2/cluster_runs"

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


def run_label(row: dict) -> str:
    return (
        f"`{row['run_dir']}`: `{row['domain']}` / `{row['condition']}` "
        f"(actual vesicles `{row['vesicle_count']}`)"
    )


def write_readme(rows: list[dict]) -> None:
    counts = Counter(row["category"] for row in rows)
    canonical_rows = [
        row
        for row in rows
        if row["category"] == "motor_xlink"
        and row["xlink_ratio_per_filament"] == base.PILOT_XLINK_RATIO
        and row["nominal_vesicle_count"] == 50
        and row["motors_per_vesicle"] == 32
    ]
    high_load_rows = [
        row
        for row in rows
        if row["category"] == "motor_xlink"
        and row["xlink_ratio_per_filament"] == 16
        and row["nominal_vesicle_count"] == 100
        and row["motors_per_vesicle"] == 64
    ]

    readme = f"""# Vesicle Paper Initial Single-Trajectory Sweep

Full one-trajectory exploratory campaign for the vesicle-bound motor-cluster paper.

## Design
- Domains: annulus, thin density-matched disc, and matched-volume native spherical shell.
- Crosslink ratios: `1:4`, `1:8`, `1:16`.
- Vesicle counts: `10`, `20`, `50`, `100`.
- Motors per vesicle: `8`, `16`, `32`, `64`.
- Disc actual counts are density-scaled by `{base.DISC_DENSITY_SCALE:.5f}` because it is radius `{base.DISC_RADIUS_UM:.3f} um` and height `{base.DISC_HEIGHT_UM:.3f} um`, not volume-matched.
- Full motorized and crosslinked grid: `{counts['motor_xlink']}` configs.
- Fair-control configs: `{len(rows) - counts['motor_xlink']}` configs.
- Total configs: `{len(rows)}`.

## Control Groups
- `actin_only`: `{counts['actin_only']}` configs.
- `xlink_only`: `{counts['xlink_only']}` configs.
- `inert_vesicles`: `{counts['inert_vesicles']}` configs.
- `xlink_inert_vesicles`: `{counts['xlink_inert_vesicles']}` configs.
- `motor_no_xlink`: `{counts['motor_no_xlink']}` configs.
- `motor_xlink`: `{counts['motor_xlink']}` configs.

## First Movies To Inspect
- Canonical full condition, ratio `1:8`, vesicles `50`, motors/vesicle `32`: {', '.join(run_label(row) for row in canonical_rows)}.
- High-load stress condition, ratio `1:16`, vesicles `100`, motors/vesicle `64`: {', '.join(run_label(row) for row in high_load_rows)}.

## Submission
Run from repo root:

```bash
vesicles/paper_initial_sweep_v2/submit_cluster.sh
```

This submits `cluster_runs/run????/config.cym`.
Set `SIM_EXE=/path/to/rebuilt/sim` if the cluster should use a specific rebuilt executable.

## Notes
- This is intentionally single-trajectory triage. Use it to decide whether the final replicated paper matrix should keep all four motor-per-vesicle values or downselect to fixed `32`.
- The spherical shell uses native hard inner/outer confinement and no visible core. Disc and shell spaces are intentionally uncolored/clear for visualization.
- The native `spherical_shell` shape requires a rebuilt Cytosim binary containing `SpaceShell`; older `sim`/`play` binaries will not parse shell configs.
- The full run map is in `manifest.csv`.
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")


def main() -> None:
    RUNS.mkdir(parents=True, exist_ok=True)
    rows = []
    idx = 0
    conditions = make_conditions()
    for domain in base.DOMAINS:
        for condition in conditions:
            run_dir = RUNS / f"run{idx:04d}"
            run_dir.mkdir(parents=True, exist_ok=True)
            (run_dir / "config.cym").write_text(base.render_config(domain, condition), encoding="utf-8")
            rows.append(make_manifest_row(idx, run_dir, domain, condition))
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
