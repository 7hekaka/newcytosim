#!/usr/bin/env python3
from __future__ import annotations

import csv
import importlib.util
import math
import random
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "turnover_continuous_supply_full_comparable"
BASE_SCRIPT = ROOT / "turnover_continuous_supply_scaled" / "generate_continuous_supply_scaled.py"


def load_base():
    spec = importlib.util.spec_from_file_location("continuous_supply_base", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


base = load_base()


@dataclass(frozen=True)
class Counts:
    initial_filaments: int = 256
    pulse_filaments: int = 64
    initial_crosslinkers: int = 1024
    pulse_crosslinkers: int = 256
    chewers_base: int = 1024
    chewers: int = 2048
    clusters: int = 80
    total_polymer: float = 19008.0


SIZE = base.SizeSpec(
    name="full_comparable",
    label="full comparable annulus matching earlier c80/m12 simulations",
    inner_radius=10.5,
    outer_radius=11.0,
    bottom_z=-20.0,
    top_z=20.0,
    filament_length=15.0,
    growth_speed=0.30,
    patch_radius=0.70,
    min_cluster_sep=0.80,
)
COUNTS = Counts()

REPLICATES = 5
N_PULSES = 13
RUN_SECONDS = 60.0
RUN_STEPS = int(RUN_SECONDS / base.TIME_STEP)
FRAMES_PER_BLOCK = 20
MOTORS_PER_CLUSTER = 12

CONDITIONS = [
    ("nomotor_xlink", "no motors, with crosslinkers", False),
    ("rotatable_xlink", "rotatable wall motors, with crosslinkers", True),
]

FIELDS = [
    "kind",
    "size",
    "condition",
    "replicate",
    "path",
    "inner_radius",
    "outer_radius",
    "bottom_z",
    "top_z",
    "chewer_zone_bottom",
    "chewer_zone_top",
    "initial_filaments",
    "pulse_filaments",
    "n_pulses",
    "initial_crosslinkers",
    "pulse_crosslinkers",
    "n_chewers",
    "chewer_count_multiplier",
    "chewing_speed",
    "max_chewing_speed",
    "n_clusters",
    "motors_per_cluster",
    "total_motors",
    "total_polymer",
    "run_seconds",
    "total_seconds",
    "notes",
]


def render_common_blocks() -> str:
    text = base.render_common_blocks(SIZE, COUNTS)
    text = text.replace("max_chewing_speed = 1.5", "max_chewing_speed = 3.0", 1)
    text = text.replace("chewing_speed = 0.8", "chewing_speed = 1.6", 1)
    return text


def run_block(label: str) -> list[str]:
    return [
        f"% {label}",
        f"run {RUN_STEPS} system",
        "{",
        f"    nb_frames = {FRAMES_PER_BLOCK}",
        "}",
        "",
    ]


def render_motor_clusters(rng: random.Random) -> list[str]:
    centers = base.sample_cluster_centers(rng, SIZE, COUNTS.clusters)
    lines = ["% --- BEGIN rotatable membrane-bound motor clusters ---"]
    for cluster, (theta, z0) in enumerate(centers):
        lines.append(f"% cluster {cluster:02d}: theta={theta:.4f} z={z0:.3f}")
        for _ in range(MOTORS_PER_CLUSTER):
            dtheta, dz = base.sample_patch_offset(rng, SIZE)
            theta_i = theta + dtheta
            z = max(SIZE.bottom_z + 0.2, min(SIZE.top_z - 0.2, z0 + dz))
            radius = SIZE.inner_radius + base.WALL_OFFSET
            x = radius * math.cos(theta_i)
            y = radius * math.sin(theta_i)
            lines.append(f"new 1 myosin1 {{ position = {x:.3f} {y:.3f} {z:.3f} }}")
    lines.append("% --- END rotatable membrane-bound motor clusters ---")
    return lines


def render_config(condition: str, has_motors: bool, replicate: int) -> str:
    rng = random.Random(20260625 + 1000 * replicate + (17 if has_motors else 0))
    lines: list[str] = [
        f"% Continuous supply full-comparable case: {SIZE.name} / {condition} / r{replicate:04d}",
        "% Matches the older c80/m12 annulus scale: inner/outer 10.5/11, z -20 to 20.",
        "% Aggressive depolymerization variant: 2x diffuse bottom chewers and 2x per-chewer chewing speed.",
        "% Initial aligned network plus pulsed top actin supply and bottom minus-end chewing.",
        "% Endpoints are minus first, plus second; plus end is green/growing.",
        render_common_blocks(),
        "",
        "% Initial full-height aligned filaments.",
    ]

    initial_low = SIZE.bottom_z + SIZE.filament_length + 0.2
    initial_high = SIZE.top_z - 0.2
    for _ in range(COUNTS.initial_filaments):
        lines.append(base.actin_line(rng, SIZE, initial_low, initial_high))

    lines.extend(
        [
            "",
            f"new {COUNTS.chewers} cutter1 {{ position = inside, chewer_zone; }}",
            f"new {COUNTS.initial_crosslinkers} crosslinker {{ position = inside; }}",
            "",
        ]
    )

    if has_motors:
        lines.extend(render_motor_clusters(rng))
        lines.append("")

    lines.extend(run_block("initial 60 s"))

    pulse_low = SIZE.top_z - SIZE.filament_length - 0.25
    pulse_high = SIZE.top_z - 0.05
    for pulse in range(1, N_PULSES + 1):
        lines.append(f"% --- top actin supply pulse {pulse:02d}: +{COUNTS.pulse_filaments} filaments ---")
        for _ in range(COUNTS.pulse_filaments):
            lines.append(base.actin_line(rng, SIZE, pulse_low, pulse_high))
        lines.append(f"new {COUNTS.pulse_crosslinkers} crosslinker {{ position = inside; }}")
        lines.append("")
        lines.extend(run_block(f"post-pulse {pulse:02d}: 60 s"))

    return "\n".join(lines) + "\n"


def write_manifest(rows: list[dict[str, object]]) -> None:
    with (OUT / "manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_submit_script() -> None:
    src = ROOT / "turnover_continuous_supply_scaled" / "submit_cluster.sh"
    dst = OUT / "submit_cluster.sh"
    shutil.copy2(src, dst)
    text = dst.read_text(encoding="utf-8")
    text = text.replace("continuous-supply scaled jobs", "full-comparable aggressive turnover jobs")
    text = text.replace("generate_continuous_supply_scaled.py", "python3 analysis/generate_turnover_full_comparable.py")
    dst.write_text(text, encoding="utf-8")


def write_readme(rows: list[dict[str, object]]) -> None:
    explicit_filament_length = (COUNTS.initial_filaments + N_PULSES * COUNTS.pulse_filaments) * SIZE.filament_length
    lines = [
        "# Continuous Top-Supply Full-Comparable Campaign",
        "",
        "Purpose: rerun the aggressive bottom-turnover model at the older annulus scale so it can be compared directly to the earlier initialized-minus-z c80/m12 simulations.",
        "",
        "## Geometry and density",
        f"- Annulus: inner radius `{SIZE.inner_radius}`, outer radius `{SIZE.outer_radius}`, z range `[{SIZE.bottom_z}, {SIZE.top_z}]`.",
        f"- Initial actin: `{COUNTS.initial_filaments}` filaments, length `{SIZE.filament_length}`.",
        f"- Top supply: `{COUNTS.pulse_filaments}` new aligned filaments every 60 s for `{N_PULSES}` pulses.",
        f"- Explicit filament length introduced by config: `{explicit_filament_length:.1f} um`.",
        f"- Actin polymer pool: `{COUNTS.total_polymer:.1f} um`.",
        f"- Crosslinkers: `{COUNTS.initial_crosslinkers}` initially, `{COUNTS.pulse_crosslinkers}` added per pulse.",
        "",
        "## Turnover and motor model",
        f"- Bottom chewer zone: lower 30% of the annulus, z `[{SIZE.chewer_bottom}, {SIZE.chewer_top}]`.",
        f"- Chewers: `{COUNTS.chewers}` diffuse bottom chewers, equal to 2x the full-size base chewer count.",
        "- `chewing_speed = 1.6` and `max_chewing_speed = 3.0`.",
        f"- Motor case: `{COUNTS.clusters}` inner-wall clusters x `{MOTORS_PER_CLUSTER}` motors/cluster = `{COUNTS.clusters * MOTORS_PER_CLUSTER}` motors.",
        "- Conditions: no motors with crosslinkers; rotatable wall motors with crosslinkers.",
        "",
        "## Runtime",
        f"- Initial block: 60 s.",
        f"- Post-onset supply blocks: `{N_PULSES}` x 60 s = `{N_PULSES * RUN_SECONDS:.0f}` s.",
        f"- Total simulated time: `{RUN_SECONDS * (N_PULSES + 1):.0f}` s.",
        f"- Replicates: `{REPLICATES}` per condition.",
        f"- Total configs: `{len(rows)}`.",
        "",
        "## Submit on cluster",
        "Compile first on a Bergamo node, then:",
        "```bash",
        "cd /lustre/isaac24/scratch/kacheamp/newcytosim/turnover_continuous_supply_full_comparable",
        "./submit_cluster.sh",
        "```",
        "",
        "The submit script defaults to `NODELIST=ber1528,ber1529` to avoid the Milan illegal-instruction issue.",
    ]
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    cluster_dir = OUT / "cluster_runs"
    if cluster_dir.exists():
        shutil.rmtree(cluster_dir)
    OUT.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for condition, _label, has_motors in CONDITIONS:
        for replicate in range(1, REPLICATES + 1):
            rel = Path("cluster_runs") / SIZE.name / condition / f"r{replicate:04d}" / "config.cym"
            path = OUT / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(render_config(condition, has_motors, replicate), encoding="utf-8")
            rows.append(
                {
                    "kind": "cluster",
                    "size": SIZE.name,
                    "condition": condition,
                    "replicate": replicate,
                    "path": rel.as_posix(),
                    "inner_radius": SIZE.inner_radius,
                    "outer_radius": SIZE.outer_radius,
                    "bottom_z": SIZE.bottom_z,
                    "top_z": SIZE.top_z,
                    "chewer_zone_bottom": SIZE.chewer_bottom,
                    "chewer_zone_top": SIZE.chewer_top,
                    "initial_filaments": COUNTS.initial_filaments,
                    "pulse_filaments": COUNTS.pulse_filaments,
                    "n_pulses": N_PULSES,
                    "initial_crosslinkers": COUNTS.initial_crosslinkers,
                    "pulse_crosslinkers": COUNTS.pulse_crosslinkers,
                    "n_chewers": COUNTS.chewers,
                    "chewer_count_multiplier": 2,
                    "chewing_speed": 1.6,
                    "max_chewing_speed": 3.0,
                    "n_clusters": COUNTS.clusters if has_motors else 0,
                    "motors_per_cluster": MOTORS_PER_CLUSTER if has_motors else 0,
                    "total_motors": COUNTS.clusters * MOTORS_PER_CLUSTER if has_motors else 0,
                    "total_polymer": COUNTS.total_polymer,
                    "run_seconds": RUN_SECONDS,
                    "total_seconds": RUN_SECONDS * (N_PULSES + 1),
                    "notes": "full comparable aggressive turnover config",
                }
            )

    write_manifest(rows)
    write_submit_script()
    write_readme(rows)
    shutil.copy2(__file__, OUT / "generate_turnover_full_comparable.py")
    print(f"Wrote {len(rows)} configs under {cluster_dir}")


if __name__ == "__main__":
    main()
