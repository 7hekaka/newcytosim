#!/usr/bin/env python3
"""Generate aggressive-depolymerization variants of the continuous-supply campaign."""

from __future__ import annotations

import csv
import importlib.util
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
BASE_DIR = ROOT / "turnover_continuous_supply_scaled"
BASE_SCRIPT = BASE_DIR / "generate_continuous_supply_scaled.py"


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
class Variant:
    name: str
    out_dir: str
    chewer_count_multiplier: int
    chewing_speed: float
    max_chewing_speed: float
    description: str
    severer_count_multiplier: int = 0
    severer_cutting_rate: float = 0.0
    severer_cut_width: float = 0.18
    severer_line_diffusion: float = 0.2


VARIANTS = [
    Variant(
        name="chew_rate2x",
        out_dir="turnover_continuous_supply_scaled_chew_rate2x",
        chewer_count_multiplier=1,
        chewing_speed=1.6,
        max_chewing_speed=3.0,
        description="same chewer count, 2x hand chewing speed with a matching fiber max-chewing cap",
    ),
    Variant(
        name="chewer_count2x",
        out_dir="turnover_continuous_supply_scaled_chewer_count2x",
        chewer_count_multiplier=2,
        chewing_speed=0.8,
        max_chewing_speed=1.5,
        description="2x diffuse bottom chewers, same per-chewer chewing speed",
    ),
    Variant(
        name="chewer_count2x_rate2x",
        out_dir="turnover_continuous_supply_scaled_chewer_count2x_rate2x",
        chewer_count_multiplier=2,
        chewing_speed=1.6,
        max_chewing_speed=3.0,
        description="2x diffuse bottom chewers and 2x per-chewer chewing speed",
    ),
    Variant(
        name="sever_chew",
        out_dir="turnover_continuous_supply_scaled_sever_chew",
        chewer_count_multiplier=1,
        chewing_speed=0.8,
        max_chewing_speed=1.5,
        description="baseline minus-end chewers plus true severing cutters in the bottom zone",
        severer_count_multiplier=1,
        severer_cutting_rate=0.25,
        severer_cut_width=0.18,
        severer_line_diffusion=0.2,
    ),
]


BASE_FIELDS = [
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
    "n_clusters",
    "motors_per_cluster",
    "total_motors",
    "total_polymer",
    "run_seconds",
    "total_seconds",
    "notes",
]


EXTRA_FIELDS = [
    "aggressive_variant",
    "chewer_count_multiplier",
    "chewing_speed",
    "max_chewing_speed",
    "n_severers",
    "severer_count_multiplier",
    "severer_cutting_rate",
    "severer_cut_width",
    "severer_line_diffusion",
]


def replace_once(text: str, pattern: str, replacement: str) -> str:
    new, count = re.subn(pattern, replacement, text, count=1)
    if count != 1:
        raise RuntimeError(f"Expected one replacement for pattern: {pattern}")
    return new


def severer_blocks(variant: Variant) -> str:
    return f"""
set hand actin_severer
{{
    binding_rate = 10
    binding_range = 0.150
    unbinding_rate = 0.2
    activity = cut
    diffusion = {variant.severer_line_diffusion}
    cutting_rate = {variant.severer_cutting_rate}
    cut_width = {variant.severer_cut_width}
    new_end_state = white, white
    display = ( color = orange; size = 4; )
}}

set single severer1
{{
    hand = actin_severer
    stiffness = 1
    activity = diffuse
    diffusion = 10.0
    confine = inside, 0, chewer_zone
}}
""".strip()


def patch_config(text: str, size, variant: Variant) -> str:
    counts = base.scaled_counts(size)
    text = text.replace(
        "% Initial aligned network plus pulsed top actin supply and bottom minus-end chewing.",
        (
            "% Initial aligned network plus pulsed top actin supply and bottom minus-end chewing.\n"
            f"% Aggressive depolymerization variant: {variant.name} ({variant.description})."
        ),
        1,
    )
    text = replace_once(
        text,
        r"max_chewing_speed\s*=\s*1\.5",
        f"max_chewing_speed = {variant.max_chewing_speed}",
    )
    text = replace_once(
        text,
        r"chewing_speed\s*=\s*0\.8",
        f"chewing_speed = {variant.chewing_speed}",
    )
    if variant.chewer_count_multiplier != 1:
        original = counts.chewers
        updated = original * variant.chewer_count_multiplier
        text = replace_once(
            text,
            rf"new\s+{original}\s+cutter1\s*\{{\s*position\s*=\s*inside,\s*chewer_zone;\s*\}}",
            f"new {updated} cutter1 {{ position = inside, chewer_zone; }}",
        )
    if variant.severer_count_multiplier:
        severers = counts.chewers * variant.severer_count_multiplier
        text = text.replace("set hand motor1\n{", f"{severer_blocks(variant)}\n\nset hand motor1\n{{", 1)
        chewer_count = counts.chewers * variant.chewer_count_multiplier
        text = replace_once(
            text,
            rf"new\s+{chewer_count}\s+cutter1\s*\{{\s*position\s*=\s*inside,\s*chewer_zone;\s*\}}",
            (
                f"new {chewer_count} cutter1 {{ position = inside, chewer_zone; }}\n"
                f"new {severers} severer1 {{ position = inside, chewer_zone; }}"
            ),
        )
    return text


def write_submit_script(out_root: Path) -> None:
    src = BASE_DIR / "submit_cluster.sh"
    dst = out_root / "submit_cluster.sh"
    shutil.copy2(src, dst)


def write_manifest(out_root: Path, rows: list[dict]) -> None:
    with (out_root / "manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=BASE_FIELDS + EXTRA_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_readme(out_root: Path, variant: Variant, rows: list[dict]) -> None:
    lines = [
        f"# Continuous Top-Supply Scaled Campaign: {variant.name}",
        "",
        f"Aggressive depolymerization follow-up to `turnover_continuous_supply_scaled`.",
        "",
        "## Variant",
        f"- `{variant.name}`: {variant.description}.",
        f"- `chewing_speed = {variant.chewing_speed}`.",
        f"- `max_chewing_speed = {variant.max_chewing_speed}`.",
        f"- chewer count multiplier = `{variant.chewer_count_multiplier}`.",
    ]
    if variant.severer_count_multiplier:
        lines.extend(
            [
                f"- severer count multiplier = `{variant.severer_count_multiplier}`.",
                f"- `severer cutting_rate = {variant.severer_cutting_rate}`.",
                f"- `severer cut_width = {variant.severer_cut_width}`.",
                f"- `severer diffusion = {variant.severer_line_diffusion}` along bound filaments.",
            ]
        )
    lines.extend(["", "## Cases"])
    for size in base.SIZES:
        counts = base.scaled_counts(size)
        chewers = counts.chewers * variant.chewer_count_multiplier
        severers = counts.chewers * variant.severer_count_multiplier
        severer_text = f", severers = {severers}" if severers else ""
        lines.append(
            f"- `{size.name}`: initial filaments = {counts.initial_filaments}, "
            f"pulse = {counts.pulse_filaments}/min, chewers = {chewers}, "
            f"clusters = {counts.clusters}{severer_text}, polymer pool = {counts.total_polymer} um."
        )
    lines.extend(
        [
            "",
            f"Replicates: {base.REPLICATES} per size/condition.",
            f"Total configs: {len(rows)}.",
            "",
            "Submit on the cluster:",
            "```bash",
            f"cd /lustre/isaac24/scratch/kacheamp/newcytosim/{out_root.name}",
            "./submit_cluster.sh",
            "```",
        ]
    )
    (out_root / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_variant(variant: Variant) -> None:
    out_root = ROOT / variant.out_dir
    cluster_dir = out_root / "cluster_runs"
    if cluster_dir.exists():
        shutil.rmtree(cluster_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []

    for size in base.SIZES:
        counts = base.scaled_counts(size)
        for condition, _label, has_motors in base.CONDITIONS:
            for replicate in range(1, base.REPLICATES + 1):
                rel = Path("cluster_runs") / size.name / condition / f"r{replicate:04d}" / "config.cym"
                path = out_root / rel
                path.parent.mkdir(parents=True, exist_ok=True)
                text = base.render_config(size, condition, has_motors, replicate)
                path.write_text(patch_config(text, size, variant), encoding="utf-8")
                chewers = counts.chewers * variant.chewer_count_multiplier
                severers = counts.chewers * variant.severer_count_multiplier
                rows.append(
                    {
                        "kind": "cluster",
                        "size": size.name,
                        "condition": condition,
                        "replicate": replicate,
                        "path": rel.as_posix(),
                        "inner_radius": size.inner_radius,
                        "outer_radius": size.outer_radius,
                        "bottom_z": size.bottom_z,
                        "top_z": size.top_z,
                        "chewer_zone_bottom": size.chewer_bottom,
                        "chewer_zone_top": size.chewer_top,
                        "initial_filaments": counts.initial_filaments,
                        "pulse_filaments": counts.pulse_filaments,
                        "n_pulses": base.N_PULSES,
                        "initial_crosslinkers": counts.initial_crosslinkers,
                        "pulse_crosslinkers": counts.pulse_crosslinkers,
                        "n_chewers": chewers,
                        "n_clusters": counts.clusters if has_motors else 0,
                        "motors_per_cluster": base.MOTORS_PER_CLUSTER if has_motors else 0,
                        "total_motors": counts.clusters * base.MOTORS_PER_CLUSTER if has_motors else 0,
                        "total_polymer": counts.total_polymer,
                        "run_seconds": base.RUN_SECONDS,
                        "total_seconds": base.RUN_SECONDS * (base.N_PULSES + 1),
                        "notes": f"aggressive depolymerization: {variant.description}",
                        "aggressive_variant": variant.name,
                        "chewer_count_multiplier": variant.chewer_count_multiplier,
                        "chewing_speed": variant.chewing_speed,
                        "max_chewing_speed": variant.max_chewing_speed,
                        "n_severers": severers,
                        "severer_count_multiplier": variant.severer_count_multiplier,
                        "severer_cutting_rate": variant.severer_cutting_rate,
                        "severer_cut_width": variant.severer_cut_width,
                        "severer_line_diffusion": variant.severer_line_diffusion,
                    }
                )

    write_submit_script(out_root)
    write_manifest(out_root, rows)
    write_readme(out_root, variant, rows)
    print(f"Wrote {len(rows)} configs under {cluster_dir}")


def main() -> None:
    for variant in VARIANTS:
        generate_variant(variant)


if __name__ == "__main__":
    main()
