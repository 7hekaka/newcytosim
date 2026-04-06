#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import random
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SEEDS = ROOT / "clu_fixed_global" / "seed_templates"

ROT_OUT = ROOT / "clu_init_minusz"
FIX_OUT = ROOT / "clu_fixed_global_init_minusz"

BEGIN_MARK = "% --- BEGIN injected motor clusters ---"
END_MARK = "% --- END injected motor clusters ---"

INNER_RADIUS = 10.5
OUTER_RADIUS = 11.5
HEIGHT = 40.0
ATTACH_OFFSET = 0.01
R_PLACE = INNER_RADIUS + ATTACH_OFFSET
CLUSTER_RADIUS = 0.800
OVERLAP_FACTOR = 1.10
MIN_CENTER_SPACING = 0.2
MASTER_HEADS = 12
GLOBAL_SEED = 12345
EXPECTED_RUN_NAMES = [f"r{idx:04d}" for idx in range(1, 31)]

RUN_DIR_RE = re.compile(r"^r\d{4}$")
MOTOR1_BLOCK_RE = re.compile(r"set hand motor1\s*\{.*?\n\}", re.S)
ACTIN_INIT_RE = re.compile(
    r"new\s+256\s+actin1\s*\{\s*length\s*=\s*15\s*;\s*position\s*=\s*inside\s*;?\s*\}",
    re.S,
)
CLUSTER_RE = re.compile(r"^% cluster\s+(\d+): center\(theta,z\)=\(([-0-9.]+),([-0-9.]+)\), heads=(\d+)$")
POSITION_RE = re.compile(r"^new 1 myosin1 \{ position = ([^ ]+) ([^ ]+) ([^ ]+) \}$")

ROTATABLE_MOTOR1_BLOCK = """set hand motor1
{
    binding_rate = 10
    binding_range = 0.02
    unbinding_rate = 0.3
    activity = move
    unloaded_speed = 2.0
    stall_force = 4
    unbinding_force = 2.5

    display = ( color=blue; size=4; )
}"""

FIXED_MOTOR1_BLOCK = """set hand motor1
{
    binding_rate = 10
    binding_range = 0.02
    unbinding_rate = 0.3
    activity = move
    unloaded_speed = 2.0
    stall_force = 4
    unbinding_force = 2.5
    orientation_mode = fixed_global
    preferred_direction = 0, 0, -1
    angular_tolerance = 0.35
    polarity_sensitive = 1

    display = ( color=blue; size=4; )
}"""

MINUSZ_ACTIN_INIT = "new 256 actin1 { length = 15; position = inside; direction = 0 0 -1; }"


def shortest_dtheta(a: float, b: float) -> float:
    return (a - b + math.pi) % (2.0 * math.pi) - math.pi


def surface_distance(theta1: float, z1: float, theta2: float, z2: float, radius: float) -> float:
    return math.hypot(radius * abs(shortest_dtheta(theta1, theta2)), z1 - z2)


def place_cluster_centers(n_clusters: int, seed: int) -> list[tuple[float, float]]:
    min_dist = max(MIN_CENTER_SPACING, OVERLAP_FACTOR * (2.0 * CLUSTER_RADIUS))
    z_min = -0.5 * HEIGHT
    z_max = 0.5 * HEIGHT
    for restart in range(64):
        rng = random.Random(seed + restart * 1000003)
        centers: list[tuple[float, float]] = []
        ok = True
        for _ in range(n_clusters):
            placed = False
            for _ in range(50000):
                theta = rng.random() * (2.0 * math.pi)
                z = rng.uniform(z_min, z_max)
                if all(
                    surface_distance(theta, z, old_theta, old_z, INNER_RADIUS) >= min_dist
                    for old_theta, old_z in centers
                ):
                    centers.append((theta, z))
                    placed = True
                    break
            if not placed:
                ok = False
                break
        if ok:
            return centers
    raise RuntimeError(f"could not place {n_clusters} clusters with seed {seed}")


def sample_patch_positions(theta0: float, z0: float, count: int, rng: random.Random) -> list[tuple[float, float, float]]:
    positions: list[tuple[float, float, float]] = []
    z_min = -0.5 * HEIGHT
    z_max = 0.5 * HEIGHT
    for _ in range(count):
        u = rng.random()
        s = CLUSTER_RADIUS * math.sqrt(u)
        ang = 2.0 * math.pi * rng.random()
        arc = s * math.cos(ang)
        dz = s * math.sin(ang)
        theta = (theta0 + arc / INNER_RADIUS) % (2.0 * math.pi)
        z = min(max(z0 + dz, z_min), z_max)
        x = R_PLACE * math.cos(theta)
        y = R_PLACE * math.sin(theta)
        positions.append((x, y, z))
    return positions


def generate_master_clusters(run_name: str) -> list[dict[str, object]]:
    run_number = int(run_name[1:])
    seed = GLOBAL_SEED + run_number
    rng = random.Random(seed)
    centers = place_cluster_centers(80, seed)
    clusters: list[dict[str, object]] = []
    for idx, (theta, z) in enumerate(centers):
        clusters.append(
            {
                "index": idx,
                "theta": theta,
                "z": z,
                "positions": sample_patch_positions(theta, z, MASTER_HEADS, rng),
            }
        )
    return clusters


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def list_run_dirs(path: Path) -> list[str]:
    if not path.is_dir():
        raise RuntimeError(f"missing seed directory: {path}")
    return sorted(child.name for child in path.iterdir() if child.is_dir() and RUN_DIR_RE.match(child.name))


def validate_run_catalog(label: str, names: list[str]) -> None:
    missing = [name for name in EXPECTED_RUN_NAMES if name not in names]
    extra = [name for name in names if name not in EXPECTED_RUN_NAMES]
    if missing or extra:
        raise RuntimeError(f"{label} run catalog mismatch: missing={missing}, extra={extra}")


def seed_config_path(template_family: str, run_name: str) -> Path:
    path = SEEDS / template_family / run_name / "config.cym"
    if not path.is_file():
        raise RuntimeError(f"missing seed config: {path}")
    return path


def apply_motor_settings(text: str, *, fixed_global: bool) -> str:
    replacement = FIXED_MOTOR1_BLOCK if fixed_global else ROTATABLE_MOTOR1_BLOCK
    new_text, count = MOTOR1_BLOCK_RE.subn(replacement, text, count=1)
    if count != 1:
        raise RuntimeError("could not replace set hand motor1 block")
    return new_text


def apply_minusz_filament_init(text: str) -> str:
    new_text, count = ACTIN_INIT_RE.subn(MINUSZ_ACTIN_INIT, text, count=1)
    if count != 1:
        raise RuntimeError("could not replace actin initialization block")
    return new_text


def extract_cluster_region(text: str) -> tuple[str, str, str]:
    start = text.index(BEGIN_MARK)
    end = text.index(END_MARK) + len(END_MARK)
    return text[:start], text[start:end], text[end:]


def parse_clusters(block: str) -> list[dict[str, object]]:
    clusters: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    for raw in block.splitlines():
        line = raw.strip()
        match = CLUSTER_RE.match(line)
        if match:
            if current is not None:
                clusters.append(current)
            current = {
                "index": int(match.group(1)),
                "theta": float(match.group(2)),
                "z": float(match.group(3)),
                "positions": [],
            }
            continue
        match = POSITION_RE.match(line)
        if match:
            if current is None:
                raise RuntimeError("position line encountered before cluster header")
            positions = current["positions"]
            assert isinstance(positions, list)
            positions.append((float(match.group(1)), float(match.group(2)), float(match.group(3))))
    if current is not None:
        clusters.append(current)
    return clusters


def build_cluster_block(clusters: list[dict[str, object]], motors_per_cluster: int) -> str:
    heads = [motors_per_cluster] * len(clusters)
    total = len(clusters) * motors_per_cluster
    lines = [
        BEGIN_MARK,
        f"% inner_radius={INNER_RADIUS:.1f} outer_radius={OUTER_RADIUS:.1f} height={HEIGHT:.1f} attach_offset={ATTACH_OFFSET:.2f}",
        "% stationary surface clusters on INNER wall",
        f"% n_clusters={len(clusters)} (requested {len(clusters)}), cluster_radius={CLUSTER_RADIUS:.3f} um",
        f"% M_total={total} -> heads_per_cluster={heads}",
        f"% r_place={R_PLACE:.4f} um, overlap_factor~{OVERLAP_FACTOR:.2f}, min_center_spacing={MIN_CENTER_SPACING:.1f} um",
    ]
    for idx, cluster in enumerate(clusters):
        theta = float(cluster["theta"])
        z = float(cluster["z"])
        positions = cluster["positions"]
        assert isinstance(positions, list)
        if len(positions) < motors_per_cluster:
            raise RuntimeError(f"cluster {idx} has only {len(positions)} positions")
        lines.append(f"% cluster {idx:02d}: center(theta,z)=({theta:.3f},{z:.3f}), heads={motors_per_cluster}")
        for x, y, zz in positions[:motors_per_cluster]:
            lines.append(f"new 1 myosin1 {{ position = {x:.3f} {y:.3f} {zz:.3f} }}")
    lines.append(END_MARK)
    return "\n".join(lines)


def replace_cluster_block(text: str, cluster_block: str) -> str:
    prefix, _, suffix = extract_cluster_region(text)
    prefix = prefix.rstrip() + "\n\n"
    suffix = "\n\n" + suffix.lstrip("\n")
    return prefix + cluster_block + suffix


def trim_old_clusters(config_path: Path, cluster_count: int, motors_per_cluster: int) -> str:
    _, block, _ = extract_cluster_region(read_text(config_path))
    clusters = parse_clusters(block)
    if len(clusters) < cluster_count:
        raise RuntimeError(f"{config_path} has only {len(clusters)} clusters")
    trimmed: list[dict[str, object]] = []
    for cluster in clusters[:cluster_count]:
        positions = cluster["positions"]
        assert isinstance(positions, list)
        trimmed.append(
            {
                "index": cluster["index"],
                "theta": cluster["theta"],
                "z": cluster["z"],
                "positions": positions[:motors_per_cluster],
            }
        )
    return build_cluster_block(trimmed, motors_per_cluster)


def generated_cluster_block(run_name: str, cluster_count: int, motors_per_cluster: int) -> str:
    master = generate_master_clusters(run_name)
    trimmed: list[dict[str, object]] = []
    for cluster in master[:cluster_count]:
        positions = cluster["positions"]
        assert isinstance(positions, list)
        trimmed.append(
            {
                "index": cluster["index"],
                "theta": cluster["theta"],
                "z": cluster["z"],
                "positions": positions[:motors_per_cluster],
            }
        )
    return build_cluster_block(trimmed, motors_per_cluster)


def run_names() -> list[str]:
    names = list_run_dirs(SEEDS / "40")
    validate_run_catalog("seed_templates/40", names)
    for template_family in ("10", "20"):
        other = list_run_dirs(SEEDS / template_family)
        validate_run_catalog(f"seed_templates/{template_family}", other)
        if other != names:
            raise RuntimeError(
                f"seed run catalogs differ between family 40 and {template_family}: "
                f"40={names}, {template_family}={other}"
            )
    return names


def build_cases() -> list[dict[str, object]]:
    return [
        {
            "group": "total480",
            "case": "c10_m48",
            "template_dir": "10",
            "cluster_mode": "keep",
            "cluster_source": "10",
            "cluster_count": 10,
            "motors_per_cluster": 48,
            "source_note": "seed_templates/10 with minus-z filament initialization",
        },
        {
            "group": "total480",
            "case": "c20_m24",
            "template_dir": "20",
            "cluster_mode": "keep",
            "cluster_source": "20",
            "cluster_count": 20,
            "motors_per_cluster": 24,
            "source_note": "seed_templates/20 with minus-z filament initialization",
        },
        {
            "group": "total480",
            "case": "c40_m12",
            "template_dir": "40",
            "cluster_mode": "keep",
            "cluster_source": "40",
            "cluster_count": 40,
            "motors_per_cluster": 12,
            "source_note": "seed_templates/40 with minus-z filament initialization",
        },
        {
            "group": "total480",
            "case": "c60_m8",
            "template_dir": "40",
            "cluster_mode": "generated",
            "cluster_source": "generated80",
            "cluster_count": 60,
            "motors_per_cluster": 8,
            "source_note": "seed_templates/40 with procedurally generated 60-cluster block and minus-z filament initialization",
        },
        {
            "group": "total480",
            "case": "c80_m6",
            "template_dir": "40",
            "cluster_mode": "generated",
            "cluster_source": "generated80",
            "cluster_count": 80,
            "motors_per_cluster": 6,
            "source_note": "seed_templates/40 with procedurally generated 80-cluster block and minus-z filament initialization",
        },
        {
            "group": "mpc12",
            "case": "c10_m12",
            "template_dir": "10",
            "cluster_mode": "trim_old",
            "cluster_source": "10",
            "cluster_count": 10,
            "motors_per_cluster": 12,
            "source_note": "seed_templates/10 trimmed to 12 heads per cluster with minus-z filament initialization",
        },
        {
            "group": "mpc12",
            "case": "c20_m12",
            "template_dir": "20",
            "cluster_mode": "trim_old",
            "cluster_source": "20",
            "cluster_count": 20,
            "motors_per_cluster": 12,
            "source_note": "seed_templates/20 trimmed to 12 heads per cluster with minus-z filament initialization",
        },
        {
            "group": "mpc12",
            "case": "c40_m12",
            "template_dir": "40",
            "cluster_mode": "keep",
            "cluster_source": "40",
            "cluster_count": 40,
            "motors_per_cluster": 12,
            "source_note": "seed_templates/40 with minus-z filament initialization",
        },
        {
            "group": "mpc12",
            "case": "c60_m12",
            "template_dir": "40",
            "cluster_mode": "generated",
            "cluster_source": "generated80",
            "cluster_count": 60,
            "motors_per_cluster": 12,
            "source_note": "seed_templates/40 with procedurally generated 60-cluster block and minus-z filament initialization",
        },
        {
            "group": "mpc12",
            "case": "c80_m12",
            "template_dir": "40",
            "cluster_mode": "generated",
            "cluster_source": "generated80",
            "cluster_count": 80,
            "motors_per_cluster": 12,
            "source_note": "seed_templates/40 with procedurally generated 80-cluster block and minus-z filament initialization",
        },
    ]


def render_case(run_name: str, case: dict[str, object], *, fixed_global: bool) -> str:
    template_dir = str(case["template_dir"])
    cluster_mode = str(case["cluster_mode"])
    cluster_count = int(case["cluster_count"])
    motors_per_cluster = int(case["motors_per_cluster"])

    template_path = seed_config_path(template_dir, run_name)
    text = read_text(template_path)
    text = apply_motor_settings(text, fixed_global=fixed_global)
    text = apply_minusz_filament_init(text)

    if cluster_mode == "trim_old":
        cluster_source = str(case["cluster_source"])
        cluster_block = trim_old_clusters(seed_config_path(cluster_source, run_name), cluster_count, motors_per_cluster)
        text = replace_cluster_block(text, cluster_block)
    elif cluster_mode == "generated":
        cluster_block = generated_cluster_block(run_name, cluster_count, motors_per_cluster)
        text = replace_cluster_block(text, cluster_block)
    elif cluster_mode == "keep":
        pass
    else:
        raise ValueError(f"unknown cluster mode: {cluster_mode}")

    return text


def regime_from_run(run_name: str) -> str:
    idx = int(run_name[1:])
    if idx <= 10:
        return "1:4"
    if idx <= 20:
        return "1:8"
    return "1:16"


def write_manifest(out_root: Path, rows: list[dict[str, object]]) -> None:
    manifest = out_root / "manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "model",
                "group",
                "case",
                "run",
                "cluster_count",
                "motors_per_cluster",
                "total_motors",
                "template_dir",
                "cluster_source",
                "initial_filament_direction",
                "notes",
                "path",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def write_xlink_map(out_root: Path, model_name: str, rows: list[dict[str, object]]) -> None:
    path = out_root / "xlink_regime_map.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["model", "group", "case", "xlink_regime", "run_dir", "run_path", "has_outputs"],
        )
        writer.writeheader()
        for row in rows:
            run_name = str(row["run"])
            group = str(row["group"])
            case = str(row["case"])
            run_dir = out_root / group / case / run_name
            writer.writerow(
                {
                    "model": model_name,
                    "group": group,
                    "case": case,
                    "xlink_regime": regime_from_run(run_name),
                    "run_dir": run_name,
                    "run_path": str(run_dir),
                    "has_outputs": int((run_dir / "log.txt").exists() and (run_dir / "objects.cmo").exists()),
                }
            )


def generate_campaign(out_root: Path, *, fixed_global: bool) -> list[dict[str, object]]:
    names = run_names()
    cases = build_cases()
    out_root.mkdir(exist_ok=True)
    rows: list[dict[str, object]] = []
    model_name = "fixed_global" if fixed_global else "rotatable"

    for case in cases:
        group = str(case["group"])
        name = str(case["case"])
        for run_name in names:
            config_text = render_case(run_name, case, fixed_global=fixed_global)
            out_path = out_root / group / name / run_name / "config.cym"
            write_text(out_path, config_text)
            rows.append(
                {
                    "model": model_name,
                    "group": group,
                    "case": name,
                    "run": run_name,
                    "cluster_count": int(case["cluster_count"]),
                    "motors_per_cluster": int(case["motors_per_cluster"]),
                    "total_motors": int(case["cluster_count"]) * int(case["motors_per_cluster"]),
                    "template_dir": case["template_dir"],
                    "cluster_source": case["cluster_source"],
                    "initial_filament_direction": "0 0 -1",
                    "notes": case["source_note"],
                    "path": str(out_path.relative_to(ROOT)),
                }
            )

    write_manifest(out_root, rows)
    write_xlink_map(out_root, model_name, rows)
    return rows


def main() -> None:
    rot_rows = generate_campaign(ROT_OUT, fixed_global=False)
    fix_rows = generate_campaign(FIX_OUT, fixed_global=True)
    print(f"Wrote {len(rot_rows)} configs under {ROT_OUT}")
    print(f"Wrote {len(fix_rows)} configs under {FIX_OUT}")
    print(f"Manifests: {ROT_OUT / 'manifest.csv'} and {FIX_OUT / 'manifest.csv'}")
    print(f"Maps: {ROT_OUT / 'xlink_regime_map.csv'} and {FIX_OUT / 'xlink_regime_map.csv'}")


if __name__ == "__main__":
    main()
