#!/usr/bin/env python3
"""Analyze returned replicated tethered-MT/free-actin runs copied into a save folder."""

from __future__ import annotations

import csv
import argparse
import math
import re
import subprocess
from collections import defaultdict
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = Path(__file__).resolve().parent
DEFAULT_SAVE = ROOT / "ast" / "tethered_mt_free_actin_pilot" / "save"
OUT = CAMPAIGN / "analysis" / "returned_save_first_pass"
REPORT = Path("/home/thekaka/project/cytosim/build_mini_turnover/bin/report")
FINAL_FRAME = 149


def run_report(run_dir: Path, what: str, frame: str) -> str:
    return subprocess.check_output(
        [str(REPORT), what, f"frame={frame}"],
        cwd=str(run_dir),
        text=True,
        stderr=subprocess.STDOUT,
    )


def parse_config_metadata(config: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    for line in config.read_text(errors="ignore").splitlines()[:40]:
        if line.startswith("% System:"):
            raw = line.split(":", 1)[1].strip()
            metadata["system_label"] = raw
            metadata["system"] = raw.split()[0]
        elif line.startswith("% Case:"):
            metadata["case"] = line.split(":", 1)[1].strip()
        elif line.startswith("% Replicate:"):
            metadata["replicate"] = line.split(":", 1)[1].strip()
        elif line.startswith("% Description:"):
            metadata["description"] = line.split(":", 1)[1].strip()
        elif "Actin total_polymer" in line:
            match = re.search(r"Actin total_polymer = ([0-9.]+)", line)
            if match:
                metadata["actin_total_polymer"] = match.group(1)
        elif "Fixed MT length" in line:
            match = re.search(r"Fixed MT length = ([0-9.]+)", line)
            if match:
                metadata["fixed_mt_length"] = match.group(1)
    return metadata


def parse_fiber_length(text: str) -> dict[int, dict[str, dict[str, float]]]:
    by_frame: dict[int, dict[str, dict[str, float]]] = defaultdict(dict)
    frame = None
    for line in text.splitlines():
        if line.startswith("% frame"):
            frame = int(line.split()[2])
            continue
        if not line or line.startswith("%") or frame is None:
            continue
        parts = line.split()
        if len(parts) < 8:
            continue
        cls = parts[0]
        by_frame[frame][cls] = {
            "count": float(parts[1]),
            "avg_len": float(parts[2]),
            "var_len": float(parts[3]),
            "min_len": float(parts[4]),
            "max_len": float(parts[5]),
            "total_len": float(parts[6]),
            "off_len": float(parts[7]),
        }
    return by_frame


def parse_positions(text: str) -> list[tuple[float, float, float]]:
    positions: list[tuple[float, float, float]] = []
    for line in text.splitlines():
        if not line or line.startswith("%"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            # solid: class identity cenX cenY cenZ ...
            # bead:  class identity posX posY posZ
            positions.append((float(parts[2]), float(parts[3]), float(parts[4])))
        except ValueError:
            continue
    return positions


def pairwise_distances(points: list[tuple[float, float, float]]) -> list[float]:
    distances: list[float] = []
    for a, b in combinations(points, 2):
        distances.append(math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3))))
    return distances


def centroid(points: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    n = len(points)
    return tuple(sum(p[i] for p in points) / n for i in range(3))


def position_metrics(points: list[tuple[float, float, float]]) -> dict[str, float]:
    if not points:
        return {
            "body_count": 0,
            "mean_pairwise_dist": math.nan,
            "min_pairwise_dist": math.nan,
            "max_pairwise_dist": math.nan,
            "centroid_radius_xy": math.nan,
            "mean_radius_from_centroid": math.nan,
        }
    dists = pairwise_distances(points)
    cen = centroid(points)
    radii = [math.sqrt(sum((p[i] - cen[i]) ** 2 for i in range(3))) for p in points]
    return {
        "body_count": len(points),
        "mean_pairwise_dist": sum(dists) / len(dists) if dists else math.nan,
        "min_pairwise_dist": min(dists) if dists else math.nan,
        "max_pairwise_dist": max(dists) if dists else math.nan,
        "centroid_radius_xy": math.sqrt(cen[0] ** 2 + cen[1] ** 2),
        "mean_radius_from_centroid": sum(radii) / len(radii),
    }


def get_body_positions(run_dir: Path, frame: int) -> list[tuple[float, float, float]]:
    positions = parse_positions(run_report(run_dir, "solid:position", str(frame)))
    if positions:
        return positions
    return parse_positions(run_report(run_dir, "bead:position", str(frame)))


def analyze_run(run_dir: Path) -> dict[str, object]:
    meta = parse_config_metadata(run_dir / "config.cym")
    row: dict[str, object] = {"run_id": run_dir.name, "run_dir": str(run_dir)}
    row.update(meta)

    length_report = run_report(run_dir, "fiber:length", f"0,{FINAL_FRAME}")
    lengths = parse_fiber_length(length_report)
    for frame_label, frame in [("early", 0), ("final", FINAL_FRAME)]:
        for cls in ["ACTIN", "MT"]:
            vals = lengths.get(frame, {}).get(cls, {})
            prefix = f"{frame_label}_{cls.lower()}"
            for key in ["count", "avg_len", "var_len", "min_len", "max_len", "total_len", "off_len"]:
                row[f"{prefix}_{key}"] = vals.get(key, math.nan)

    early_actin = float(row.get("early_actin_total_len", math.nan))
    final_actin = float(row.get("final_actin_total_len", math.nan))
    row["actin_total_len_change"] = final_actin - early_actin
    row["actin_total_len_fold_change"] = final_actin / early_actin if early_actin else math.nan

    early_pos = get_body_positions(run_dir, 0)
    final_pos = get_body_positions(run_dir, FINAL_FRAME)
    for label, points in [("early", early_pos), ("final", final_pos)]:
        for key, val in position_metrics(points).items():
            row[f"{label}_{key}"] = val
    if early_pos and final_pos and len(early_pos) == len(final_pos):
        displacements = [
            math.sqrt(sum((final_pos[i][j] - early_pos[i][j]) ** 2 for j in range(3)))
            for i in range(len(early_pos))
        ]
        row["mean_body_displacement"] = sum(displacements) / len(displacements)
        row["max_body_displacement"] = max(displacements)
    else:
        row["mean_body_displacement"] = math.nan
        row["max_body_displacement"] = math.nan

    return row


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[(str(row.get("system", "")), str(row.get("case", "")))].append(row)

    numeric_cols = [
        "final_actin_count",
        "final_actin_avg_len",
        "final_actin_max_len",
        "final_actin_total_len",
        "final_actin_off_len",
        "actin_total_len_change",
        "actin_total_len_fold_change",
        "final_mt_count",
        "final_mt_avg_len",
        "final_mt_max_len",
        "final_mt_total_len",
        "final_body_count",
        "final_mean_pairwise_dist",
        "final_min_pairwise_dist",
        "final_max_pairwise_dist",
        "final_mean_radius_from_centroid",
        "mean_body_displacement",
        "max_body_displacement",
    ]
    out: list[dict[str, object]] = []
    for (system, case), vals in sorted(groups.items()):
        summary: dict[str, object] = {
            "system": system,
            "case": case,
            "n": len(vals),
            "run_ids": ";".join(str(v["run_id"]) for v in vals),
            "replicates": ";".join(str(v.get("replicate", "")) for v in vals),
        }
        for col in numeric_cols:
            xs = [float(v[col]) for v in vals if col in v and not math.isnan(float(v[col]))]
            if xs:
                mean = sum(xs) / len(xs)
                if len(xs) > 1:
                    sd = math.sqrt(sum((x - mean) ** 2 for x in xs) / (len(xs) - 1))
                    sem = sd / math.sqrt(len(xs))
                else:
                    sd = math.nan
                    sem = math.nan
                summary[f"{col}_mean"] = mean
                summary[f"{col}_sd"] = sd
                summary[f"{col}_sem"] = sem
            else:
                summary[f"{col}_mean"] = math.nan
                summary[f"{col}_sd"] = math.nan
                summary[f"{col}_sem"] = math.nan
        out.append(summary)
    return out


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def maybe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def write_representative_movies(rows: list[dict[str, object]], condition_rows: list[dict[str, object]]) -> None:
    by_condition: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_condition[(str(row.get("system")), str(row.get("case")))].append(row)

    metric = "final_mean_pairwise_dist"
    lines = [
        "Representative runs chosen as the replicate closest to the condition mean final nuclear/body pairwise distance.",
        "Use these first for visual inspection; then inspect outliers if the condition looks heterogeneous.",
        "",
    ]
    for summary in condition_rows:
        key = (str(summary["system"]), str(summary["case"]))
        vals = by_condition[key]
        target = maybe_float(summary.get(f"{metric}_mean"))
        valid = [v for v in vals if not math.isnan(maybe_float(v.get(metric)))]
        if not valid:
            continue
        best = min(valid, key=lambda v: abs(maybe_float(v.get(metric)) - target))
        lines.append(
            f"{key[0]},{key[1]}: {best['run_id']} ({best.get('replicate')}) "
            f"path={best['run_dir']} final_pairwise={maybe_float(best.get(metric)):.3f}"
        )
    (OUT / "representative_movies.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_plots(condition_rows: list[dict[str, object]]) -> None:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as exc:  # pragma: no cover - optional plotting dependency
        (OUT / "plot_error.txt").write_text(f"matplotlib unavailable: {exc}\n", encoding="utf-8")
        return

    plots = OUT / "plots"
    plots.mkdir(exist_ok=True)

    systems = ["small_system", "large_system"]
    cases = [
        "01_actin_only_perinuclear_long",
        "02_fixed_mt_uniform_actin_long",
        "03_growing_mt_uniform_actin_long",
        "04_fixed_mt_perinuclear_actin_long",
        "05_growing_mt_perinuclear_actin_long",
        "06_division_fixed_mt_perinuclear_actin_long",
        "07_division_growing_mt_perinuclear_actin_long",
    ]
    case_labels = ["Actin only", "Fixed MT\nuniform", "Growing MT\nuniform", "Fixed MT\nperi", "Growing MT\nperi", "Div fixed\nperi", "Div growing\nperi"]
    complete = {(r["system"], r["case"]): int(r["n"]) for r in condition_rows}
    matrix = np.array([[complete.get((system, case), 0) for case in cases] for system in systems], dtype=float)

    fig, ax = plt.subplots(figsize=(10, 2.7))
    im = ax.imshow(matrix, vmin=0, vmax=5, cmap="Greys")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, str(int(matrix[i, j])), ha="center", va="center", fontsize=11)
    ax.set_yticks(range(len(systems)), ["Small", "Large"], fontsize=11)
    ax.set_xticks(range(len(cases)), case_labels, fontsize=9, rotation=35, ha="right")
    ax.set_ylabel("System", fontsize=12)
    ax.set_title("Returned replicate count", fontsize=13)
    fig.colorbar(im, ax=ax, label="completed replicates")
    fig.tight_layout()
    fig.savefig(plots / "completion_matrix.png", dpi=300)
    fig.savefig(plots / "completion_matrix.pdf")
    plt.close(fig)

    metrics = [
        ("final_actin_total_len_mean", "Final actin total length (um)"),
        ("final_mt_total_len_mean", "Final MT total length (um)"),
        ("final_mean_pairwise_dist_mean", "Final mean body pairwise distance (um)"),
        ("mean_body_displacement_mean", "Mean body displacement (um)"),
    ]
    for system in systems:
        rows = [r for r in condition_rows if r["system"] == system]
        if not rows:
            continue
        fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
        axes = axes.ravel()
        x = np.arange(len(cases))
        row_by_case = {r["case"]: r for r in rows}
        for ax, (col, ylabel) in zip(axes, metrics):
            means = [maybe_float(row_by_case.get(case, {}).get(col)) for case in cases]
            sems = [maybe_float(row_by_case.get(case, {}).get(col.replace("_mean", "_sem"))) for case in cases]
            ax.bar(x, means, yerr=[0 if math.isnan(s) else s for s in sems], color="#6f8faf", edgecolor="black", linewidth=0.6, capsize=3)
            ax.set_ylabel(ylabel, fontsize=11)
            ax.grid(axis="y", alpha=0.18)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        for ax in axes:
            ax.set_xticks(x, case_labels, rotation=35, ha="right", fontsize=8)
        fig.suptitle(system.replace("_", " "), fontsize=14)
        fig.tight_layout()
        fig.savefig(plots / f"{system}_condition_summary.png", dpi=300)
        fig.savefig(plots / f"{system}_condition_summary.pdf")
        plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--save-dir",
        type=Path,
        default=DEFAULT_SAVE,
        help=f"Returned save directory to analyze. Default: {DEFAULT_SAVE}",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=OUT,
        help=f"Output directory. Default: {OUT}",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    global OUT
    OUT = args.out_dir
    OUT.mkdir(parents=True, exist_ok=True)
    save_dir = args.save_dir
    runs = sorted(p for p in save_dir.glob("r*") if (p / "config.cym").exists() and (p / "objects.cmo").exists())
    rows = [analyze_run(p) for p in runs]
    condition_rows = summarize(rows)

    write_csv(OUT / "returned_run_metrics.csv", rows)
    write_csv(OUT / "returned_condition_summary.csv", condition_rows)
    write_representative_movies(rows, condition_rows)
    make_plots(condition_rows)

    with (OUT / "completion_summary.txt").open("w", encoding="utf-8") as handle:
        handle.write(f"save_dir: {save_dir}\n")
        handle.write(f"completed_runs: {len(rows)}\n\n")
        for row in condition_rows:
            handle.write(
                f"{row['system']},{row['case']}: n={row['n']} runs={row['run_ids']} reps={row['replicates']}\n"
            )
    print(OUT)


if __name__ == "__main__":
    main()
