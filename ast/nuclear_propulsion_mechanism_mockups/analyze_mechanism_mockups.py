#!/usr/bin/env python3
"""Summarize nuclear displacement in mechanism mockups."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
from collections import defaultdict
from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


ROOT = Path(__file__).resolve().parent
REPORT = Path("/home/thekaka/project/cytosim/build_mini_turnover/bin/report")
SECONDS_PER_FRAME = 1.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", default="production")
    parser.add_argument("--output-dir", default="analysis")
    return parser.parse_args()


def resolve_path(path: str) -> Path:
    out = Path(path)
    return out if out.is_absolute() else ROOT / out


def discover_runs(input_root: Path) -> list[Path]:
    return sorted(
        p for p in input_root.rglob("*")
        if p.is_dir() and (p / "objects.cmo").exists() and (p / "config.cym").exists()
    )


def run_report(run_dir: Path, what: str, frame: str | None = None) -> str:
    cmd = [str(REPORT), what]
    if frame is not None:
        cmd.append(f"frame={frame}")
    return subprocess.check_output(cmd, cwd=str(run_dir), text=True)


def parse_positions_by_frame(text: str) -> dict[int, tuple[float, float, float]]:
    positions: dict[int, tuple[float, float, float]] = {}
    frame: int | None = None
    for line in text.splitlines():
        if line.startswith("% frame"):
            parts = line.split()
            if len(parts) >= 3:
                try:
                    frame = int(parts[2])
                except ValueError:
                    frame = None
            continue
        if not line or line.startswith("%") or frame is None:
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            positions[frame] = (float(parts[2]), float(parts[3]), float(parts[4]))
        except ValueError:
            continue
    return positions


def get_body_positions(run_dir: Path) -> dict[int, tuple[float, float, float]]:
    for what in ("solid:position", "bead:position"):
        try:
            positions = parse_positions_by_frame(run_report(run_dir, what))
        except subprocess.CalledProcessError:
            positions = {}
        if positions:
            return positions
    return {}


def mean_sem(vals: list[float]) -> tuple[float, float]:
    xs = [x for x in vals if not math.isnan(x)]
    if not xs:
        return math.nan, math.nan
    mean = sum(xs) / len(xs)
    if len(xs) < 2:
        return mean, math.nan
    sd = math.sqrt(sum((x - mean) ** 2 for x in xs) / (len(xs) - 1))
    return mean, sd / math.sqrt(len(xs))


def metadata(config: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in config.read_text(errors="ignore").splitlines()[:24]:
        if line.startswith("% Case:"):
            out["case"] = line.split(":", 1)[1].strip()
        elif line.startswith("% Replicate:"):
            out["replicate"] = line.split(":", 1)[1].strip()
        elif line.startswith("% Description:"):
            out["description"] = line.split(":", 1)[1].strip()
    return out


def parse_fiber_length(text: str) -> dict[str, float]:
    vals = {
        "actin_count": math.nan,
        "actin_avg_len": math.nan,
        "actin_total_len": math.nan,
        "actin_max_len": math.nan,
    }
    for line in text.splitlines():
        if line.startswith("ACTIN"):
            parts = line.split()
            vals["actin_count"] = float(parts[1])
            vals["actin_avg_len"] = float(parts[2])
            vals["actin_max_len"] = float(parts[5])
            vals["actin_total_len"] = float(parts[6])
    return vals


def parse_fiber_end_geometry(text: str) -> dict[str, float]:
    plus_x = []
    minus_x = []
    for line in text.splitlines():
        if not line or line.startswith("%"):
            continue
        parts = line.split()
        if len(parts) < 17:
            continue
        try:
            plus_x.append(float(parts[4]))
            minus_x.append(float(parts[11]))
        except ValueError:
            pass
    px_mean, _ = mean_sem(plus_x)
    mx_mean, _ = mean_sem(minus_x)
    return {
        "plus_end_x_mean": px_mean,
        "minus_end_x_mean": mx_mean,
        "plus_minus_x_delta": px_mean - mx_mean if not math.isnan(px_mean) and not math.isnan(mx_mean) else math.nan,
    }


def analyze_run(run_dir: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    positions = get_body_positions(run_dir)
    row: dict[str, object] = {"run_dir": str(run_dir)}
    row.update(metadata(run_dir / "config.cym"))
    if not positions:
        row["final_frame"] = 0
        return row, []
    final_frame = max(positions)
    p0 = positions.get(0)
    p1 = positions.get(final_frame)
    row["final_frame"] = final_frame
    row["final_time_s"] = final_frame * SECONDS_PER_FRAME
    if p0 and p1:
        dx = p1[0] - p0[0]
        dy = p1[1] - p0[1]
        dz = p1[2] - p0[2]
        row.update(
            {
                "x0": p0[0],
                "y0": p0[1],
                "z0": p0[2],
                "xf": p1[0],
                "yf": p1[1],
                "zf": p1[2],
                "dx": dx,
                "dy": dy,
                "dz": dz,
                "net_displacement": math.sqrt(dx * dx + dy * dy + dz * dz),
                "x_velocity": dx / max(final_frame * SECONDS_PER_FRAME, 1.0),
            }
        )
    row.update(parse_fiber_length(run_report(run_dir, "fiber:length", str(final_frame))))
    row.update(parse_fiber_end_geometry(run_report(run_dir, "fiber:end", str(final_frame))))

    time_rows = []
    if p0:
        for frame, pos in sorted(positions.items()):
            dx = pos[0] - p0[0]
            dy = pos[1] - p0[1]
            dz = pos[2] - p0[2]
            time_rows.append(
                {
                    "run_dir": str(run_dir),
                    "case": row.get("case", ""),
                    "replicate": row.get("replicate", ""),
                    "frame": frame,
                    "time_s": frame * SECONDS_PER_FRAME,
                    "time_min": frame * SECONDS_PER_FRAME / 60,
                    "dx": dx,
                    "dy": dy,
                    "dz": dz,
                    "net_displacement": math.sqrt(dx * dx + dy * dy + dz * dz),
                }
            )
    return row, time_rows


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[str(row.get("case", ""))].append(row)
    cols = [
        "dx",
        "dy",
        "dz",
        "net_displacement",
        "x_velocity",
        "actin_count",
        "actin_total_len",
        "actin_avg_len",
        "actin_max_len",
        "plus_minus_x_delta",
    ]
    out = []
    for case, vals in sorted(groups.items()):
        summary: dict[str, object] = {"case": case, "n": len(vals)}
        for col in cols:
            mean, sem = mean_sem([float(v.get(col, math.nan)) for v in vals])
            summary[f"{col}_mean"] = mean
            summary[f"{col}_sem"] = sem
        out.append(summary)
    return out


def summarize_timecourse(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[(str(row["case"]), int(row["frame"]))].append(row)
    out = []
    for (case, frame), vals in sorted(groups.items()):
        summary: dict[str, object] = {
            "case": case,
            "frame": frame,
            "time_s": vals[0]["time_s"],
            "time_min": vals[0]["time_min"],
            "n": len(vals),
        }
        for col in ("dx", "dy", "dz", "net_displacement"):
            mean, sem = mean_sem([float(v[col]) for v in vals])
            summary[f"{col}_mean"] = mean
            summary[f"{col}_sem"] = sem
        out.append(summary)
    return out


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_final_summary(summary: list[dict[str, object]], out_dir: Path) -> None:
    if plt is None or not summary:
        return
    labels = [str(row["case"]).replace("_", "\n", 2) for row in summary]
    xs = range(len(summary))
    dx = [float(row.get("dx_mean", math.nan)) for row in summary]
    dx_sem = [0.0 if math.isnan(float(row.get("dx_sem", math.nan))) else float(row.get("dx_sem", 0)) for row in summary]
    net = [float(row.get("net_displacement_mean", math.nan)) for row in summary]
    net_sem = [
        0.0 if math.isnan(float(row.get("net_displacement_sem", math.nan))) else float(row.get("net_displacement_sem", 0))
        for row in summary
    ]

    fig, axes = plt.subplots(2, 1, figsize=(10.8, 7.0), sharex=True)
    axes[0].bar(xs, dx, yerr=dx_sem, color="#4C78A8", alpha=0.45, capsize=4)
    axes[0].axhline(0, color="0.4", lw=1, ls="--")
    axes[0].set_ylabel("Final x displacement (um)", fontsize=12)
    axes[1].bar(xs, net, yerr=net_sem, color="#F58518", alpha=0.45, capsize=4)
    axes[1].set_ylabel("Net displacement (um)", fontsize=12)
    axes[1].set_xticks(list(xs))
    axes[1].set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color="0.9")
        ax.tick_params(axis="y", labelsize=10)
    fig.tight_layout()
    fig.savefig(out_dir / "mechanism_final_displacement.png", dpi=300)
    fig.savefig(out_dir / "mechanism_final_displacement.svg")
    plt.close(fig)


def plot_timecourse(time_summary: list[dict[str, object]], out_dir: Path) -> None:
    if plt is None or not time_summary:
        return
    groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in time_summary:
        groups[str(row["case"])].append(row)
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    for case, vals in sorted(groups.items()):
        vals = sorted(vals, key=lambda r: int(r["frame"]))
        t = [float(v["time_min"]) for v in vals]
        y = [float(v["dx_mean"]) for v in vals]
        ax.plot(t, y, lw=1.9, label=case.replace("_", " "))
    ax.axhline(0, color="0.5", lw=1, ls="--")
    ax.set_xlabel("Time (min)", fontsize=12)
    ax.set_ylabel("Nucleus displacement along x (um)", fontsize=12)
    ax.tick_params(labelsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(out_dir / "mechanism_dx_timecourse.png", dpi=300)
    fig.savefig(out_dir / "mechanism_dx_timecourse.svg")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_root = resolve_path(args.input_root)
    out_dir = resolve_path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    runs = discover_runs(input_root)
    if not runs:
        raise SystemExit(f"No completed runs found under {input_root}")

    run_rows = []
    time_rows = []
    for run in runs:
        row, trows = analyze_run(run)
        run_rows.append(row)
        time_rows.extend(trows)

    summary = summarize(run_rows)
    time_summary = summarize_timecourse(time_rows)
    write_csv(out_dir / "run_metrics.csv", run_rows)
    write_csv(out_dir / "condition_summary.csv", summary)
    write_csv(out_dir / "timecourse_run_metrics.csv", time_rows)
    write_csv(out_dir / "timecourse_condition_summary.csv", time_summary)
    plot_final_summary(summary, out_dir)
    plot_timecourse(time_summary, out_dir)
    print(out_dir / "condition_summary.csv")


if __name__ == "__main__":
    main()
