#!/usr/bin/env python3
"""Analyze long treadmilling-only nuclear-aster pilot outputs."""

from __future__ import annotations

import csv
import argparse
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
FINAL_FRAME = 300
N_FRAMES = 301
SECONDS_PER_FRAME = 1.0


def run_report(run_dir: Path, what: str, frame: str | None = None) -> str:
    cmd = [str(REPORT), what]
    if frame is not None:
        cmd.append(f"frame={frame}")
    return subprocess.check_output(cmd, cwd=str(run_dir), text=True)


def parse_position(text: str) -> tuple[float, float, float] | None:
    for line in text.splitlines():
        if not line or line.startswith("%"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            return (float(parts[2]), float(parts[3]), float(parts[4]))
        except ValueError:
            continue
    return None


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


def get_body_position(run_dir: Path, frame: int) -> tuple[float, float, float] | None:
    for what in ("solid:position", "bead:position"):
        try:
            pos = parse_position(run_report(run_dir, what, str(frame)))
        except subprocess.CalledProcessError:
            pos = None
        if pos is not None:
            return pos
    return None


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


def parse_fiber_length(text: str) -> dict[str, float]:
    vals = {
        "actin_count": math.nan,
        "actin_total_len": math.nan,
        "actin_avg_len": math.nan,
        "actin_min_len": math.nan,
        "actin_max_len": math.nan,
    }
    for line in text.splitlines():
        if line.startswith("ACTIN"):
            parts = line.split()
            vals["actin_count"] = float(parts[1])
            vals["actin_avg_len"] = float(parts[2])
            vals["actin_min_len"] = float(parts[4])
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


def metadata(config: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in config.read_text(errors="ignore").splitlines()[:36]:
        if line.startswith("% Case:"):
            out["case"] = line.split(":", 1)[1].strip()
        elif line.startswith("% Replicate:"):
            out["replicate"] = line.split(":", 1)[1].strip()
        elif line.startswith("% Description:"):
            out["description"] = line.split(":", 1)[1].strip()
    return out


def available_final_frame(run_dir: Path) -> int:
    for frame in range(FINAL_FRAME, -1, -1):
        if get_body_position(run_dir, frame) is not None:
            return frame
    return 0


def analyze_run(run_dir: Path, positions: dict[int, tuple[float, float, float]] | None = None) -> dict[str, object]:
    row: dict[str, object] = {"run_dir": str(run_dir)}
    row.update(metadata(run_dir / "config.cym"))
    if positions is None:
        positions = get_body_positions(run_dir)
    final_frame = max(positions) if positions else available_final_frame(run_dir)
    row["final_frame"] = final_frame
    row["final_time_s"] = final_frame * SECONDS_PER_FRAME
    p0 = positions.get(0) if positions else get_body_position(run_dir, 0)
    p1 = positions.get(final_frame) if positions else get_body_position(run_dir, final_frame)
    if p0 and p1:
        dx = p1[0] - p0[0]
        dy = p1[1] - p0[1]
        dz = p1[2] - p0[2]
        elapsed = max(final_frame * SECONDS_PER_FRAME, 1.0)
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
                "x_velocity": dx / elapsed,
            }
        )
    row.update(parse_fiber_length(run_report(run_dir, "fiber:length", str(final_frame))))
    row.update(parse_fiber_end_geometry(run_report(run_dir, "fiber:end", str(final_frame))))
    return row


def analyze_timecourse(run_dir: Path, positions: dict[int, tuple[float, float, float]] | None = None) -> list[dict[str, object]]:
    meta = metadata(run_dir / "config.cym")
    if positions is None:
        positions = get_body_positions(run_dir)
    if not positions:
        return []
    final_frame = max(positions)
    p0 = positions.get(0)
    if p0 is None:
        return []
    rows = []
    for frame in range(final_frame + 1):
        pos = positions.get(frame)
        if pos is None:
            continue
        dx = pos[0] - p0[0]
        dy = pos[1] - p0[1]
        dz = pos[2] - p0[2]
        rows.append(
            {
                "run_dir": str(run_dir),
                "case": meta.get("case", ""),
                "replicate": meta.get("replicate", ""),
                "frame": frame,
                "time_s": frame * SECONDS_PER_FRAME,
                "time_min": frame * SECONDS_PER_FRAME / 60,
                "dx": dx,
                "dy": dy,
                "dz": dz,
                "net_displacement": math.sqrt(dx * dx + dy * dy + dz * dz),
            }
        )
    return rows


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[str(row["case"])].append(row)
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
        "plus_end_x_mean",
        "minus_end_x_mean",
        "plus_minus_x_delta",
    ]
    out = []
    for case, vals in sorted(groups.items()):
        summary: dict[str, object] = {"case": case, "n": len(vals)}
        summary["completed_5min"] = sum(int(float(v.get("final_frame", 0))) >= FINAL_FRAME for v in vals)
        for col in cols:
            mean, sem = mean_sem([float(v[col]) for v in vals if col in v])
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


def plot_timecourse(rows: list[dict[str, object]], out_dir: Path) -> None:
    if plt is None or not rows:
        return
    colors = {
        "01_symmetric_dense_treadmill": "#666666",
        "02_minus_x_cap_dense_treadmill": "#0072B2",
        "03_plus_x_cap_dense_treadmill": "#D55E00",
        "04_no_nucleator_control": "#CC79A7",
    }
    labels = {
        "01_symmetric_dense_treadmill": "symmetric dense aster",
        "02_minus_x_cap_dense_treadmill": "minus-x cap aster",
        "03_plus_x_cap_dense_treadmill": "plus-x cap aster",
        "04_no_nucleator_control": "no nucleator control",
    }
    by_case: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_case[str(row["case"])].append(row)

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    for case, vals in sorted(by_case.items()):
        vals = sorted(vals, key=lambda r: int(r["frame"]))
        t = [float(v["time_min"]) for v in vals]
        y = [float(v["dx_mean"]) for v in vals]
        sem = [0.0 if math.isnan(float(v["dx_sem"])) else float(v["dx_sem"]) for v in vals]
        ax.plot(t, y, lw=2.2, color=colors.get(case, "black"), label=labels.get(case, case))
        ax.fill_between(t, [a - b for a, b in zip(y, sem)], [a + b for a, b in zip(y, sem)],
                        color=colors.get(case, "black"), alpha=0.16, linewidth=0)
    ax.axhline(0, color="0.55", lw=1.0, ls="--")
    ax.set_xlabel("Time (min)", fontsize=13)
    ax.set_ylabel("Nucleus displacement along x, dx (um)", fontsize=13)
    ax.tick_params(labelsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="best")
    fig.tight_layout()
    fig.savefig(out_dir / "dx_timecourse.png", dpi=300)
    fig.savefig(out_dir / "dx_timecourse.svg")
    plt.close(fig)


def plot_final_summary(
    run_rows: list[dict[str, object]],
    summary_rows: list[dict[str, object]],
    out_dir: Path,
) -> None:
    if plt is None or not run_rows or not summary_rows:
        return
    labels = {
        "01_symmetric_dense_treadmill": "symmetric\naster",
        "02_minus_x_cap_dense_treadmill": "minus-x\ncap",
        "03_plus_x_cap_dense_treadmill": "plus-x\ncap",
        "04_no_nucleator_control": "no\nnucleator",
    }
    colors = {
        "01_symmetric_dense_treadmill": "#666666",
        "02_minus_x_cap_dense_treadmill": "#0072B2",
        "03_plus_x_cap_dense_treadmill": "#D55E00",
        "04_no_nucleator_control": "#CC79A7",
    }
    metrics = [
        ("dx", "Final x displacement (um)", "A"),
        ("net_displacement", "Net displacement (um)", "B"),
        ("actin_avg_len", "Final actin length (um)", "C"),
        ("plus_minus_x_delta", "Plus-minus end x-separation (um)", "D"),
    ]
    cases = [str(row["case"]) for row in summary_rows]
    summary_by_case = {str(row["case"]): row for row in summary_rows}
    by_case: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in run_rows:
        by_case[str(row.get("case", ""))].append(row)

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.8))
    jitter_offsets = [-0.09, 0.0, 0.09, -0.05, 0.05]
    for ax, (metric, ylabel, panel) in zip(axes.flat, metrics):
        means = [float(summary_by_case[c].get(f"{metric}_mean", math.nan)) for c in cases]
        sems = [
            0.0 if math.isnan(float(summary_by_case[c].get(f"{metric}_sem", math.nan)))
            else float(summary_by_case[c].get(f"{metric}_sem", math.nan))
            for c in cases
        ]
        xs = list(range(len(cases)))
        ax.bar(
            xs,
            means,
            yerr=sems,
            color=[colors.get(c, "0.7") for c in cases],
            alpha=0.24,
            edgecolor=[colors.get(c, "0.3") for c in cases],
            linewidth=1.2,
            capsize=4,
        )
        for xpos, case in zip(xs, cases):
            vals = []
            for row in by_case.get(case, []):
                try:
                    val = float(row.get(metric, math.nan))
                except (TypeError, ValueError):
                    val = math.nan
                if not math.isnan(val):
                    vals.append(val)
            for idx, val in enumerate(vals):
                ax.plot(
                    xpos + jitter_offsets[idx % len(jitter_offsets)],
                    val,
                    "o",
                    ms=4.8,
                    color=colors.get(case, "0.2"),
                    markeredgecolor="white",
                    markeredgewidth=0.5,
                )
        ax.axhline(0, color="0.55", lw=1.0, ls="--")
        ax.set_xticks(xs)
        ax.set_xticklabels([labels.get(c, c) for c in cases], fontsize=9)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.tick_params(axis="y", labelsize=9)
        ax.text(-0.14, 1.04, panel, transform=ax.transAxes, fontsize=13, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color="0.9", linewidth=0.8)
    fig.tight_layout()
    fig.savefig(out_dir / "final_summary.png", dpi=300)
    fig.savefig(out_dir / "final_summary.svg")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-root",
        default="production",
        help="Run folder root relative to this pilot folder, e.g. production or job03/save.",
    )
    parser.add_argument(
        "--output-dir",
        default="analysis",
        help="Output folder relative to this pilot folder.",
    )
    return parser.parse_args()


def resolve_path(path: str) -> Path:
    out = Path(path)
    if not out.is_absolute():
        out = ROOT / out
    return out


def discover_runs(input_root: Path) -> list[Path]:
    return sorted(
        p for p in input_root.rglob("*")
        if p.is_dir() and (p / "objects.cmo").exists() and (p / "config.cym").exists()
    )


def write_representative_runs(path: Path, rows: list[dict[str, object]]) -> None:
    groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[str(row.get("case", ""))].append(row)

    lines = [
        "Representative run suggestions",
        "",
        "Chosen as the replicate closest to the condition mean dx, so these are good first movies to inspect.",
        "",
        "Run mapping:",
    ]
    for case, vals in sorted(groups.items()):
        vals = sorted(vals, key=lambda r: str(r.get("replicate", "")))
        ids = ", ".join(Path(str(v["run_dir"])).name for v in vals)
        lines.append(f"- {case}: {ids}")

    lines.extend(["", "Suggested movies:"])
    for case, vals in sorted(groups.items()):
        dxs = [float(v.get("dx", math.nan)) for v in vals]
        mean, _ = mean_sem(dxs)
        if math.isnan(mean):
            continue
        chosen = min(vals, key=lambda v: abs(float(v.get("dx", math.nan)) - mean))
        run_dir = Path(str(chosen["run_dir"]))
        lines.append(
            f"- {case}: {run_dir.name} "
            f"(dx={float(chosen['dx']):+.3f} um; path={run_dir})"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    input_root = resolve_path(args.input_root)
    out_dir = resolve_path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    runs = discover_runs(input_root)
    if not runs:
        raise SystemExit(f"No completed run folders found under {input_root}")

    run_rows = []
    time_rows = []
    for run in runs:
        positions = get_body_positions(run)
        run_rows.append(analyze_run(run, positions))
        time_rows.extend(analyze_timecourse(run, positions))
    summary = summarize(run_rows)
    time_summary = summarize_timecourse(time_rows)
    write_csv(out_dir / "run_metrics.csv", run_rows)
    write_csv(out_dir / "condition_summary.csv", summary)
    write_csv(out_dir / "timecourse_run_metrics.csv", time_rows)
    write_csv(out_dir / "timecourse_condition_summary.csv", time_summary)
    write_representative_runs(out_dir / "representative_runs.txt", run_rows)
    plot_timecourse(time_summary, out_dir)
    plot_final_summary(run_rows, summary, out_dir)
    print(out_dir / "condition_summary.csv")


if __name__ == "__main__":
    main()
