#!/usr/bin/env python3
"""Analyze displacement and actin length in the one-sided turnover pilot."""

from __future__ import annotations

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
FINAL_FRAME = 60
N_FRAMES = 61
SECONDS_PER_FRAME = 1.0
TURNOVER_ONSET_FRAME = 4


def run_report(run_dir: Path, what: str, frame: str) -> str:
    return subprocess.check_output([str(REPORT), what, f"frame={frame}"], cwd=str(run_dir), text=True)


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


def get_body_position(run_dir: Path, frame: int) -> tuple[float, float, float] | None:
    for what in ("solid:position", "bead:position"):
        try:
            pos = parse_position(run_report(run_dir, what, str(frame)))
        except subprocess.CalledProcessError:
            pos = None
        if pos is not None:
            return pos
    return None


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


def parse_single_state(text: str) -> dict[str, float]:
    total = 0
    bound = 0
    for line in text.splitlines():
        if not line or line.startswith("%"):
            continue
        parts = line.split()
        if len(parts) < 10:
            continue
        total += 1
        try:
            if int(parts[8]) > 0:
                bound += 1
        except ValueError:
            pass
    return {
        "chewer_count": float(total),
        "bound_chewer_count": float(bound),
        "bound_chewer_fraction": float(bound) / total if total else math.nan,
    }


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


def analyze_run(run_dir: Path) -> dict[str, object]:
    row: dict[str, object] = {"run_dir": str(run_dir)}
    row.update(metadata(run_dir / "config.cym"))
    p0 = get_body_position(run_dir, 0)
    p_onset = get_body_position(run_dir, TURNOVER_ONSET_FRAME)
    p1 = get_body_position(run_dir, FINAL_FRAME)
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
                "x_velocity": dx / 60.0,
            }
        )
    if p_onset and p1:
        row.update(
            {
                "x_onset": p_onset[0],
                "dx_4to60": p1[0] - p_onset[0],
                "x_velocity_4to60": (p1[0] - p_onset[0]) / (FINAL_FRAME - TURNOVER_ONSET_FRAME),
            }
        )
    row.update(parse_fiber_length(run_report(run_dir, "fiber:length", str(FINAL_FRAME))))
    try:
        row.update(parse_single_state(run_report(run_dir, "single:state", str(FINAL_FRAME))))
    except subprocess.CalledProcessError:
        row.update({"chewer_count": 0.0, "bound_chewer_count": 0.0, "bound_chewer_fraction": math.nan})
    try:
        row.update(parse_fiber_end_geometry(run_report(run_dir, "fiber:end", str(FINAL_FRAME))))
    except subprocess.CalledProcessError:
        row.update({"plus_end_x_mean": math.nan, "minus_end_x_mean": math.nan, "plus_minus_x_delta": math.nan})
    return row


def analyze_timecourse(run_dir: Path) -> list[dict[str, object]]:
    meta = metadata(run_dir / "config.cym")
    p0 = get_body_position(run_dir, 0)
    if p0 is None:
        return []
    rows = []
    for frame in range(N_FRAMES):
        pos = get_body_position(run_dir, frame)
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
                "time_min": frame * SECONDS_PER_FRAME / 60.0,
                "dx": dx,
                "dy": dy,
                "dz": dz,
                "net_displacement": math.sqrt(dx * dx + dy * dy + dz * dz),
            }
        )
    return rows


def mean_sem(vals: list[float]) -> tuple[float, float]:
    xs = [x for x in vals if not math.isnan(x)]
    if not xs:
        return math.nan, math.nan
    mean = sum(xs) / len(xs)
    if len(xs) < 2:
        return mean, math.nan
    sd = math.sqrt(sum((x - mean) ** 2 for x in xs) / (len(xs) - 1))
    return mean, sd / math.sqrt(len(xs))


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[str(row["case"])].append(row)
    cols = [
        "dx",
        "dx_4to60",
        "net_displacement",
        "x_velocity",
        "x_velocity_4to60",
        "actin_count",
        "actin_total_len",
        "actin_avg_len",
        "actin_max_len",
        "chewer_count",
        "bound_chewer_count",
        "bound_chewer_fraction",
        "plus_end_x_mean",
        "minus_end_x_mean",
        "plus_minus_x_delta",
    ]
    out = []
    for case, vals in sorted(groups.items()):
        summary: dict[str, object] = {"case": case, "n": len(vals)}
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
        "01_symmetric_short_no_chewer": "#666666",
        "02_symmetric_short_one_sided_chewer": "#0072B2",
        "03_cap_short_no_chewer": "#D55E00",
        "04_cap_short_one_sided_chewer": "#009E73",
        "05_no_nucleator_fixed_chewer_control": "#CC79A7",
    }
    labels = {
        "01_symmetric_short_no_chewer": "symmetric short aster",
        "02_symmetric_short_one_sided_chewer": "symmetric + one-sided turnover",
        "03_cap_short_no_chewer": "one-sided cap aster",
        "04_cap_short_one_sided_chewer": "one-sided cap + turnover",
        "05_no_nucleator_fixed_chewer_control": "no nucleator control",
        "06_symmetric_short_delayed_one_sided_chewer": "symmetric + delayed turnover",
        "07_cap_short_delayed_one_sided_chewer": "one-sided cap + delayed turnover",
    }
    by_case: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_case[str(row["case"])].append(row)

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    for case, vals in sorted(by_case.items()):
        vals = sorted(vals, key=lambda r: int(r["frame"]))
        t = [float(v["time_s"]) for v in vals]
        y = [float(v["dx_mean"]) for v in vals]
        sem = [0.0 if math.isnan(float(v["dx_sem"])) else float(v["dx_sem"]) for v in vals]
        ax.plot(t, y, lw=2.2, color=colors.get(case, "black"), label=labels.get(case, case))
        ax.fill_between(t, [a - b for a, b in zip(y, sem)], [a + b for a, b in zip(y, sem)],
                        color=colors.get(case, "black"), alpha=0.16, linewidth=0)
    ax.axhline(0, color="0.55", lw=1.0, ls="--")
    ax.set_xlabel("Time (s)", fontsize=13)
    ax.set_ylabel("Nucleus displacement along x, dx (um)", fontsize=13)
    ax.tick_params(labelsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=8.5, loc="best")
    fig.tight_layout()
    fig.savefig(out_dir / "dx_timecourse.png", dpi=300)
    fig.savefig(out_dir / "dx_timecourse.svg")
    plt.close(fig)


def main() -> None:
    out_dir = ROOT / "analysis"
    out_dir.mkdir(exist_ok=True)
    runs = sorted(p for p in ROOT.glob("*/r*") if (p / "objects.cmo").exists())
    run_rows = [analyze_run(run) for run in runs]
    time_rows = []
    for run in runs:
        time_rows.extend(analyze_timecourse(run))
    summary = summarize(run_rows)
    time_summary = summarize_timecourse(time_rows)
    write_csv(out_dir / "run_metrics.csv", run_rows)
    write_csv(out_dir / "condition_summary.csv", summary)
    write_csv(out_dir / "timecourse_run_metrics.csv", time_rows)
    write_csv(out_dir / "timecourse_condition_summary.csv", time_summary)
    plot_timecourse(time_summary, out_dir)
    print(out_dir / "condition_summary.csv")


if __name__ == "__main__":
    main()
