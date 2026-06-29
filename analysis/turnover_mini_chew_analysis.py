#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS_DIR = ROOT / "analysis"
for path in (str(ROOT), str(ANALYSIS_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)

from analysis import turnover_campaign_analysis as core


SCENARIO_LABELS = {
    "full_height_aligned": "Full height\naligned",
    "top_two_thirds_aligned": "Top 2/3\naligned",
    "top_two_thirds_mixed_polarity": "Top 2/3\nmixed polarity",
}

SCENARIO_Z_BOUNDS = {
    "full_height_aligned": (-10.0, 10.0),
    "top_two_thirds_aligned": (-10.0 + 20.0 / 3.0, 10.0),
    "top_two_thirds_mixed_polarity": (-10.0 + 20.0 / 3.0, 10.0),
}

CONDITION_LABELS = {
    "nomotor_xlink": "No motors\n+ xlinks",
    "rotatable_xlink": "Rotatable\n+ xlinks",
    "fixed_global_xlink": "Fixed global\n+ xlinks",
}

CONDITION_COLORS = {
    "nomotor_xlink": "#4d4d4d",
    "rotatable_xlink": "#1b9e77",
    "fixed_global_xlink": "#d95f02",
}

CONFIG_CASE_RE = re.compile(r"Mini turnover case:\s*([^/]+?)\s*/\s*(\S+)")

NUMERIC_FIELDS = [
    "mean_vz_abs",
    "mean_vz_signed",
    "peak_vz_abs",
    "time_to_peak_min",
    "auc_vz_abs",
    "excess_auc_vz_abs",
    "directionality_index",
    "active_transport_fraction",
    "matched_active_transport_fraction",
    "longest_active_streak_min",
    "onset_time_after_motor_min",
    "net_z_displacement",
    "abs_net_z_displacement",
    "mass_normalized_auc",
    "mass_retention",
    "mean_swirl_rate_abs",
    "swirl_penalty",
    "seed_zone_mass_fraction_motor_onset",
    "seed_zone_mass_fraction_final",
    "seed_zone_retention",
    "seed_zone_mass_fraction_mean",
    "chewer_band_mass_fraction_motor_onset",
    "chewer_band_mass_fraction_final",
    "chewer_band_mass_fraction_gain",
    "chewer_band_mass_fraction_mean",
    "bottom_half_mass_fraction_motor_onset",
    "bottom_half_mass_fraction_final",
    "bottom_half_mass_fraction_gain",
    "bottom_half_mass_fraction_mean",
    "axial_spread_entropy_motor_onset",
    "axial_spread_entropy_final",
    "axial_spread_entropy_gain",
    "axial_spread_entropy_mean",
    "speckle_length_proxy_motor_onset_um",
    "speckle_length_proxy_final_um",
    "top_bias_final",
    "nematic_xy_final",
]

TIMECOURSE_FIELDS = [
    ("vz_abs", "Mean |v_z| (um/s)", True),
    ("speckle_length_proxy_um", "Filament mass proxy (um)", True),
    ("chewer_band_mass_fraction", "Chewer-band mass fraction", True),
    ("bottom_half_mass_fraction", "Bottom-half mass fraction", True),
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Analyze the replicated mini minus-end-chew turnover campaign.")
    ap.add_argument("--campaign-root", type=Path, default=ROOT / "turnover_mini_minus_chew")
    ap.add_argument("--save-dir", type=Path, default=ROOT / "turnover_mini_minus_chew" / "job00" / "save")
    ap.add_argument("--outdir", type=Path, default=ROOT / "analysis" / "results" / "turnover_mini_minus_chew")
    ap.add_argument("--report-bin", type=Path, default=ROOT / "build_mini_turnover" / "bin" / "report")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--interval", type=float, default=2.0)
    ap.add_argument("--speckle-name", default="speckle_i2.txt")
    return ap.parse_args()


def load_mini_manifest(campaign_root: Path, save_dir: Path) -> list[dict]:
    rows = []
    manifest_path = campaign_root / "manifest.csv"
    cluster_rows = [row for row in csv.DictReader(manifest_path.open(encoding="utf-8")) if row["kind"] == "cluster"]
    manifest_by_case: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in cluster_rows:
        manifest_by_case[(row["scenario"], row["condition"])].append(row)

    saved_dirs = sorted(p for p in save_dir.glob("r[0-9][0-9][0-9][0-9]") if p.is_dir())
    if saved_dirs:
        source_rows = []
        case_counts: dict[tuple[str, str], int] = defaultdict(int)
        for index, run_path in enumerate(saved_dirs):
            config_path = run_path / "config.cym"
            if not config_path.exists():
                continue
            scenario = condition = None
            for line in config_path.read_text(encoding="utf-8", errors="replace").splitlines()[:20]:
                match = CONFIG_CASE_RE.search(line)
                if match:
                    scenario = match.group(1).strip()
                    condition = match.group(2).strip()
                    break
            if not scenario or not condition:
                continue
            key = (scenario, condition)
            case_counts[key] += 1
            meta_rows = manifest_by_case.get(key, [])
            meta = meta_rows[case_counts[key] - 1] if case_counts[key] <= len(meta_rows) else {}
            source_rows.append((index, run_path, case_counts[key], scenario, condition, meta))
    else:
        source_rows = []
        for index, row in enumerate(cluster_rows):
            source_rows.append((index, save_dir / f"r{index + 1:04d}", int(row["replicate"]), row["scenario"], row["condition"], row))

    for index, run_path, replicate, scenario, condition, row in source_rows:
        z_low, z_high = SCENARIO_Z_BOUNDS[scenario]
        rows.append(
            {
                "index": index,
                "scenario": scenario,
                "scenario_label": SCENARIO_LABELS[scenario].replace("\n", ", "),
                "scenario_short": SCENARIO_LABELS[scenario],
                "condition": condition,
                "condition_label": CONDITION_LABELS[condition].replace("\n", ", "),
                "condition_short": CONDITION_LABELS[condition],
                "motor_mode": row.get("motor_mode", "none" if condition == "nomotor_xlink" else condition.replace("_xlink", "")),
                "include_xlinks": 1,
                "n_clusters": int(row.get("n_clusters", 0) or 0),
                "motors_per_cluster": int(row.get("motors_per_cluster", 0) or 0),
                "total_motors": int(row.get("total_motors", 0) or 0),
                "replicate": int(replicate),
                "run_dir": run_path.name,
                "run_path": str(run_path),
                "scenario_z_low": z_low,
                "scenario_z_high": z_high,
                "chewer_zone_bottom": -10.0,
                "chewer_zone_top": -4.0,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_run_status(path: Path, run_rows: list[dict]) -> None:
    rows = []
    for row in run_rows:
        rows.append(
            {
                "scenario": row.get("scenario"),
                "condition": row.get("condition"),
                "run_dir": row.get("run_dir"),
                "run_path": row.get("run_path"),
                "status": row.get("status"),
                "frame_mode": row.get("frame_mode", ""),
                "growth_end_min": row.get("growth_end_min", ""),
                "motor_onset_min": row.get("motor_onset_min", ""),
            }
        )
    write_csv(path, rows)


def parse_fiber_length_output(text: str, row: dict) -> list[dict]:
    out = []
    current_frame: int | None = None
    for line in text.splitlines():
        if line.startswith("% frame"):
            parts = line.split()
            if len(parts) >= 3:
                current_frame = int(parts[2])
            continue
        if current_frame is None or line.startswith("%") or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 8 or parts[0] == "class":
            continue
        try:
            out.append(
                {
                    "scenario": row["scenario"],
                    "condition": row["condition"],
                    "replicate": row["replicate"],
                    "run_dir": row["run_dir"],
                    "frame": current_frame,
                    "time_s": float(current_frame),
                    "time_min": float(current_frame) / 60.0,
                    "fiber_class": parts[0],
                    "count": int(parts[1]),
                    "avg_len": float(parts[2]),
                    "var_len": float(parts[3]),
                    "min_len": float(parts[4]),
                    "max_len": float(parts[5]),
                    "total_len": float(parts[6]),
                    "off_len": float(parts[7]),
                }
            )
        except ValueError:
            continue
    return out


def compute_fiber_length_rows(manifest_rows: list[dict], report_bin: Path) -> list[dict]:
    rows = []
    for row in manifest_rows:
        result = subprocess.run(
            [str(report_bin), "fiber:length", "frame=0,120,180,780"],
            cwd=row["run_path"],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=True,
        )
        rows.extend(parse_fiber_length_output(result.stdout, row))
    return rows


def aggregate_fiber_length(fiber_rows: list[dict]) -> list[dict]:
    grouped: dict[tuple[str, str, int], list[dict]] = defaultdict(list)
    for row in fiber_rows:
        grouped[(row["scenario"], row["condition"], int(row["frame"]))].append(row)

    out = []
    for (scenario, condition, frame), rows in sorted(grouped.items()):
        item = {
            "scenario": scenario,
            "condition": condition,
            "frame": frame,
            "time_min": frame / 60.0,
            "n_replicates": len(rows),
        }
        for field in ["total_len", "off_len", "avg_len", "min_len", "max_len"]:
            mean, sd, sem, n = finite_mean_sd([float(row[field]) for row in rows])
            item[field] = mean
            item[f"{field}_sd"] = sd
            item[f"{field}_sem"] = sem
            item[f"{field}_n"] = n
        out.append(item)
    return out


def finite_mean_sd(values: list[float]) -> tuple[float, float, float, int]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan"), float("nan"), float("nan"), 0
    mean = float(np.mean(arr))
    sd = float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0
    sem = float(sd / math.sqrt(arr.size)) if arr.size > 0 else float("nan")
    return mean, sd, sem, int(arr.size)


def aggregate_metrics(metric_rows: list[dict]) -> list[dict]:
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in metric_rows:
        grouped[(row["scenario"], row["condition"])].append(row)

    out = []
    for (_scenario, _condition), rows in grouped.items():
        first = rows[0]
        agg = {
            "scenario": first["scenario"],
            "scenario_label": first["scenario_label"],
            "scenario_short": first["scenario_short"],
            "condition": first["condition"],
            "condition_label": first["condition_label"],
            "condition_short": first["condition_short"],
            "motor_mode": first["motor_mode"],
            "include_xlinks": first["include_xlinks"],
            "n_clusters": first["n_clusters"],
            "motors_per_cluster": first["motors_per_cluster"],
            "total_motors": first["total_motors"],
            "n_replicates": len(rows),
            "run_dir": "mean",
            "run_path": "",
            "growth_end_min": first["growth_end_min"],
            "motor_onset_min": first["motor_onset_min"],
            "threshold_vz_abs": first["threshold_vz_abs"],
            "threshold_vz_abs_matched": first["threshold_vz_abs_matched"],
            "zmin": first["zmin"],
            "zmax": first["zmax"],
            "scenario_z_low": first["scenario_z_low"],
            "scenario_z_high": first["scenario_z_high"],
            "chewer_zone_bottom": first["chewer_zone_bottom"],
            "chewer_zone_top": first["chewer_zone_top"],
        }
        for field in NUMERIC_FIELDS:
            mean, sd, sem, n = finite_mean_sd([float(row.get(field, float("nan"))) for row in rows])
            agg[field] = mean
            agg[f"{field}_sd"] = sd
            agg[f"{field}_sem"] = sem
            agg[f"{field}_n"] = n
        out.append(agg)
    return out


def pick_medoid_rows(run_rows: list[dict], metric_rows: list[dict], condition_metrics: list[dict]) -> list[dict]:
    metric_by_run = {(row["scenario"], row["condition"], row["run_dir"]): row for row in metric_rows}
    run_by_key = {(row["scenario"], row["condition"], row["run_dir"]): row for row in run_rows}
    out = []
    for cond in condition_metrics:
        candidates = [
            row for row in metric_rows
            if row["scenario"] == cond["scenario"] and row["condition"] == cond["condition"]
        ]
        if not candidates:
            continue
        target = float(cond.get("auc_vz_abs", float("nan")))
        if not np.isfinite(target):
            chosen = candidates[0]
        else:
            chosen = min(candidates, key=lambda row: abs(float(row.get("auc_vz_abs", float("nan"))) - target))
        key = (chosen["scenario"], chosen["condition"], chosen["run_dir"])
        medoid_run = dict(run_by_key[key])
        medoid_run["run_dir"] = chosen["run_dir"]
        medoid_run["medoid_for_auc_vz_abs"] = metric_by_run[key]["auc_vz_abs"]
        out.append(medoid_run)
    return out


def aggregate_timecourses(run_rows: list[dict]) -> list[dict]:
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in run_rows:
        if row.get("status") == "ok":
            grouped[(row["scenario"], row["condition"])].append(row)

    rows = []
    for (scenario, condition), group in grouped.items():
        first = group[0]
        time = np.asarray(first["time_min"], dtype=float)
        phase = first["phase_name"]
        for idx, tval in enumerate(time):
            out = {
                "scenario": scenario,
                "condition": condition,
                "time_min": float(tval),
                "phase_name": phase[idx],
                "n_replicates": len(group),
            }
            for field, _label, _clip in TIMECOURSE_FIELDS:
                vals = [float(np.asarray(row[field], dtype=float)[idx]) for row in group]
                mean, sd, sem, n = finite_mean_sd(vals)
                out[f"{field}_mean"] = mean
                out[f"{field}_sd"] = sd
                out[f"{field}_sem"] = sem
                out[f"{field}_n"] = n
            rows.append(out)
    return rows


def plot_condition_timecourses(condition_timecourses: list[dict], scenario_order: list[str], condition_order: list[str], outdir: Path) -> list[Path]:
    by_key: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in condition_timecourses:
        by_key[(row["scenario"], row["condition"])].append(row)

    outpaths = []
    for scenario in scenario_order:
        fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
        axes = axes.ravel()
        for ax, (field, ylabel, clip_zero) in zip(axes, TIMECOURSE_FIELDS):
            ymin = float("inf")
            ymax = float("-inf")
            for condition in condition_order:
                rows = sorted(by_key.get((scenario, condition), []), key=lambda row: row["time_min"])
                if not rows:
                    continue
                x = np.asarray([row["time_min"] for row in rows], dtype=float)
                y = np.asarray([row[f"{field}_mean"] for row in rows], dtype=float)
                sem = np.asarray([row[f"{field}_sem"] for row in rows], dtype=float)
                color = CONDITION_COLORS.get(condition, "#333333")
                ax.plot(x, y, color=color, lw=2.0, label=CONDITION_LABELS[condition].replace("\n", " "))
                ax.fill_between(x, y - sem, y + sem, color=color, alpha=0.15, linewidth=0)
                finite = np.concatenate([y[np.isfinite(y)], (y - sem)[np.isfinite(y - sem)], (y + sem)[np.isfinite(y + sem)]])
                if finite.size:
                    ymin = min(ymin, float(np.nanmin(finite)))
                    ymax = max(ymax, float(np.nanmax(finite)))
            ax.axvline(2.0, color="#666666", linestyle="--", lw=1.1)
            ax.axvline(3.0, color="#111111", linestyle=":", lw=1.3)
            if np.isfinite(ymin) and np.isfinite(ymax):
                span = ymax - ymin if ymax > ymin else max(abs(ymax), 1.0)
                ymin = max(0.0, ymin - 0.08 * span) if clip_zero else ymin - 0.08 * span
                ymax = ymax + 0.08 * span
                ax.set_ylim(ymin, ymax)
            ax.set_ylabel(ylabel, fontsize=12)
            core.style(ax)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=11)
        for ax in axes[-2:]:
            ax.set_xlabel("Time (min)", fontsize=12)
        fig.suptitle(f"{SCENARIO_LABELS[scenario].replace(chr(10), ' ')}: replicate mean +/- SEM", fontsize=17, y=0.98)
        fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.14, wspace=0.24, hspace=0.26)
        outpath = outdir / f"{scenario}_condition_mean_sem.png"
        outpath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outpath, dpi=300, bbox_inches="tight")
        plt.close(fig)
        outpaths.append(outpath)
    return outpaths


def write_summary(path: Path, condition_metrics: list[dict], medoid_rows: list[dict], run_rows: list[dict], fiber_condition_rows: list[dict]) -> None:
    ranked_auc = sorted(condition_metrics, key=lambda row: row["auc_vz_abs"], reverse=True)
    ranked_excess = sorted(condition_metrics, key=lambda row: row["excess_auc_vz_abs"], reverse=True)
    ranked_chewer = sorted(condition_metrics, key=lambda row: row["chewer_band_mass_fraction_final"])
    failed = [row for row in run_rows if row.get("status") != "ok"]
    medoids = {(row["scenario"], row["condition"]): row["run_dir"] for row in medoid_rows}
    final_fiber = {
        (row["scenario"], row["condition"]): row
        for row in fiber_condition_rows
        if int(row["frame"]) == 780
    }
    lines = [
        "# Mini minus-end-chew turnover analysis",
        "",
        "## Dataset",
        f"- Runs parsed successfully: `{len(run_rows) - len(failed)}` / `{len(run_rows)}`.",
        "- Replicates: `5` per scenario/condition.",
        "- Growth-only phase ends at `2.0 min`; crosslink conditioning ends and motor readout starts at `3.0 min`.",
        "- Movement metrics are computed only from the post-motor window.",
        "- Matrix plots use replicate means. Kymographs and unwrapped line snapshots use the replicate whose AUC is closest to that condition mean.",
        "",
        "## Top conditions by mean AUC(|v_z|)",
    ]
    for row in ranked_auc[:6]:
        medoid = medoids.get((row["scenario"], row["condition"]), "")
        lines.append(
            f"- `{row['scenario']}` + `{row['condition']}`: AUC `{row['auc_vz_abs']:.3f} +/- {row['auc_vz_abs_sem']:.3f} um SEM, "
            f"mean |v_z| `{row['mean_vz_abs']:.5f} um/s`, medoid `{medoid}`."
        )
    lines.extend(["", "## Top conditions by excess AUC above matched no-motor control"])
    for row in ranked_excess[:6]:
        lines.append(
            f"- `{row['scenario']}` + `{row['condition']}`: excess AUC `{row['excess_auc_vz_abs']:.3f} +/- {row['excess_auc_vz_abs_sem']:.3f} um SEM, "
            f"matched-active fraction `{row['matched_active_transport_fraction']:.3f}`."
        )
    lines.extend(["", "## Lowest final chewer-band mass fraction"])
    for row in ranked_chewer[:6]:
        lines.append(
            f"- `{row['scenario']}` + `{row['condition']}`: final chewer-band fraction `{row['chewer_band_mass_fraction_final']:.3f}`, "
            f"bottom-half fraction `{row['bottom_half_mass_fraction_final']:.3f}`."
        )
    lines.extend(["", "## Exact fiber:length at final frame"])
    for row in ranked_auc[:6]:
        exact = final_fiber.get((row["scenario"], row["condition"]))
        if exact is None:
            continue
        lines.append(
            f"- `{row['scenario']}` + `{row['condition']}`: total length `{exact['total_len']:.1f} um`, "
            f"off/chewed length `{exact['off_len']:.1f} um`."
        )
    if failed:
        lines.extend(["", "## Failed runs"])
        for row in failed:
            lines.append(f"- `{row.get('run_dir')}`: {row.get('status')}")
    lines.extend(
        [
            "",
            "## Files",
            "- `run_metrics.csv`: one row per replicate.",
            "- `condition_metrics.csv`: mean/SD/SEM across five replicates.",
            "- `condition_timecourses.csv`: replicate-mean time traces.",
            "- `fiber_length_frames.csv` and `fiber_length_condition_summary.csv`: exact `report fiber:length` summaries.",
            "- `movement_matrices/`: condition-level metric maps.",
            "- `timecourses/by_scenario/`: mean +/- SEM traces.",
            "- `kymographs/full_timeline_kymograph_grid_medoids.png`: representative axial-density trajectories.",
            "- `unwrapped_annulus/lines/`: representative line snapshots.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

    core.SNAPSHOT_TIMES_MIN = [
        ("t000", 0.0, "0.0 min"),
        ("t120", 120.0 / 60.0, "2.0 min"),
        ("t180", 180.0 / 60.0, "3.0 min"),
        ("t480", 480.0 / 60.0, "8.0 min"),
        ("t780", 780.0 / 60.0, "13.0 min"),
    ]
    core.CONDITION_COLORS.update(CONDITION_COLORS)
    plt.rcParams.update({"font.size": 11, "axes.linewidth": 1.1, "savefig.facecolor": "white", "figure.facecolor": "white"})

    campaign_root = args.campaign_root.resolve()
    save_dir = args.save_dir.resolve()
    outdir = args.outdir.resolve()
    report_bin = args.report_bin.resolve()
    if not report_bin.exists():
        raise FileNotFoundError(report_bin)

    manifest_rows = load_mini_manifest(campaign_root, save_dir)
    tasks = [(row, str(report_bin), args.interval, args.speckle_name) for row in manifest_rows]
    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            run_rows = list(executor.map(core.analyze_run, tasks))
    else:
        run_rows = [core.analyze_run(task) for task in tasks]

    thresholds = core.compute_active_thresholds(run_rows)
    matched_thresholds = core.compute_scenario_matched_thresholds(run_rows)
    run_metrics = core.compute_run_metrics(run_rows, thresholds, matched_thresholds)
    condition_metrics = aggregate_metrics(run_metrics)
    medoid_rows = pick_medoid_rows(run_rows, run_metrics, condition_metrics)
    condition_timecourses = aggregate_timecourses(run_rows)
    fiber_rows = compute_fiber_length_rows(manifest_rows, report_bin)
    fiber_condition_rows = aggregate_fiber_length(fiber_rows)

    scenario_order = list(SCENARIO_LABELS)
    condition_order = list(CONDITION_LABELS)
    scenario_labels = SCENARIO_LABELS
    condition_labels = CONDITION_LABELS

    outdir.mkdir(parents=True, exist_ok=True)
    write_run_status(outdir / "run_status.csv", run_rows)
    write_csv(outdir / "run_metrics.csv", run_metrics)
    write_csv(outdir / "condition_metrics.csv", condition_metrics)
    write_csv(outdir / "condition_timecourses.csv", condition_timecourses)
    write_csv(outdir / "fiber_length_frames.csv", fiber_rows)
    write_csv(outdir / "fiber_length_condition_summary.csv", fiber_condition_rows)
    write_csv(outdir / "active_thresholds.csv", [{"include_xlinks": key, "threshold_vz_abs": value} for key, value in sorted(thresholds.items())])
    write_csv(
        outdir / "scenario_matched_active_thresholds.csv",
        [
            {"scenario": scenario, "include_xlinks": include_xlinks, "threshold_vz_abs_matched": value}
            for (scenario, include_xlinks), value in sorted(matched_thresholds.items())
        ],
    )
    core.write_timecourses_long(outdir / "run_timecourses_long.csv", run_rows)

    metric_paths = []
    for spec in core.SUMMARY_METRICS:
        metric_paths.append(
            core.plot_metric_matrix(
                condition_metrics,
                scenario_order,
                condition_order,
                scenario_labels,
                condition_labels,
                spec,
                outdir / "movement_matrices" / f"{spec['slug']}.png",
            )
        )
    overview_path = core.plot_metric_overview(
        condition_metrics,
        scenario_order,
        condition_order,
        scenario_labels,
        condition_labels,
        outdir / "movement_matrices" / "movement_metric_overview.png",
    )
    timecourse_paths = plot_condition_timecourses(condition_timecourses, scenario_order, condition_order, outdir / "timecourses" / "by_scenario")
    kymograph_path = core.plot_kymograph_grid(
        medoid_rows,
        scenario_order,
        scenario_labels,
        condition_order,
        condition_labels,
        outdir / "kymographs" / "full_timeline_kymograph_grid_medoids.png",
    )
    heatmap_paths = core.plot_snapshot_heatmaps(medoid_rows, scenario_order, scenario_labels, condition_order, condition_labels, outdir / "unwrapped_annulus" / "heatmaps")
    line_paths = core.plot_snapshot_lines(medoid_rows, scenario_order, scenario_labels, condition_order, condition_labels, outdir / "unwrapped_annulus" / "lines")
    write_summary(outdir / "README.md", condition_metrics, medoid_rows, run_rows, fiber_condition_rows)

    print(outdir / "README.md")
    print(outdir / "run_status.csv")
    print(outdir / "run_metrics.csv")
    print(outdir / "condition_metrics.csv")
    print(outdir / "condition_timecourses.csv")
    print(outdir / "fiber_length_frames.csv")
    print(outdir / "fiber_length_condition_summary.csv")
    print(overview_path)
    print(kymograph_path)
    for path in metric_paths + timecourse_paths + heatmap_paths + line_paths:
        print(path)


if __name__ == "__main__":
    main()
