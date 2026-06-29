#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
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
for path in (str(ROOT), str(ROOT / "analysis")):
    if path not in sys.path:
        sys.path.insert(0, path)

from analysis import turnover_campaign_analysis as core
from analysis import turnover_mini_chew_analysis as mini


SCENARIO_LABELS = {
    "small_long": "Small elongated",
    "large_long": "Large elongated",
}

SCENARIO_SHORT = {
    "small_long": "Small\nelongated",
    "large_long": "Large\nelongated",
}

CONDITION_LABELS = {
    "nomotor_xlink": "No motors, + xlinks",
    "rotatable_xlink": "Rotatable, + xlinks",
}

CONDITION_SHORT = {
    "nomotor_xlink": "No motors\n+ xlinks",
    "rotatable_xlink": "Rotatable\n+ xlinks",
}

CONDITION_COLORS = {
    "nomotor_xlink": "#4d4d4d",
    "rotatable_xlink": "#1b9e77",
}

TIMECOURSE_FIELDS = [
    ("vz_abs", "Mean |v_z| (um/s)", True),
    ("speckle_length_proxy_um", "Filament mass proxy (um)", True),
    ("chewer_band_mass_fraction", "Chewer-band mass fraction", True),
    ("bottom_half_mass_fraction", "Bottom-half mass fraction", True),
    ("axial_spread_entropy", "Axial spread entropy", True),
]

FIBER_LENGTH_FRAMES = (0, 20, 40, 80, 120, 160)
FRAME_TO_MIN = 1.0 / 20.0
ANALYSIS_START_MIN = 1.0


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Analyze scaled continuous top-supply minus-end-chew campaign outputs.")
    ap.add_argument("--campaign-root", type=Path, default=ROOT / "turnover_continuous_supply_scaled")
    ap.add_argument("--save-dir", type=Path, default=ROOT / "turnover_continuous_supply_scaled" / "job00" / "save")
    ap.add_argument("--outdir", type=Path, default=ROOT / "analysis" / "results" / "turnover_continuous_supply_scaled")
    ap.add_argument("--report-bin", type=Path, default=ROOT / "build_mini_turnover" / "bin" / "report")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--interval", type=float, default=2.0)
    ap.add_argument("--speckle-name", default="speckle_i2.txt")
    ap.add_argument("--analysis-start-min", type=float, default=ANALYSIS_START_MIN)
    return ap.parse_args()


def as_float(row: dict, key: str, default: float = float("nan")) -> float:
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return default


def source_zone_low(row: dict) -> float:
    bottom = as_float(row, "bottom_z")
    top = as_float(row, "top_z")
    if not np.isfinite(bottom) or not np.isfinite(top):
        return bottom
    return top - 0.30 * (top - bottom)


def load_manifest(campaign_root: Path, save_dir: Path) -> list[dict]:
    manifest_path = campaign_root / "manifest.csv"
    rows = []
    manifest_rows = [
        row
        for row in csv.DictReader(manifest_path.open(encoding="utf-8"))
        if row.get("kind") == "cluster"
    ]
    # The cluster submit helper runs `find cluster_runs -path '*/config.cym' | sort`,
    # so save/r0001 follows sorted config-path order, not manifest insertion order.
    for index, row in enumerate(sorted(manifest_rows, key=lambda item: item.get("path", ""))):
        size = row["size"]
        condition = row["condition"]
        run_dir = f"r{index + 1:04d}"
        rows.append(
            {
                "index": index,
                "scenario": size,
                "scenario_label": SCENARIO_LABELS.get(size, size),
                "scenario_short": SCENARIO_SHORT.get(size, size),
                "condition": condition,
                "condition_label": CONDITION_LABELS.get(condition, condition),
                "condition_short": CONDITION_SHORT.get(condition, condition),
                "motor_mode": "none" if condition == "nomotor_xlink" else "rotatable",
                "include_xlinks": 1,
                "n_clusters": int(float(row.get("n_clusters", 0) or 0)),
                "motors_per_cluster": int(float(row.get("motors_per_cluster", 0) or 0)),
                "total_motors": int(float(row.get("total_motors", 0) or 0)),
                "replicate": int(float(row["replicate"])),
                "run_dir": run_dir,
                "run_path": str(save_dir / run_dir),
                "scenario_z_low": source_zone_low(row),
                "scenario_z_high": as_float(row, "top_z"),
                "chewer_zone_bottom": as_float(row, "chewer_zone_bottom"),
                "chewer_zone_top": as_float(row, "chewer_zone_top"),
                "bottom_z": as_float(row, "bottom_z"),
                "top_z": as_float(row, "top_z"),
                "pulse_filaments": int(float(row.get("pulse_filaments", 0) or 0)),
                "n_pulses": int(float(row.get("n_pulses", 0) or 0)),
            }
        )
    return rows


def nearest_index(values: np.ndarray, target: float) -> int:
    if values.size == 0:
        return 0
    return int(np.nanargmin(np.abs(values - target)))


def reset_analysis_start(run_rows: list[dict], start_min: float) -> None:
    for row in run_rows:
        if row.get("status") != "ok":
            continue
        time_min = np.asarray(row["time_min"], dtype=float)
        idx = nearest_index(time_min, start_min)
        row["growth_end_min"] = start_min
        row["motor_onset_min"] = start_min
        row["speckle_length_proxy_motor_onset_um"] = float(np.asarray(row["speckle_length_proxy_um"], dtype=float)[idx])
        row["seed_zone_mass_fraction_motor_onset"] = float(np.asarray(row["seed_zone_mass_fraction"], dtype=float)[idx])
        row["chewer_band_mass_fraction_motor_onset"] = float(np.asarray(row["chewer_band_mass_fraction"], dtype=float)[idx])
        row["bottom_half_mass_fraction_motor_onset"] = float(np.asarray(row["bottom_half_mass_fraction"], dtype=float)[idx])
        row["axial_spread_entropy_motor_onset"] = float(np.asarray(row["axial_spread_entropy"], dtype=float)[idx])


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
                    "time_min": current_frame * FRAME_TO_MIN,
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
    frame_arg = "frame=" + ",".join(str(frame) for frame in FIBER_LENGTH_FRAMES)
    for row in manifest_rows:
        result = subprocess.run(
            [str(report_bin), "fiber:length", frame_arg],
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
            "time_min": frame * FRAME_TO_MIN,
            "n_replicates": len(rows),
        }
        for field in ["total_len", "off_len", "avg_len", "min_len", "max_len"]:
            mean, sd, sem, n = mini.finite_mean_sd([float(row[field]) for row in rows])
            item[field] = mean
            item[f"{field}_sd"] = sd
            item[f"{field}_sem"] = sem
            item[f"{field}_n"] = n
        out.append(item)
    return out


def write_summary(path: Path, condition_metrics: list[dict], fiber_condition_rows: list[dict], run_rows: list[dict], analysis_start_min: float) -> None:
    failed = [row for row in run_rows if row.get("status") != "ok"]
    ranked_mass = sorted(condition_metrics, key=lambda row: row.get("mass_retention", float("nan")))
    ranked_chewer = sorted(condition_metrics, key=lambda row: row.get("chewer_band_mass_fraction_final", float("nan")))
    final_fiber = {
        (row["scenario"], row["condition"]): row
        for row in fiber_condition_rows
        if int(row["frame"]) == max(FIBER_LENGTH_FRAMES)
    }
    lines = [
        "# Continuous top-supply scaled turnover analysis",
        "",
        "## Dataset",
        f"- Runs parsed successfully: `{len(run_rows) - len(failed)}` / `{len(run_rows)}`.",
        "- Cases: `small_long` and `large_long`, each with no-motor and rotatable-motor conditions.",
        f"- Movement and retention metrics are computed from `{analysis_start_min:.2f} min` onward.",
        "- Top/source-zone mass fraction uses the upper 30% of the annulus.",
        "- Chewer-band mass fraction uses the bottom 30% chewer zone from the campaign manifest.",
        "",
        "## Lowest final chewer-band mass fraction",
    ]
    for row in ranked_chewer:
        lines.append(
            f"- `{row['scenario']}` + `{row['condition']}`: final chewer-band fraction "
            f"`{row['chewer_band_mass_fraction_final']:.3f}`, bottom-half fraction `{row['bottom_half_mass_fraction_final']:.3f}`."
        )
    lines.extend(["", "## Lowest mass retention"])
    for row in ranked_mass:
        exact = final_fiber.get((row["scenario"], row["condition"]))
        extra = ""
        if exact is not None:
            extra = f" final fiber length `{exact['total_len']:.1f} um`, off length `{exact['off_len']:.1f} um`."
        lines.append(
            f"- `{row['scenario']}` + `{row['condition']}`: mass retention `{row['mass_retention']:.3f}`, "
            f"AUC(|v_z|) `{row['auc_vz_abs']:.3f} um`.{extra}"
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
            "- `condition_timecourses.csv`: mean time traces.",
            "- `fiber_length_frames.csv` and `fiber_length_condition_summary.csv`: exact `report fiber:length` summaries.",
            "- `movement_matrices/`: metric maps.",
            "- `timecourses/by_size/`: condition mean +/- SEM traces.",
            "- `kymographs/full_timeline_kymograph_grid_medoids.png`: representative axial-density trajectories.",
            "- `unwrapped_annulus/`: representative unwrapped heatmaps and line snapshots.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

    core.PHASE_NAMES = ["initial", "pulse_01", "pulse_02", "pulse_03", "pulse_04", "pulse_05", "pulse_06", "pulse_07"]
    core.SNAPSHOT_TIMES_MIN = [
        ("t000", 0.0, "0.0 min"),
        ("t060", 1.0, "1.0 min"),
        ("t120", 2.0, "2.0 min"),
        ("t240", 4.0, "4.0 min"),
        ("t360", 6.0, "6.0 min"),
        ("t480", 8.0, "8.0 min"),
    ]
    core.CONDITION_COLORS.update(CONDITION_COLORS)
    mini.SCENARIO_LABELS = SCENARIO_SHORT
    mini.CONDITION_LABELS = CONDITION_SHORT
    mini.CONDITION_COLORS = CONDITION_COLORS
    mini.TIMECOURSE_FIELDS = TIMECOURSE_FIELDS
    plt.rcParams.update({"font.size": 11, "axes.linewidth": 1.1, "savefig.facecolor": "white", "figure.facecolor": "white"})

    campaign_root = args.campaign_root.resolve()
    save_dir = args.save_dir.resolve()
    outdir = args.outdir.resolve()
    report_bin = args.report_bin.resolve()
    if not report_bin.exists():
        raise FileNotFoundError(report_bin)

    manifest_rows = load_manifest(campaign_root, save_dir)
    tasks = [(row, str(report_bin), args.interval, args.speckle_name) for row in manifest_rows]
    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            run_rows = list(executor.map(core.analyze_run, tasks))
    else:
        run_rows = [core.analyze_run(task) for task in tasks]

    reset_analysis_start(run_rows, args.analysis_start_min)
    thresholds = core.compute_active_thresholds(run_rows)
    matched_thresholds = core.compute_scenario_matched_thresholds(run_rows)
    run_metrics = core.compute_run_metrics(run_rows, thresholds, matched_thresholds)
    condition_metrics = mini.aggregate_metrics(run_metrics)
    medoid_rows = mini.pick_medoid_rows(run_rows, run_metrics, condition_metrics)
    condition_timecourses = mini.aggregate_timecourses(run_rows)
    fiber_rows = compute_fiber_length_rows(manifest_rows, report_bin)
    fiber_condition_rows = aggregate_fiber_length(fiber_rows)

    scenario_order = list(SCENARIO_LABELS)
    condition_order = list(CONDITION_LABELS)
    outdir.mkdir(parents=True, exist_ok=True)
    mini.write_run_status(outdir / "run_status.csv", run_rows)
    mini.write_csv(outdir / "run_metrics.csv", run_metrics)
    mini.write_csv(outdir / "condition_metrics.csv", condition_metrics)
    mini.write_csv(outdir / "condition_timecourses.csv", condition_timecourses)
    mini.write_csv(outdir / "fiber_length_frames.csv", fiber_rows)
    mini.write_csv(outdir / "fiber_length_condition_summary.csv", fiber_condition_rows)
    mini.write_csv(outdir / "active_thresholds.csv", [{"include_xlinks": key, "threshold_vz_abs": value} for key, value in sorted(thresholds.items())])
    mini.write_csv(
        outdir / "scenario_matched_active_thresholds.csv",
        [
            {"scenario": scenario, "include_xlinks": include_xlinks, "threshold_vz_abs_matched": value}
            for (scenario, include_xlinks), value in sorted(matched_thresholds.items())
        ],
    )
    core.write_timecourses_long(outdir / "run_timecourses_long.csv", run_rows)

    for spec in core.SUMMARY_METRICS:
        core.plot_metric_matrix(
            condition_metrics,
            scenario_order,
            condition_order,
            SCENARIO_SHORT,
            CONDITION_SHORT,
            spec,
            outdir / "movement_matrices" / f"{spec['slug']}.png",
        )
    core.plot_metric_overview(
        condition_metrics,
        scenario_order,
        condition_order,
        SCENARIO_SHORT,
        CONDITION_SHORT,
        outdir / "movement_matrices" / "movement_metric_overview.png",
    )
    mini.plot_condition_timecourses(condition_timecourses, scenario_order, condition_order, outdir / "timecourses" / "by_size")
    core.plot_kymograph_grid(
        medoid_rows,
        scenario_order,
        SCENARIO_SHORT,
        condition_order,
        CONDITION_SHORT,
        outdir / "kymographs" / "full_timeline_kymograph_grid_medoids.png",
    )
    core.plot_snapshot_heatmaps(medoid_rows, scenario_order, SCENARIO_SHORT, condition_order, CONDITION_SHORT, outdir / "unwrapped_annulus" / "heatmaps")
    core.plot_snapshot_lines(medoid_rows, scenario_order, SCENARIO_SHORT, condition_order, CONDITION_SHORT, outdir / "unwrapped_annulus" / "lines")
    write_summary(outdir / "README.md", condition_metrics, fiber_condition_rows, run_rows, args.analysis_start_min)

    print(outdir / "README.md")
    print(outdir / "run_status.csv")
    print(outdir / "condition_metrics.csv")
    print(outdir / "fiber_length_condition_summary.csv")


if __name__ == "__main__":
    main()
