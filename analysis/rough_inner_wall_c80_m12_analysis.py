#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(ROOT / "analysis"))

import plot_timecourse_panels as base
import run_speckle_kymograph_batch as kymo


REPORT = ROOT / "build_rough_tools" / "bin" / "report"
OUT = ROOT / "analysis" / "results" / "rough_inner_wall_c80_m12"
PLOTS = OUT / "plots"
SMOOTH_ROT = ROOT / "clu_init_minusz" / "mpc12" / "c80_m12"
SMOOTH_FIX = ROOT / "clu_fixed_global_init_minusz" / "mpc12" / "c80_m12"
ROUGH_FIX = ROOT / "clu_rough_inner_wall_c80_m12_init_minusz" / "job11" / "save"
ROUGH_ROT = ROOT / "clu_rough_inner_wall_c80_m12_init_minusz" / "job12" / "save"
THRESHOLDS = ROOT / "analysis" / "results" / "init_minusz_comparison" / "transport_regime" / "control_thresholds.csv"

TARGET_FINAL_MIN = 800.0 / 60.0
REGIMES = ["1:4", "1:8", "1:16"]
WINDOWS = [
    ("early", 0.0, 4.0),
    ("mid", 4.0, 9.0),
    ("late", 9.0, TARGET_FINAL_MIN + 1e-9),
]
MODEL_SPECS = {
    "smooth_rotatable": ("Smooth rotatable", "#222222"),
    "rough_rotatable": ("Rough rotatable", "#00a6b4"),
    "smooth_fixed_global": ("Smooth fixed-global", "#d95f02"),
    "rough_fixed_global": ("Rough fixed-global", "#7b61c8"),
}
MODEL_ORDER = ["smooth_rotatable", "rough_rotatable", "smooth_fixed_global", "rough_fixed_global"]
METRICS = [
    ("mean_vz_abs", "Mean |v_z| (um/s)"),
    ("late_mean_vz_abs", "Late mean |v_z| (um/s)"),
    ("active_transport_fraction", "Active transport fraction"),
    ("abs_net_z_displacement", "|Net z shift| (um)"),
    ("net_z_displacement", "Net z shift (um)"),
]


def regime_for_run(idx: int) -> str:
    if idx <= 10:
        return "1:4"
    if idx <= 20:
        return "1:8"
    return "1:16"


def metadata() -> list[dict]:
    specs = [
        ("smooth_rotatable", SMOOTH_ROT),
        ("rough_rotatable", ROUGH_ROT),
        ("smooth_fixed_global", SMOOTH_FIX),
        ("rough_fixed_global", ROUGH_FIX),
    ]
    rows = []
    for model, root in specs:
        for idx in range(1, 31):
            run_dir = f"r{idx:04d}"
            run_path = root / run_dir
            rows.append(
                {
                    "model": model,
                    "group": "mpc12",
                    "case": "c80_m12",
                    "xlink_regime": regime_for_run(idx),
                    "run_dir": run_dir,
                    "run_path": str(run_path),
                    "has_outputs": int((run_path / "objects.cmo").exists() and (run_path / "config.cym").exists()),
                }
            )
    return [row for row in rows if row["has_outputs"]]


def read_thresholds() -> dict[str, float]:
    rows = csv.DictReader(THRESHOLDS.open(encoding="utf-8"))
    return {row["xlink_regime"]: float(row["threshold_vz_abs"]) for row in rows}


def trim_to_target(out: dict) -> dict:
    if out.get("status") != "ok":
        return out
    time = np.asarray(out["time_min"], dtype=float)
    keep = int(np.count_nonzero(time <= TARGET_FINAL_MIN + 1e-9))
    if keep <= 0:
        return out
    for key in ["time_min", "z_com", "vz_abs", "rho_z"]:
        if key in out:
            out[key] = out[key][:keep]
    return out


def analyze_row(task: tuple[dict, str]) -> dict:
    row, report_s = task
    out = kymo.analyze_run((row, 2.0, "speckle_i2.txt", report_s, 80))
    if out.get("status") != "ok":
        return out
    cfg = base.parse_config_info(Path(row["run_path"]) / "config.cym")
    rho = np.asarray(out["rho_z"], dtype=float)
    z_edges = np.asarray(out["z_edges"], dtype=float)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    z_com = rho @ z_centers / np.clip(rho.sum(axis=1), 1e-12, None)
    time = np.asarray(out["time_min"], dtype=float)
    if time.size > 1:
        dt_s = float((time[1] - time[0]) * 60.0)
    else:
        dt_s = float(cfg["runs"][-1]["frame_dt"])
    out["dt_s"] = dt_s
    out["z_com"] = z_com.tolist()
    out["vz_abs"] = np.abs(np.gradient(z_com, dt_s)).tolist() if z_com.size > 1 else [0.0]
    return trim_to_target(out)


def first_sustained_time(time_min: np.ndarray, signal: np.ndarray, threshold: float, consecutive: int = 3) -> float:
    if not np.isfinite(threshold):
        return float("nan")
    count = 0
    for idx, flag in enumerate(signal > threshold):
        count = count + 1 if flag else 0
        if count >= consecutive:
            return float(time_min[idx - consecutive + 1])
    return float("nan")


def window_mean(time_min: np.ndarray, values: np.ndarray, start: float, stop: float) -> float:
    mask = (time_min >= start) & (time_min < stop)
    return float(np.nanmean(values[mask])) if np.any(mask) else float("nan")


def compute_run_metrics(run_rows: list[dict], thresholds: dict[str, float]) -> list[dict]:
    out = []
    for row in run_rows:
        if row.get("status") != "ok":
            continue
        time_min = np.asarray(row["time_min"], dtype=float)
        time_s = time_min * 60.0
        vz_abs = np.asarray(row["vz_abs"], dtype=float)
        z_com = np.asarray(row["z_com"], dtype=float)
        threshold = thresholds.get(row["xlink_regime"], float("nan"))
        peak_idx = int(np.nanargmax(vz_abs)) if vz_abs.size else 0
        metric_row = {
            "model": row["model"],
            "group": row["group"],
            "case": row["case"],
            "xlink_regime": row["xlink_regime"],
            "run_dir": row["run_dir"],
            "run_path": row["run_path"],
            "n_frames_used": int(vz_abs.size),
            "threshold_vz_abs": threshold,
            "mean_vz_abs": float(np.nanmean(vz_abs)),
            "peak_vz_abs": float(np.nanmax(vz_abs)),
            "time_to_peak_min": float(time_min[peak_idx]) if time_min.size else float("nan"),
            "auc_vz_abs": float(np.trapezoid(vz_abs, time_s)) if time_s.size else float("nan"),
            "active_transport_fraction": float(np.mean(vz_abs > threshold)) if np.isfinite(threshold) else float("nan"),
            "onset_time_min": first_sustained_time(time_min, vz_abs, threshold),
            "net_z_displacement": float(z_com[-1] - z_com[0]) if z_com.size else float("nan"),
            "abs_net_z_displacement": float(abs(z_com[-1] - z_com[0])) if z_com.size else float("nan"),
        }
        for slug, start, stop in WINDOWS:
            metric_row[f"{slug}_mean_vz_abs"] = window_mean(time_min, vz_abs, start, stop)
        out.append(metric_row)
    return out


def aggregate_condition_metrics(rows: list[dict]) -> list[dict]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["model"], row["xlink_regime"])].append(row)
    metric_fields = [
        "mean_vz_abs",
        "late_mean_vz_abs",
        "active_transport_fraction",
        "abs_net_z_displacement",
        "net_z_displacement",
        "peak_vz_abs",
        "auc_vz_abs",
        "early_mean_vz_abs",
        "mid_mean_vz_abs",
    ]
    out = []
    for (model, regime), items in sorted(grouped.items()):
        row = {"model": model, "xlink_regime": regime, "n_runs_used": len(items)}
        for field in metric_fields:
            vals = np.asarray([float(item[field]) for item in items], dtype=float)
            vals = vals[np.isfinite(vals)]
            row[f"{field}_mean"] = float(np.nanmean(vals)) if vals.size else float("nan")
            row[f"{field}_sem"] = float(np.nanstd(vals, ddof=1) / math.sqrt(vals.size)) if vals.size > 1 else 0.0
        out.append(row)
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def effect_rows(condition_rows: list[dict]) -> list[dict]:
    lookup = {(row["model"], row["xlink_regime"]): row for row in condition_rows}
    out = []
    for regime in REGIMES:
        for orientation in ["rotatable", "fixed_global"]:
            smooth = lookup.get((f"smooth_{orientation}", regime))
            rough = lookup.get((f"rough_{orientation}", regime))
            if smooth is None or rough is None:
                continue
            for metric, _ in METRICS:
                s = float(smooth[f"{metric}_mean"])
                r = float(rough[f"{metric}_mean"])
                out.append(
                    {
                        "orientation": orientation,
                        "xlink_regime": regime,
                        "metric": metric,
                        "smooth_mean": s,
                        "rough_mean": r,
                        "rough_minus_smooth": r - s,
                        "rough_over_smooth": r / s if s != 0 and np.isfinite(s) else float("nan"),
                    }
                )
    return out


def grouped_timecourses(run_rows: list[dict]) -> dict[tuple[str, str], dict]:
    grouped = defaultdict(list)
    for row in run_rows:
        if row.get("status") == "ok":
            grouped[(row["model"], row["xlink_regime"])].append(row)
    out = {}
    for key, items in grouped.items():
        n = min(len(item["time_min"]) for item in items)
        time = np.asarray(items[0]["time_min"][:n], dtype=float)
        z_stack = np.stack([np.asarray(item["z_com"][:n], dtype=float) for item in items], axis=0)
        rho_stack = np.stack([np.asarray(item["rho_z"][:n], dtype=float) for item in items], axis=0)
        out[key] = {
            "time_min": time,
            "z_edges": np.asarray(items[0]["z_edges"], dtype=float),
            "z_mean": np.nanmean(z_stack, axis=0),
            "z_sem": np.nanstd(z_stack, axis=0, ddof=1) / math.sqrt(z_stack.shape[0]) if z_stack.shape[0] > 1 else np.zeros(n),
            "rho_mean": np.nanmean(rho_stack, axis=0),
            "n": z_stack.shape[0],
        }
    return out


def plot_metric_summary(condition_rows: list[dict]) -> Path:
    lookup = {(row["model"], row["xlink_regime"]): row for row in condition_rows}
    fig, axes = plt.subplots(len(METRICS), len(REGIMES), figsize=(15.5, 14.0), sharex=True)
    for ridx, (metric, ylabel) in enumerate(METRICS):
        for cidx, regime in enumerate(REGIMES):
            ax = axes[ridx, cidx]
            vals = []
            errs = []
            labels = []
            colors = []
            for model in MODEL_ORDER:
                row = lookup.get((model, regime))
                vals.append(float(row[f"{metric}_mean"]) if row else float("nan"))
                errs.append(float(row[f"{metric}_sem"]) if row else float("nan"))
                labels.append(MODEL_SPECS[model][0].replace(" ", "\n"))
                colors.append(MODEL_SPECS[model][1])
            x = np.arange(len(vals))
            ax.bar(x, vals, yerr=errs, capsize=4, color=colors, edgecolor="#222222", linewidth=0.8)
            ax.axhline(0, color="#555555", lw=0.8)
            if ridx == 0:
                ax.set_title(regime, fontsize=16)
            if cidx == 0:
                ax.set_ylabel(ylabel, fontsize=13)
            ax.tick_params(axis="both", labelsize=11)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.grid(axis="y", color="#dddddd", linewidth=0.6, alpha=0.75)
            if ridx == len(METRICS) - 1:
                ax.set_xticks(x)
                ax.set_xticklabels(labels, rotation=0, ha="center", fontsize=10)
            else:
                ax.set_xticks([])
    fig.suptitle("c80/m12 initially minus-z: smooth vs rough inner wall", fontsize=18, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    out = PLOTS / "metric_summary_by_regime.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=260)
    plt.close(fig)
    return out


def plot_z_com_1to8(timecourses: dict[tuple[str, str], dict]) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), sharey=True)
    for ax, orientation in zip(axes, ["rotatable", "fixed_global"]):
        for surface in ["smooth", "rough"]:
            model = f"{surface}_{orientation}"
            tc = timecourses.get((model, "1:8"))
            if tc is None:
                continue
            label, color = MODEL_SPECS[model]
            t = tc["time_min"]
            z = tc["z_mean"]
            sem = tc["z_sem"]
            ax.plot(t, z, color=color, lw=2.4, label=f"{label} (n={tc['n']})")
            ax.fill_between(t, z - sem, z + sem, color=color, alpha=0.18, linewidth=0)
        ax.axhline(0, color="#777777", lw=0.8)
        ax.set_title(orientation.replace("_", "-"), fontsize=15)
        ax.set_xlabel("Time (min)", fontsize=13)
        ax.grid(axis="y", color="#dddddd", linewidth=0.7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False, fontsize=10)
    axes[0].set_ylabel("Axial center of actin density, z (um)", fontsize=13)
    fig.tight_layout()
    out = PLOTS / "z_com_timecourse_1to8.png"
    fig.savefig(out, dpi=260)
    plt.close(fig)
    return out


def plot_mean_kymographs(timecourses: dict[tuple[str, str], dict]) -> Path:
    vals = []
    for tc in timecourses.values():
        vals.append(tc["rho_mean"].ravel())
    vmax = float(np.nanpercentile(np.concatenate(vals), 99.0)) if vals else 1.0
    fig, axes = plt.subplots(len(MODEL_ORDER), len(REGIMES), figsize=(13.5, 10.8), sharex=True, sharey=True)
    for ridx, model in enumerate(MODEL_ORDER):
        for cidx, regime in enumerate(REGIMES):
            ax = axes[ridx, cidx]
            tc = timecourses.get((model, regime))
            if tc is None:
                ax.text(0.5, 0.5, "n/a", transform=ax.transAxes, ha="center", va="center")
                continue
            t = tc["time_min"]
            z_edges = tc["z_edges"]
            rho = tc["rho_mean"]
            mids = 0.5 * (t[1:] + t[:-1])
            t_edges = np.concatenate(([t[0] - (mids[0] - t[0])], mids, [t[-1] + (t[-1] - mids[-1])]))
            im = ax.pcolormesh(t_edges, z_edges, rho.T, shading="auto", cmap="viridis", vmin=0, vmax=vmax)
            ax.set_ylim(z_edges[-1], z_edges[0])
            if ridx == 0:
                ax.set_title(regime, fontsize=14)
            if cidx == 0:
                ax.set_ylabel(MODEL_SPECS[model][0] + "\nz (um)", fontsize=11)
            if ridx == len(MODEL_ORDER) - 1:
                ax.set_xlabel("Time (min)", fontsize=12)
            ax.tick_params(labelsize=10)
    fig.subplots_adjust(right=0.90, wspace=0.12, hspace=0.18)
    cax = fig.add_axes((0.92, 0.15, 0.018, 0.70))
    fig.colorbar(im, cax=cax, label="Normalized axial density")
    fig.suptitle("Mean axial density kymographs, c80/m12", fontsize=17, y=0.98)
    out = PLOTS / "mean_axial_density_kymographs.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=260)
    plt.close(fig)
    return out


def main() -> None:
    if not REPORT.exists():
        raise FileNotFoundError(f"Missing rough-aware report binary: {REPORT}")
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)

    rows = metadata()
    write_csv(OUT / "input_runs.csv", rows)
    tasks = [(row, str(REPORT)) for row in rows]
    with ProcessPoolExecutor(max_workers=4) as pool:
        run_rows = list(pool.map(analyze_row, tasks))

    status_rows = [{k: v for k, v in row.items() if k not in {"time_min", "rho_z", "z_edges", "z_com", "vz_abs"}} for row in run_rows]
    write_csv(OUT / "run_status.csv", status_rows)

    thresholds = read_thresholds()
    run_metric_rows = compute_run_metrics(run_rows, thresholds)
    condition_rows = aggregate_condition_metrics(run_metric_rows)
    effects = effect_rows(condition_rows)
    write_csv(OUT / "run_metrics.csv", run_metric_rows)
    write_csv(OUT / "condition_metrics.csv", condition_rows)
    write_csv(OUT / "rough_vs_smooth_effects.csv", effects)

    timecourses = grouped_timecourses(run_rows)
    plots = [
        plot_metric_summary(condition_rows),
        plot_z_com_1to8(timecourses),
        plot_mean_kymographs(timecourses),
    ]
    print(OUT / "run_status.csv")
    print(OUT / "run_metrics.csv")
    print(OUT / "condition_metrics.csv")
    print(OUT / "rough_vs_smooth_effects.csv")
    for plot in plots:
        print(plot)


if __name__ == "__main__":
    main()
