#!/usr/bin/env python3
"""Publication-style 1:16 transport plots for random/aligned annulus runs.

This is a focused replacement for the older slide plots that varied
crosslinker ratio.  It keeps only the 1:16 crosslinker condition and compares:

- random rotatable
- random fixed
- aligned rotatable
- aligned fixed

Mean plots use the existing transport-regime CSVs.  The timecourse plot uses
speckle kymograph data to compute signed axial velocity over time.
"""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "analysis") not in sys.path:
    sys.path.insert(0, str(ROOT / "analysis"))

from analysis import plot_timecourse_panels as base
from analysis import run_speckle_kymograph_batch as kymo
from analysis import run_speckle_timecourse_batch as speck


OUT = ROOT / "analysis" / "results" / "clu_minusz_slides_2026-06-02" / "transport_1to16_publication"
MEAN_OUT = OUT / "mean_plots"
TIME_OUT = OUT / "timecourse_plots"
TABLE_OUT = OUT / "tables"

RANDOM_COND = ROOT / "analysis" / "results" / "fixed_global_story" / "transport_regime" / "condition_metrics.csv"
RANDOM_RUN = ROOT / "analysis" / "results" / "fixed_global_story" / "transport_regime" / "run_metrics.csv"
ALIGNED_COND = ROOT / "analysis" / "results" / "init_minusz_comparison" / "transport_regime" / "condition_metrics.csv"
ALIGNED_RUN = ROOT / "analysis" / "results" / "init_minusz_comparison" / "transport_regime" / "run_metrics.csv"

TIMECOURSE_MAPS = {
    "random_rotatable": (ROOT / "clu" / "xlink_regime_map.csv", "Random rotatable"),
    "random_fixed": (ROOT / "clu_fixed_global" / "xlink_regime_map.csv", "Random fixed"),
    "aligned_rotatable": (ROOT / "clu_init_minusz" / "xlink_regime_map.csv", "Aligned rotatable"),
    "aligned_fixed": (ROOT / "clu_fixed_global_init_minusz" / "xlink_regime_map.csv", "Aligned fixed"),
}

REGIME = "1:16"
TARGET_FINAL_MIN = 800.0 / 60.0
TARGET_FINAL_S = 800.0
SPECKLE_INTERVAL = 2.0
SPECKLE_NAME = "speckle_i2.txt"
NBINS_Z = 80

FAMILIES = [
    ("total480", "Total motors = 480", "Cluster count"),
    ("mpc12", "12 motors per cluster", "Cluster count"),
]
CASE_ORDER = {
    "total480": [
        ("c10_m48", 10),
        ("c20_m24", 20),
        ("c40_m12", 40),
        ("c60_m8", 60),
        ("c80_m6", 80),
    ],
    "mpc12": [
        ("c10_m12", 10),
        ("c20_m12", 20),
        ("c40_m12", 40),
        ("c60_m12", 60),
        ("c80_m12", 80),
    ],
}

SERIES = [
    ("random", "rotatable", "Random rotatable", "#222222", "s", "-", "white"),
    ("random", "fixed_global", "Random fixed", "#D55E00", "o", "-", "#D55E00"),
    ("aligned", "rotatable", "Aligned rotatable", "#009E73", "^", "--", "white"),
    ("aligned", "fixed_global", "Aligned fixed", "#7570B3", "D", "--", "#7570B3"),
]

TIMECOURSE_PANEL_ORDER = [
    ("random_rotatable", "Random rotatable"),
    ("aligned_rotatable", "Aligned rotatable"),
    ("random_fixed", "Random fixed"),
    ("aligned_fixed", "Aligned fixed"),
]

COUNT_COLORS = {
    10: "#0072B2",
    20: "#D55E00",
    40: "#009E73",
    60: "#CC79A7",
    80: "#000000",
}

MEAN_METRICS = [
    {
        "field": "mean_vz_signed",
        "ylabel": r"Mean axial velocity, $v_z$ ($\mu$m/s)",
        "slug": "mean_axial_velocity",
        "clip_zero": False,
        "zero_line": True,
    },
    {
        "field": "active_transport_fraction",
        "ylabel": "Active transport fraction",
        "slug": "active_transport_fraction",
        "clip_zero": True,
    },
    {
        "field": "directionality_index",
        "ylabel": r"Directionality index, $|\Delta z|/\mathrm{AUC}(|v_z|)$",
        "slug": "directionality_index",
        "clip_zero": True,
        "force_upper": 1.05,
    },
    {
        "field": "abs_net_z_displacement",
        "ylabel": r"$|\Delta z|$ ($\mu$m)",
        "slug": "net_z_shift",
        "clip_zero": True,
    },
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def finite_float(value: object) -> float:
    if value in (None, ""):
        return float("nan")
    try:
        return float(value)
    except Exception:
        return float("nan")


def sem(values: list[float]) -> float:
    arr = np.asarray([value for value in values if np.isfinite(value)], dtype=float)
    if arr.size == 0:
        return float("nan")
    if arr.size == 1:
        return 0.0
    return float(np.nanstd(arr, ddof=1) / math.sqrt(arr.size))


def condition_only_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (row["model"], row["group"], row["case"], row["xlink_regime"])


def enrich_condition_rows(cond_path: Path, run_path: Path, orientation: str) -> list[dict[str, object]]:
    cond_rows = [dict(row) for row in read_csv(cond_path) if row.get("xlink_regime") == REGIME]
    run_rows = [dict(row) for row in read_csv(run_path) if row.get("xlink_regime") == REGIME]

    grouped_runs: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in run_rows:
        grouped_runs[condition_only_key(row)].append(row)

    enriched: list[dict[str, object]] = []
    for row in cond_rows:
        row = dict(row)
        row["orientation"] = orientation

        net_mean = finite_float(row.get("net_z_displacement_mean"))
        net_sem = finite_float(row.get("net_z_displacement_sem"))
        row["mean_vz_signed_mean"] = net_mean / TARGET_FINAL_S if np.isfinite(net_mean) else float("nan")
        row["mean_vz_signed_sem"] = net_sem / TARGET_FINAL_S if np.isfinite(net_sem) else float("nan")

        directionality_values: list[float] = []
        for run in grouped_runs.get(condition_only_key(row), []):
            auc = finite_float(run.get("auc_vz_abs"))
            disp = finite_float(run.get("abs_net_z_displacement"))
            if np.isfinite(auc) and auc > 0 and np.isfinite(disp):
                directionality_values.append(disp / auc)
        row["directionality_index_mean"] = float(np.nanmean(directionality_values)) if directionality_values else float("nan")
        row["directionality_index_sem"] = sem(directionality_values)
        enriched.append(row)
    return enriched


def mean_lookup(rows: list[dict[str, object]]) -> dict[tuple[str, str, str, str, str], dict[str, object]]:
    return {
        (str(row["orientation"]), str(row["model"]), str(row["group"]), str(row["case"]), str(row["xlink_regime"])): row
        for row in rows
    }


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.tick_params(direction="out", width=1.4, length=5.5, labelsize=14)
    ax.grid(axis="y", color="#D8D8D8", linewidth=0.8, alpha=0.75)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        0.02,
        0.96,
        label,
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        ha="left",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5, "alpha": 0.92},
        zorder=10,
    )


def metric_bounds(lookup: dict, metric: dict) -> tuple[float, float]:
    mean_key = f"{metric['field']}_mean"
    sem_key = f"{metric['field']}_sem"
    values: list[float] = []
    for family, _title, _xlabel in FAMILIES:
        for case, _x in CASE_ORDER[family]:
            for orientation, model, *_rest in SERIES:
                row = lookup.get((orientation, model, family, case, REGIME))
                if row is None:
                    continue
                mean = finite_float(row.get(mean_key))
                err = finite_float(row.get(sem_key))
                if np.isfinite(mean):
                    values.extend([mean - err if np.isfinite(err) else mean, mean + err if np.isfinite(err) else mean])
    if not values:
        return 0.0, 1.0
    lo = min(values)
    hi = max(values)
    pad = 0.11 * (hi - lo) if hi > lo else 0.10 * max(abs(hi), 1.0)
    lo -= pad
    hi += pad
    if metric.get("clip_zero"):
        lo = max(0.0, lo)
    if metric.get("zero_line"):
        lo = min(lo, -0.0007)
        hi = max(hi, 0.0007)
    if metric.get("force_upper") is not None:
        hi = min(max(hi, 1.0), float(metric["force_upper"]))
    return lo, hi


def draw_family_metric(ax: plt.Axes, lookup: dict, family: str, metric: dict) -> None:
    mean_key = f"{metric['field']}_mean"
    sem_key = f"{metric['field']}_sem"
    for orientation, model, _label, color, marker, linestyle, face in SERIES:
        xs: list[int] = []
        ys: list[float] = []
        es: list[float] = []
        for case, xpos in CASE_ORDER[family]:
            row = lookup.get((orientation, model, family, case, REGIME))
            if row is None:
                continue
            xs.append(xpos)
            ys.append(finite_float(row.get(mean_key)))
            es.append(finite_float(row.get(sem_key)))
        ax.errorbar(
            xs,
            ys,
            yerr=es,
            color=color,
            marker=marker,
            linestyle=linestyle,
            linewidth=2.8,
            markersize=8.5,
            markerfacecolor=face,
            markeredgewidth=1.8,
            capsize=4.5,
            elinewidth=1.5,
        )
    if metric.get("zero_line"):
        ax.axhline(0, color="#777777", linewidth=1.0, linestyle="--", zorder=0)
    ax.set_xticks([x for _case, x in CASE_ORDER[family]])
    ax.set_xticklabels([str(x) for _case, x in CASE_ORDER[family]])
    style_axis(ax)


def legend_handles() -> list[plt.Line2D]:
    return [
        plt.Line2D(
            [0],
            [0],
            color=color,
            marker=marker,
            linestyle=linestyle,
            linewidth=2.8,
            markersize=8.5,
            markerfacecolor=face,
            markeredgewidth=1.8,
            markeredgecolor=color,
            label=label,
        )
        for _orientation, _model, label, color, marker, linestyle, face in SERIES
    ]


def save_individual_mean_plot(lookup: dict, metric: dict) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.2), sharey=True)
    ylo, yhi = metric_bounds(lookup, metric)
    for idx, (family, title, xlabel) in enumerate(FAMILIES):
        ax = axes[idx]
        draw_family_metric(ax, lookup, family, metric)
        ax.set_ylim(ylo, yhi)
        ax.set_title(title, fontsize=18, pad=10)
        ax.set_xlabel(xlabel, fontsize=16)
        if idx == 0:
            ax.set_ylabel(metric["ylabel"], fontsize=17)
        add_panel_label(ax, chr(ord("A") + idx))
    fig.legend(handles=legend_handles(), loc="lower center", ncol=4, frameon=False, fontsize=13, bbox_to_anchor=(0.5, -0.03))
    fig.subplots_adjust(left=0.09, right=0.99, top=0.88, bottom=0.24, wspace=0.10)
    MEAN_OUT.mkdir(parents=True, exist_ok=True)
    out = MEAN_OUT / f"{metric['slug']}_1to16.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    return out


def save_combined_mean_plot(lookup: dict) -> Path:
    fig, axes = plt.subplots(len(MEAN_METRICS), 2, figsize=(13.0, 18.0), sharex="col")
    for ridx, metric in enumerate(MEAN_METRICS):
        ylo, yhi = metric_bounds(lookup, metric)
        for cidx, (family, title, xlabel) in enumerate(FAMILIES):
            ax = axes[ridx, cidx]
            draw_family_metric(ax, lookup, family, metric)
            ax.set_ylim(ylo, yhi)
            if ridx == 0:
                ax.set_title(title, fontsize=18, pad=10)
            if ridx == len(MEAN_METRICS) - 1:
                ax.set_xlabel(xlabel, fontsize=16)
            if cidx == 0:
                ax.set_ylabel(metric["ylabel"], fontsize=16)
            add_panel_label(ax, chr(ord("A") + ridx * 2 + cidx))
    fig.legend(handles=legend_handles(), loc="lower center", ncol=4, frameon=False, fontsize=13, bbox_to_anchor=(0.5, 0.006))
    fig.subplots_adjust(left=0.10, right=0.99, top=0.98, bottom=0.075, wspace=0.10, hspace=0.30)
    MEAN_OUT.mkdir(parents=True, exist_ok=True)
    out = MEAN_OUT / "mean_transport_metrics_1to16_combined.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    return out


def load_timecourse_metadata(path: Path, model_key: str, family: str) -> list[dict[str, str]]:
    rows = base.load_metadata(path, model_key)
    out = []
    allowed_cases = {case for case, _count in CASE_ORDER[family]}
    for row in rows:
        if row.get("xlink_regime") != REGIME:
            continue
        if row.get("group") != family:
            continue
        if row.get("case") not in allowed_cases:
            continue
        out.append(row)
    return out


def analyze_timecourse_task(task: tuple[dict[str, str], str]) -> dict[str, object]:
    row, report_bin = task
    out = kymo.analyze_run((row, SPECKLE_INTERVAL, SPECKLE_NAME, report_bin, NBINS_Z))
    out["panel_key"] = row["model"]
    if out.get("status") != "ok":
        return out

    time = np.asarray(out["time_min"], dtype=float)
    keep = time <= TARGET_FINAL_MIN + 1e-9
    if np.any(keep):
        time = time[keep]
        rho = np.asarray(out["rho_z"], dtype=float)[keep]
    else:
        rho = np.asarray(out["rho_z"], dtype=float)
    z_edges = np.asarray(out["z_edges"], dtype=float)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    z_com = rho @ z_centers / np.clip(rho.sum(axis=1), 1e-12, None)
    if time.size > 1:
        dt_s = float((time[1] - time[0]) * 60.0)
        vz_signed = np.gradient(z_com, dt_s)
    else:
        vz_signed = np.zeros_like(z_com)

    out["time_min"] = time.tolist()
    out["z_com"] = z_com.tolist()
    out["vz_signed"] = vz_signed.tolist()
    out.pop("rho_z", None)
    out.pop("z_edges", None)
    return out


def aggregate_timecourse(run_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in run_rows:
        if row.get("status") == "ok":
            grouped[(str(row["panel_key"]), str(row["case"]), str(row["xlink_regime"]))].append(row)

    rows: list[dict[str, object]] = []
    for (panel_key, case, regime), items in sorted(grouped.items()):
        n_frames = min(len(item["time_min"]) for item in items)
        time = np.asarray(items[0]["time_min"][:n_frames], dtype=float)
        stack = np.vstack([np.asarray(item["vz_signed"][:n_frames], dtype=float) for item in items])
        mean = np.nanmean(stack, axis=0)
        err = np.nanstd(stack, axis=0, ddof=1) / math.sqrt(stack.shape[0]) if stack.shape[0] > 1 else np.zeros_like(mean)
        for idx, t in enumerate(time):
            rows.append(
                {
                    "panel_key": panel_key,
                    "case": case,
                    "xlink_regime": regime,
                    "frame_index": idx + 1,
                    "time_min": float(t),
                    "n_runs_used": len(items),
                    "vz_signed_mean": float(mean[idx]),
                    "vz_signed_sem": float(err[idx]),
                }
            )
    return rows


def timecourse_lookup(rows: list[dict[str, object]]) -> dict[tuple[str, str], list[dict[str, object]]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["panel_key"]), str(row["case"]))].append(row)
    for key in grouped:
        grouped[key].sort(key=lambda row: int(row["frame_index"]))
    return grouped


def save_timecourse_plot(rows: list[dict[str, object]], family: str) -> Path:
    lookup = timecourse_lookup(rows)
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.2), sharex=True, sharey=True)
    axes = axes.ravel()
    values: list[float] = []
    for key_rows in lookup.values():
        for row in key_rows:
            mean = finite_float(row.get("vz_signed_mean"))
            err = finite_float(row.get("vz_signed_sem"))
            if np.isfinite(mean):
                values.extend([mean - err if np.isfinite(err) else mean, mean + err if np.isfinite(err) else mean])
    ylo = min(values) if values else -0.001
    yhi = max(values) if values else 0.001
    pad = 0.12 * (yhi - ylo) if yhi > ylo else 0.001
    ylo -= pad
    yhi += pad
    ylo = min(ylo, -0.0008)
    yhi = max(yhi, 0.0008)

    for ax, (panel_key, panel_label) in zip(axes, TIMECOURSE_PANEL_ORDER):
        style_axis(ax)
        ax.axhline(0, color="#777777", linewidth=1.0, linestyle="--", zorder=0)
        ax.set_title(panel_label, fontsize=17, pad=8)
        for case, count in CASE_ORDER[family]:
            series_rows = lookup.get((panel_key, case), [])
            if not series_rows:
                continue
            t = np.asarray([finite_float(row["time_min"]) for row in series_rows], dtype=float)
            y = np.asarray([finite_float(row["vz_signed_mean"]) for row in series_rows], dtype=float)
            e = np.asarray([finite_float(row["vz_signed_sem"]) for row in series_rows], dtype=float)
            color = COUNT_COLORS[count]
            ax.plot(t, y, color=color, linewidth=2.6, label=f"{count} clusters")
            ax.fill_between(t, y - e, y + e, color=color, alpha=0.14, linewidth=0)
        ax.set_xlim(0, TARGET_FINAL_MIN)
        ax.set_ylim(ylo, yhi)
    for idx, ax in enumerate(axes):
        add_panel_label(ax, chr(ord("A") + idx))
    axes[0].set_ylabel(r"Axial velocity, $v_z$ ($\mu$m/s)", fontsize=16)
    axes[2].set_ylabel(r"Axial velocity, $v_z$ ($\mu$m/s)", fontsize=16)
    axes[2].set_xlabel("Time after motor introduction (min)", fontsize=16)
    axes[3].set_xlabel("Time after motor introduction (min)", fontsize=16)

    handles = [
        plt.Line2D([0], [0], color=COUNT_COLORS[count], lw=2.8, label=f"{count} clusters")
        for _case, count in CASE_ORDER[family]
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=13, bbox_to_anchor=(0.5, 0.005))
    fig.subplots_adjust(left=0.10, right=0.99, top=0.93, bottom=0.13, wspace=0.10, hspace=0.30)

    TIME_OUT.mkdir(parents=True, exist_ok=True)
    out = TIME_OUT / f"{family}_mean_axial_velocity_timecourse_1to16.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    return out


def build_timecourse(family: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    cache = TABLE_OUT / f"{family}_timecourse_1to16.csv"
    status_cache = TABLE_OUT / f"{family}_timecourse_run_status_1to16.csv"
    if cache.exists() and status_cache.exists():
        return read_csv(cache), read_csv(status_cache)

    # The selected runs already contain speckle_i2.txt.  A real report binary
    # is only needed if a speckle file has to be regenerated.
    report_bin = "/bin/true"
    tasks = []
    for panel_key, (path, _label) in TIMECOURSE_MAPS.items():
        for row in load_timecourse_metadata(path, panel_key, family):
            tasks.append((row, report_bin))

    with ProcessPoolExecutor(max_workers=8) as executor:
        run_rows = list(executor.map(analyze_timecourse_task, tasks))

    status_rows = [
        {key: value for key, value in row.items() if key not in {"time_min", "z_com", "vz_signed"}}
        for row in run_rows
    ]
    agg_rows = aggregate_timecourse(run_rows)
    write_csv(status_cache, status_rows)
    write_csv(cache, agg_rows)
    return agg_rows, status_rows


def write_readme(paths: list[Path]) -> None:
    lines = [
        "# 1:16 Transport Plots",
        "",
        "Focused publication-style plots for the random/aligned annulus transport story.",
        "",
        "Inputs:",
        f"- Random condition metrics: `{RANDOM_COND.relative_to(ROOT)}`",
        f"- Random run metrics: `{RANDOM_RUN.relative_to(ROOT)}`",
        f"- Aligned condition metrics: `{ALIGNED_COND.relative_to(ROOT)}`",
        f"- Aligned run metrics: `{ALIGNED_RUN.relative_to(ROOT)}`",
        "",
        "Definitions:",
        "- Mean axial velocity is computed as mean signed net z displacement divided by the 800 s motorized analysis window.",
        "- Directionality is computed at replicate level as `abs_net_z_displacement / auc_vz_abs`, then averaged across replicates.",
        "- Error bars and shaded timecourse bands are SEM across replicates.",
        "- All plots use crosslinker ratio `1:16` only.",
        "",
        "Generated outputs:",
    ]
    for path in paths:
        lines.append(f"- `{path.relative_to(ROOT)}`")
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.linewidth": 1.4,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "font.family": "DejaVu Sans",
        }
    )

    OUT.mkdir(parents=True, exist_ok=True)
    TABLE_OUT.mkdir(parents=True, exist_ok=True)

    enriched_rows: list[dict[str, object]] = []
    enriched_rows.extend(enrich_condition_rows(RANDOM_COND, RANDOM_RUN, "random"))
    enriched_rows.extend(enrich_condition_rows(ALIGNED_COND, ALIGNED_RUN, "aligned"))
    write_csv(TABLE_OUT / "mean_transport_metrics_1to16_source.csv", enriched_rows)
    lookup = mean_lookup(enriched_rows)

    outputs: list[Path] = []
    outputs.append(save_combined_mean_plot(lookup))
    for metric in MEAN_METRICS:
        outputs.append(save_individual_mean_plot(lookup, metric))

    for family in ("mpc12", "total480"):
        rows, _status = build_timecourse(family)
        outputs.append(save_timecourse_plot(rows, family))

    write_readme(outputs)
    print(OUT)
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
