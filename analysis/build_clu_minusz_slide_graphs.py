#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
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

from analysis import comparison_unwrapped_annulus_story as story
from analysis import plot_timecourse_panels as base


OUT = ROOT / "analysis" / "results" / "clu_minusz_slides_2026-06-02"
TABLES = OUT / "tables"
METRIC_OUT = OUT / "transport_metrics"
FOUR_WAY_OUT = OUT / "four_way_comparison"
HEATMAP_OUT = OUT / "unwrapped_heatmaps"
LINE_OUT = OUT / "unwrapped_lines"

ALIGNED_COND = ROOT / "analysis" / "results" / "init_minusz_comparison" / "transport_regime" / "condition_metrics.csv"
ALIGNED_RUN = ROOT / "analysis" / "results" / "init_minusz_comparison" / "transport_regime" / "run_metrics.csv"
RANDOM_COND = ROOT / "analysis" / "results" / "fixed_global_story" / "transport_regime" / "condition_metrics.csv"
RANDOM_RUN = ROOT / "analysis" / "results" / "fixed_global_story" / "transport_regime" / "run_metrics.csv"

XLINK_ORDER = ["1:4", "1:8", "1:16"]
XLINK_SLUG = {"1:4": "1to4", "1:8": "1to8", "1:16": "1to16"}
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
TARGET_FINAL_S = 800.0

ALIGNED_SERIES = [
    ("aligned", "rotatable", "Rotatable", "#111111", "s", "-", "white"),
    ("aligned", "fixed_global", "Fixed orientation", "#d95f02", "o", "-", "#d95f02"),
]
FOUR_WAY_SERIES = [
    ("random", "rotatable", "Random rotatable", "#111111", "s", "-", "white"),
    ("random", "fixed_global", "Random fixed orientation", "#d95f02", "o", "-", "#d95f02"),
    ("aligned", "rotatable", "Aligned rotatable", "#1b9e77", "^", "--", "white"),
    ("aligned", "fixed_global", "Aligned fixed orientation", "#7570b3", "D", "--", "#7570b3"),
]
METRICS = [
    {
        "field": "mean_vz_abs",
        "ylabel": "Mean |v_z| (um/s)",
        "filename": "mean_vz_abs.png",
        "title": "Mean axial transport speed",
        "clip_zero": True,
    },
    {
        "field": "active_transport_fraction",
        "ylabel": "Active transport fraction",
        "filename": "active_transport_fraction.png",
        "title": "Active transport fraction",
        "clip_zero": True,
    },
    {
        "field": "directionality_index",
        "ylabel": "Directionality index",
        "filename": "directionality_index.png",
        "title": "Coherent axial directionality |Delta z| / AUC(|v_z|)",
        "clip_zero": True,
        "force_upper": 1.05,
    },
    {
        "field": "auc_vz_abs",
        "ylabel": "AUC(|v_z|) (um)",
        "filename": "auc_vz_abs.png",
        "title": "Total axial motion over motor-on window",
        "clip_zero": True,
    },
    {
        "field": "abs_net_z_displacement",
        "ylabel": "|Net z shift| (um)",
        "filename": "abs_net_z_displacement.png",
        "title": "Net axial displacement",
        "clip_zero": True,
    },
    {
        "field": "mean_vz_signed",
        "ylabel": "Mean signed v_z (um/s)",
        "filename": "mean_vz_signed.png",
        "title": "Signed axial drift",
        "clip_zero": False,
        "zero_line": True,
    },
]
PRIMARY_FOUR_WAY = {"mean_vz_abs", "active_transport_fraction", "directionality_index", "auc_vz_abs"}


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
    arr = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    if arr.size <= 1:
        return 0.0 if arr.size == 1 else float("nan")
    return float(np.nanstd(arr, ddof=1) / math.sqrt(arr.size))


def condition_key(row: dict[str, str], orientation: str) -> tuple[str, str, str, str, str]:
    return (orientation, row["model"], row["group"], row["case"], row["xlink_regime"])


def condition_only_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (row["model"], row["group"], row["case"], row["xlink_regime"])


def enrich_condition_rows(cond_path: Path, run_path: Path, orientation: str) -> list[dict[str, object]]:
    cond_rows = [dict(row) for row in read_csv(cond_path)]
    run_rows = [dict(row) for row in read_csv(run_path)]

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
        if directionality_values:
            row["directionality_index_mean"] = float(np.nanmean(directionality_values))
            row["directionality_index_sem"] = sem(directionality_values)
        else:
            row["directionality_index_mean"] = float("nan")
            row["directionality_index_sem"] = float("nan")
        enriched.append(row)
    return enriched


def build_lookup(rows: list[dict[str, object]]) -> dict[tuple[str, str, str, str, str], dict[str, object]]:
    lookup: dict[tuple[str, str, str, str, str], dict[str, object]] = {}
    for row in rows:
        lookup[(str(row["orientation"]), str(row["model"]), str(row["group"]), str(row["case"]), str(row["xlink_regime"]))] = row
    return lookup


def metric_bounds(lookup: dict, metric: dict, series: list[tuple], include_control: bool) -> tuple[float, float]:
    mean_key = f"{metric['field']}_mean"
    sem_key = f"{metric['field']}_sem"
    values: list[float] = []
    for regime in XLINK_ORDER:
        if include_control:
            row = lookup.get(("aligned", "control", "control", "c0_m0", regime))
            if row:
                m = finite_float(row.get(mean_key))
                s = finite_float(row.get(sem_key))
                if np.isfinite(m):
                    values.extend([m - s, m + s])
        for family, _title, _xlabel in FAMILIES:
            for case, _x in CASE_ORDER[family]:
                for orientation, model, *_rest in series:
                    row = lookup.get((orientation, model, family, case, regime))
                    if row is None:
                        continue
                    m = finite_float(row.get(mean_key))
                    s = finite_float(row.get(sem_key))
                    if np.isfinite(m):
                        if np.isfinite(s):
                            values.extend([m - s, m + s])
                        else:
                            values.append(m)
    if not values:
        return 0.0, 1.0
    lo = min(values)
    hi = max(values)
    pad = 0.10 * (hi - lo) if hi > lo else 0.10 * max(abs(hi), 1.0)
    lo -= pad
    hi += pad
    if metric.get("clip_zero"):
        lo = max(0.0, lo)
    if metric.get("force_upper") is not None:
        hi = min(max(hi, 1.0), float(metric["force_upper"]))
    return lo, hi


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.tick_params(direction="out", width=1.3, labelsize=15, length=5)
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.8, alpha=0.9)


def draw_aligned_metric(lookup: dict, metric: dict) -> Path:
    plt.rcParams.update(
        {
            "font.size": 16,
            "axes.linewidth": 1.4,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )
    fig, axes = plt.subplots(3, 2, figsize=(16.0, 14.0), sharey=True)
    mean_key = f"{metric['field']}_mean"
    sem_key = f"{metric['field']}_sem"
    ylo, yhi = metric_bounds(lookup, metric, ALIGNED_SERIES, include_control=False)

    for ridx, regime in enumerate(XLINK_ORDER):
        for cidx, (family, family_title, xlabel) in enumerate(FAMILIES):
            ax = axes[ridx, cidx]
            for orientation, model, _label, color, marker, linestyle, face in ALIGNED_SERIES:
                xs: list[int] = []
                ys: list[float] = []
                es: list[float] = []
                for case, xpos in CASE_ORDER[family]:
                    row = lookup.get((orientation, model, family, case, regime))
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
                    linewidth=3.0,
                    markersize=9,
                    markerfacecolor=face,
                    markeredgewidth=1.8,
                    capsize=4,
                )
            if metric.get("zero_line"):
                ax.axhline(0.0, color="#8f8f8f", linewidth=1.0, linestyle="--", zorder=0)
            ax.text(0.03, 0.93, regime, transform=ax.transAxes, ha="left", va="top", fontsize=18, color="#333333")
            ax.set_ylim(ylo, yhi)
            ax.set_xticks([x for _case, x in CASE_ORDER[family]])
            ax.set_xticklabels([str(x) for _case, x in CASE_ORDER[family]])
            if ridx == 0:
                ax.set_title(family_title, fontsize=22, pad=12)
            if ridx == 2:
                ax.set_xlabel(xlabel, fontsize=19)
            style_axis(ax)

    handles = [
        plt.Line2D([0], [0], color=spec[3], lw=3.0, marker=spec[4],
                   markersize=9, markeredgewidth=1.8, markeredgecolor=spec[3], markerfacecolor=spec[6],
                   linestyle=spec[5], label=spec[2])
        for spec in ALIGNED_SERIES
    ]
    fig.suptitle(f"Aligned minus-z filaments: {metric['title']}", fontsize=26, y=0.985)
    fig.text(0.035, 0.53, metric["ylabel"], ha="center", va="center", rotation="vertical", fontsize=22)
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.018), fontsize=17)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.93, bottom=0.13, wspace=0.14, hspace=0.24)
    out = METRIC_OUT / metric["filename"]
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def draw_four_way_metric(lookup: dict, metric: dict) -> Path:
    plt.rcParams.update(
        {
            "font.size": 16,
            "axes.linewidth": 1.4,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )
    fig, axes = plt.subplots(3, 2, figsize=(16.0, 14.0), sharey=True)
    mean_key = f"{metric['field']}_mean"
    sem_key = f"{metric['field']}_sem"
    ylo, yhi = metric_bounds(lookup, metric, FOUR_WAY_SERIES, include_control=False)

    for ridx, regime in enumerate(XLINK_ORDER):
        for cidx, (family, family_title, xlabel) in enumerate(FAMILIES):
            ax = axes[ridx, cidx]
            for orientation, model, _label, color, marker, linestyle, face in FOUR_WAY_SERIES:
                xs: list[int] = []
                ys: list[float] = []
                es: list[float] = []
                for case, xpos in CASE_ORDER[family]:
                    row = lookup.get((orientation, model, family, case, regime))
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
                    linewidth=3.0,
                    markersize=9,
                    markerfacecolor=face,
                    markeredgewidth=1.8,
                    capsize=4,
                )
            if metric.get("zero_line"):
                ax.axhline(0.0, color="#8f8f8f", linewidth=1.0, linestyle="--", zorder=0)
            ax.text(0.03, 0.93, regime, transform=ax.transAxes, ha="left", va="top", fontsize=18, color="#333333")
            ax.set_ylim(ylo, yhi)
            ax.set_xticks([x for _case, x in CASE_ORDER[family]])
            ax.set_xticklabels([str(x) for _case, x in CASE_ORDER[family]])
            if ridx == 0:
                ax.set_title(family_title, fontsize=22, pad=12)
            if ridx == 2:
                ax.set_xlabel(xlabel, fontsize=19)
            style_axis(ax)

    handles = [
        plt.Line2D([0], [0], color=spec[3], lw=3.0, marker=spec[4], markersize=9,
                   markeredgewidth=1.8, markeredgecolor=spec[3], markerfacecolor=spec[6],
                   linestyle=spec[5], label=spec[2])
        for spec in FOUR_WAY_SERIES
    ]
    fig.suptitle(f"Random vs aligned filaments: {metric['title']}", fontsize=26, y=0.985)
    fig.text(0.035, 0.53, metric["ylabel"], ha="center", va="center", rotation="vertical", fontsize=22)
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.013), fontsize=16)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.93, bottom=0.15, wspace=0.14, hspace=0.24)
    out = FOUR_WAY_OUT / metric["filename"]
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def configure_story_module() -> None:
    story.ROT_MAP = ROOT / "clu_init_minusz" / "xlink_regime_map.csv"
    story.FIX_MAP = ROOT / "clu_fixed_global_init_minusz" / "xlink_regime_map.csv"
    story.MEDOID = ALIGNED_COND.parent / "medoid_selection_summary.csv"
    story.OUT = OUT / "unwrapped_annulus_source"
    story.FAMILIES = {"total480", "mpc12"}
    # The representative runs already have speckle_i2.txt files. The shared
    # renderer asks for a report binary before checking those files; this keeps
    # the builder usable on machines where Cytosim binaries are not built.
    story.speck.pick_report_bin = lambda: Path("/bin/true")


def family_specs() -> list[dict]:
    return [spec for spec in base.FAMILY_SPECS if spec["group"] in {"total480", "mpc12"}]


def split_cases(family_spec: dict, part: str) -> list[tuple[str, str]]:
    order = list(family_spec["order"])
    if part == "low_counts":
        return order[:3]
    return order[3:]


def draw_snapshot_axis(ax: plt.Axes, panel: story.Panel | None, *, mode: str, vmax: float | None) -> object | None:
    if mode == "line":
        story.draw_line_panel(ax, panel)
        for line in ax.lines:
            line.set_linewidth(0.65)
            line.set_alpha(0.45)
        return None
    im = story.draw_density_panel(ax, panel, vmax=float(vmax))
    return im


def draw_split_timeline(
    family_spec: dict,
    regime: str,
    panels: list[story.Panel],
    *,
    mode: str,
    vmax: float | None,
    part: str,
) -> tuple[plt.Figure, object | None]:
    chosen = split_cases(family_spec, part)
    rows = []
    for case, label in chosen:
        rows.append({"model": "rotatable", "case": case, "label": f"{label}\nrotatable"})
        rows.append({"model": "fixed_global", "case": case, "label": f"{label}\nfixed orientation"})

    lookup = story.panel_lookup(panels)
    nrows = len(rows)
    ncols = len(story.SNAPSHOTS)
    fig, axes = plt.subplots(nrows, ncols, figsize=(16.0, 2.15 * nrows), sharex=True, sharey=True)
    if nrows == 1:
        axes = np.array([axes])

    im = None
    for ridx, row in enumerate(rows):
        for cidx, (snapshot_slug, _fraction, snapshot_label) in enumerate(story.SNAPSHOTS):
            ax = axes[ridx, cidx]
            panel = lookup.get((regime, row["model"], row["case"], snapshot_slug))
            im = draw_snapshot_axis(ax, panel, mode=mode, vmax=vmax) or im
            ax.tick_params(labelsize=11, width=1.0, length=4)
            if ridx == 0:
                title = snapshot_label
                if panel is not None:
                    title = f"{snapshot_label}\n{panel.time_min:.1f} min"
                ax.set_title(title, fontsize=16, pad=8)
            if cidx == 0:
                ax.set_ylabel("z (um)", fontsize=15)
                text_color = "white" if mode == "density" else "black"
                ax.text(
                    0.015,
                    0.93,
                    row["label"],
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=12,
                    color=text_color,
                    weight="bold",
                )
            else:
                ax.set_ylabel("")
                ax.set_yticklabels([])
            if ridx == nrows - 1:
                ax.set_xlabel("Unwrapped circumference (um)", fontsize=14)
            else:
                ax.set_xticklabels([])

    count_label = "10-40 clusters" if part == "low_counts" else "60-80 clusters"
    title_family = "Total motors = 480" if family_spec["group"] == "total480" else "12 motors per cluster"
    fig.suptitle(f"{title_family} | xlink ratio {regime} | {count_label}", fontsize=23, y=0.99)
    right = 0.96 if mode == "line" else 0.91
    fig.subplots_adjust(left=0.10, right=right, bottom=0.07, top=0.91, wspace=0.06, hspace=0.14)
    return fig, im


def build_unwrapped_figures() -> list[Path]:
    configure_story_module()
    HEATMAP_OUT.mkdir(parents=True, exist_ok=True)
    LINE_OUT.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.size": 13,
            "axes.linewidth": 1.0,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )

    family_panels, family_vmax = story.collect_panels()
    outputs: list[Path] = []

    for family_spec in family_specs():
        family = family_spec["group"]
        panels = family_panels[family]
        vmax = family_vmax[family]
        for regime in XLINK_ORDER:
            for part in ("low_counts", "high_counts"):
                line_fig, _ = draw_split_timeline(family_spec, regime, panels, mode="line", vmax=None, part=part)
                line_path = LINE_OUT / f"{family}_{XLINK_SLUG[regime]}_{part}_timeline_lines.png"
                line_fig.savefig(line_path, dpi=260, bbox_inches="tight")
                plt.close(line_fig)
                outputs.append(line_path)

                heat_fig, im = draw_split_timeline(family_spec, regime, panels, mode="density", vmax=vmax, part=part)
                if im is not None:
                    cax = heat_fig.add_axes([0.925, 0.16, 0.014, 0.68])
                    cbar = heat_fig.colorbar(im, cax=cax)
                    cbar.set_label("Normalized areal density", fontsize=14)
                    cbar.ax.tick_params(labelsize=11, width=1.0)
                    cbar.outline.set_linewidth(1.0)
                heat_path = HEATMAP_OUT / f"{family}_{XLINK_SLUG[regime]}_{part}_timeline_heatmap.png"
                heat_fig.savefig(heat_path, dpi=260, bbox_inches="tight")
                plt.close(heat_fig)
                outputs.append(heat_path)
    return outputs


def write_readme(paths: list[Path]) -> Path:
    readme = OUT / "README.md"
    lines = [
        "# CLU Minus-Z Slide Graphs",
        "",
        "Generated slide-ready assets for initially minus-z-aligned filament simulations.",
        "",
        "## Inputs",
        f"- Aligned metrics: `{ALIGNED_COND.relative_to(ROOT)}`",
        f"- Aligned run metrics: `{ALIGNED_RUN.relative_to(ROOT)}`",
        f"- Random-orientation metrics: `{RANDOM_COND.relative_to(ROOT)}`",
        f"- Random-orientation run metrics: `{RANDOM_RUN.relative_to(ROOT)}`",
        "",
        "## Metric Notes",
        "- Error bars are SEM across runs.",
        "- Directionality index is computed at run level as `abs_net_z_displacement / auc_vz_abs` and then averaged.",
        "- Signed mean velocity is computed as `net_z_displacement / 800 s`.",
        "- No-motor controls are omitted from the initialized-minus-z slide figures because the available control runs use random initial filament orientations.",
        "",
        "## Generated Assets",
    ]
    for path in paths:
        lines.append(f"- `{path.relative_to(ROOT)}`")
    readme.parent.mkdir(parents=True, exist_ok=True)
    readme.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return readme


def main() -> None:
    aligned_rows = enrich_condition_rows(ALIGNED_COND, ALIGNED_RUN, "aligned")
    random_rows = enrich_condition_rows(RANDOM_COND, RANDOM_RUN, "random")
    all_rows = aligned_rows + random_rows
    lookup = build_lookup(all_rows)

    TABLES.mkdir(parents=True, exist_ok=True)
    write_csv(TABLES / "aligned_condition_metrics_enriched.csv", aligned_rows)
    write_csv(TABLES / "four_way_condition_metrics_enriched.csv", all_rows)

    outputs: list[Path] = []
    for metric in METRICS:
        outputs.append(draw_aligned_metric(lookup, metric))
    for metric in METRICS:
        if metric["field"] in PRIMARY_FOUR_WAY:
            outputs.append(draw_four_way_metric(lookup, metric))
    outputs.extend(build_unwrapped_figures())
    outputs.append(TABLES / "aligned_condition_metrics_enriched.csv")
    outputs.append(TABLES / "four_way_condition_metrics_enriched.csv")
    outputs.append(write_readme(outputs))

    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
