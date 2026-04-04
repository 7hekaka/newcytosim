#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_ROTATABLE = ROOT / "clu" / "analysis_full" / "condition_summary.csv"
DEFAULT_FIXED = ROOT / "clu_fixed_global" / "analysis_full" / "condition_summary.csv"
DEFAULT_OUTDIR = ROOT / "comparison_plots"

XLINK_ORDER = ["1:4", "1:8", "1:16"]

FAMILY_SPECS = [
    {
        "group": "total480",
        "kind": "numeric",
        "xlabel": "Cluster count (total motors = 480)",
        "order": [
            ("c10_m48", 10),
            ("c20_m24", 20),
            ("c40_m12", 40),
            ("c60_m8", 60),
            ("c80_m6", 80),
        ],
    },
    {
        "group": "mpc12",
        "kind": "numeric",
        "xlabel": "Cluster count (12 motors per cluster)",
        "order": [
            ("c10_m12", 10),
            ("c20_m12", 20),
            ("c40_m12", 40),
            ("c60_m12", 60),
            ("c80_m12", 80),
        ],
    },
    {
        "group": "focus40x12",
        "kind": "categorical",
        "xlabel": "Crosslinker condition (40 clusters, 12 motors/cluster)",
        "order": [
            ("parallel_only", "parallel"),
            ("parallel_antiparallel_50_50", "50:50"),
            ("offrate_0p5", "off-rate 0.5"),
            ("offrate_1p0", "off-rate 1.0"),
        ],
    },
]

MODEL_SPECS = [
    {
        "key": "rotatable",
        "label": "Rotatable",
        "color": "#222222",
        "marker": "s",
        "markerfacecolor": "white",
    },
    {
        "key": "fixed_global",
        "label": "Fixed global",
        "color": "#d95f02",
        "marker": "o",
        "markerfacecolor": "#d95f02",
    },
]

METRICS = [
    {
        "field": "vz_abs_mean",
        "ylabel": "Mean |v_z| (um/s)",
        "slug": "vz_abs_mean",
        "clip_zero": True,
    },
    {
        "field": "Dz_mean",
        "ylabel": "Mean depletion index",
        "slug": "depletion_mean",
        "clip_zero": True,
    },
    {
        "field": "shell_occ_mean",
        "ylabel": "Mean occupied shell fraction",
        "slug": "shell_occupancy_mean",
        "clip_zero": True,
    },
    {
        "field": "nematic_xy_mean",
        "ylabel": "Mean XY nematic order",
        "slug": "nematic_xy_mean",
        "clip_zero": True,
    },
    {
        "field": "top_bias_mean",
        "ylabel": "Mean top-bottom bias",
        "slug": "top_bias_mean",
        "clip_zero": False,
    },
    {
        "field": "swirl_rate_abs_mean",
        "ylabel": "Mean |swirl rate| (rad/s)",
        "slug": "swirl_rate_abs_mean",
        "clip_zero": True,
    },
]


def load_rows(path: Path, model: str) -> list[dict]:
    with path.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        row["model"] = model
    return rows


def as_float(row: dict | None, key: str) -> float:
    if not row:
        return float("nan")
    value = row.get(key, "")
    return float(value) if value not in ("", None) else float("nan")


def build_lookup(rows: list[dict]) -> dict[tuple[str, str, str, str], dict]:
    out = {}
    for row in rows:
        out[(row["model"], row["group"], row["case"], row["xlink_regime"])] = row
    return out


def style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.1)
    ax.spines["bottom"].set_linewidth(1.1)
    ax.tick_params(direction="out", width=1.1, labelsize=10)
    ax.grid(axis="y", color="#d7d7d7", linewidth=0.6, alpha=0.8)


def compute_ylim(lookup: dict, metric: dict) -> tuple[float, float]:
    mean_key = f'{metric["field"]}_mean'
    sem_key = f'{metric["field"]}_sem'
    vals = []
    for row in lookup.values():
        mean = as_float(row, mean_key)
        sem = as_float(row, sem_key)
        if np.isfinite(mean):
            if np.isfinite(sem):
                vals.extend([mean - sem, mean + sem])
            else:
                vals.append(mean)
    if not vals:
        return 0.0, 1.0
    lo = min(vals)
    hi = max(vals)
    pad = 0.1 * max(abs(hi), 1.0) if hi <= lo else 0.08 * (hi - lo)
    lo -= pad
    hi += pad
    if metric["clip_zero"]:
        lo = max(0.0, lo)
    return lo, hi


def series_for_family(
    lookup: dict,
    *,
    model: str,
    group: str,
    regime: str,
    family_spec: dict,
    metric: dict,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    mean_key = f'{metric["field"]}_mean'
    sem_key = f'{metric["field"]}_sem'

    if family_spec["kind"] == "numeric":
        xs = np.array([x for _, x in family_spec["order"]], dtype=float)
        labels = [str(x) for _, x in family_spec["order"]]
    else:
        xs = np.arange(len(family_spec["order"]), dtype=float)
        labels = [label for _, label in family_spec["order"]]

    ys = []
    es = []
    for case, _ in family_spec["order"]:
        row = lookup.get((model, group, case, regime))
        ys.append(as_float(row, mean_key))
        es.append(as_float(row, sem_key))

    return xs, np.array(ys, dtype=float), np.array(es, dtype=float), labels


def plot_model_series(ax, xs: np.ndarray, ys: np.ndarray, es: np.ndarray, model_spec: dict) -> None:
    ax.errorbar(
        xs,
        ys,
        yerr=es,
        color=model_spec["color"],
        lw=1.7,
        elinewidth=1.1,
        capsize=2.8,
        marker=model_spec["marker"],
        markersize=5.4,
        markeredgewidth=1.2,
        markeredgecolor=model_spec["color"],
        markerfacecolor=model_spec["markerfacecolor"],
        linestyle="-",
        zorder=3,
    )


def save_figure(fig, outbase: Path) -> None:
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")


def make_metric_figure(lookup: dict, metric: dict, outdir: Path) -> list[Path]:
    fig, axes = plt.subplots(3, 3, figsize=(12.2, 10.2), sharey=True)
    ylo, yhi = compute_ylim(lookup, metric)

    for ridx, regime in enumerate(XLINK_ORDER):
        for cidx, family_spec in enumerate(FAMILY_SPECS):
            ax = axes[ridx, cidx]
            style_axes(ax)

            for model_spec in MODEL_SPECS:
                xs, ys, es, labels = series_for_family(
                    lookup,
                    model=model_spec["key"],
                    group=family_spec["group"],
                    regime=regime,
                    family_spec=family_spec,
                    metric=metric,
                )
                plot_model_series(ax, xs, ys, es, model_spec)

            ax.set_ylim(ylo, yhi)

            if family_spec["kind"] == "numeric":
                xticks = [x for _, x in family_spec["order"]]
                ax.set_xticks(xticks)
                ax.set_xticklabels([str(x) for x in xticks])
            else:
                xticks = np.arange(len(family_spec["order"]))
                ax.set_xticks(xticks)
                ax.set_xticklabels(labels, rotation=18, ha="right")

            if cidx == 0:
                ax.set_ylabel(metric["ylabel"], fontsize=11)
            if ridx == 2:
                ax.set_xlabel(family_spec["xlabel"], fontsize=11)

            ax.text(
                0.03,
                0.94,
                regime,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=10,
                color="#444444",
            )

    handles = []
    labels = []
    for model_spec in MODEL_SPECS:
        handle = plt.Line2D(
            [0],
            [0],
            color=model_spec["color"],
            lw=1.7,
            marker=model_spec["marker"],
            markersize=5.4,
            markeredgewidth=1.2,
            markeredgecolor=model_spec["color"],
            markerfacecolor=model_spec["markerfacecolor"],
        )
        handles.append(handle)
        labels.append(model_spec["label"])

    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, -0.01),
        fontsize=11,
    )
    fig.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.12, wspace=0.20, hspace=0.22)

    outbase = outdir / metric["slug"]
    save_figure(fig, outbase)
    plt.close(fig)
    return [outbase.with_suffix(".png"), outbase.with_suffix(".pdf")]


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Plot paper-style fixed-global vs rotatable comparison panels with one metric per figure."
    )
    ap.add_argument("--rotatable", type=Path, default=DEFAULT_ROTATABLE)
    ap.add_argument("--fixed", type=Path, default=DEFAULT_FIXED)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = ap.parse_args()

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.linewidth": 1.1,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )

    rot_rows = load_rows(args.rotatable, "rotatable")
    fixed_rows = load_rows(args.fixed, "fixed_global")
    lookup = build_lookup(rot_rows + fixed_rows)

    args.outdir.mkdir(parents=True, exist_ok=True)
    outfiles = []
    for metric in METRICS:
        outfiles.extend(make_metric_figure(lookup, metric, args.outdir))

    for path in outfiles:
        print(path)


if __name__ == "__main__":
    main()
