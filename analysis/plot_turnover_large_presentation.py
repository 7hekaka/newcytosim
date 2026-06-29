#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "analysis" / "results"
OUT = RESULTS / "turnover_continuous_supply_scaled_large_presentation"

INPUTS = [
    (
        "baseline",
        "Baseline chewing",
        RESULTS / "turnover_continuous_supply_scaled_job09",
    ),
    (
        "ultra",
        "2x chewers + 2x rate",
        RESULTS / "turnover_continuous_supply_scaled_job12_chewer_count2x_rate2x",
    ),
]

CONDITION_ORDER = ["nomotor_xlink", "rotatable_xlink"]
CONDITION_LABELS = {
    "nomotor_xlink": "No motors",
    "rotatable_xlink": "Motors",
}
CONDITION_COLORS = {
    "nomotor_xlink": "#4d4d4d",
    "rotatable_xlink": "#0a9f7a",
}

ENDPOINT_VELOCITY_METRICS = [
    ("mean_vz_abs", "Mean |v_z| (um/s)"),
    ("auc_vz_abs", "AUC |v_z| (um)"),
    ("peak_vz_abs", "Peak |v_z| (um/s)"),
    ("active_transport_fraction", "Active transport fraction"),
    ("matched_active_transport_fraction", "Matched active fraction"),
    ("mean_vz_signed", "Mean v_z (um/s)"),
]

ENDPOINT_TURNOVER_METRICS = [
    ("chewer_band_mass_fraction_final", "Final chewer-band fraction"),
    ("bottom_half_mass_fraction_final", "Final bottom-half fraction"),
    ("mass_retention", "Mass retention"),
    ("speckle_length_proxy_final_um", "Final filament proxy (um)"),
]

STORY_METRICS = [
    ("auc_vz_abs", "AUC |v_z| (um)"),
    ("active_transport_fraction", "Active transport fraction"),
    ("chewer_band_mass_fraction_final", "Final chewer-band fraction"),
    ("bottom_half_mass_fraction_final", "Final bottom-half fraction"),
]


def f(value: object, default: float = float("nan")) -> float:
    try:
        if value in (None, ""):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def fmt(value: float, digits: int = 3) -> str:
    if not math.isfinite(value):
        return "n/a"
    return f"{value:.{digits}f}"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
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


def load_metrics() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for variant, variant_label, folder in INPUTS:
        for row in read_csv(folder / "condition_metrics.csv"):
            if row.get("scenario") != "large_long":
                continue
            if row.get("condition") not in CONDITION_ORDER:
                continue
            out: dict[str, object] = {
                "variant": variant,
                "variant_label": variant_label,
                **row,
            }
            rows.append(out)
    return rows


def load_timecourses() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for variant, variant_label, folder in INPUTS:
        for row in read_csv(folder / "condition_timecourses.csv"):
            if row.get("scenario") != "large_long":
                continue
            if row.get("condition") not in CONDITION_ORDER:
                continue
            out: dict[str, object] = {
                "variant": variant,
                "variant_label": variant_label,
                **row,
            }
            rows.append(out)
    return rows


def metric_row_index(rows: list[dict[str, object]]) -> dict[tuple[str, str], dict[str, object]]:
    return {
        (str(row["variant"]), str(row["condition"])): row
        for row in rows
    }


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.labelsize": 14,
            "axes.titlesize": 14,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 11,
            "axes.linewidth": 1.1,
            "figure.dpi": 160,
            "savefig.dpi": 300,
            "svg.fonttype": "none",
        }
    )


def save_figure(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.png", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        0.01,
        0.98,
        label,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=15,
        fontweight="bold",
    )


def bar_positions() -> list[tuple[str, str, str]]:
    return [
        ("baseline", "nomotor_xlink", "Baseline\nno motors"),
        ("baseline", "rotatable_xlink", "Baseline\nmotors"),
        ("ultra", "nomotor_xlink", "Ultra\nno motors"),
        ("ultra", "rotatable_xlink", "Ultra\nmotors"),
    ]


def plot_bar_metric(
    ax: plt.Axes,
    indexed: dict[tuple[str, str], dict[str, object]],
    metric: str,
    ylabel: str,
) -> None:
    positions = bar_positions()
    x = np.arange(len(positions))
    values = []
    errors = []
    colors = []
    hatches = []
    for variant, condition, _label in positions:
        row = indexed[(variant, condition)]
        values.append(f(row.get(metric)))
        errors.append(f(row.get(f"{metric}_sem"), 0.0))
        colors.append(CONDITION_COLORS[condition])
        hatches.append("" if variant == "baseline" else "///")
    bars = ax.bar(
        x,
        values,
        yerr=errors,
        color=colors,
        edgecolor="#222222",
        linewidth=0.9,
        capsize=3.0,
        width=0.68,
    )
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x, [label for *_rest, label in positions])
    ax.tick_params(axis="x", rotation=0)
    ax.grid(axis="y", color="#d8d8d8", linewidth=0.8, alpha=0.8)
    ax.set_axisbelow(True)
    ax.margins(x=0.04)


def plot_endpoint_panels(rows: list[dict[str, object]]) -> None:
    indexed = metric_row_index(rows)

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4))
    for ax, (metric, ylabel), label in zip(
        axes.flat,
        ENDPOINT_VELOCITY_METRICS,
        ["A", "B", "C", "D", "E", "F"],
    ):
        plot_bar_metric(ax, indexed, metric, ylabel)
        panel_label(ax, label)
    fig.tight_layout(w_pad=2.0, h_pad=2.0)
    save_figure(fig, "large_velocity_endpoint_metrics")

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4))
    for ax, (metric, ylabel), label in zip(
        axes.flat,
        ENDPOINT_TURNOVER_METRICS,
        ["A", "B", "C", "D"],
    ):
        plot_bar_metric(ax, indexed, metric, ylabel)
        panel_label(ax, label)
    fig.tight_layout(w_pad=2.0, h_pad=2.0)
    save_figure(fig, "large_turnover_endpoint_metrics")

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4))
    for ax, (metric, ylabel), label in zip(axes.flat, STORY_METRICS, ["A", "B", "C", "D"]):
        plot_bar_metric(ax, indexed, metric, ylabel)
        panel_label(ax, label)
    fig.tight_layout(w_pad=2.0, h_pad=2.0)
    save_figure(fig, "large_story_summary")


def group_timecourses(rows: list[dict[str, object]]) -> dict[tuple[str, str], list[dict[str, object]]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in rows:
        key = (str(row["variant"]), str(row["condition"]))
        grouped.setdefault(key, []).append(row)
    for key in grouped:
        grouped[key].sort(key=lambda row: f(row.get("time_min")))
    return grouped


def cumulative_auc(time_min: np.ndarray, speed_um_s: np.ndarray) -> np.ndarray:
    if time_min.size == 0:
        return np.array([])
    auc = np.zeros_like(time_min, dtype=float)
    for i in range(1, time_min.size):
        dt_s = (time_min[i] - time_min[i - 1]) * 60.0
        auc[i] = auc[i - 1] + 0.5 * (speed_um_s[i] + speed_um_s[i - 1]) * dt_s
    return auc


def plot_timecourse(
    ax: plt.Axes,
    grouped: dict[tuple[str, str], list[dict[str, object]]],
    variant: str,
    field: str,
    ylabel: str,
    *,
    cumulative: bool = False,
) -> None:
    for condition in CONDITION_ORDER:
        rows = grouped[(variant, condition)]
        time = np.array([f(row.get("time_min")) for row in rows], dtype=float)
        mean = np.array([f(row.get(f"{field}_mean")) for row in rows], dtype=float)
        sem = np.array([f(row.get(f"{field}_sem"), 0.0) for row in rows], dtype=float)
        if cumulative:
            post_onset = time >= 1.0
            time = time[post_onset]
            mean = mean[post_onset]
            mean = cumulative_auc(time, mean)
            sem = np.zeros_like(mean)
        color = CONDITION_COLORS[condition]
        ax.plot(time, mean, color=color, linewidth=2.3, label=CONDITION_LABELS[condition])
        if not cumulative:
            ax.fill_between(time, mean - sem, mean + sem, color=color, alpha=0.18, linewidth=0)
    ax.axvline(1.0, color="#555555", linestyle=":", linewidth=1.1)
    ax.set_ylabel(ylabel)
    ax.grid(axis="both", color="#d8d8d8", linewidth=0.8, alpha=0.75)
    ax.set_axisbelow(True)
    ax.margins(x=0.02)


def plot_velocity_timecourses(rows: list[dict[str, object]]) -> None:
    grouped = group_timecourses(rows)
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.6), sharex=True)
    for row_index, (variant, label, _folder) in enumerate(INPUTS):
        speed_ax = axes[row_index, 0]
        auc_ax = axes[row_index, 1]
        plot_timecourse(speed_ax, grouped, variant, "vz_abs", "Mean |v_z| (um/s)")
        plot_timecourse(
            auc_ax,
            grouped,
            variant,
            "vz_abs",
            "Post-onset AUC |v_z| (um)",
            cumulative=True,
        )
        speed_ax.text(0.03, 0.92, label, transform=speed_ax.transAxes, ha="left", va="top", fontsize=13)
        auc_ax.text(0.03, 0.92, label, transform=auc_ax.transAxes, ha="left", va="top", fontsize=13)
    for ax, label in zip(axes.flat, ["A", "B", "C", "D"]):
        panel_label(ax, label)
    for ax in axes[-1, :]:
        ax.set_xlabel("Time (min)")
    axes[0, 1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout(w_pad=2.4, h_pad=1.8)
    save_figure(fig, "large_velocity_timecourses")


def plot_turnover_timecourses(rows: list[dict[str, object]]) -> None:
    grouped = group_timecourses(rows)
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.6), sharex=True)
    for row_index, (variant, label, _folder) in enumerate(INPUTS):
        chewer_ax = axes[row_index, 0]
        bottom_ax = axes[row_index, 1]
        plot_timecourse(
            chewer_ax,
            grouped,
            variant,
            "chewer_band_mass_fraction",
            "Chewer-band mass fraction",
        )
        plot_timecourse(
            bottom_ax,
            grouped,
            variant,
            "bottom_half_mass_fraction",
            "Bottom-half mass fraction",
        )
        chewer_ax.text(0.03, 0.92, label, transform=chewer_ax.transAxes, ha="left", va="top", fontsize=13)
        bottom_ax.text(0.03, 0.92, label, transform=bottom_ax.transAxes, ha="left", va="top", fontsize=13)
    for ax, label in zip(axes.flat, ["A", "B", "C", "D"]):
        panel_label(ax, label)
    for ax in axes[-1, :]:
        ax.set_xlabel("Time (min)")
    axes[0, 1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout(w_pad=2.4, h_pad=1.8)
    save_figure(fig, "large_turnover_timecourses")


def write_readme(metrics: list[dict[str, object]]) -> None:
    indexed = metric_row_index(metrics)
    base_motor = indexed[("baseline", "rotatable_xlink")]
    base_nomotor = indexed[("baseline", "nomotor_xlink")]
    ultra_motor = indexed[("ultra", "rotatable_xlink")]
    ultra_nomotor = indexed[("ultra", "nomotor_xlink")]

    base_auc_gain = f(base_motor["auc_vz_abs"]) / f(base_nomotor["auc_vz_abs"])
    ultra_auc_gain = f(ultra_motor["auc_vz_abs"]) / f(ultra_nomotor["auc_vz_abs"])
    chewer_drop = f(base_motor["chewer_band_mass_fraction_final"]) - f(
        ultra_motor["chewer_band_mass_fraction_final"]
    )

    text = f"""# Large Continuous-Supply Turnover Presentation Plots

Inputs:
- Baseline: `analysis/results/turnover_continuous_supply_scaled_job09`
- Ultra aggressive chewing: `analysis/results/turnover_continuous_supply_scaled_job12_chewer_count2x_rate2x`

Included conditions:
- `large_long/nomotor_xlink`
- `large_long/rotatable_xlink`

Excluded:
- Small-system tests are intentionally omitted.
- `job13_sever_chew` is not plotted because all 20 runs were partial/incomplete in `run_status.csv`; representative run `r0001` stopped with Cytosim status 11 and `Segmentation fault`.

Generated figures:
- `large_story_summary.svg/.png`: compact comparison for the main slide.
- `large_velocity_endpoint_metrics.svg/.png`: endpoint velocity and directionality metrics.
- `large_velocity_timecourses.svg/.png`: speed and cumulative speed timecourses.
- `large_turnover_endpoint_metrics.svg/.png`: final bottom-clearing and mass readouts.
- `large_turnover_timecourses.svg/.png`: bottom-region mass timecourses.

Key large-system readout:
- Baseline motors increase AUC |v_z| by {fmt(base_auc_gain, 2)}x over no motors.
- Ultra-aggressive chewing motors increase AUC |v_z| by {fmt(ultra_auc_gain, 2)}x over no motors.
- Ultra-aggressive chewing drops the motor-condition final chewer-band mass fraction by {fmt(chewer_drop, 3)} versus baseline.
- Motor AUC remains comparable after aggressive chewing: baseline {fmt(f(base_motor['auc_vz_abs']), 2)} um, ultra {fmt(f(ultra_motor['auc_vz_abs']), 2)} um.
"""
    (OUT / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup_style()
    OUT.mkdir(parents=True, exist_ok=True)
    metrics = load_metrics()
    timecourses = load_timecourses()
    write_csv(OUT / "large_condition_summary.csv", metrics)
    write_csv(OUT / "large_timecourse_summary.csv", timecourses)
    plot_endpoint_panels(metrics)
    plot_velocity_timecourses(timecourses)
    plot_turnover_timecourses(timecourses)
    write_readme(metrics)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
