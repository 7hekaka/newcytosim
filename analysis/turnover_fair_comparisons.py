#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "analysis" / "results"
OUT = RESULTS / "turnover_fair_comparisons"

CAMPAIGN_SPECS = [
    {"name": "turnover_growth_sweep", "label": "Growth sweep"},
    {"name": "turnover_growth_pool_bias_sweep", "label": "Pool-bias sweep"},
    {"name": "turnover_growth_pool_bias_bottom_chewers", "label": "Bottom chewers"},
    {"name": "turnover_growth_pool_bias_bottom_chewers_diffuse", "label": "Bottom chewers diffuse"},
    {"name": "turnover_growth_pool_bias_motor_escalation", "label": "Motor escalation"},
]
CAMPAIGN_LABELS = {spec["name"]: spec["label"] for spec in CAMPAIGN_SPECS}
CAMPAIGN_ORDER = [spec["name"] for spec in CAMPAIGN_SPECS]

SCENARIO_ORDER = [
    "sc00_full_height_bipolar",
    "sc01_full_height_aligned",
    "sc02_top_two_thirds_aligned",
    "sc03_top_one_third_aligned",
    "sc04_top_cap_aligned",
]
SCENARIO_LABELS = {
    "sc00_full_height_bipolar": "Full height mixed",
    "sc01_full_height_aligned": "Full height aligned",
    "sc02_top_two_thirds_aligned": "Top 2/3 aligned",
    "sc03_top_one_third_aligned": "Top 1/3 aligned",
    "sc04_top_cap_aligned": "Top cap aligned",
}

CONDITION_ORDER = [
    "c00_nomotor_noxlink",
    "c01_nomotor_xlink",
    "c02_rotatable_noxlink",
    "c03_rotatable_xlink",
    "c04_fixed_global_noxlink",
    "c05_fixed_global_xlink",
]
CONDITION_LABELS = {
    "c00_nomotor_noxlink": "No motors / no xlinks",
    "c01_nomotor_xlink": "No motors / + xlinks",
    "c02_rotatable_noxlink": "Rotatable / no xlinks",
    "c03_rotatable_xlink": "Rotatable / + xlinks",
    "c04_fixed_global_noxlink": "Fixed-global / no xlinks",
    "c05_fixed_global_xlink": "Fixed-global / + xlinks",
}

RAW_METRICS = [
    "mean_vz_abs",
    "auc_vz_abs",
    "excess_auc_vz_abs",
    "longest_active_streak_min",
    "abs_net_z_displacement",
    "mean_vz_signed",
    "directionality_index",
    "mass_normalized_auc",
    "mass_retention",
    "swirl_penalty",
    "chewer_band_mass_fraction_final",
    "bottom_half_mass_fraction_final",
    "seed_zone_retention",
    "axial_spread_entropy_mean",
]
DELTA_METRICS = RAW_METRICS + ["matched_active_transport_fraction"]
PLOTTED_DELTA_METRICS = [
    "auc_vz_abs",
    "matched_active_transport_fraction",
    "excess_auc_vz_abs",
    "longest_active_streak_min",
    "mass_normalized_auc",
    "mass_retention",
    "swirl_penalty",
    "chewer_band_mass_fraction_final",
    "bottom_half_mass_fraction_final",
    "seed_zone_retention",
    "axial_spread_entropy_mean",
]
METRIC_LABELS = {
    "mean_vz_abs": "Mean |v_z| delta (um/s)",
    "auc_vz_abs": "AUC(|v_z|) delta (um)",
    "excess_auc_vz_abs": "Excess AUC delta (um)",
    "longest_active_streak_min": "Longest active streak delta (min)",
    "abs_net_z_displacement": "|Net z shift| delta (um)",
    "mean_vz_signed": "Mean v_z delta (um/s)",
    "directionality_index": "Directionality delta",
    "mass_normalized_auc": "Mass-normalized AUC delta",
    "mass_retention": "Mass-retention delta",
    "swirl_penalty": "Swirl-penalty delta",
    "chewer_band_mass_fraction_final": "Final chewer-band mass fraction delta",
    "bottom_half_mass_fraction_final": "Final bottom-half mass fraction delta",
    "seed_zone_retention": "Seed-zone retention delta",
    "axial_spread_entropy_mean": "Mean axial entropy delta",
    "matched_active_transport_fraction": "Matched active fraction delta",
}
METRIC_FORMATS = {
    "mean_vz_abs": "{:.3f}",
    "auc_vz_abs": "{:.2f}",
    "excess_auc_vz_abs": "{:.2f}",
    "longest_active_streak_min": "{:.2f}",
    "abs_net_z_displacement": "{:.2f}",
    "mean_vz_signed": "{:.3f}",
    "directionality_index": "{:.2f}",
    "mass_normalized_auc": "{:.3f}",
    "mass_retention": "{:.2f}",
    "swirl_penalty": "{:.2f}",
    "chewer_band_mass_fraction_final": "{:.2f}",
    "bottom_half_mass_fraction_final": "{:.2f}",
    "seed_zone_retention": "{:.2f}",
    "axial_spread_entropy_mean": "{:.2f}",
    "matched_active_transport_fraction": "{:.2f}",
}

EFFECT_SPECS = [
    {
        "slug": "motor_gain_rotatable_noxlink",
        "label": "Rotatable motor gain without xlinks",
        "numerator": "c02_rotatable_noxlink",
        "denominator": "c00_nomotor_noxlink",
    },
    {
        "slug": "motor_gain_rotatable_xlink",
        "label": "Rotatable motor gain with xlinks",
        "numerator": "c03_rotatable_xlink",
        "denominator": "c01_nomotor_xlink",
    },
    {
        "slug": "motor_gain_fixed_global_noxlink",
        "label": "Fixed-global motor gain without xlinks",
        "numerator": "c04_fixed_global_noxlink",
        "denominator": "c00_nomotor_noxlink",
    },
    {
        "slug": "motor_gain_fixed_global_xlink",
        "label": "Fixed-global motor gain with xlinks",
        "numerator": "c05_fixed_global_xlink",
        "denominator": "c01_nomotor_xlink",
    },
    {
        "slug": "xlink_gain_nomotor",
        "label": "Xlink gain in no-motor controls",
        "numerator": "c01_nomotor_xlink",
        "denominator": "c00_nomotor_noxlink",
    },
    {
        "slug": "xlink_gain_rotatable",
        "label": "Xlink gain in rotatable runs",
        "numerator": "c03_rotatable_xlink",
        "denominator": "c02_rotatable_noxlink",
    },
    {
        "slug": "xlink_gain_fixed_global",
        "label": "Xlink gain in fixed-global runs",
        "numerator": "c05_fixed_global_xlink",
        "denominator": "c04_fixed_global_noxlink",
    },
    {
        "slug": "orientation_gain_noxlink",
        "label": "Rotatable minus fixed-global without xlinks",
        "numerator": "c02_rotatable_noxlink",
        "denominator": "c04_fixed_global_noxlink",
    },
    {
        "slug": "orientation_gain_xlink",
        "label": "Rotatable minus fixed-global with xlinks",
        "numerator": "c03_rotatable_xlink",
        "denominator": "c05_fixed_global_xlink",
    },
]

REFERENCE_COMPARISONS = [
    {
        "target": "turnover_growth_pool_bias_sweep",
        "reference": "turnover_growth_sweep",
        "label": "Pool-bias minus growth sweep",
    },
    {
        "target": "turnover_growth_pool_bias_bottom_chewers",
        "reference": "turnover_growth_pool_bias_sweep",
        "label": "Bottom chewers minus pool-bias sweep",
    },
    {
        "target": "turnover_growth_pool_bias_bottom_chewers_diffuse",
        "reference": "turnover_growth_pool_bias_bottom_chewers",
        "label": "Diffuse bottom chewers minus fixed bottom chewers",
    },
    {
        "target": "turnover_growth_pool_bias_bottom_chewers_diffuse",
        "reference": "turnover_growth_pool_bias_sweep",
        "label": "Diffuse bottom chewers minus pool-bias sweep",
    },
    {
        "target": "turnover_growth_pool_bias_motor_escalation",
        "reference": "turnover_growth_pool_bias_sweep",
        "label": "Motor escalation minus pool-bias sweep",
    },
]


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def as_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value in {"", None}:
        return float("nan")
    try:
        return float(value)
    except ValueError:
        return float("nan")


def load_campaign_runs(campaign_name: str) -> tuple[list[dict], dict[tuple[str, str], dict]]:
    path = RESULTS / campaign_name / "run_metrics.csv"
    rows = read_csv_rows(path)
    for row in rows:
        row["campaign"] = campaign_name
        for key in RAW_METRICS:
            row[f"_{key}"] = as_float(row, key)
    lookup = {(row["scenario"], row["condition"]): row for row in rows}
    return rows, lookup


def load_timecourse_rows(campaign_name: str) -> list[dict[str, str]]:
    path = RESULTS / campaign_name / "run_timecourses_long.csv"
    if not path.is_file():
        return []
    return read_csv_rows(path)


def compute_matched_thresholds(
    campaign_name: str,
    run_rows: list[dict],
) -> tuple[list[dict], dict[tuple[str, int], float], dict[tuple[str, str, str], float]]:
    timecourse_rows = load_timecourse_rows(campaign_name)
    meta = {(row["scenario"], row["condition"], row["run_dir"]): row for row in run_rows}

    grouped: dict[tuple[str, int], list[float]] = defaultdict(list)
    per_run_values: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in timecourse_rows:
        if row.get("phase_name") != "motor":
            continue
        key = (row["scenario"], row["condition"], row["run_dir"])
        info = meta.get(key)
        if info is None:
            continue
        value = as_float(row, "vz_abs")
        per_run_values[key].append(value)
        if info["condition"] in {"c00_nomotor_noxlink", "c01_nomotor_xlink"} and math.isfinite(value):
            grouped[(info["scenario"], int(info["include_xlinks"]))].append(value)

    threshold_rows: list[dict] = []
    thresholds: dict[tuple[str, int], float] = {}
    for scenario in SCENARIO_ORDER:
        for include_xlinks in (0, 1):
            values = np.asarray(grouped.get((scenario, include_xlinks), []), dtype=float)
            threshold = float(np.nanpercentile(values, 95.0)) if values.size else float("nan")
            thresholds[(scenario, include_xlinks)] = threshold
            threshold_rows.append(
                {
                    "campaign": campaign_name,
                    "campaign_label": CAMPAIGN_LABELS[campaign_name],
                    "scenario": scenario,
                    "scenario_label": SCENARIO_LABELS[scenario],
                    "include_xlinks": include_xlinks,
                    "threshold_vz_abs_matched": threshold,
                    "n_control_frame_values": int(values.size),
                }
            )

    matched_active: dict[tuple[str, str, str], float] = {}
    for key, values in per_run_values.items():
        scenario, _condition, _run_dir = key
        info = meta[key]
        threshold = thresholds.get((scenario, int(info["include_xlinks"])), float("nan"))
        arr = np.asarray(values, dtype=float)
        if arr.size and math.isfinite(threshold):
            matched_active[key] = float(np.mean(arr > threshold))
        else:
            matched_active[key] = float("nan")

    return threshold_rows, thresholds, matched_active


def enrich_runs_with_matched_active(
    run_rows: list[dict],
    matched_active: dict[tuple[str, str, str], float],
) -> None:
    for row in run_rows:
        key = (row["scenario"], row["condition"], row["run_dir"])
        row["_matched_active_transport_fraction"] = matched_active.get(key, float("nan"))


def compute_effect_rows(campaign_runs: dict[str, dict[tuple[str, str], dict]]) -> list[dict]:
    rows: list[dict] = []
    for campaign_name in CAMPAIGN_ORDER:
        lookup = campaign_runs[campaign_name]
        for effect in EFFECT_SPECS:
            numerator = effect["numerator"]
            denominator = effect["denominator"]
            for scenario in SCENARIO_ORDER:
                num_row = lookup.get((scenario, numerator))
                den_row = lookup.get((scenario, denominator))
                out = {
                    "campaign": campaign_name,
                    "campaign_label": CAMPAIGN_LABELS[campaign_name],
                    "scenario": scenario,
                    "scenario_label": SCENARIO_LABELS[scenario],
                    "effect": effect["slug"],
                    "effect_label": effect["label"],
                    "numerator_condition": numerator,
                    "numerator_label": CONDITION_LABELS[numerator],
                    "denominator_condition": denominator,
                    "denominator_label": CONDITION_LABELS[denominator],
                }
                if num_row is None or den_row is None:
                    for metric in DELTA_METRICS:
                        out[f"delta_{metric}"] = float("nan")
                    rows.append(out)
                    continue
                for metric in RAW_METRICS:
                    out[f"numerator_{metric}"] = num_row[f"_{metric}"]
                    out[f"denominator_{metric}"] = den_row[f"_{metric}"]
                    out[f"delta_{metric}"] = num_row[f"_{metric}"] - den_row[f"_{metric}"]
                num_active = num_row.get("_matched_active_transport_fraction", float("nan"))
                den_active = den_row.get("_matched_active_transport_fraction", float("nan"))
                out["numerator_matched_active_transport_fraction"] = num_active
                out["denominator_matched_active_transport_fraction"] = den_active
                out["delta_matched_active_transport_fraction"] = num_active - den_active
                rows.append(out)
    return rows


def compute_reference_rows(campaign_runs: dict[str, dict[tuple[str, str], dict]]) -> list[dict]:
    rows: list[dict] = []
    for spec in REFERENCE_COMPARISONS:
        target = spec["target"]
        reference = spec["reference"]
        target_lookup = campaign_runs[target]
        reference_lookup = campaign_runs[reference]
        all_conditions = sorted(set(condition for _scenario, condition in target_lookup.keys()) | set(condition for _scenario, condition in reference_lookup.keys()))
        for scenario in SCENARIO_ORDER:
            for condition in all_conditions:
                target_row = target_lookup.get((scenario, condition))
                reference_row = reference_lookup.get((scenario, condition))
                out = {
                    "target_campaign": target,
                    "target_label": CAMPAIGN_LABELS[target],
                    "reference_campaign": reference,
                    "reference_label": CAMPAIGN_LABELS[reference],
                    "comparison_label": spec["label"],
                    "scenario": scenario,
                    "scenario_label": SCENARIO_LABELS[scenario],
                    "condition": condition,
                    "condition_label": CONDITION_LABELS.get(condition, condition),
                }
                if target_row is None or reference_row is None:
                    for metric in DELTA_METRICS:
                        out[f"delta_{metric}"] = float("nan")
                    rows.append(out)
                    continue
                for metric in RAW_METRICS:
                    out[f"target_{metric}"] = target_row[f"_{metric}"]
                    out[f"reference_{metric}"] = reference_row[f"_{metric}"]
                    out[f"delta_{metric}"] = target_row[f"_{metric}"] - reference_row[f"_{metric}"]
                target_active = target_row.get("_matched_active_transport_fraction", float("nan"))
                reference_active = reference_row.get("_matched_active_transport_fraction", float("nan"))
                out["target_matched_active_transport_fraction"] = target_active
                out["reference_matched_active_transport_fraction"] = reference_active
                out["delta_matched_active_transport_fraction"] = target_active - reference_active
                rows.append(out)
    return rows


def value_limits(arr: np.ndarray) -> tuple[float, float]:
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return -1.0, 1.0
    bound = float(np.nanmax(np.abs(finite)))
    if bound < 1e-12:
        bound = 1.0
    return -bound, bound


def draw_matrix(
    arr: np.ndarray,
    row_labels: list[str],
    col_labels: list[str],
    title: str,
    outpath: Path,
    metric: str,
) -> Path:
    vmin, vmax = value_limits(arr)
    fig, ax = plt.subplots(figsize=(1.9 * len(col_labels) + 2.8, 0.75 * len(row_labels) + 2.0))
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad("#efefef")
    image = ax.imshow(arr, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=25, ha="right", fontsize=10)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=10)
    ax.set_title(title, fontsize=15, pad=10)
    fmt = METRIC_FORMATS[metric]
    midpoint = 0.5 * (vmin + vmax)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            value = arr[i, j]
            text = "n/a" if not np.isfinite(value) else fmt.format(value)
            color = "white" if np.isfinite(value) and value > midpoint + 0.15 * (vmax - vmin) else "black"
            ax.text(j, i, text, ha="center", va="center", fontsize=8.5, color=color)
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label(METRIC_LABELS[metric], fontsize=10)
    fig.subplots_adjust(left=0.22, right=0.96, top=0.90, bottom=0.22)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outpath


def plot_effect_matrices(effect_rows: list[dict]) -> list[Path]:
    outpaths: list[Path] = []
    for metric in PLOTTED_DELTA_METRICS:
        for effect in EFFECT_SPECS:
            arr = np.full((len(CAMPAIGN_ORDER), len(SCENARIO_ORDER)), np.nan, dtype=float)
            row_labels = [CAMPAIGN_LABELS[name] for name in CAMPAIGN_ORDER]
            col_labels = [SCENARIO_LABELS[name] for name in SCENARIO_ORDER]
            for i, campaign_name in enumerate(CAMPAIGN_ORDER):
                for j, scenario in enumerate(SCENARIO_ORDER):
                    value = float("nan")
                    for row in effect_rows:
                        if row["campaign"] == campaign_name and row["scenario"] == scenario and row["effect"] == effect["slug"]:
                            value = row[f"delta_{metric}"]
                            break
                    arr[i, j] = value
            outpaths.append(
                draw_matrix(
                    arr,
                    row_labels,
                    col_labels,
                    f"{effect['label']} across campaigns",
                    OUT / "within_campaign" / metric / f"{effect['slug']}.png",
                    metric,
                )
            )
    return outpaths


def plot_reference_matrices(reference_rows: list[dict]) -> list[Path]:
    outpaths: list[Path] = []
    for metric in PLOTTED_DELTA_METRICS:
        for spec in REFERENCE_COMPARISONS:
            condition_set = sorted(
                {
                    row["condition"]
                    for row in reference_rows
                    if row["target_campaign"] == spec["target"] and row["reference_campaign"] == spec["reference"]
                },
                key=lambda value: CONDITION_ORDER.index(value) if value in CONDITION_ORDER else value,
            )
            arr = np.full((len(SCENARIO_ORDER), len(condition_set)), np.nan, dtype=float)
            row_labels = [SCENARIO_LABELS[name] for name in SCENARIO_ORDER]
            col_labels = [CONDITION_LABELS.get(name, name) for name in condition_set]
            for i, scenario in enumerate(SCENARIO_ORDER):
                for j, condition in enumerate(condition_set):
                    value = float("nan")
                    for row in reference_rows:
                        if (
                            row["target_campaign"] == spec["target"]
                            and row["reference_campaign"] == spec["reference"]
                            and row["scenario"] == scenario
                            and row["condition"] == condition
                        ):
                            value = row[f"delta_{metric}"]
                            break
                    arr[i, j] = value
            slug = f"{spec['target']}__minus__{spec['reference']}"
            outpaths.append(
                draw_matrix(
                    arr,
                    row_labels,
                    col_labels,
                    spec["label"],
                    OUT / "cross_campaign" / metric / f"{slug}.png",
                    metric,
                )
            )
    return outpaths


def top_rows(rows: list[dict], key: str, n: int = 5) -> list[dict]:
    finite = [row for row in rows if math.isfinite(row.get(key, float("nan")))]
    return sorted(finite, key=lambda row: row[key], reverse=True)[:n]


def write_readme(
    threshold_rows: list[dict],
    effect_rows: list[dict],
    reference_rows: list[dict],
) -> None:
    top_motor_gain = top_rows(
        [row for row in effect_rows if row["effect"] == "motor_gain_rotatable_xlink"],
        "delta_auc_vz_abs",
    )
    top_excess_motor_gain = top_rows(
        [row for row in effect_rows if row["effect"] == "motor_gain_rotatable_xlink"],
        "delta_excess_auc_vz_abs",
    )
    top_orientation = top_rows(
        [row for row in effect_rows if row["effect"] == "orientation_gain_xlink"],
        "delta_auc_vz_abs",
    )
    top_diffuse_gain = top_rows(
        [
            row
            for row in reference_rows
            if row["target_campaign"] == "turnover_growth_pool_bias_bottom_chewers_diffuse"
            and row["reference_campaign"] == "turnover_growth_pool_bias_bottom_chewers"
        ],
        "delta_auc_vz_abs",
    )
    top_escalation_gain = top_rows(
        [
            row
            for row in reference_rows
            if row["target_campaign"] == "turnover_growth_pool_bias_motor_escalation"
            and row["reference_campaign"] == "turnover_growth_pool_bias_sweep"
        ],
        "delta_auc_vz_abs",
    )
    diffuse_bottom_half_reductions = sorted(
        [
            row
            for row in reference_rows
            if row["target_campaign"] == "turnover_growth_pool_bias_bottom_chewers_diffuse"
            and row["reference_campaign"] == "turnover_growth_pool_bias_bottom_chewers"
            and math.isfinite(row.get("delta_bottom_half_mass_fraction_final", float("nan")))
        ],
        key=lambda row: row["delta_bottom_half_mass_fraction_final"],
    )[:5]

    lines = [
        "# Turnover Fair Comparisons",
        "",
        "## Fairness rules",
        "- Scenario-to-scenario raw transport is not treated as a causal ranking by itself.",
        "- Motor effects are evaluated against matched no-motor baselines from the same campaign, scenario, and xlink state.",
        "- Xlink effects are evaluated within the same campaign, scenario, and motor mode.",
        "- Rotatable vs fixed-global is evaluated within the same campaign, scenario, and xlink state.",
        "- Active transport is recomputed with scenario-matched thresholds from the motor-phase `vz_abs` distribution of the no-motor controls.",
        "- Cross-campaign differences are reported only for matched scenario-condition pairs.",
        "",
        "## Scenario-matched thresholds",
    ]
    for row in threshold_rows:
        if row["campaign"] != "turnover_growth_pool_bias_bottom_chewers_diffuse":
            continue
        lines.append(
            f"- Diffuse bottom chewers, `{row['scenario']}`, xlinks={row['include_xlinks']}: threshold = `{row['threshold_vz_abs_matched']:.6f} um/s` from `{row['n_control_frame_values']}` control frames."
        )

    lines.extend(["", "## Largest rotatable-plus-xlink motor gains by delta AUC"])
    for row in top_motor_gain:
        lines.append(
            f"- `{row['campaign']}` / `{row['scenario']}`: `{row['delta_auc_vz_abs']:.3f} um` (`{row['numerator_condition']}` minus `{row['denominator_condition']}`)."
        )

    lines.extend(["", "## Largest rotatable-plus-xlink motor gains by excess AUC"])
    for row in top_excess_motor_gain:
        lines.append(
            f"- `{row['campaign']}` / `{row['scenario']}`: `{row['delta_excess_auc_vz_abs']:.3f} um` above the matched no-motor baseline."
        )

    lines.extend(["", "## Largest rotatable-minus-fixed-global gains with xlinks"])
    for row in top_orientation:
        lines.append(
            f"- `{row['campaign']}` / `{row['scenario']}`: `{row['delta_auc_vz_abs']:.3f} um` (`c03_rotatable_xlink - c05_fixed_global_xlink`)."
        )

    lines.extend(["", "## Largest diffuse-minus-fixed bottom-chewer gains"])
    for row in top_diffuse_gain:
        lines.append(
            f"- `{row['scenario']}` / `{row['condition']}`: `{row['delta_auc_vz_abs']:.3f} um` in AUC and `{row['delta_matched_active_transport_fraction']:.3f}` in matched active fraction."
        )

    lines.extend(["", "## Largest motor-escalation gains over the pool-bias sweep"])
    for row in top_escalation_gain:
        lines.append(
            f"- `{row['scenario']}` / `{row['condition']}`: `{row['delta_auc_vz_abs']:.3f} um` in AUC."
        )

    if diffuse_bottom_half_reductions:
        lines.extend(["", "## Largest diffuse-chewer reductions in final common bottom-half mass"])
        for row in diffuse_bottom_half_reductions:
            lines.append(
                f"- `{row['scenario']}` / `{row['condition']}`: `{row['delta_bottom_half_mass_fraction_final']:.3f}` final bottom-half fraction vs fixed bottom chewers."
            )

    lines.extend(
        [
            "",
            "## Files",
            "- `scenario_matched_thresholds.csv`: scenario-specific control thresholds used for fair active-transport comparisons.",
            "- `effect_deltas.csv`: within-campaign matched comparisons.",
            "- `reference_deltas.csv`: matched cross-campaign comparisons.",
            "- `within_campaign/<metric>/`: heatmaps of fair within-campaign deltas for AUC plus the eight derived metrics.",
            "- `cross_campaign/<metric>/`: heatmaps of matched condition deltas across campaign transitions.",
        ]
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.linewidth": 1.0,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )

    all_run_rows: dict[str, list[dict]] = {}
    campaign_runs: dict[str, dict[tuple[str, str], dict]] = {}
    threshold_rows: list[dict] = []

    for campaign_name in CAMPAIGN_ORDER:
        run_rows, lookup = load_campaign_runs(campaign_name)
        local_threshold_rows, _thresholds, matched_active = compute_matched_thresholds(campaign_name, run_rows)
        enrich_runs_with_matched_active(run_rows, matched_active)
        all_run_rows[campaign_name] = run_rows
        campaign_runs[campaign_name] = lookup
        threshold_rows.extend(local_threshold_rows)

    effect_rows = compute_effect_rows(campaign_runs)
    reference_rows = compute_reference_rows(campaign_runs)

    write_csv(OUT / "scenario_matched_thresholds.csv", threshold_rows)
    write_csv(OUT / "effect_deltas.csv", effect_rows)
    write_csv(OUT / "reference_deltas.csv", reference_rows)

    plot_effect_matrices(effect_rows)
    plot_reference_matrices(reference_rows)
    write_readme(threshold_rows, effect_rows, reference_rows)

    print(OUT / "README.md")
    print(OUT / "scenario_matched_thresholds.csv")
    print(OUT / "effect_deltas.csv")
    print(OUT / "reference_deltas.csv")


if __name__ == "__main__":
    main()
