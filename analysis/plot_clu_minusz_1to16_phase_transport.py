#!/usr/bin/env python3
"""Phase-resolved 1:16 transport plots for random/aligned annulus runs.

This companion script keeps the original full timecourse plots intact and adds
early/middle/late summaries so transient axial transport is not visually washed
out by the late relaxed regime.
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

from analysis import plot_clu_minusz_1to16_publication_transport as transport


OUT = transport.OUT / "phase_regime_plots"
TABLE_OUT = OUT / "tables"

PHASES = [
    ("early", "Early", 0.0, 2.0, "#F3D6A4"),
    ("middle", "Middle", 2.0, 6.0, "#D7E7F4"),
    ("late", "Late", 6.0, transport.TARGET_FINAL_MIN, "#E4E4E4"),
]

PHASE_COLORS = {
    "early": "#D55E00",
    "middle": "#0072B2",
    "late": "#666666",
}


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
        zorder=20,
    )


def save_all_formats(fig: plt.Figure, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=600, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    return path


def analyze_timecourse_task(task: tuple[dict[str, str], str]) -> dict[str, object]:
    return transport.analyze_timecourse_task(task)


def load_run_timecourses(family: str) -> list[dict[str, object]]:
    cache = TABLE_OUT / f"{family}_run_timecourses_1to16.csv"
    if cache.exists():
        rows = read_csv(cache)
        return [
            {
                **row,
                "time_min": [finite_float(value) for value in row["time_min"].split(";") if value],
                "vz_signed": [finite_float(value) for value in row["vz_signed"].split(";") if value],
            }
            for row in rows
        ]

    report_bin = "/bin/true"
    tasks: list[tuple[dict[str, str], str]] = []
    for panel_key, (path, _label) in transport.TIMECOURSE_MAPS.items():
        for row in transport.load_timecourse_metadata(path, panel_key, family):
            tasks.append((row, report_bin))

    with ProcessPoolExecutor(max_workers=8) as executor:
        raw_rows = list(executor.map(analyze_timecourse_task, tasks))

    rows: list[dict[str, object]] = []
    status_rows: list[dict[str, object]] = []
    for row in raw_rows:
        status_rows.append({key: value for key, value in row.items() if key not in {"time_min", "z_com", "vz_signed"}})
        if row.get("status") != "ok":
            continue
        time = np.asarray(row["time_min"], dtype=float)
        keep = time <= transport.TARGET_FINAL_MIN + 1e-9
        time = time[keep]
        vz = np.asarray(row["vz_signed"], dtype=float)[keep]
        rows.append(
            {
                "panel_key": row["panel_key"],
                "case": row["case"],
                "xlink_regime": row["xlink_regime"],
                "run_dir": row["run_dir"],
                "run_path": row["run_path"],
                "time_min": time.tolist(),
                "vz_signed": vz.tolist(),
            }
        )

    serial_rows = []
    for row in rows:
        serial_rows.append(
            {
                **{key: value for key, value in row.items() if key not in {"time_min", "vz_signed"}},
                "time_min": ";".join(f"{value:.10g}" for value in row["time_min"]),
                "vz_signed": ";".join(f"{value:.10g}" for value in row["vz_signed"]),
            }
        )
    write_csv(cache, serial_rows)
    write_csv(TABLE_OUT / f"{family}_run_status_1to16.csv", status_rows)
    return rows


def integrate_phase(time_min: np.ndarray, velocity: np.ndarray, start_min: float, end_min: float) -> tuple[float, float]:
    if time_min.size == 0 or velocity.size == 0 or end_min <= start_min:
        return float("nan"), float("nan")
    keep = (time_min > start_min) & (time_min < end_min)
    t = np.concatenate([[start_min], time_min[keep], [end_min]])
    v = np.interp(t, time_min, velocity)
    displacement = float(np.trapezoid(v, t * 60.0))
    mean_velocity = displacement / ((end_min - start_min) * 60.0)
    return mean_velocity, displacement


def cumulative_displacement(time_min: np.ndarray, velocity: np.ndarray) -> np.ndarray:
    if time_min.size == 0:
        return np.asarray([], dtype=float)
    disp = np.zeros_like(time_min, dtype=float)
    if time_min.size > 1:
        dt = np.diff(time_min) * 60.0
        disp[1:] = np.cumsum(0.5 * (velocity[1:] + velocity[:-1]) * dt)
    return disp


def aggregate_timecourses(rows: list[dict[str, object]], value_key: str = "vz_signed") -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["panel_key"]), str(row["case"]))].append(row)

    out: list[dict[str, object]] = []
    for (panel_key, case), items in sorted(grouped.items()):
        n_frames = min(len(item["time_min"]) for item in items)
        time = np.asarray(items[0]["time_min"][:n_frames], dtype=float)
        stack = []
        for item in items:
            values = np.asarray(item[value_key][:n_frames], dtype=float)
            stack.append(values)
        arr = np.vstack(stack)
        mean = np.nanmean(arr, axis=0)
        err = np.nanstd(arr, axis=0, ddof=1) / math.sqrt(arr.shape[0]) if arr.shape[0] > 1 else np.zeros_like(mean)
        for idx, t in enumerate(time):
            out.append(
                {
                    "panel_key": panel_key,
                    "case": case,
                    "frame_index": idx + 1,
                    "time_min": float(t),
                    f"{value_key}_mean": float(mean[idx]),
                    f"{value_key}_sem": float(err[idx]),
                    "n_runs_used": arr.shape[0],
                }
            )
    return out


def add_cumulative_to_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in rows:
        time = np.asarray(row["time_min"], dtype=float)
        velocity = np.asarray(row["vz_signed"], dtype=float)
        out.append({**row, "delta_z": cumulative_displacement(time, velocity).tolist()})
    return out


def phase_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in rows:
        time = np.asarray(row["time_min"], dtype=float)
        velocity = np.asarray(row["vz_signed"], dtype=float)
        for phase_key, phase_label, start, end, _color in PHASES:
            mean_v, disp = integrate_phase(time, velocity, start, end)
            out.append(
                {
                    "panel_key": row["panel_key"],
                    "case": row["case"],
                    "run_dir": row["run_dir"],
                    "phase_key": phase_key,
                    "phase_label": phase_label,
                    "phase_start_min": start,
                    "phase_end_min": end,
                    "phase_mean_vz": mean_v,
                    "phase_delta_z": disp,
                }
            )
    return out


def aggregate_phase(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["panel_key"]), str(row["case"]), str(row["phase_key"]))].append(row)

    out: list[dict[str, object]] = []
    for (panel_key, case, phase_key), items in sorted(grouped.items()):
        phase = next(phase for phase in PHASES if phase[0] == phase_key)
        mean_v = [finite_float(item["phase_mean_vz"]) for item in items]
        disp = [finite_float(item["phase_delta_z"]) for item in items]
        out.append(
            {
                "panel_key": panel_key,
                "case": case,
                "phase_key": phase_key,
                "phase_label": phase[1],
                "phase_start_min": phase[2],
                "phase_end_min": phase[3],
                "n_runs_used": len(items),
                "phase_mean_vz_mean": float(np.nanmean(mean_v)),
                "phase_mean_vz_sem": sem(mean_v),
                "phase_delta_z_mean": float(np.nanmean(disp)),
                "phase_delta_z_sem": sem(disp),
            }
        )
    return out


def tc_lookup(rows: list[dict[str, object]], value_key: str) -> dict[tuple[str, str], list[dict[str, object]]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["panel_key"]), str(row["case"]))].append(row)
    for key in grouped:
        grouped[key].sort(key=lambda row: int(row["frame_index"]))
    return grouped


def plot_timecourse_grid(
    rows: list[dict[str, object]],
    family: str,
    value_key: str,
    ylabel: str,
    outfile: Path,
    xlim: tuple[float, float] | None = None,
    shade_phases: bool = False,
) -> Path:
    lookup = tc_lookup(rows, value_key)
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.2), sharex=True, sharey=True)
    axes = axes.ravel()

    values: list[float] = []
    for series_rows in lookup.values():
        for row in series_rows:
            t = finite_float(row["time_min"])
            if xlim is not None and not (xlim[0] <= t <= xlim[1]):
                continue
            mean = finite_float(row[f"{value_key}_mean"])
            err = finite_float(row[f"{value_key}_sem"])
            if np.isfinite(mean):
                values.extend([mean - err if np.isfinite(err) else mean, mean + err if np.isfinite(err) else mean])
    ylo = min(values) if values else -0.001
    yhi = max(values) if values else 0.001
    pad = 0.12 * (yhi - ylo) if yhi > ylo else 0.001
    ylo -= pad
    yhi += pad
    if value_key == "vz_signed":
        ylo = min(ylo, -0.0008)
        yhi = max(yhi, 0.0008)

    for ax, (panel_key, panel_label) in zip(axes, transport.TIMECOURSE_PANEL_ORDER):
        style_axis(ax)
        if shade_phases:
            for _key, label, start, end, color in PHASES:
                if xlim is not None and (end < xlim[0] or start > xlim[1]):
                    continue
                left = max(start, xlim[0] if xlim else start)
                right = min(end, xlim[1] if xlim else end)
                ax.axvspan(left, right, color=color, alpha=0.22, zorder=0)
                if xlim is None:
                    ax.text((start + end) / 2, 0.985, label, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=10, color="#444444")
        if value_key == "vz_signed":
            ax.axhline(0, color="#777777", linewidth=1.0, linestyle="--", zorder=1)
        ax.set_title(panel_label, fontsize=17, pad=8)
        for case, count in transport.CASE_ORDER[family]:
            series_rows = lookup.get((panel_key, case), [])
            if not series_rows:
                continue
            t = np.asarray([finite_float(row["time_min"]) for row in series_rows], dtype=float)
            y = np.asarray([finite_float(row[f"{value_key}_mean"]) for row in series_rows], dtype=float)
            e = np.asarray([finite_float(row[f"{value_key}_sem"]) for row in series_rows], dtype=float)
            if xlim is not None:
                keep = (t >= xlim[0]) & (t <= xlim[1])
                t, y, e = t[keep], y[keep], e[keep]
            color = transport.COUNT_COLORS[count]
            ax.plot(t, y, color=color, linewidth=2.6, label=f"{count} clusters", zorder=3)
            ax.fill_between(t, y - e, y + e, color=color, alpha=0.14, linewidth=0, zorder=2)
        ax.set_ylim(ylo, yhi)
        ax.set_xlim(*(xlim if xlim is not None else (0.0, transport.TARGET_FINAL_MIN)))

    for idx, ax in enumerate(axes):
        add_panel_label(ax, chr(ord("A") + idx))
    axes[0].set_ylabel(ylabel, fontsize=16)
    axes[2].set_ylabel(ylabel, fontsize=16)
    axes[2].set_xlabel("Time after motor introduction (min)", fontsize=16)
    axes[3].set_xlabel("Time after motor introduction (min)", fontsize=16)

    handles = [
        plt.Line2D([0], [0], color=transport.COUNT_COLORS[count], lw=2.8, label=f"{count} clusters")
        for _case, count in transport.CASE_ORDER[family]
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=13, bbox_to_anchor=(0.5, 0.005))
    fig.subplots_adjust(left=0.10, right=0.99, top=0.93, bottom=0.13, wspace=0.10, hspace=0.30)
    return save_all_formats(fig, outfile)


def plot_phase_grid(rows: list[dict[str, object]], family: str, value_key: str, ylabel: str, outfile: Path) -> Path:
    grouped: dict[tuple[str, str], dict[str, dict[str, object]]] = defaultdict(dict)
    for row in rows:
        grouped[(str(row["panel_key"]), str(row["case"]))][str(row["phase_key"])] = row

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.2), sharex=True, sharey=True)
    axes = axes.ravel()

    values: list[float] = []
    for row in rows:
        mean = finite_float(row[f"{value_key}_mean"])
        err = finite_float(row[f"{value_key}_sem"])
        if np.isfinite(mean):
            values.extend([mean - err if np.isfinite(err) else mean, mean + err if np.isfinite(err) else mean])
    ylo = min(values) if values else -0.001
    yhi = max(values) if values else 0.001
    pad = 0.14 * (yhi - ylo) if yhi > ylo else 0.001
    ylo -= pad
    yhi += pad
    ylo = min(ylo, -0.0005 if value_key == "phase_mean_vz" else -0.05)

    for ax, (panel_key, panel_label) in zip(axes, transport.TIMECOURSE_PANEL_ORDER):
        style_axis(ax)
        ax.axhline(0, color="#777777", linewidth=1.0, linestyle="--", zorder=1)
        ax.set_title(panel_label, fontsize=17, pad=8)
        for phase_key, phase_label, _start, _end, _shade in PHASES:
            xs: list[int] = []
            ys: list[float] = []
            es: list[float] = []
            for case, count in transport.CASE_ORDER[family]:
                row = grouped.get((panel_key, case), {}).get(phase_key)
                if row is None:
                    continue
                xs.append(count)
                ys.append(finite_float(row[f"{value_key}_mean"]))
                es.append(finite_float(row[f"{value_key}_sem"]))
            ax.errorbar(
                xs,
                ys,
                yerr=es,
                color=PHASE_COLORS[phase_key],
                marker="o",
                linestyle="-",
                linewidth=2.6,
                markersize=7.0,
                capsize=4.2,
                elinewidth=1.4,
                label=phase_label,
            )
        ax.set_xticks([count for _case, count in transport.CASE_ORDER[family]])
        ax.set_ylim(ylo, yhi)

    for idx, ax in enumerate(axes):
        add_panel_label(ax, chr(ord("A") + idx))
    axes[0].set_ylabel(ylabel, fontsize=16)
    axes[2].set_ylabel(ylabel, fontsize=16)
    axes[2].set_xlabel("Cluster count", fontsize=16)
    axes[3].set_xlabel("Cluster count", fontsize=16)

    handles = [
        plt.Line2D([0], [0], color=PHASE_COLORS[key], lw=2.8, marker="o", markersize=7, label=label)
        for key, label, _start, _end, _shade in PHASES
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=13, bbox_to_anchor=(0.5, 0.005))
    fig.subplots_adjust(left=0.10, right=0.99, top=0.93, bottom=0.13, wspace=0.10, hspace=0.30)
    return save_all_formats(fig, outfile)


def write_readme(outputs: list[Path]) -> None:
    lines = [
        "# Phase-Resolved 1:16 Transport Plots",
        "",
        "These figures split the post-motor-introduction transport dynamics into early, middle, and late windows.",
        "",
        "Phase definitions:",
    ]
    for key, label, start, end, _color in PHASES:
        lines.append(f"- {label}: {start:g}-{end:.3g} min (`{key}`)")
    lines.extend(
        [
            "",
            "Interpretation note:",
            "- The aligned cases show a transient early axial-velocity pulse that relaxes toward baseline later in the simulation.",
            "- Phase averages and phase displacements are computed per replicate first, then summarized as mean +/- SEM.",
            "",
            "Generated outputs:",
        ]
    )
    for path in outputs:
        lines.append(f"- `{path.relative_to(ROOT)}`")
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def process_family(family: str) -> list[Path]:
    run_rows = load_run_timecourses(family)
    run_rows_with_disp = add_cumulative_to_rows(run_rows)

    phase_rows = phase_summary(run_rows)
    phase_agg = aggregate_phase(phase_rows)
    write_csv(TABLE_OUT / f"{family}_phase_run_metrics_1to16.csv", phase_rows)
    write_csv(TABLE_OUT / f"{family}_phase_condition_metrics_1to16.csv", phase_agg)

    vz_agg = aggregate_timecourses(run_rows, "vz_signed")
    disp_agg = aggregate_timecourses(run_rows_with_disp, "delta_z")
    write_csv(TABLE_OUT / f"{family}_velocity_timecourse_1to16.csv", vz_agg)
    write_csv(TABLE_OUT / f"{family}_displacement_timecourse_1to16.csv", disp_agg)

    outputs = [
        plot_timecourse_grid(
            vz_agg,
            family,
            "vz_signed",
            r"Axial velocity, $v_z$ ($\mu$m/s)",
            OUT / f"{family}_velocity_timecourse_shaded_phases_1to16.png",
            shade_phases=True,
        ),
        plot_timecourse_grid(
            vz_agg,
            family,
            "vz_signed",
            r"Axial velocity, $v_z$ ($\mu$m/s)",
            OUT / f"{family}_velocity_timecourse_early_zoom_1to16.png",
            xlim=(0.0, 2.0),
            shade_phases=True,
        ),
        plot_phase_grid(
            phase_agg,
            family,
            "phase_mean_vz",
            r"Phase-averaged $v_z$ ($\mu$m/s)",
            OUT / f"{family}_phase_mean_velocity_1to16.png",
        ),
        plot_phase_grid(
            phase_agg,
            family,
            "phase_delta_z",
            r"Phase displacement, $\Delta z$ ($\mu$m)",
            OUT / f"{family}_phase_displacement_1to16.png",
        ),
        plot_timecourse_grid(
            disp_agg,
            family,
            "delta_z",
            r"Cumulative displacement, $\Delta z$ ($\mu$m)",
            OUT / f"{family}_cumulative_displacement_timecourse_1to16.png",
            shade_phases=True,
        ),
    ]
    return outputs


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

    outputs: list[Path] = []
    for family in ("mpc12", "total480"):
        outputs.extend(process_family(family))
    write_readme(outputs)

    print(OUT)
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
