#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path('/home/thekaka/project/cytosim')
RANDOM_COND = ROOT / 'analysis/results/fixed_global_story/transport_regime/condition_metrics.csv'
RANDOM_RUN = ROOT / 'analysis/results/fixed_global_story/transport_regime/run_metrics.csv'
ALIGNED_COND = ROOT / 'analysis/results/init_minusz_comparison/transport_regime/condition_metrics.csv'
ALIGNED_RUN = ROOT / 'analysis/results/init_minusz_comparison/transport_regime/run_metrics.csv'
OUT = ROOT / 'analysis/results/full_story_2026-04-12/four_way_comparison'

XLINK_ORDER = ['1:4', '1:8', '1:16']
FAMILIES = [
    ('total480', 'Cluster count (total motors = 480)'),
    ('mpc12', 'Cluster count (12 motors per cluster)'),
]
CASE_ORDER = {
    'total480': [('c10_m48', 10), ('c20_m24', 20), ('c40_m12', 40), ('c60_m8', 60), ('c80_m6', 80)],
    'mpc12': [('c10_m12', 10), ('c20_m12', 20), ('c40_m12', 40), ('c60_m12', 60), ('c80_m12', 80)],
}
SERIES = [
    ('random', 'rotatable', 'Random rotatable', '#222222', 's', '-', 'white'),
    ('random', 'fixed_global', 'Random fixed global', '#d95f02', 'o', '-', '#d95f02'),
    ('aligned', 'rotatable', 'Aligned rotatable', '#1b9e77', '^', '--', 'white'),
    ('aligned', 'fixed_global', 'Aligned fixed global', '#7570b3', 'D', '--', '#7570b3'),
]
METRICS = [
    {'field': 'mean_vz_abs', 'ylabel': 'Mean |v_z| (um/s)', 'filename': 'mean_vz_abs.png', 'clip_zero': True, 'zero_line': False},
    {'field': 'auc_vz_abs', 'ylabel': 'AUC(|v_z|) (um)', 'filename': 'auc_vz_abs.png', 'clip_zero': True, 'zero_line': False},
    {'field': 'active_transport_fraction', 'ylabel': 'Active transport fraction', 'filename': 'active_transport_fraction.png', 'clip_zero': True, 'zero_line': False},
    {'field': 'abs_net_z_displacement', 'ylabel': '|Net z shift| (um)', 'filename': 'abs_net_z_displacement.png', 'clip_zero': True, 'zero_line': False},
    {'field': 'directionality_index', 'ylabel': 'Directionality index |Δz| / AUC(|v_z|)', 'filename': 'directionality_index.png', 'clip_zero': True, 'zero_line': False},
    {'field': 'mean_vz_signed', 'ylabel': 'Mean v_z (um/s)', 'filename': 'mean_vz_signed.png', 'clip_zero': False, 'zero_line': True},
]
TARGET_FINAL_S = 800.0


def load_rows(path: Path, orientation: str) -> list[dict]:
    rows = [dict(row) for row in csv.DictReader(path.open(encoding='utf-8'))]
    for row in rows:
        row['orientation'] = orientation
    return rows


def key(row: dict) -> tuple[str, str, str, str, str]:
    return (row['orientation'], row['model'], row['group'], row['case'], row['xlink_regime'])


def finite_float(value: str | None) -> float:
    if value in (None, ''):
        return float('nan')
    try:
        return float(value)
    except Exception:
        return float('nan')


def build_condition_lookup(cond_rows: list[dict], run_rows: list[dict]) -> dict[tuple[str, str, str, str, str], dict]:
    lookup = {key(row): dict(row) for row in cond_rows}

    grouped: dict[tuple[str, str, str, str, str], list[dict]] = defaultdict(list)
    for row in run_rows:
        grouped[key(row)].append(row)

    for k, row in lookup.items():
        net_mean = finite_float(row.get('net_z_displacement_mean'))
        net_sem = finite_float(row.get('net_z_displacement_sem'))
        row['mean_vz_signed_mean'] = net_mean / TARGET_FINAL_S if np.isfinite(net_mean) else float('nan')
        row['mean_vz_signed_sem'] = net_sem / TARGET_FINAL_S if np.isfinite(net_sem) else float('nan')

        values = []
        for run in grouped.get(k, []):
            auc = finite_float(run.get('auc_vz_abs'))
            disp = finite_float(run.get('abs_net_z_displacement'))
            if np.isfinite(auc) and auc > 0 and np.isfinite(disp):
                values.append(disp / auc)
        if values:
            arr = np.asarray(values, dtype=float)
            row['directionality_index_mean'] = float(np.nanmean(arr))
            row['directionality_index_sem'] = float(np.nanstd(arr, ddof=1) / math.sqrt(arr.size)) if arr.size > 1 else 0.0
        else:
            row['directionality_index_mean'] = float('nan')
            row['directionality_index_sem'] = float('nan')
    return lookup


def finite_bounds(values: list[float]) -> tuple[float, float]:
    vals = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    if vals.size == 0:
        return 0.0, 1.0
    lo = float(vals.min())
    hi = float(vals.max())
    pad = 0.10 * (hi - lo) if hi > lo else 0.10 * max(abs(hi), 1.0)
    return lo - pad, hi + pad


def plot_metric(lookup: dict, metric: dict) -> Path:
    plt.rcParams.update({
        'font.size': 12,
        'axes.linewidth': 1.2,
        'savefig.facecolor': 'white',
        'figure.facecolor': 'white',
    })
    fig, axes = plt.subplots(3, 2, figsize=(15.5, 14.0), sharex='col', constrained_layout=False)
    mean_key = f"{metric['field']}_mean"
    sem_key = f"{metric['field']}_sem"

    all_vals = []
    for orientation, model, _, _, _, _, _ in SERIES:
        for family, _ in FAMILIES:
            for regime in XLINK_ORDER:
                for case, _x in CASE_ORDER[family]:
                    row = lookup.get((orientation, model, family, case, regime))
                    if not row:
                        continue
                    m = finite_float(row.get(mean_key))
                    s = finite_float(row.get(sem_key))
                    if np.isfinite(m):
                        all_vals.append(m - s if np.isfinite(s) else m)
                        all_vals.append(m + s if np.isfinite(s) else m)
    ymin, ymax = finite_bounds(all_vals)
    if metric['clip_zero']:
        ymin = max(0.0, ymin)
    if metric['field'] == 'directionality_index':
        ymax = min(max(ymax, 1.0), 1.05)

    for ridx, regime in enumerate(XLINK_ORDER):
        for cidx, (family, xlabel) in enumerate(FAMILIES):
            ax = axes[ridx, cidx]
            for orientation, model, label, color, marker, linestyle, mfc in SERIES:
                xs = []
                ys = []
                es = []
                for case, xpos in CASE_ORDER[family]:
                    row = lookup.get((orientation, model, family, case, regime))
                    if not row:
                        continue
                    xs.append(xpos)
                    ys.append(finite_float(row.get(mean_key)))
                    es.append(finite_float(row.get(sem_key)))
                if xs:
                    ax.errorbar(
                        xs,
                        ys,
                        yerr=es,
                        color=color,
                        marker=marker,
                        linestyle=linestyle,
                        linewidth=2.2,
                        markersize=7.5,
                        markerfacecolor=mfc,
                        markeredgewidth=1.8,
                        capsize=3.5,
                        label=label if ridx == 0 and cidx == 0 else None,
                    )
            if metric['zero_line']:
                ax.axhline(0.0, color='#999999', linewidth=1.0, linestyle='--', zorder=0)
            ax.text(0.03, 0.94, regime, transform=ax.transAxes, ha='left', va='top', fontsize=16, color='#333333')
            if ridx == 0:
                ax.set_title('Total motors = 480' if family == 'total480' else '12 motors per cluster', fontsize=17, pad=10)
            if cidx == 0:
                ax.set_ylabel(metric['ylabel'], fontsize=17)
            if ridx == 2:
                ax.set_xlabel(xlabel, fontsize=16)
            ax.set_ylim(ymin, ymax)
            ax.tick_params(direction='out', width=1.1, labelsize=13)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(axis='y', color='#dddddd', linewidth=0.8)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=2, loc='lower center', bbox_to_anchor=(0.5, 0.02), frameon=False, fontsize=14)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.95, bottom=0.12, wspace=0.16, hspace=0.20)
    OUT.mkdir(parents=True, exist_ok=True)
    outpath = OUT / metric['filename']
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return outpath


def main() -> None:
    cond_rows = load_rows(RANDOM_COND, 'random') + load_rows(ALIGNED_COND, 'aligned')
    run_rows = load_rows(RANDOM_RUN, 'random') + load_rows(ALIGNED_RUN, 'aligned')
    lookup = build_condition_lookup(cond_rows, run_rows)
    outputs = [plot_metric(lookup, metric) for metric in METRICS]
    for path in outputs:
        print(path)


if __name__ == '__main__':
    main()
