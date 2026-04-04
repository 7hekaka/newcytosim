#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import plot_timecourse_panels as base
import run_speckle_timecourse_batch as speck

ROOT = Path(__file__).resolve().parent.parent
ROT = ROOT / 'clu' / 'analysis_full' / 'condition_summary.csv'
FIX = ROOT / 'clu_fixed_global' / 'analysis_full' / 'condition_summary.csv'
CONTROL = ROOT / 'control' / '2'
OUT = ROOT / 'comparison_plots_presentation'

XLINK_ORDER = ['1:4', '1:8', '1:16']
TOTAL480_ORDER = [('c10_m48', 10), ('c20_m24', 20), ('c40_m12', 40), ('c60_m8', 60), ('c80_m6', 80)]
MPC12_ORDER = [('c10_m12', 10), ('c20_m12', 20), ('c40_m12', 40), ('c60_m12', 60), ('c80_m12', 80)]
MODEL_SPECS = [
    ('rotatable', 'Rotatable', '#222222', 's', 'white'),
    ('fixed_global', 'Fixed global', '#d95f02', 'o', '#d95f02'),
]
CONTROL_SPEC = ('control', 'No motors', '#7f7f7f', 'D', '#d9d9d9')
MATCH_FINAL_MIN = 800.0 / 60.0


def load_rows(path: Path, model: str):
    rows = list(csv.DictReader(path.open(encoding='utf-8')))
    for row in rows:
        row['model'] = model
    return rows


def as_float(row, key):
    if row is None:
        return float('nan')
    v = row.get(key, '')
    return float(v) if v not in ('', None) else float('nan')


def build_lookup(rows):
    return {(r['model'], r['group'], r['case'], r['xlink_regime']): r for r in rows}


def control_metadata() -> list[dict]:
    rows = []
    for idx in range(1, 31):
        if idx <= 10:
            regime = '1:4'
        elif idx <= 20:
            regime = '1:8'
        else:
            regime = '1:16'
        run_dir = f'r{idx:04d}'
        run_path = CONTROL / run_dir
        rows.append(
            {
                'model': 'control',
                'group': 'control',
                'case': 'c0',
                'xlink_regime': regime,
                'run_dir': run_dir,
                'run_path': str(run_path),
            }
        )
    return rows


def ensure_point(run_path: Path, point_path: Path, report_bin: Path) -> None:
    if point_path.exists() and point_path.stat().st_size > 0:
        return
    with point_path.open('w', encoding='utf-8') as fh:
        subprocess.run(
            [str(report_bin), 'fiber:point'],
            cwd=run_path,
            stdout=fh,
            stderr=subprocess.DEVNULL,
            check=True,
        )


def analyze_control_task(task: tuple[dict, str]) -> dict:
    row, report_bin_s = task
    run_path = Path(row['run_path'])
    ensure_point(run_path, run_path / 'point.txt', Path(report_bin_s))
    return base.analyze_run((row, 8, 72, 80))


def summarize_control(run_rows: list[dict]) -> list[dict]:
    grouped = defaultdict(list)
    for row in run_rows:
        if row.get('status') != 'ok':
            continue
        arr = np.asarray(row['vz_abs'], dtype=float)
        time = np.asarray(row['time_min'], dtype=float)
        arr = arr[time <= MATCH_FINAL_MIN + 1e-9]
        if arr.size == 0:
            continue
        grouped[row['xlink_regime']].append(float(np.nanmean(arr)))
    out = []
    for regime in XLINK_ORDER:
        vals = np.asarray(grouped.get(regime, []), dtype=float)
        if vals.size == 0:
            out.append({'xlink_regime': regime, 'n_runs_used': 0, 'vz_abs_mean_mean': float('nan'), 'vz_abs_mean_sem': float('nan')})
            continue
        sem = float(np.nanstd(vals, ddof=1) / math.sqrt(vals.size)) if vals.size > 1 else 0.0
        out.append({'xlink_regime': regime, 'n_runs_used': int(vals.size), 'vz_abs_mean_mean': float(np.nanmean(vals)), 'vz_abs_mean_sem': sem})
    return out


def compute_ylim(lookup, control_summary):
    vals = []
    for row in lookup.values():
        m = as_float(row, 'vz_abs_mean_mean')
        s = as_float(row, 'vz_abs_mean_sem')
        if np.isfinite(m):
            vals.extend([m - s, m + s])
    for row in control_summary:
        m = row['vz_abs_mean_mean']
        s = row['vz_abs_mean_sem']
        if np.isfinite(m):
            vals.extend([m - s, m + s])
    lo = max(0.0, min(vals) - 0.0003)
    hi = max(vals) + 0.00045
    return lo, hi


def style(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.4)
    ax.spines['bottom'].set_linewidth(1.4)
    ax.tick_params(direction='out', width=1.3, labelsize=15)
    ax.grid(axis='y', color='#d7d7d7', linewidth=0.7, alpha=0.85)


def draw_family(ax, family_order, lookup, control_lookup, group, regime):
    ctrl = control_lookup[regime]
    c_model, _, c_color, c_marker, c_face = CONTROL_SPEC
    ax.errorbar([0], [ctrl['vz_abs_mean_mean']], yerr=[ctrl['vz_abs_mean_sem']], color=c_color, marker=c_marker,
                markersize=8, lw=0, capsize=4, markeredgewidth=1.8, markeredgecolor=c_color, markerfacecolor=c_face, zorder=5)
    for model, _, color, marker, face in MODEL_SPECS:
        xs = [x for _, x in family_order]
        ys = []
        es = []
        for case, _ in family_order:
            row = lookup.get((model, group, case, regime))
            ys.append(as_float(row, 'vz_abs_mean_mean'))
            es.append(as_float(row, 'vz_abs_mean_sem'))
        ax.errorbar(xs, ys, yerr=es, color=color, marker=marker, markersize=8, lw=2.3, capsize=4,
                    markeredgewidth=1.8, markeredgecolor=color, markerfacecolor=face)
    ax.set_xticks([0] + [x for _, x in family_order])
    ax.set_xticklabels(['0'] + [str(x) for _, x in family_order], fontsize=15)
    style(ax)


def main():
    rows = load_rows(ROT, 'rotatable') + load_rows(FIX, 'fixed_global')
    lookup = build_lookup(rows)

    report_bin = speck.pick_report_bin()
    control_tasks = [(row, str(report_bin)) for row in control_metadata()]
    with ProcessPoolExecutor(max_workers=4) as ex:
        control_run_rows = list(ex.map(analyze_control_task, control_tasks))
    control_summary = summarize_control(control_run_rows)
    control_lookup = {row['xlink_regime']: row for row in control_summary}

    ylo, yhi = compute_ylim(lookup, control_summary)

    plt.rcParams.update({'font.size': 15, 'axes.linewidth': 1.3, 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
    fig, axes = plt.subplots(3, 2, figsize=(12.8, 13.2), sharey=True)

    for ridx, regime in enumerate(XLINK_ORDER):
        draw_family(axes[ridx, 0], TOTAL480_ORDER, lookup, control_lookup, 'total480', regime)
        draw_family(axes[ridx, 1], MPC12_ORDER, lookup, control_lookup, 'mpc12', regime)
        for cidx in range(2):
            axes[ridx, cidx].set_ylim(ylo, yhi)
            axes[ridx, cidx].text(0.03, 0.94, regime, transform=axes[ridx, cidx].transAxes, ha='left', va='top', fontsize=16, color='#444444')

    for ridx in range(3):
        axes[ridx, 0].set_ylabel('Mean |v_z| (um/s)', fontsize=18)
    axes[2, 0].set_xlabel('Cluster count (0 = no motors)\n(total motors = 480)', fontsize=18)
    axes[2, 1].set_xlabel('Cluster count (0 = no motors)\n(12 motors per cluster)', fontsize=18)

    handles = []
    labels = []
    for _, label, color, marker, face in [CONTROL_SPEC] + MODEL_SPECS:
        handles.append(plt.Line2D([0], [0], color=color, lw=2.3 if label != 'No motors' else 0, marker=marker, markersize=8, markeredgewidth=1.8, markeredgecolor=color, markerfacecolor=face))
        labels.append(label)
    fig.legend(handles, labels, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=17)
    fig.subplots_adjust(left=0.10, right=0.99, top=0.98, bottom=0.11, wspace=0.18, hspace=0.22)

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / 'vz_abs_mean_with_control.png'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)

    csv_path = OUT / 'control_velocity_summary.csv'
    with csv_path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=['xlink_regime', 'n_runs_used', 'vz_abs_mean_mean', 'vz_abs_mean_sem'])
        writer.writeheader()
        writer.writerows(control_summary)
    print(csv_path)
    print(out)


if __name__ == '__main__':
    main()

