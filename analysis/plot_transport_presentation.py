#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
ROT = ROOT / 'clu' / 'analysis_full' / 'condition_summary.csv'
FIX = ROOT / 'clu_fixed_global' / 'analysis_full' / 'condition_summary.csv'
OUT = ROOT / 'comparison_plots_presentation'

XLINK_ORDER = ['1:4', '1:8', '1:16']
TOTAL480_ORDER = [('c10_m48', 10), ('c20_m24', 20), ('c40_m12', 40), ('c60_m8', 60), ('c80_m6', 80)]
MPC12_ORDER = [('c10_m12', 10), ('c20_m12', 20), ('c40_m12', 40), ('c60_m12', 60), ('c80_m12', 80)]
FOCUS_ORDER = [('parallel_only', 'parallel'), ('parallel_antiparallel_50_50', '50:50'), ('offrate_0p5', 'off-rate 0.5'), ('offrate_1p0', 'off-rate 1.0')]
MODEL_SPECS = [
    ('rotatable', 'Rotatable', '#222222', 's', 'white'),
    ('fixed_global', 'Fixed global', '#d95f02', 'o', '#d95f02'),
]

METRICS = [
    {
        'field': 'vz_std',
        'ylabel': 'Axial speed variability (um/s)',
        'slug': 'vz_std_presentation',
        'clip_zero': True,
    },
    {
        'field': 'Dz_max',
        'ylabel': 'Peak depletion index',
        'slug': 'depletion_peak_presentation',
        'clip_zero': True,
    },
    {
        'field': 'swirl_rate_abs_mean',
        'ylabel': 'Mean |swirl rate| (rad/s)',
        'slug': 'swirl_abs_presentation',
        'clip_zero': True,
    },
    {
        'field': 'nematic_xy_final',
        'ylabel': 'Final XY nematic order',
        'slug': 'nematic_final_presentation',
        'clip_zero': True,
    },
]


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


def compute_ylim(lookup, field, clip_zero):
    vals = []
    for row in lookup.values():
        m = as_float(row, f'{field}_mean')
        s = as_float(row, f'{field}_sem')
        if np.isfinite(m):
            if np.isfinite(s):
                vals.extend([m - s, m + s])
            else:
                vals.append(m)
    if not vals:
        return 0.0, 1.0
    lo = min(vals)
    hi = max(vals)
    pad = 0.10 * (hi - lo) if hi > lo else 0.10 * max(abs(hi), 1.0)
    lo -= pad
    hi += pad
    if clip_zero:
        lo = max(0.0, lo)
    return lo, hi


def style(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.4)
    ax.spines['bottom'].set_linewidth(1.4)
    ax.tick_params(direction='out', width=1.3, labelsize=15)
    ax.grid(axis='y', color='#d7d7d7', linewidth=0.7, alpha=0.85)


def draw(ax, xs, case_order, lookup, group, regime, field, categorical=False):
    for model, _, color, marker, face in MODEL_SPECS:
        ys = []
        es = []
        for case in case_order:
            row = lookup.get((model, group, case, regime))
            ys.append(as_float(row, f'{field}_mean'))
            es.append(as_float(row, f'{field}_sem'))
        ax.errorbar(
            xs,
            ys,
            yerr=es,
            color=color,
            marker=marker,
            markersize=8,
            lw=2.3,
            capsize=4,
            markeredgewidth=1.8,
            markeredgecolor=color,
            markerfacecolor=face,
        )
    if categorical:
        ax.set_xticks(xs)
        ax.set_xticklabels([lab for _, lab in FOCUS_ORDER], rotation=18, ha='right', fontsize=15)
    else:
        ax.set_xticks(xs)
        ax.set_xticklabels([str(x) for x in xs], fontsize=15)
    style(ax)


def make_metric_figure(lookup, metric):
    plt.rcParams.update({'font.size': 15, 'axes.linewidth': 1.3, 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
    fig, axes = plt.subplots(3, 3, figsize=(16.8, 13.2), sharey=True)
    ylo, yhi = compute_ylim(lookup, metric['field'], metric['clip_zero'])

    for ridx, regime in enumerate(XLINK_ORDER):
        xs_num = [x for _, x in TOTAL480_ORDER]
        draw(axes[ridx, 0], xs_num, [c for c, _ in TOTAL480_ORDER], lookup, 'total480', regime, metric['field'], categorical=False)
        draw(axes[ridx, 1], xs_num, [c for c, _ in MPC12_ORDER], lookup, 'mpc12', regime, metric['field'], categorical=False)
        draw(axes[ridx, 2], np.arange(len(FOCUS_ORDER)), [c for c, _ in FOCUS_ORDER], lookup, 'focus40x12', regime, metric['field'], categorical=True)
        for cidx in range(3):
            axes[ridx, cidx].set_ylim(ylo, yhi)
            axes[ridx, cidx].text(0.03, 0.94, regime, transform=axes[ridx, cidx].transAxes, ha='left', va='top', fontsize=16, color='#444444')

    for ridx in range(3):
        axes[ridx, 0].set_ylabel(metric['ylabel'], fontsize=18)
    axes[2, 0].set_xlabel('Cluster count\n(total motors = 480)', fontsize=18)
    axes[2, 1].set_xlabel('Cluster count\n(12 motors per cluster)', fontsize=18)
    axes[2, 2].set_xlabel('Crosslinker condition\n(40 clusters, 12 motors/cluster)', fontsize=18)

    handles = []
    labels = []
    for _, label, color, marker, face in MODEL_SPECS:
        handles.append(plt.Line2D([0], [0], color=color, lw=2.3, marker=marker, markersize=8, markeredgewidth=1.8, markeredgecolor=color, markerfacecolor=face))
        labels.append(label)
    fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=17)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.11, wspace=0.22, hspace=0.22)
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"{metric['slug']}.png"
    fig.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out


def main():
    rows = load_rows(ROT, 'rotatable') + load_rows(FIX, 'fixed_global')
    lookup = build_lookup(rows)
    outs = [make_metric_figure(lookup, metric) for metric in METRICS]
    for out in outs:
        print(out)


if __name__ == '__main__':
    main()
