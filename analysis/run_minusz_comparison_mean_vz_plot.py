#!/usr/bin/env python3
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

root = Path('/home/thekaka/project/cytosim')
cond = root / 'analysis/results/init_minusz_comparison/transport_regime/condition_metrics.csv'
out = root / 'analysis/results/init_minusz_comparison/transport_regime/plots/mean_vz_abs.png'

XLINK_ORDER = ['1:4', '1:8', '1:16']
TOTAL480_ORDER = [('c0_m0', 0), ('c10_m48', 10), ('c20_m24', 20), ('c40_m12', 40), ('c60_m8', 60), ('c80_m6', 80)]
MPC12_ORDER = [('c0_m0', 0), ('c10_m12', 10), ('c20_m12', 20), ('c40_m12', 40), ('c60_m12', 60), ('c80_m12', 80)]
MODEL_SPECS = [
    ('rotatable', 'Rotatable', '#222222', 's', 'white'),
    ('fixed_global', 'Fixed global', '#d95f02', 'o', '#d95f02'),
]
CONTROL_SPEC = ('control', 'No motors', '#7f7f7f', 'D', '#d9d9d9')

rows = list(csv.DictReader(cond.open(encoding='utf-8')))
lookup = {(r['model'], r['group'], r['case'], r['xlink_regime']): r for r in rows}


def val(model, group, case, regime, key):
    row = lookup.get((model, group, case, regime))
    if row is None:
        return np.nan
    raw = row.get(key, '')
    return float(raw) if raw not in ('', None) else np.nan


vals = []
for regime in XLINK_ORDER:
    for model, *_ in MODEL_SPECS:
        for case, _ in TOTAL480_ORDER[1:]:
            mean = val(model, 'total480', case, regime, 'mean_vz_abs_mean')
            sem = val(model, 'total480', case, regime, 'mean_vz_abs_sem')
            if np.isfinite(mean):
                vals += [mean - (sem if np.isfinite(sem) else 0.0), mean + (sem if np.isfinite(sem) else 0.0)]
        for case, _ in MPC12_ORDER[1:]:
            mean = val(model, 'mpc12', case, regime, 'mean_vz_abs_mean')
            sem = val(model, 'mpc12', case, regime, 'mean_vz_abs_sem')
            if np.isfinite(mean):
                vals += [mean - (sem if np.isfinite(sem) else 0.0), mean + (sem if np.isfinite(sem) else 0.0)]
    cmean = val('control', 'control', 'c0_m0', regime, 'mean_vz_abs_mean')
    csem = val('control', 'control', 'c0_m0', regime, 'mean_vz_abs_sem')
    if np.isfinite(cmean):
        vals += [cmean - (csem if np.isfinite(csem) else 0.0), cmean + (csem if np.isfinite(csem) else 0.0)]

ylo = max(0.0, min(vals) - 0.0003) if vals else 0.0
yhi = max(vals) + 0.0004 if vals else 0.01

plt.rcParams.update({'font.size': 14, 'axes.linewidth': 1.3, 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
fig, axes = plt.subplots(3, 2, figsize=(12.5, 12.0), sharey=True)

for i, regime in enumerate(XLINK_ORDER):
    for j, (group, order, xlabel) in enumerate([
        ('total480', TOTAL480_ORDER, 'Cluster count (0 = no motors; total motors = 480)'),
        ('mpc12', MPC12_ORDER, 'Cluster count (0 = no motors; 12 motors per cluster)'),
    ]):
        ax = axes[i, j]
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.3)
        ax.spines['bottom'].set_linewidth(1.3)
        ax.tick_params(direction='out', width=1.2, labelsize=13)
        ax.grid(axis='y', color='#d7d7d7', linewidth=0.7, alpha=0.85)

        cmean = val('control', 'control', 'c0_m0', regime, 'mean_vz_abs_mean')
        csem = val('control', 'control', 'c0_m0', regime, 'mean_vz_abs_sem')
        if np.isfinite(cmean):
            ax.errorbar([0], [cmean], yerr=[0.0 if not np.isfinite(csem) else csem], color=CONTROL_SPEC[2], marker=CONTROL_SPEC[3], markersize=7, lw=0, capsize=4, markerfacecolor=CONTROL_SPEC[4], markeredgewidth=1.6)

        xs = [x for _, x in order if x > 0]
        for model, label, color, marker, face in MODEL_SPECS:
            ys = [val(model, group, case, regime, 'mean_vz_abs_mean') for case, x in order if x > 0]
            es = [val(model, group, case, regime, 'mean_vz_abs_sem') for case, x in order if x > 0]
            ax.errorbar(xs, ys, yerr=es, color=color, marker=marker, markersize=7, lw=2.2, capsize=4, markerfacecolor=face, markeredgewidth=1.6)

        ax.set_xticks([x for _, x in order])
        ax.set_xticklabels([str(x) for _, x in order])
        ax.set_ylim(ylo, yhi)
        ax.text(0.03, 0.94, regime, transform=ax.transAxes, ha='left', va='top', fontsize=15, color='#444444')
        if j == 0:
            ax.set_ylabel('Mean |v_z| (um/s)', fontsize=16)
        if i == 2:
            ax.set_xlabel(xlabel, fontsize=15)

handles = [
    plt.Line2D([0], [0], color=CONTROL_SPEC[2], lw=0, marker=CONTROL_SPEC[3], markersize=7, markerfacecolor=CONTROL_SPEC[4], markeredgewidth=1.6, label='No motors'),
    plt.Line2D([0], [0], color='#222222', lw=2.2, marker='s', markersize=7, markerfacecolor='white', markeredgewidth=1.6, label='Rotatable'),
    plt.Line2D([0], [0], color='#d95f02', lw=2.2, marker='o', markersize=7, markerfacecolor='#d95f02', markeredgewidth=1.6, label='Fixed global'),
]
fig.legend(handles=handles, labels=['No motors', 'Rotatable', 'Fixed global'], loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=15)
fig.subplots_adjust(left=0.11, right=0.99, top=0.98, bottom=0.11, wspace=0.18, hspace=0.22)
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=300, bbox_inches='tight')
print(out)
