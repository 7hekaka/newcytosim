#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
BASE = ROOT / 'analysis' / 'results' / 'fixed_global_story' / 'transport_regime' / 'condition_metrics.csv'
MOTOR = ROOT / 'analysis' / 'results' / 'motor_force_sweep' / 'condition_metrics.csv'
OUT = ROOT / 'analysis' / 'results' / 'fixed_global_story' / 'total_motor_sweep'

XLINK_ORDER = ['1:4', '1:8', '1:16']
SELECTION = [
    {'source': 'base', 'model': 'control', 'group': 'control', 'case': 'c0_m0', 'total_motors': 0, 'label': '0', 'marker': 'D', 'edge': '#7f7f7f', 'face': '#d9d9d9'},
    {'source': 'base', 'model': 'fixed_global', 'group': 'mpc12', 'case': 'c10_m12', 'total_motors': 120, 'label': '120', 'marker': 'o', 'edge': '#d95f02', 'face': '#d95f02'},
    {'source': 'base', 'model': 'fixed_global', 'group': 'mpc12', 'case': 'c20_m12', 'total_motors': 240, 'label': '240', 'marker': 'o', 'edge': '#d95f02', 'face': '#d95f02'},
    {'source': 'base', 'model': 'fixed_global', 'group': 'mpc12', 'case': 'c40_m12', 'total_motors': 480, 'label': '480', 'marker': 'o', 'edge': '#d95f02', 'face': '#d95f02'},
    {'source': 'base', 'model': 'fixed_global', 'group': 'mpc12', 'case': 'c80_m12', 'total_motors': 960, 'label': '960', 'marker': 'o', 'edge': '#d95f02', 'face': '#d95f02'},
    {'source': 'motor', 'model': 'fixed_global_motor_force', 'group': 'c40_m50', 'case': 'uf2p5', 'total_motors': 2000, 'label': '2000', 'marker': 's', 'edge': '#a63d00', 'face': '#f28e2b'},
    {'source': 'motor', 'model': 'fixed_global_motor_force', 'group': 'c40_m100', 'case': 'uf2p5', 'total_motors': 4000, 'label': '4000', 'marker': '^', 'edge': '#7f2704', 'face': '#fdae6b'},
]
METRICS = [
    {'field': 'mean_vz_abs', 'ylabel': 'Mean |v_z| (um/s)', 'slug': 'mean_vz_abs'},
    {'field': 'peak_vz_abs', 'ylabel': 'Peak |v_z| (um/s)', 'slug': 'peak_vz_abs'},
    {'field': 'auc_vz_abs', 'ylabel': 'AUC(|v_z|) (um)', 'slug': 'auc_vz_abs'},
    {'field': 'active_transport_fraction', 'ylabel': 'Active transport fraction', 'slug': 'active_transport_fraction'},
]


def load_lookup(path: Path) -> dict[tuple[str, str, str, str], dict]:
    rows = [dict(row) for row in csv.DictReader(path.open(encoding='utf-8'))]
    return {(r['model'], r['group'], r['case'], r['xlink_regime']): r for r in rows}


def style(ax) -> None:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.3)
    ax.spines['bottom'].set_linewidth(1.3)
    ax.tick_params(direction='out', width=1.2, labelsize=12)
    ax.grid(axis='y', color='#d9d9d9', linewidth=0.7, alpha=0.85)


def build_rows(base_lookup: dict, motor_lookup: dict) -> list[dict]:
    rows = []
    for regime in XLINK_ORDER:
        for spec in SELECTION:
            lookup = base_lookup if spec['source'] == 'base' else motor_lookup
            row = lookup.get((spec['model'], spec['group'], spec['case'], regime))
            if row is None:
                continue
            rows.append({**spec, 'xlink_regime': regime, **row})
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = []
    for row in rows:
        for k in row:
            if k not in keys:
                keys.append(k)
    with path.open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def plot_metric(rows: list[dict], metric: dict) -> Path:
    fig, axes = plt.subplots(3, 1, figsize=(9.2, 11.5), sharex=True)
    for ax, regime in zip(axes, XLINK_ORDER):
        sub = [r for r in rows if r['xlink_regime'] == regime]
        sub.sort(key=lambda r: float(r['total_motors']))
        xs = [float(r['total_motors']) for r in sub]
        ys = [float(r[f"{metric['field']}_mean"]) for r in sub]
        es = [float(r[f"{metric['field']}_sem"]) for r in sub]
        ax.plot(xs, ys, color='#d95f02', lw=2.0, zorder=2)
        for r, x, y, e in zip(sub, xs, ys, es):
            ax.errorbar([x], [y], yerr=[e], color=r['edge'], marker=r['marker'], markersize=8, lw=0, capsize=3.5, markeredgewidth=1.5, markeredgecolor=r['edge'], markerfacecolor=r['face'], zorder=3)
        ax.text(0.02, 0.92, regime, transform=ax.transAxes, ha='left', va='top', fontsize=13)
        ax.set_ylabel(metric['ylabel'], fontsize=14)
        style(ax)
    axes[-1].set_xlabel('Total motors in the system', fontsize=15)
    axes[-1].set_xticks([0, 120, 240, 480, 960, 2000, 4000])
    axes[-1].set_xticklabels(['0', '120', '240', '480', '960', '2000', '4000'], fontsize=12)
    fig.subplots_adjust(left=0.16, right=0.98, top=0.98, bottom=0.09, hspace=0.16)
    path = OUT / 'plots' / f"{metric['slug']}.png"
    fig.savefig(path, dpi=240, bbox_inches='tight')
    plt.close(fig)
    return path


def save_legend() -> Path:
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker='D', color='none', markeredgecolor='#7f7f7f', markerfacecolor='#d9d9d9', markeredgewidth=1.4, markersize=8, label='No motors'),
        Line2D([0], [0], marker='o', color='none', markeredgecolor='#d95f02', markerfacecolor='#d95f02', markeredgewidth=1.4, markersize=8, label='Baseline fixed-global (12 motors/cluster)'),
        Line2D([0], [0], marker='s', color='none', markeredgecolor='#a63d00', markerfacecolor='#f28e2b', markeredgewidth=1.4, markersize=8, label='40 clusters, 50 motors/cluster'),
        Line2D([0], [0], marker='^', color='none', markeredgecolor='#7f2704', markerfacecolor='#fdae6b', markeredgewidth=1.4, markersize=8, label='40 clusters, 100 motors/cluster'),
    ]
    fig = plt.figure(figsize=(8.5, 1.2), facecolor='white')
    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
    ax.axis('off')
    ax.legend(handles=handles, loc='center', ncol=2, frameon=False, fontsize=12, columnspacing=1.8, handletextpad=0.6)
    path = OUT / 'plots' / 'legend_only.png'
    fig.savefig(path, dpi=220, bbox_inches='tight')
    plt.close(fig)
    return path


def main() -> None:
    plt.rcParams.update({'font.size': 12, 'axes.linewidth': 1.2, 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
    (OUT / 'plots').mkdir(parents=True, exist_ok=True)
    base_lookup = load_lookup(BASE)
    motor_lookup = load_lookup(MOTOR)
    rows = build_rows(base_lookup, motor_lookup)
    written = []
    write_csv(OUT / 'selected_conditions.csv', rows)
    written.append(OUT / 'selected_conditions.csv')
    for metric in METRICS:
        written.append(plot_metric(rows, metric))
    written.append(save_legend())
    for path in written:
        print(path)


if __name__ == '__main__':
    main()
