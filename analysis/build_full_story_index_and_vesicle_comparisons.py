#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import shutil
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path('/home/thekaka/project/cytosim')
BUNDLE = ROOT / 'analysis/results/full_story_2026-04-12'
VESICLE_COND = ROOT / 'vesicles/vesicle_transport_out_postprocessed/condition_metrics.csv'
VESICLE_PLOTS = ROOT / 'vesicles/vesicle_transport_out_postprocessed/plots'
RANDOM_COND = ROOT / 'analysis/results/fixed_global_story/transport_regime/condition_metrics.csv'
RANDOM_RUN = ROOT / 'analysis/results/fixed_global_story/transport_regime/run_metrics.csv'
ALIGNED_COND = ROOT / 'analysis/results/init_minusz_comparison/transport_regime/condition_metrics.csv'
ALIGNED_RUN = ROOT / 'analysis/results/init_minusz_comparison/transport_regime/run_metrics.csv'
TARGET_FINAL_S = 800.0
VES_LIST = [20, 50, 100]
VIS_LIST = [1, 2, 4]
XLINKS = [4, 8, 16]
METRICS = [
    {'field': 'mean_vz_abs', 'ylabel': 'Mean |v_z| (um/s)', 'filename': 'mean_vz_abs.png', 'clip_zero': True, 'zero_line': False},
    {'field': 'auc_vz_abs', 'ylabel': 'AUC(|v_z|) (um)', 'filename': 'auc_vz_abs.png', 'clip_zero': True, 'zero_line': False},
    {'field': 'active_transport_fraction', 'ylabel': 'Active transport fraction', 'filename': 'active_transport_fraction.png', 'clip_zero': True, 'zero_line': False},
    {'field': 'mean_vz_signed', 'ylabel': 'Mean v_z (um/s)', 'filename': 'mean_vz_signed.png', 'clip_zero': False, 'zero_line': True},
    {'field': 'directionality_index', 'ylabel': 'Directionality index |Δz| / AUC(|v_z|)', 'filename': 'directionality_index.png', 'clip_zero': True, 'zero_line': False},
]
THEMES = {
    'vesicle_vs_rotatable': {
        'title': 'Vesicle vs rotatable wall-bound reference',
        'refs': [
            ('random', 'rotatable', 'Random rotatable', '#222222', 's', '-'),
            ('aligned', 'rotatable', 'Aligned rotatable', '#1b9e77', '^', '--'),
        ],
    },
    'vesicle_vs_fixed': {
        'title': 'Vesicle vs fixed wall-bound reference',
        'refs': [
            ('random', 'fixed_global', 'Random fixed global', '#d95f02', 'o', '-'),
            ('aligned', 'fixed_global', 'Aligned fixed global', '#7570b3', 'D', '--'),
        ],
    },
    'vesicle_vs_aligned': {
        'title': 'Vesicle vs aligned wall-bound reference',
        'refs': [
            ('aligned', 'rotatable', 'Aligned rotatable', '#1b9e77', '^', '--'),
            ('aligned', 'fixed_global', 'Aligned fixed global', '#7570b3', 'D', '--'),
        ],
    },
    'vesicle_vs_random': {
        'title': 'Vesicle vs random wall-bound reference',
        'refs': [
            ('random', 'rotatable', 'Random rotatable', '#222222', 's', '-'),
            ('random', 'fixed_global', 'Random fixed global', '#d95f02', 'o', '-'),
        ],
    },
}


def copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding='utf-8')


def read_csv(path: Path) -> list[dict]:
    return [dict(row) for row in csv.DictReader(path.open(encoding='utf-8', newline=''))]


def finite_float(value):
    if value in (None, ''):
        return float('nan')
    try:
        return float(value)
    except Exception:
        return float('nan')


def wall_key(row: dict) -> tuple[str, str, str, str]:
    return (row['orientation'], row['model'], row['group'], row['xlink_regime'])


def load_wall_benchmarks() -> dict[tuple[str, str, str, str], dict]:
    cond_rows = read_csv(RANDOM_COND) + read_csv(ALIGNED_COND)
    run_rows = read_csv(RANDOM_RUN) + read_csv(ALIGNED_RUN)
    for row in cond_rows[:len(read_csv(RANDOM_COND))]:
        row['orientation'] = 'random'
    for row in cond_rows[len(read_csv(RANDOM_COND)):]:
        row['orientation'] = 'aligned'
    for row in run_rows[:len(read_csv(RANDOM_RUN))]:
        row['orientation'] = 'random'
    for row in run_rows[len(read_csv(RANDOM_RUN)):]:
        row['orientation'] = 'aligned'

    grouped = defaultdict(list)
    for row in run_rows:
        if row['group'] == 'mpc12' and row['case'] == 'c40_m12':
            grouped[(row['orientation'], row['model'], row['group'], row['xlink_regime'])].append(row)

    lookup = {}
    for row in cond_rows:
        if row['group'] != 'mpc12' or row['case'] != 'c40_m12':
            continue
        row = dict(row)
        net_mean = finite_float(row.get('net_z_displacement_mean'))
        net_sem = finite_float(row.get('net_z_displacement_sem'))
        row['mean_vz_signed_mean'] = net_mean / TARGET_FINAL_S if np.isfinite(net_mean) else float('nan')
        row['mean_vz_signed_sem'] = net_sem / TARGET_FINAL_S if np.isfinite(net_sem) else float('nan')
        vals = []
        for run in grouped.get((row['orientation'], row['model'], row['group'], row['xlink_regime']), []):
            auc = finite_float(run.get('auc_vz_abs'))
            disp = finite_float(run.get('abs_net_z_displacement'))
            if np.isfinite(auc) and auc > 0 and np.isfinite(disp):
                vals.append(disp / auc)
        if vals:
            arr = np.asarray(vals, dtype=float)
            row['directionality_index_mean'] = float(np.nanmean(arr))
            row['directionality_index_sem'] = float(np.nanstd(arr, ddof=1) / math.sqrt(arr.size)) if arr.size > 1 else 0.0
        else:
            row['directionality_index_mean'] = float('nan')
            row['directionality_index_sem'] = float('nan')
        lookup[wall_key(row)] = row
    return lookup


def load_vesicle_lookup() -> dict[tuple[int, int, int], dict]:
    rows = read_csv(VESICLE_COND)
    return {(int(r['vis']), int(r['xlink']), int(r['ves'])): r for r in rows}


def compute_bounds(ves_lookup, wall_lookup, refs, metric):
    mean_key = f"{metric['field']}_mean"
    sem_key = f"{metric['field']}_sem"
    vals = []
    for vis in VIS_LIST:
        for xlink in XLINKS:
            for ves in VES_LIST:
                row = ves_lookup.get((vis, xlink, ves))
                if row:
                    m = finite_float(row.get(mean_key))
                    s = finite_float(row.get(sem_key))
                    if np.isfinite(m):
                        vals.extend([m - (s if np.isfinite(s) else 0.0), m + (s if np.isfinite(s) else 0.0)])
            regime = f'1:{xlink}'
            for orientation, model, *_ in refs:
                row = wall_lookup.get((orientation, model, 'mpc12', regime))
                if row:
                    m = finite_float(row.get(mean_key))
                    s = finite_float(row.get(sem_key))
                    if np.isfinite(m):
                        vals.extend([m - (s if np.isfinite(s) else 0.0), m + (s if np.isfinite(s) else 0.0)])
    vals = np.asarray([v for v in vals if np.isfinite(v)], dtype=float)
    if vals.size == 0:
        return 0.0, 1.0
    lo = float(vals.min())
    hi = float(vals.max())
    pad = 0.10 * (hi - lo) if hi > lo else 0.10 * max(abs(hi), 1.0)
    lo -= pad
    hi += pad
    if metric['clip_zero']:
        lo = max(0.0, lo)
    if metric['field'] == 'directionality_index':
        hi = min(max(hi, 1.0), 1.05)
    return lo, hi


def plot_theme(theme_key: str, theme: dict, ves_lookup: dict, wall_lookup: dict, metric: dict, outdir: Path) -> Path:
    plt.rcParams.update({'font.size': 12, 'axes.linewidth': 1.2, 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
    fig, axes = plt.subplots(3, 3, figsize=(16.0, 13.5), sharex=True, constrained_layout=False)
    mean_key = f"{metric['field']}_mean"
    sem_key = f"{metric['field']}_sem"
    ylo, yhi = compute_bounds(ves_lookup, wall_lookup, theme['refs'], metric)

    for ridx, xlink in enumerate(XLINKS):
        regime = f'1:{xlink}'
        for cidx, vis in enumerate(VIS_LIST):
            ax = axes[ridx, cidx]
            ves_y = []
            ves_e = []
            for ves in VES_LIST:
                row = ves_lookup.get((vis, xlink, ves))
                ves_y.append(finite_float(row.get(mean_key)) if row else float('nan'))
                ves_e.append(finite_float(row.get(sem_key)) if row else float('nan'))
            ax.errorbar(VES_LIST, ves_y, yerr=ves_e, color='#c05a00', marker='o', markersize=6.5, linewidth=2.2, capsize=3.5, markerfacecolor='#c05a00', markeredgewidth=1.4, label='Vesicle' if ridx == 0 and cidx == 0 else None)

            xspan = np.asarray([min(VES_LIST), max(VES_LIST)], dtype=float)
            for orientation, model, label, color, marker, linestyle in theme['refs']:
                row = wall_lookup.get((orientation, model, 'mpc12', regime))
                if not row:
                    continue
                mean = finite_float(row.get(mean_key))
                sem = finite_float(row.get(sem_key))
                if not np.isfinite(mean):
                    continue
                lo = mean - sem if np.isfinite(sem) else mean
                hi = mean + sem if np.isfinite(sem) else mean
                ax.fill_between(xspan, [lo, lo], [hi, hi], color=color, alpha=0.10)
                ax.plot(xspan, [mean, mean], color=color, linestyle=linestyle, linewidth=2.0, marker=marker, markersize=6.0, markevery=[0, 1], markerfacecolor='white' if model == 'rotatable' else color, markeredgewidth=1.4, label=label if ridx == 0 and cidx == 0 else None)

            if metric['zero_line']:
                ax.axhline(0.0, color='#999999', linewidth=1.0, linestyle='--', zorder=0)
            ax.text(0.03, 0.94, regime, transform=ax.transAxes, ha='left', va='top', fontsize=15, color='#444444')
            if ridx == 0:
                ax.set_title(f'vis={vis}', fontsize=16, pad=10)
            if cidx == 0:
                ax.set_ylabel(metric['ylabel'], fontsize=16)
            if ridx == 2:
                ax.set_xlabel('Vesicle count', fontsize=15)
            ax.set_ylim(ylo, yhi)
            ax.set_xticks(VES_LIST)
            ax.tick_params(direction='out', width=1.1, labelsize=12)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(axis='y', color='#dddddd', linewidth=0.8)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.02), ncol=3, frameon=False, fontsize=13)
    fig.suptitle(theme['title'], fontsize=18, y=0.98)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.93, bottom=0.11, wspace=0.18, hspace=0.22)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / metric['filename']
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return outpath


def main() -> None:
    BUNDLE.mkdir(parents=True, exist_ok=True)

    # top-level README
    write_text(BUNDLE / 'README.md', '''# Full Story Comparison Bundle

This folder is the presentation-facing index for the transport story.

## Folder map
- `random_orientation/`: wall-bound random-initialization comparisons between fixed-global and rotatable clusters.
- `aligned_minusz/`: wall-bound minus-z-aligned comparisons between fixed-global and rotatable clusters.
- `four_way_comparison/`: combined wall-bound comparison plots for random vs aligned and fixed vs rotatable.
- `vesicle_only/`: compact vesicle-only transport plots.
- `vesicle_cross_comparisons/`: vesicle curves compared against wall-bound reference conditions.
- `notes/`: short markdown summaries for slide building.

## Benchmark convention used in vesicle cross-comparisons
To avoid forcing incompatible x-axes together, the vesicle comparisons plot vesicle curves versus vesicle count and use wall-bound `c40_m12` as a horizontal benchmark reference. This is the canonical 40-cluster, 12-motors/cluster wall-bound condition and is identical between the `mpc12` and `total480` families.
''')

    # add README to existing subfolders
    write_text(BUNDLE / 'random_orientation/README.md', '''# Random Orientation

Wall-bound cases with random initial orientation. Use these plots for the original fixed-vs-rotatable story before polarity is imposed.
''')
    write_text(BUNDLE / 'aligned_minusz/README.md', '''# Aligned Minus-Z

Wall-bound cases initialized along minus-z. Use these plots to show what transport looks like once polarity is supplied at the start.
''')
    write_text(BUNDLE / 'four_way_comparison/README.md', '''# Four-Way Comparison

These plots compare random fixed, random rotatable, aligned fixed, and aligned rotatable on the same axes.
''')

    # vesicle-only folder
    ves_only = BUNDLE / 'vesicle_only' / 'transport_metrics'
    ves_only.mkdir(parents=True, exist_ok=True)
    for name in ['mean_vz_abs.png', 'auc_vz_abs.png', 'active_transport_fraction.png', 'mean_vz_signed.png', 'directionality_index.png']:
        src = VESICLE_PLOTS / name
        if src.exists():
            copy(src, ves_only / name)
    write_text(BUNDLE / 'vesicle_only/README.md', '''# Vesicle Only

These are the compact vesicle-only transport metrics generated from the postprocessed speckle bundle.
- `mean_vz_abs.png`: magnitude of axial transport.
- `auc_vz_abs.png`: accumulated axial transport activity.
- `active_transport_fraction.png`: fraction of time above the control-derived threshold.
- `mean_vz_signed.png`: signed axial velocity; values near zero indicate cancellation.
- `directionality_index.png`: one-way transport score defined as `|Δz| / AUC(|v_z|)`.
''')

    ves_lookup = load_vesicle_lookup()
    wall_lookup = load_wall_benchmarks()

    cross_root = BUNDLE / 'vesicle_cross_comparisons'
    write_text(cross_root / 'README.md', '''# Vesicle Cross-Comparisons

Each subfolder compares vesicle trajectories against wall-bound benchmark references taken from the canonical wall-bound condition `c40_m12`.

## Layout convention
- Rows: crosslinker ratio (`1:4`, `1:8`, `1:16`)
- Columns: vesicle viscosity-like condition (`vis=1`, `vis=2`, `vis=4`)
- X-axis: vesicle count (`20`, `50`, `100`)
- Orange curve: vesicle system
- Horizontal benchmark lines: selected wall-bound reference branches
''')

    for theme_key, theme in THEMES.items():
        outdir = cross_root / theme_key
        outdir.mkdir(parents=True, exist_ok=True)
        for metric in METRICS:
            plot_theme(theme_key, theme, ves_lookup, wall_lookup, metric, outdir)
        refs_text = '\n'.join([f'- `{label}`' for _, _, label, *_ in theme['refs']])
        write_text(outdir / 'README.md', f'''# {theme_key.replace('_', ' ').title()}

{theme['title']}.

## What is shown
- Vesicle curves vary with vesicle count.
- Wall-bound references are horizontal benchmark lines from `c40_m12`.
- This lets us compare systems with different native control parameters without faking a one-to-one x-axis match.

## Reference branches in this folder
{refs_text}
''')

    # extend notes summary
    summary = BUNDLE / 'notes/story_summary.md'
    existing = summary.read_text(encoding='utf-8') if summary.exists() else '# Story Summary\n'
    if '## Vesicle comparisons' not in existing:
        existing += '\n\n## Vesicle comparisons\n- Vesicle-only transport metrics are now copied into `vesicle_only/transport_metrics`.\n- Vesicle cross-comparison folders benchmark vesicle curves against wall-bound `c40_m12` references instead of forcing incompatible x-axis mappings.\n- The four vesicle benchmark themes are `vesicle_vs_rotatable`, `vesicle_vs_fixed`, `vesicle_vs_aligned`, and `vesicle_vs_random`.\n'
        summary.write_text(existing, encoding='utf-8')

    print(BUNDLE)


if __name__ == '__main__':
    main()
