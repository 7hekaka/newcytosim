#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis import plot_timecourse_panels as base
from analysis import run_speckle_kymograph_batch as kymo
from analysis import run_speckle_timecourse_batch as speck
from analysis import transport_regime_batch as treg

MAP = ROOT / 'clu_fixed_global_xlink_force_sweep' / 'xlink_regime_map.csv'
CONTROL = ROOT / 'control' / '2'
OUT = ROOT / 'analysis' / 'results' / 'xlink_force_sweep'
PLOTS = OUT / 'plots'

XLINK_ORDER = ['1:4', '1:8', '1:16']
FORCE_ORDER = [('c0_m0', 0.0), ('xlink_uf0p5', 0.5), ('xlink_uf1p0', 1.0), ('xlink_uf2p5', 2.5), ('xlink_uf5p0', 5.0)]
GROUP_SPECS = [
    ('c40_m12', 'c40_m12', '40 clusters, 12 motors/cluster'),
    ('c40_m50', 'c40_m50', '40 clusters, 50 motors/cluster'),
    ('c40_m100', 'c40_m100', '40 clusters, 100 motors/cluster'),
]
MODEL_SPEC = ('fixed_global_xlink_force', 'Fixed global', '#d95f02', 'o', '#d95f02')
CONTROL_SPEC = ('control', 'No motors', '#7f7f7f', 'D', '#d9d9d9')
TARGET_FINAL_MIN = 800.0 / 60.0
WINDOWS = [
    ('early', 0.0, 4.0),
    ('mid', 4.0, 9.0),
    ('late', 9.0, TARGET_FINAL_MIN + 1e-9),
]
METRICS = [
    {'field': 'mean_vz_abs', 'ylabel': 'Mean |v_z| (um/s)', 'slug': 'mean_vz_abs', 'clip_zero': True},
    {'field': 'std_vz_abs', 'ylabel': 'SD(|v_z|) (um/s)', 'slug': 'std_vz_abs', 'clip_zero': True},
    {'field': 'peak_vz_abs', 'ylabel': 'Peak |v_z| (um/s)', 'slug': 'peak_vz_abs', 'clip_zero': True},
    {'field': 'auc_vz_abs', 'ylabel': 'AUC(|v_z|) (um)', 'slug': 'auc_vz_abs', 'clip_zero': True},
    {'field': 'active_transport_fraction', 'ylabel': 'Active transport fraction', 'slug': 'active_transport_fraction', 'clip_zero': True},
    {'field': 'time_to_peak_min', 'ylabel': 'Time to peak |v_z| (min)', 'slug': 'time_to_peak_min', 'clip_zero': True},
    {'field': 'abs_net_z_displacement', 'ylabel': '|Net z shift| (um)', 'slug': 'abs_net_z_displacement', 'clip_zero': True},
]


def load_campaign_rows() -> list[dict]:
    rows = list(csv.DictReader(MAP.open(encoding='utf-8')))
    keep = []
    for row in rows:
        row = dict(row)
        if row.get('has_outputs') != '1':
            continue
        keep.append(row)
    return keep


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
        rows.append(
            {
                'model': 'control',
                'group': 'control',
                'case': 'c0_m0',
                'xlink_regime': regime,
                'run_dir': run_dir,
                'run_path': str(CONTROL / run_dir),
            }
        )
    return rows


def trim_to_target(out: dict) -> dict:
    if out.get('status') != 'ok':
        return out
    time = np.asarray(out['time_min'], dtype=float)
    keep = int(np.count_nonzero(time <= TARGET_FINAL_MIN + 1e-9))
    if keep <= 0:
        return out
    for key in ['time_min', 'vz_abs', 'z_com', 'rho_z']:
        if key in out:
            out[key] = out[key][:keep]
    return out


def analyze_transport_run(task: tuple[dict, float, str, str, int]) -> dict:
    row, interval, speckle_name, report_bin_s, nbins_z = task
    out = kymo.analyze_run((row, interval, speckle_name, report_bin_s, nbins_z))
    if out.get('status') != 'ok':
        return out
    run_path = Path(row['run_path'])
    cfg = base.parse_config_info(run_path / 'config.cym')
    dt = cfg['runs'][-1]['frame_dt']
    rho = np.asarray(out['rho_z'], dtype=float)
    z_edges = np.asarray(out['z_edges'], dtype=float)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    z_com = rho @ z_centers / np.clip(rho.sum(axis=1), 1e-12, None)
    out['dt_s'] = float(dt)
    out['z_com'] = z_com.tolist()
    out = trim_to_target(out)
    time = np.asarray(out['time_min'], dtype=float)
    z_com = np.asarray(out['z_com'], dtype=float)
    if time.size > 1:
        dt_local = (time[1] - time[0]) * 60.0
    else:
        dt_local = dt
    out['dt_s'] = float(dt_local)
    out['vz_abs'] = np.abs(np.gradient(z_com, dt_local)).tolist() if z_com.size > 1 else [0.0]
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def compute_run_metrics(run_rows: list[dict], thresholds: dict[str, float]) -> list[dict]:
    out = []
    for row in run_rows:
        if row.get('status') != 'ok':
            continue
        time_min = np.asarray(row['time_min'], dtype=float)
        time_s = time_min * 60.0
        vz_abs = np.asarray(row['vz_abs'], dtype=float)
        z_com = np.asarray(row['z_com'], dtype=float)
        threshold = thresholds.get(row['xlink_regime'], float('nan'))
        peak_idx = int(np.nanargmax(vz_abs)) if vz_abs.size else 0
        metric_row = {
            'model': row['model'],
            'group': row['group'],
            'case': row['case'],
            'xlink_regime': row['xlink_regime'],
            'run_dir': row['run_dir'],
            'run_path': row['run_path'],
            'n_frames_used': int(vz_abs.size),
            'threshold_vz_abs': threshold,
            'mean_vz_abs': float(np.nanmean(vz_abs)),
            'std_vz_abs': float(np.nanstd(vz_abs, ddof=1)) if vz_abs.size > 1 else 0.0,
            'peak_vz_abs': float(np.nanmax(vz_abs)),
            'time_to_peak_min': float(time_min[peak_idx]) if time_min.size else float('nan'),
            'auc_vz_abs': float(np.trapezoid(vz_abs, time_s)) if time_s.size else float('nan'),
            'active_transport_fraction': float(np.mean(vz_abs > threshold)) if np.isfinite(threshold) else float('nan'),
            'onset_time_min': treg.first_sustained_time(time_min, vz_abs, threshold),
            'net_z_displacement': float(z_com[-1] - z_com[0]) if z_com.size else float('nan'),
            'abs_net_z_displacement': float(abs(z_com[-1] - z_com[0])) if z_com.size else float('nan'),
        }
        for slug, start, stop in WINDOWS:
            metric_row[f'{slug}_mean_vz_abs'] = treg.window_mean(time_min, vz_abs, start, stop)
        out.append(metric_row)
    return out


def aggregate_condition_metrics(run_metric_rows: list[dict]) -> list[dict]:
    metric_fields = [
        'mean_vz_abs', 'std_vz_abs', 'peak_vz_abs', 'time_to_peak_min', 'auc_vz_abs',
        'active_transport_fraction', 'onset_time_min', 'net_z_displacement', 'abs_net_z_displacement',
        'early_mean_vz_abs', 'mid_mean_vz_abs', 'late_mean_vz_abs'
    ]
    grouped = defaultdict(list)
    for row in run_metric_rows:
        grouped[(row['model'], row['group'], row['case'], row['xlink_regime'])].append(row)
    out = []
    for (model, group, case, regime), items in sorted(grouped.items()):
        row = {'model': model, 'group': group, 'case': case, 'xlink_regime': regime, 'n_runs_used': len(items)}
        for field in metric_fields:
            vals = np.asarray([float(item[field]) for item in items], dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                row[f'{field}_mean'] = float('nan')
                row[f'{field}_sem'] = float('nan')
            else:
                row[f'{field}_mean'] = float(np.nanmean(vals))
                row[f'{field}_sem'] = float(np.nanstd(vals, ddof=1) / math.sqrt(vals.size)) if vals.size > 1 else 0.0
        out.append(row)
    return out


def build_lookup(rows: list[dict]) -> dict[tuple[str, str, str, str], dict]:
    return {(r['model'], r['group'], r['case'], r['xlink_regime']): r for r in rows}


def compute_ylim(rows_lookup: dict, metric: dict) -> tuple[float, float]:
    vals = []
    mean_key = f"{metric['field']}_mean"
    sem_key = f"{metric['field']}_sem"
    for row in rows_lookup.values():
        try:
            m = float(row.get(mean_key, float('nan')))
            s = float(row.get(sem_key, float('nan')))
        except Exception:
            continue
        if np.isfinite(m):
            vals.extend([m - s if np.isfinite(s) else m, m + s if np.isfinite(s) else m])
    if not vals:
        return 0.0, 1.0
    lo = min(vals)
    hi = max(vals)
    pad = 0.10 * (hi - lo) if hi > lo else 0.10 * max(abs(hi), 1.0)
    lo -= pad
    hi += pad
    if metric['clip_zero']:
        lo = max(0.0, lo)
    return lo, hi


def style(ax) -> None:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.4)
    ax.spines['bottom'].set_linewidth(1.4)
    ax.tick_params(direction='out', width=1.3, labelsize=14)
    ax.grid(axis='y', color='#d7d7d7', linewidth=0.7, alpha=0.85)


def draw_force_family(ax, lookup, group: str, regime: str, metric: dict) -> None:
    mean_key = f"{metric['field']}_mean"
    sem_key = f"{metric['field']}_sem"
    control_row = lookup.get(('control', 'control', 'c0_m0', regime))
    if control_row is not None:
        ax.errorbar([0.0], [float(control_row[mean_key])], yerr=[float(control_row[sem_key])], color=CONTROL_SPEC[2], marker=CONTROL_SPEC[3], markersize=8, lw=0, capsize=4, markeredgewidth=1.8, markeredgecolor=CONTROL_SPEC[2], markerfacecolor=CONTROL_SPEC[4], zorder=5)
    xs = [x for _, x in FORCE_ORDER if x > 0]
    ys = []
    es = []
    for case, x in FORCE_ORDER:
        if x == 0:
            continue
        row = lookup.get((MODEL_SPEC[0], group, case, regime))
        ys.append(float(row[mean_key]) if row else float('nan'))
        es.append(float(row[sem_key]) if row else float('nan'))
    ax.errorbar(xs, ys, yerr=es, color=MODEL_SPEC[2], marker=MODEL_SPEC[3], markersize=8, lw=2.3, capsize=4, markeredgewidth=1.8, markeredgecolor=MODEL_SPEC[2], markerfacecolor=MODEL_SPEC[4])
    ax.set_xticks([x for _, x in FORCE_ORDER])
    ax.set_xticklabels(['0', '0.5', '1.0', '2.5', '5.0'], fontsize=14)
    style(ax)


def save_metric_figure(lookup: dict, metric: dict) -> Path:
    plt.rcParams.update({'font.size': 14, 'axes.linewidth': 1.3, 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
    ncols = len(GROUP_SPECS)
    fig, axes = plt.subplots(3, ncols, figsize=(6.0 * ncols + 0.8, 13.2), sharey=True)
    if ncols == 1:
        axes = np.array(axes).reshape(3, 1)
    ylo, yhi = compute_ylim(lookup, metric)
    for ridx, regime in enumerate(XLINK_ORDER):
        for cidx, (group, _, xlabel) in enumerate(GROUP_SPECS):
            draw_force_family(axes[ridx, cidx], lookup, group, regime, metric)
            axes[ridx, cidx].set_ylim(ylo, yhi)
            axes[ridx, cidx].text(0.03, 0.94, regime, transform=axes[ridx, cidx].transAxes, ha='left', va='top', fontsize=16, color='#444444')
    for ridx in range(3):
        axes[ridx, 0].set_ylabel(metric['ylabel'], fontsize=18)
    for cidx, (_, _, xlabel) in enumerate(GROUP_SPECS):
        axes[2, cidx].set_xlabel(f'Crosslinker unbinding force (0 = no motors)\n({xlabel})', fontsize=17)
    handles = [
        plt.Line2D([0], [0], color=CONTROL_SPEC[2], lw=0, marker=CONTROL_SPEC[3], markersize=8, markeredgewidth=1.8, markeredgecolor=CONTROL_SPEC[2], markerfacecolor=CONTROL_SPEC[4]),
        plt.Line2D([0], [0], color=MODEL_SPEC[2], lw=2.3, marker=MODEL_SPEC[3], markersize=8, markeredgewidth=1.8, markeredgecolor=MODEL_SPEC[2], markerfacecolor=MODEL_SPEC[4]),
    ]
    labels = [CONTROL_SPEC[1], MODEL_SPEC[1]]
    fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=16)
    fig.subplots_adjust(left=0.10, right=0.99, top=0.98, bottom=0.12, wspace=0.18, hspace=0.22)
    out = PLOTS / f"{metric['slug']}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out


def save_window_figure(lookup: dict) -> Path:
    plt.rcParams.update({'font.size': 14, 'axes.linewidth': 1.3, 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
    fig, axes = plt.subplots(3, 3, figsize=(16.2, 12.6), sharey='row')
    vals = []
    for row in lookup.values():
        for wslug, _, _ in WINDOWS:
            key = f'{wslug}_mean_vz_abs_mean'
            sem_key = f'{wslug}_mean_vz_abs_sem'
            try:
                m = float(row[key]); s = float(row[sem_key])
            except Exception:
                continue
            if np.isfinite(m):
                vals.extend([m - s if np.isfinite(s) else m, m + s if np.isfinite(s) else m])
    ylo = max(0.0, min(vals) - 0.0003) if vals else 0.0
    yhi = max(vals) + 0.00045 if vals else 1.0
    colors = {'c40_m12': '#1b9e77', 'c40_m50': '#d95f02', 'c40_m100': '#7570b3'}
    markers = {'c40_m12': 'o', 'c40_m50': 's', 'c40_m100': '^'}
    faces = {'c40_m12': 'white', 'c40_m50': '#d95f02', 'c40_m100': 'white'}
    for ridx, regime in enumerate(XLINK_ORDER):
        for cidx, (wslug, start, stop) in enumerate(WINDOWS):
            ax = axes[ridx, cidx]
            control_row = lookup.get(('control', 'control', 'c0_m0', regime))
            if control_row is not None:
                ax.errorbar([0.0], [float(control_row[f'{wslug}_mean_vz_abs_mean'])], yerr=[float(control_row[f'{wslug}_mean_vz_abs_sem'])], color=CONTROL_SPEC[2], marker=CONTROL_SPEC[3], markersize=8, lw=0, capsize=4, markeredgewidth=1.8, markeredgecolor=CONTROL_SPEC[2], markerfacecolor=CONTROL_SPEC[4], zorder=5)
            for group, _, _ in GROUP_SPECS:
                xs = [x for _, x in FORCE_ORDER if x > 0]
                ys = []
                es = []
                for case, x in FORCE_ORDER:
                    if x == 0:
                        continue
                    row = lookup.get((MODEL_SPEC[0], group, case, regime))
                    ys.append(float(row[f'{wslug}_mean_vz_abs_mean']) if row else float('nan'))
                    es.append(float(row[f'{wslug}_mean_vz_abs_sem']) if row else float('nan'))
                ax.errorbar(xs, ys, yerr=es, color=colors[group], marker=markers[group], markersize=8, lw=2.3, capsize=4, markeredgewidth=1.8, markeredgecolor=colors[group], markerfacecolor=faces[group])
            ax.set_ylim(ylo, yhi)
            ax.set_xticks([x for _, x in FORCE_ORDER])
            ax.set_xticklabels(['0', '0.5', '1.0', '2.5', '5.0'], fontsize=14)
            style(ax)
            if ridx == 0:
                title = f'{int(start)}-{int(stop if stop < TARGET_FINAL_MIN else round(TARGET_FINAL_MIN))} min'
                ax.set_title(title, fontsize=16, pad=8)
            if cidx == 0:
                ax.text(0.03, 0.94, regime, transform=ax.transAxes, ha='left', va='top', fontsize=16, color='#444444')
    for ridx in range(3):
        axes[ridx, 0].set_ylabel('Window mean |v_z| (um/s)', fontsize=18)
    for cidx in range(3):
        axes[2, cidx].set_xlabel('Crosslinker unbinding force (0 = no motors)', fontsize=17)
    handles = [
        plt.Line2D([0], [0], color=CONTROL_SPEC[2], lw=0, marker=CONTROL_SPEC[3], markersize=8, markeredgewidth=1.8, markeredgecolor=CONTROL_SPEC[2], markerfacecolor=CONTROL_SPEC[4])
    ]
    labels = ['No motors']
    for group, _, xlabel in GROUP_SPECS:
        handles.append(plt.Line2D([0], [0], color=colors[group], lw=2.3, marker=markers[group], markersize=8, markeredgewidth=1.8, markeredgecolor=colors[group], markerfacecolor=faces[group]))
        labels.append(xlabel)
    fig.legend(handles, labels, loc='lower center', ncol=len(labels), frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=15)
    fig.subplots_adjust(left=0.10, right=0.99, top=0.95, bottom=0.14, wspace=0.18, hspace=0.22)
    out = PLOTS / 'window_mean_vz_abs.png'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out


def duplicate_control_rows_for_groups(rows: list[dict]) -> list[dict]:
    out = [row for row in rows if row['model'] != 'control']
    control_rows = [row for row in rows if row['model'] == 'control']
    for row in control_rows:
        for group, _, _ in GROUP_SPECS:
            clone = dict(row)
            clone['model'] = MODEL_SPEC[0]
            clone['group'] = group
            out.append(clone)
    return out


def save_medoid_kymographs(run_rows: list[dict]) -> list[Path]:
    medoid_rows, summary_rows = treg.build_medoid_condition_rows(run_rows)
    medoid_rows = duplicate_control_rows_for_groups(medoid_rows)
    summary_rows = duplicate_control_rows_for_groups(summary_rows)
    write_csv(OUT / 'medoid_selection_summary.csv', summary_rows)
    kymo.write_condition_csv(OUT / 'medoid_condition_summary.csv', medoid_rows)
    lookup = kymo.build_lookup(medoid_rows)
    family_specs = [
        {'group': group, 'slug': f'{group}_medoid', 'order': [('c0_m0', '0'), ('uf0p5', '0.5'), ('uf1p0', '1.0'), ('uf2p5', '2.5'), ('uf5p0', '5.0')]}
        for group, _, _ in GROUP_SPECS
    ]
    outs = []
    for family in family_specs:
        outs.append(kymo.plot_family_model(lookup, family, model=MODEL_SPEC[0], outdir=PLOTS, vmin=0.0, vmax=2.0))
        outs.append(kymo.plot_family_model(lookup, family, model=MODEL_SPEC[0], outdir=PLOTS, vmin=0.0, vmax=2.0, with_colorbar=True))
    kymo.save_standalone_colorbar(PLOTS, vmin=0.0, vmax=2.0)
    return outs


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    report_bin = speck.pick_report_bin()
    tasks = []
    for row in control_metadata():
        tasks.append((row, 2.0, 'speckle_i2.txt', str(report_bin), 80))
    for row in load_campaign_rows():
        tasks.append((row, 2.0, 'speckle_i2.txt', str(report_bin), 80))

    with ProcessPoolExecutor(max_workers=4) as ex:
        run_rows = list(ex.map(analyze_transport_run, tasks))

    write_csv(OUT / 'run_status.csv', [{k: v for k, v in row.items() if k not in {'time_min', 'vz_abs', 'z_com', 'rho_z', 'z_edges'}} for row in run_rows])

    threshold_rows = treg.compute_control_thresholds(run_rows)
    write_csv(OUT / 'control_thresholds.csv', threshold_rows)
    threshold_lookup = {row['xlink_regime']: row['threshold_vz_abs'] for row in threshold_rows}

    run_metric_rows = compute_run_metrics(run_rows, threshold_lookup)
    write_csv(OUT / 'run_metrics.csv', run_metric_rows)
    condition_rows = aggregate_condition_metrics(run_metric_rows)
    write_csv(OUT / 'condition_metrics.csv', condition_rows)
    lookup = build_lookup(condition_rows)

    outs = []
    for metric in METRICS:
        outs.append(save_metric_figure(lookup, metric))
    outs.append(save_window_figure(lookup))
    outs.extend(save_medoid_kymographs(run_rows))

    print(OUT / 'run_status.csv')
    print(OUT / 'control_thresholds.csv')
    print(OUT / 'run_metrics.csv')
    print(OUT / 'condition_metrics.csv')
    for path in outs:
        print(path)


if __name__ == '__main__':
    main()
