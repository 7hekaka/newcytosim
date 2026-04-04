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

import plot_timecourse_panels as base
import run_speckle_kymograph_batch as kymo
import run_speckle_timecourse_batch as speck

ROT = base.DEFAULT_ROTATABLE
FIX = base.DEFAULT_FIXED
CONTROL = ROOT / 'control' / '2'
OUT = ROOT / 'transport_regime_analysis'
PLOTS = OUT / 'plots'

XLINK_ORDER = ['1:4', '1:8', '1:16']
TOTAL480_ORDER = [('c0_m0', 0), ('c10_m48', 10), ('c20_m24', 20), ('c40_m12', 40), ('c60_m8', 60), ('c80_m6', 80)]
MPC12_ORDER = [('c0_m0', 0), ('c10_m12', 10), ('c20_m12', 20), ('c40_m12', 40), ('c60_m12', 60), ('c80_m12', 80)]
MODEL_SPECS = [
    ('rotatable', 'Rotatable', '#222222', 's', 'white'),
    ('fixed_global', 'Fixed global', '#d95f02', 'o', '#d95f02'),
]
CONTROL_SPEC = ('control', 'No motors', '#7f7f7f', 'D', '#d9d9d9')
TARGET_FINAL_MIN = 800.0 / 60.0
WINDOWS = [
    ('early', 0.0, 4.0),
    ('mid', 4.0, 9.0),
    ('late', 9.0, TARGET_FINAL_MIN + 1e-9),
]
METRICS = [
    {'field': 'peak_vz_abs', 'ylabel': 'Peak |v_z| (um/s)', 'slug': 'peak_vz_abs', 'clip_zero': True},
    {'field': 'auc_vz_abs', 'ylabel': 'AUC(|v_z|) (um)', 'slug': 'auc_vz_abs', 'clip_zero': True},
    {'field': 'active_transport_fraction', 'ylabel': 'Active transport fraction', 'slug': 'active_transport_fraction', 'clip_zero': True},
    {'field': 'time_to_peak_min', 'ylabel': 'Time to peak |v_z| (min)', 'slug': 'time_to_peak_min', 'clip_zero': True},
    {'field': 'abs_net_z_displacement', 'ylabel': '|Net z shift| (um)', 'slug': 'abs_net_z_displacement', 'clip_zero': True},
]
WINDOW_METRIC = {'field': 'window_mean_vz_abs', 'ylabel': 'Window mean |v_z| (um/s)', 'slug': 'window_mean_vz_abs', 'clip_zero': True}


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


def load_tasks() -> list[tuple[dict, float, str, str, int]]:
    report_bin = speck.pick_report_bin()
    tasks = []
    for row in control_metadata():
        tasks.append((row, 2.0, 'speckle_i2.txt', str(report_bin), 80))
    for row in base.load_metadata(ROT, 'rotatable'):
        if row['group'] in {'total480', 'mpc12'}:
            tasks.append((row, 2.0, 'speckle_i2.txt', str(report_bin), 80))
    for row in base.load_metadata(FIX, 'fixed_global'):
        if row['group'] in {'total480', 'mpc12'}:
            tasks.append((row, 2.0, 'speckle_i2.txt', str(report_bin), 80))
    return tasks


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


def compute_control_thresholds(run_rows: list[dict]) -> list[dict]:
    grouped = defaultdict(list)
    for row in run_rows:
        if row.get('model') != 'control' or row.get('status') != 'ok':
            continue
        grouped[row['xlink_regime']].extend(float(v) for v in row['vz_abs'])
    out = []
    for regime in XLINK_ORDER:
        vals = np.asarray(grouped.get(regime, []), dtype=float)
        threshold = float(np.nanpercentile(vals, 95.0)) if vals.size else float('nan')
        out.append({'xlink_regime': regime, 'threshold_vz_abs': threshold, 'n_frame_values': int(vals.size)})
    return out


def first_sustained_time(time_min: np.ndarray, signal: np.ndarray, threshold: float, consecutive: int = 3) -> float:
    if not np.isfinite(threshold):
        return float('nan')
    mask = signal > threshold
    count = 0
    for idx, flag in enumerate(mask):
        count = count + 1 if flag else 0
        if count >= consecutive:
            return float(time_min[idx - consecutive + 1])
    return float('nan')


def window_mean(time_min: np.ndarray, values: np.ndarray, start: float, stop: float) -> float:
    mask = (time_min >= start) & (time_min < stop)
    if not np.any(mask):
        return float('nan')
    return float(np.nanmean(values[mask]))


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
            'peak_vz_abs': float(np.nanmax(vz_abs)),
            'time_to_peak_min': float(time_min[peak_idx]) if time_min.size else float('nan'),
            'auc_vz_abs': float(np.trapezoid(vz_abs, time_s)) if time_s.size else float('nan'),
            'active_transport_fraction': float(np.mean(vz_abs > threshold)) if np.isfinite(threshold) else float('nan'),
            'onset_time_min': first_sustained_time(time_min, vz_abs, threshold),
            'net_z_displacement': float(z_com[-1] - z_com[0]) if z_com.size else float('nan'),
            'abs_net_z_displacement': float(abs(z_com[-1] - z_com[0])) if z_com.size else float('nan'),
        }
        for slug, start, stop in WINDOWS:
            metric_row[f'{slug}_mean_vz_abs'] = window_mean(time_min, vz_abs, start, stop)
        out.append(metric_row)
    return out


def aggregate_condition_metrics(run_metric_rows: list[dict]) -> list[dict]:
    metric_fields = [
        'mean_vz_abs', 'peak_vz_abs', 'time_to_peak_min', 'auc_vz_abs',
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
        m = row.get(mean_key, float('nan'))
        s = row.get(sem_key, float('nan'))
        try:
            m = float(m)
            s = float(s)
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


def draw_count_family(ax, family_order, lookup, group: str, regime: str, metric: dict) -> None:
    control_row = lookup.get(('control', 'control', 'c0_m0', regime))
    mean_key = f"{metric['field']}_mean"
    sem_key = f"{metric['field']}_sem"
    if control_row is not None:
        ax.errorbar([0], [float(control_row[mean_key])], yerr=[float(control_row[sem_key])], color=CONTROL_SPEC[2], marker=CONTROL_SPEC[3], markersize=8, lw=0, capsize=4, markeredgewidth=1.8, markeredgecolor=CONTROL_SPEC[2], markerfacecolor=CONTROL_SPEC[4], zorder=5)
    for model, _, color, marker, face in MODEL_SPECS:
        xs = [x for _, x in family_order if x > 0]
        ys = []
        es = []
        for case, x in family_order:
            if x == 0:
                continue
            row = lookup.get((model, group, case, regime))
            ys.append(float(row[mean_key]) if row else float('nan'))
            es.append(float(row[sem_key]) if row else float('nan'))
        ax.errorbar(xs, ys, yerr=es, color=color, marker=marker, markersize=8, lw=2.3, capsize=4, markeredgewidth=1.8, markeredgecolor=color, markerfacecolor=face)
    ax.set_xticks([x for _, x in family_order])
    ax.set_xticklabels([str(x) for _, x in family_order], fontsize=14)
    style(ax)


def save_metric_figure(lookup: dict, metric: dict) -> Path:
    plt.rcParams.update({'font.size': 14, 'axes.linewidth': 1.3, 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
    fig, axes = plt.subplots(3, 2, figsize=(13.2, 13.2), sharey=True)
    ylo, yhi = compute_ylim(lookup, metric)
    for ridx, regime in enumerate(XLINK_ORDER):
        draw_count_family(axes[ridx, 0], TOTAL480_ORDER, lookup, 'total480', regime, metric)
        draw_count_family(axes[ridx, 1], MPC12_ORDER, lookup, 'mpc12', regime, metric)
        for cidx in range(2):
            axes[ridx, cidx].set_ylim(ylo, yhi)
            axes[ridx, cidx].text(0.03, 0.94, regime, transform=axes[ridx, cidx].transAxes, ha='left', va='top', fontsize=16, color='#444444')
    for ridx in range(3):
        axes[ridx, 0].set_ylabel(metric['ylabel'], fontsize=18)
    axes[2, 0].set_xlabel('Cluster count (0 = no motors)\n(total motors = 480)', fontsize=17)
    axes[2, 1].set_xlabel('Cluster count (0 = no motors)\n(12 motors per cluster)', fontsize=17)
    handles = []
    labels = []
    for _, label, color, marker, face in [CONTROL_SPEC] + MODEL_SPECS:
        handles.append(plt.Line2D([0], [0], color=color, lw=2.3 if label != 'No motors' else 0, marker=marker, markersize=8, markeredgewidth=1.8, markeredgecolor=color, markerfacecolor=face))
        labels.append(label)
    fig.legend(handles, labels, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=16)
    fig.subplots_adjust(left=0.11, right=0.99, top=0.98, bottom=0.12, wspace=0.18, hspace=0.22)
    out = PLOTS / f"{metric['slug']}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out


def save_window_figure(lookup: dict) -> list[Path]:
    outs = []
    for group, family_order, slug, xlabel in [
        ('total480', TOTAL480_ORDER, 'total480_window_mean_vz_abs', 'Cluster count (0 = no motors, total motors = 480)'),
        ('mpc12', MPC12_ORDER, 'mpc12_window_mean_vz_abs', 'Cluster count (0 = no motors, 12 motors per cluster)'),
    ]:
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
        for ridx, regime in enumerate(XLINK_ORDER):
            for cidx, (wslug, start, stop) in enumerate(WINDOWS):
                ax = axes[ridx, cidx]
                metric = {'field': f'{wslug}_mean_vz_abs', 'ylabel': WINDOW_METRIC['ylabel'], 'clip_zero': True}
                draw_count_family(ax, family_order, lookup, group, regime, metric)
                ax.set_ylim(ylo, yhi)
                if ridx == 0:
                    ax.set_title(f'{int(start)}-{int(stop if stop < TARGET_FINAL_MIN else round(TARGET_FINAL_MIN))} min', fontsize=16, pad=8)
                if cidx == 0:
                    ax.text(0.03, 0.94, regime, transform=ax.transAxes, ha='left', va='top', fontsize=16, color='#444444')
                if cidx > 0:
                    ax.set_yticklabels([])
        for ridx in range(3):
            axes[ridx, 0].set_ylabel(WINDOW_METRIC['ylabel'], fontsize=18)
        for cidx in range(3):
            axes[2, cidx].set_xlabel(xlabel, fontsize=17)
        handles = []
        labels = []
        for _, label, color, marker, face in [CONTROL_SPEC] + MODEL_SPECS:
            handles.append(plt.Line2D([0], [0], color=color, lw=2.3 if label != 'No motors' else 0, marker=marker, markersize=8, markeredgewidth=1.8, markeredgecolor=color, markerfacecolor=face))
            labels.append(label)
        fig.legend(handles, labels, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=16)
        fig.subplots_adjust(left=0.10, right=0.99, top=0.95, bottom=0.12, wspace=0.18, hspace=0.22)
        out = PLOTS / f'{slug}.png'
        fig.savefig(out, dpi=300, bbox_inches='tight')
        plt.close(fig)
        outs.append(out)
    return outs


def compute_medoid(rows: list[dict]) -> dict | None:
    usable = [row for row in rows if row.get('status') == 'ok']
    if not usable:
        return None
    n_frames = min(len(row['time_min']) for row in usable)
    mats = [np.asarray(row['rho_z'][:n_frames], dtype=float).reshape(-1) for row in usable]
    dists = np.zeros((len(mats), len(mats)), dtype=float)
    for i in range(len(mats)):
        for j in range(i + 1, len(mats)):
            d = float(np.linalg.norm(mats[i] - mats[j]))
            dists[i, j] = dists[j, i] = d
    idx = int(np.argmin(np.mean(dists, axis=1)))
    return usable[idx]


def build_medoid_condition_rows(run_rows: list[dict]) -> tuple[list[dict], list[dict]]:
    grouped = defaultdict(list)
    for row in run_rows:
        if row.get('status') == 'ok':
            grouped[(row['model'], row['group'], row['case'], row['xlink_regime'])].append(row)
    medoid_rows = []
    summary_rows = []
    for key, items in sorted(grouped.items()):
        medoid = compute_medoid(items)
        if medoid is None:
            continue
        model, group, case, regime = key
        summary_rows.append({'model': model, 'group': group, 'case': case, 'xlink_regime': regime, 'medoid_run_dir': medoid['run_dir'], 'n_runs_candidates': len(items)})
        medoid_rows.append({
            'model': model,
            'group': group,
            'case': case,
            'xlink_regime': regime,
            'n_runs_used': 1,
            'time_min': np.asarray(medoid['time_min'], dtype=float),
            'z_edges': np.asarray(medoid['z_edges'], dtype=float),
            'rho_z_mean': np.asarray(medoid['rho_z'], dtype=float),
        })
    return medoid_rows, summary_rows


def duplicate_control_condition_rows(rows: list[dict]) -> list[dict]:
    out = list(rows)
    control_rows = [row for row in rows if row['model'] == 'control']
    for row in control_rows:
        for model in ['rotatable', 'fixed_global']:
            for group in ['total480', 'mpc12']:
                clone = dict(row)
                clone['model'] = model
                clone['group'] = group
                out.append(clone)
    return out


def save_medoid_kymographs(run_rows: list[dict]) -> list[Path]:
    medoid_rows, summary_rows = build_medoid_condition_rows(run_rows)
    medoid_rows = duplicate_control_condition_rows(medoid_rows)
    summary_rows = duplicate_control_condition_rows(summary_rows)
    write_csv(OUT / 'medoid_selection_summary.csv', summary_rows)
    kymo.write_condition_csv(OUT / 'medoid_condition_summary.csv', medoid_rows)
    lookup = kymo.build_lookup(medoid_rows)
    family_specs = [
        {'group': 'total480', 'slug': 'total480_medoid', 'order': [('c0_m0', '0'), ('c10_m48', '10'), ('c20_m24', '20'), ('c40_m12', '40'), ('c60_m8', '60'), ('c80_m6', '80')]},
        {'group': 'mpc12', 'slug': 'mpc12_medoid', 'order': [('c0_m0', '0'), ('c10_m12', '10'), ('c20_m12', '20'), ('c40_m12', '40'), ('c60_m12', '60'), ('c80_m12', '80')]},
    ]
    outs = []
    for family in family_specs:
        outs.append(kymo.plot_family_model(lookup, family, model='rotatable', outdir=PLOTS, vmin=0.0, vmax=2.0))
        outs.append(kymo.plot_family_model(lookup, family, model='fixed_global', outdir=PLOTS, vmin=0.0, vmax=2.0))
    kymo.save_standalone_colorbar(PLOTS, vmin=0.0, vmax=2.0)
    return outs


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    tasks = load_tasks()
    with ProcessPoolExecutor(max_workers=4) as ex:
        run_rows = list(ex.map(analyze_transport_run, tasks))
    write_csv(OUT / 'run_status.csv', [{k: v for k, v in row.items() if k not in {'time_min', 'vz_abs', 'z_com', 'rho_z', 'z_edges'}} for row in run_rows])

    threshold_rows = compute_control_thresholds(run_rows)
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
    outs.extend(save_window_figure(lookup))
    outs.extend(save_medoid_kymographs(run_rows))

    print(OUT / 'run_status.csv')
    print(OUT / 'control_thresholds.csv')
    print(OUT / 'run_metrics.csv')
    print(OUT / 'condition_metrics.csv')
    for out in outs:
        print(out)


if __name__ == '__main__':
    main()

