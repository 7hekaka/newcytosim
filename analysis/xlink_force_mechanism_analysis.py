#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import re
import subprocess
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
from analysis import run_speckle_timecourse_batch as speck

MAP = ROOT / 'clu_fixed_global_xlink_force_sweep' / 'xlink_regime_map.csv'
CONTROL = ROOT / 'control' / '2'
OUT = ROOT / 'analysis' / 'results' / 'xlink_force_mechanism'
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

FRAME_RE = re.compile(r'^% frame\s+(\d+)')
COUPLE_FORCE_RE = re.compile(r'^\s*(\S+)\s+(\d+)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s*$')
COUPLE_CFG_RE = re.compile(r'^\s*(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$')
FIBER_LEN_RE = re.compile(r'^\s*(\S+)\s+(\d+)\s+([+-]?\d+(?:\.\d+)?)')
CLUSTER_COUNT_RE = re.compile(r'^\s*(\d+)\s+clusters:\s*$')
CLUSTER_LINE_RE = re.compile(r'^\s*(\d+)\s+(\d+):')

METRICS = [
    {'field': 'avg_force_mean', 'ylabel': 'Mean crosslink force (pN)', 'slug': 'couple_avg_force', 'clip_zero': True},
    {'field': 'max_force_peak', 'ylabel': 'Peak crosslink force (pN)', 'slug': 'couple_peak_force', 'clip_zero': True},
    {'field': 'couple_total_mean', 'ylabel': 'Mean bound crosslinkers', 'slug': 'couple_total', 'clip_zero': True},
    {'field': 'largest_cluster_fraction_mean', 'ylabel': 'Largest cluster / all fibers', 'slug': 'largest_cluster_fraction', 'clip_zero': True},
    {'field': 'clustered_fiber_fraction_mean', 'ylabel': 'Clustered fibers / all fibers', 'slug': 'clustered_fiber_fraction', 'clip_zero': True},
]
WINDOW_METRICS = [
    {'field': 'avg_force', 'ylabel': 'Window mean crosslink force (pN)', 'slug': 'window_avg_force', 'clip_zero': True},
    {'field': 'largest_cluster_fraction', 'ylabel': 'Window mean largest-cluster fraction', 'slug': 'window_largest_cluster_fraction', 'clip_zero': True},
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
        rows.append({'model': 'control', 'group': 'control', 'case': 'c0_m0', 'xlink_regime': regime, 'run_dir': run_dir, 'run_path': str(CONTROL / run_dir)})
    return rows


def trim_to_target(time_min: np.ndarray, arrays: dict[str, np.ndarray]) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    keep = int(np.count_nonzero(time_min <= TARGET_FINAL_MIN + 1e-9))
    if keep <= 0:
        return time_min, arrays
    trimmed = {key: np.asarray(val)[:keep] for key, val in arrays.items()}
    return time_min[:keep], trimmed


def _parse_framewise(text: str, parser) -> list[dict]:
    frames = []
    current_frame = None
    block_lines: list[str] = []

    def finalize() -> None:
        nonlocal block_lines, current_frame
        if current_frame is None:
            return
        frames.append(parser(current_frame, block_lines))
        block_lines = []

    for line in text.splitlines():
        m = FRAME_RE.match(line)
        if m:
            finalize()
            current_frame = int(m.group(1))
            continue
        block_lines.append(line.rstrip('\n'))
    finalize()
    return frames


def parse_couple_force(text: str) -> list[dict]:
    def parser(frame: int, lines: list[str]) -> dict:
        row = {'frame': frame, 'couple_name': '', 'count': 0, 'avg_force': 0.0, 'max_force': 0.0, 'max_len': 0.0}
        for line in lines:
            m = COUPLE_FORCE_RE.match(line)
            if m:
                row.update({'couple_name': m.group(1), 'count': int(m.group(2)), 'avg_force': float(m.group(3)), 'max_force': float(m.group(4)), 'max_len': float(m.group(5))})
                break
        return row
    return _parse_framewise(text, parser)


def parse_couple_configuration(text: str) -> list[dict]:
    def parser(frame: int, lines: list[str]) -> dict:
        row = {'frame': frame, 'couples_name': 'all', 'Total': 0, 'P': 0, 'A': 0, 'X': 0, 'T+': 0, 'V+': 0, 'T-': 0, 'V-': 0}
        for line in lines:
            m = COUPLE_CFG_RE.match(line)
            if m:
                row.update({'couples_name': m.group(1), 'Total': int(m.group(2)), 'P': int(m.group(3)), 'A': int(m.group(4)), 'X': int(m.group(5)), 'T+': int(m.group(6)), 'V+': int(m.group(7)), 'T-': int(m.group(8)), 'V-': int(m.group(9))})
                break
        return row
    return _parse_framewise(text, parser)


def parse_fiber_length(text: str) -> list[dict]:
    def parser(frame: int, lines: list[str]) -> dict:
        row = {'frame': frame, 'fiber_class': '', 'count': 0}
        for line in lines:
            m = FIBER_LEN_RE.match(line)
            if m:
                row.update({'fiber_class': m.group(1), 'count': int(m.group(2))})
                break
        return row
    return _parse_framewise(text, parser)


def parse_fiber_cluster(text: str) -> list[dict]:
    def parser(frame: int, lines: list[str]) -> dict:
        cluster_count = 0
        sizes: list[int] = []
        for line in lines:
            m_count = CLUSTER_COUNT_RE.match(line)
            if m_count:
                cluster_count = int(m_count.group(1))
                continue
            m_line = CLUSTER_LINE_RE.match(line)
            if m_line:
                sizes.append(int(m_line.group(2)))
        return {'frame': frame, 'n_clusters': cluster_count, 'largest_cluster_size': max(sizes) if sizes else 0, 'clustered_fibers': int(sum(sizes))}
    return _parse_framewise(text, parser)


def run_report(run_path: Path, report_bin: Path, *args: str) -> str:
    proc = subprocess.run([str(report_bin), *args], cwd=run_path, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=True, text=True, encoding='utf-8')
    return proc.stdout


def trim_frame_rows(rows: list[dict], runs: list[dict]) -> list[dict]:
    total_expected = sum(r['nb_frames'] for r in runs)
    total_expected_with_initial = sum(r['nb_frames'] + 1 for r in runs)
    final_nb = runs[-1]['nb_frames']
    count = len(rows)
    if count == total_expected:
        return rows[-final_nb:]
    if count == total_expected_with_initial:
        return rows[-(final_nb + 1):][1:]
    raise ValueError(f'expected {total_expected} or {total_expected_with_initial} frames, got {count}')


def analyze_run(task: tuple[dict, str]) -> dict:
    row, report_bin_s = task
    run_path = Path(row['run_path'])
    cfg = base.parse_config_info(run_path / 'config.cym')
    dt = cfg['runs'][-1]['frame_dt']
    report_bin = Path(report_bin_s)
    out = {'model': row['model'], 'group': row['group'], 'case': row['case'], 'xlink_regime': row['xlink_regime'], 'run_dir': row['run_dir'], 'run_path': row['run_path'], 'status': 'pending'}
    try:
        force_rows = trim_frame_rows(parse_couple_force(run_report(run_path, report_bin, 'couple:force')), cfg['runs'])
        cfg_rows = trim_frame_rows(parse_couple_configuration(run_report(run_path, report_bin, 'couple:configuration')), cfg['runs'])
        length_rows = trim_frame_rows(parse_fiber_length(run_report(run_path, report_bin, 'fiber:length')), cfg['runs'])
        cluster_rows = trim_frame_rows(parse_fiber_cluster(run_report(run_path, report_bin, 'fiber:cluster', 'couple=1')), cfg['runs'])
        n = len(force_rows)
        if not (len(cfg_rows) == len(length_rows) == len(cluster_rows) == n):
            raise ValueError('report frame counts do not match')
        time_min = (dt * np.arange(1, n + 1, dtype=float)) / 60.0
        arrays = {
            'avg_force': np.asarray([r['avg_force'] for r in force_rows], dtype=float),
            'max_force': np.asarray([r['max_force'] for r in force_rows], dtype=float),
            'couple_total': np.asarray([r['Total'] for r in cfg_rows], dtype=float),
            'couple_X': np.asarray([r['X'] for r in cfg_rows], dtype=float),
            'couple_T': np.asarray([r['T+'] + r['T-'] for r in cfg_rows], dtype=float),
            'fiber_count': np.asarray([r['count'] for r in length_rows], dtype=float),
            'n_clusters': np.asarray([r['n_clusters'] for r in cluster_rows], dtype=float),
            'largest_cluster_size': np.asarray([r['largest_cluster_size'] for r in cluster_rows], dtype=float),
            'clustered_fibers': np.asarray([r['clustered_fibers'] for r in cluster_rows], dtype=float),
        }
        time_min, arrays = trim_to_target(time_min, arrays)
        with np.errstate(divide='ignore', invalid='ignore'):
            total = np.clip(arrays['couple_total'], 1.0, None)
            fiber_total = np.clip(arrays['fiber_count'], 1.0, None)
            arrays['x_fraction'] = arrays['couple_X'] / total
            arrays['t_fraction'] = arrays['couple_T'] / total
            arrays['largest_cluster_fraction'] = arrays['largest_cluster_size'] / fiber_total
            arrays['clustered_fiber_fraction'] = arrays['clustered_fibers'] / fiber_total
        out.update({'status': 'ok', 'time_min': time_min.tolist()})
        for key, value in arrays.items():
            out[key] = np.asarray(value, dtype=float).tolist()
        return out
    except Exception as exc:
        out['status'] = f'error:{exc}'
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


def mean_sem(vals: np.ndarray) -> tuple[float, float]:
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return float('nan'), float('nan')
    if vals.size == 1:
        return float(vals[0]), 0.0
    return float(np.nanmean(vals)), float(np.nanstd(vals, ddof=1) / math.sqrt(vals.size))


def compute_run_metrics(run_rows: list[dict]) -> list[dict]:
    out = []
    for row in run_rows:
        if row.get('status') != 'ok':
            continue
        time_min = np.asarray(row['time_min'], dtype=float)
        metric_row = {'model': row['model'], 'group': row['group'], 'case': row['case'], 'xlink_regime': row['xlink_regime'], 'run_dir': row['run_dir'], 'run_path': row['run_path'], 'n_frames_used': int(time_min.size)}
        for field in ['avg_force', 'couple_total', 'x_fraction', 'largest_cluster_fraction', 'clustered_fiber_fraction', 'n_clusters']:
            arr = np.asarray(row[field], dtype=float)
            metric_row[f'{field}_mean'] = float(np.nanmean(arr)) if arr.size else float('nan')
            for slug, start, stop in WINDOWS:
                mask = (time_min >= start) & (time_min < stop)
                metric_row[f'{slug}_{field}_mean'] = float(np.nanmean(arr[mask])) if np.any(mask) else float('nan')
        max_force = np.asarray(row['max_force'], dtype=float)
        metric_row['max_force_peak'] = float(np.nanmax(max_force)) if max_force.size else float('nan')
        out.append(metric_row)
    return out


def aggregate_condition_metrics(run_metric_rows: list[dict]) -> list[dict]:
    metric_fields = ['avg_force_mean', 'couple_total_mean', 'x_fraction_mean', 'largest_cluster_fraction_mean', 'clustered_fiber_fraction_mean', 'n_clusters_mean', 'max_force_peak']
    for slug, _, _ in WINDOWS:
        metric_fields.extend([f'{slug}_avg_force_mean', f'{slug}_largest_cluster_fraction_mean'])
    grouped = defaultdict(list)
    for row in run_metric_rows:
        grouped[(row['model'], row['group'], row['case'], row['xlink_regime'])].append(row)
    out = []
    for (model, group, case, regime), items in sorted(grouped.items()):
        row = {'model': model, 'group': group, 'case': case, 'xlink_regime': regime, 'n_runs_used': len(items)}
        for field in metric_fields:
            vals = np.asarray([float(item[field]) for item in items], dtype=float)
            row[f'{field}_mean'], row[f'{field}_sem'] = mean_sem(vals)
        out.append(row)
    return out


def build_lookup(rows: list[dict]) -> dict[tuple[str, str, str, str], dict]:
    return {(r['model'], r['group'], r['case'], r['xlink_regime']): r for r in rows}


def style(ax) -> None:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.4)
    ax.spines['bottom'].set_linewidth(1.4)
    ax.tick_params(direction='out', width=1.3, labelsize=13)
    ax.grid(axis='y', color='#d7d7d7', linewidth=0.7, alpha=0.85)


def compute_ylim(rows_lookup: dict, mean_key: str, sem_key: str, clip_zero: bool) -> tuple[float, float]:
    vals = []
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
    if clip_zero:
        lo = max(0.0, lo)
    return lo, hi


def draw_force_family(ax, lookup, group: str, regime: str, mean_key: str, sem_key: str) -> None:
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
    ax.set_xticklabels(['0', '0.5', '1.0', '2.5', '5.0'], fontsize=13)
    style(ax)


def save_metric_figure(lookup: dict, metric: dict) -> Path:
    mean_key = f"{metric['field']}_mean"
    sem_key = f"{metric['field']}_sem"
    fig, axes = plt.subplots(3, len(GROUP_SPECS), figsize=(5.0 * len(GROUP_SPECS) + 0.8, 10.8), sharex='col', sharey=True)
    ymin, ymax = compute_ylim(lookup, mean_key, sem_key, metric['clip_zero'])
    for ridx, regime in enumerate(XLINK_ORDER):
        for cidx, (group, _, xlabel) in enumerate(GROUP_SPECS):
            ax = axes[ridx, cidx]
            draw_force_family(ax, lookup, group, regime, mean_key, sem_key)
            ax.set_ylim(ymin, ymax)
            ax.text(0.03, 0.95, regime, transform=ax.transAxes, ha='left', va='top', fontsize=16)
            if cidx == 0:
                ax.set_ylabel(metric['ylabel'], fontsize=16)
            if ridx == len(XLINK_ORDER) - 1:
                ax.set_xlabel(f"Crosslinker unbinding force (0 = no motors)\n({xlabel})", fontsize=15)
    legend_handles = [
        plt.Line2D([0], [0], color=CONTROL_SPEC[2], marker=CONTROL_SPEC[3], lw=0, markersize=8, markerfacecolor=CONTROL_SPEC[4], markeredgewidth=1.8),
        plt.Line2D([0], [0], color=MODEL_SPEC[2], marker=MODEL_SPEC[3], lw=2.3, markersize=8),
    ]
    fig.legend(legend_handles, ['No motors', 'Fixed global'], loc='lower center', ncol=2, frameon=False, fontsize=14)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.14, wspace=0.18, hspace=0.18)
    outpath = PLOTS / f"{metric['slug']}.png"
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return outpath


def save_group_window_figure(lookup: dict, metric: dict, group_spec: tuple[str, str, str]) -> Path:
    group, slug_group, xlabel = group_spec
    fig, axes = plt.subplots(len(XLINK_ORDER), len(WINDOWS), figsize=(13.2, 10.8), sharex=True, sharey='row')
    mean_keys = [f'{slug}_{metric["field"]}_mean_mean' for slug, _, _ in WINDOWS]
    sem_keys = [f'{slug}_{metric["field"]}_mean_sem' for slug, _, _ in WINDOWS]
    ymins = []
    ymaxs = []
    for mk, sk in zip(mean_keys, sem_keys):
        lo, hi = compute_ylim(lookup, mk, sk, metric['clip_zero'])
        ymins.append(lo)
        ymaxs.append(hi)
    ymin = min(ymins)
    ymax = max(ymaxs)
    for ridx, regime in enumerate(XLINK_ORDER):
        for cidx, (win_slug, start, stop) in enumerate(WINDOWS):
            ax = axes[ridx, cidx]
            mean_key = f'{win_slug}_{metric["field"]}_mean_mean'
            sem_key = f'{win_slug}_{metric["field"]}_mean_sem'
            draw_force_family(ax, lookup, group, regime, mean_key, sem_key)
            ax.set_ylim(ymin, ymax)
            if ridx == 0:
                ax.set_title(f'{int(start)}-{int(stop)} min', fontsize=16, pad=10)
            ax.text(0.03, 0.95, regime, transform=ax.transAxes, ha='left', va='top', fontsize=16)
            if cidx == 0:
                ax.set_ylabel(metric['ylabel'], fontsize=16)
            if ridx == len(XLINK_ORDER) - 1:
                ax.set_xlabel(f"Crosslinker unbinding force (0 = no motors)\n({xlabel})", fontsize=14)
    legend_handles = [
        plt.Line2D([0], [0], color=CONTROL_SPEC[2], marker=CONTROL_SPEC[3], lw=0, markersize=8, markerfacecolor=CONTROL_SPEC[4], markeredgewidth=1.8),
        plt.Line2D([0], [0], color=MODEL_SPEC[2], marker=MODEL_SPEC[3], lw=2.3, markersize=8),
    ]
    fig.legend(legend_handles, ['No motors', 'Fixed global'], loc='lower center', ncol=2, frameon=False, fontsize=14)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.96, bottom=0.14, wspace=0.18, hspace=0.18)
    outpath = PLOTS / f"{slug_group}_{metric['slug']}.png"
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return outpath


def main() -> None:
    plt.rcParams.update({'font.size': 12, 'axes.linewidth': 1.2, 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    report_bin = speck.pick_report_bin()
    rows = load_campaign_rows()
    tasks = [(dict(row, model=MODEL_SPEC[0]), str(report_bin)) for row in rows]
    tasks.extend((row, str(report_bin)) for row in control_metadata())
    with ProcessPoolExecutor(max_workers=4) as ex:
        run_rows = list(ex.map(analyze_run, tasks))
    write_csv(OUT / 'run_status.csv', [{k: v for k, v in row.items() if not isinstance(v, list)} for row in run_rows])
    run_metric_rows = compute_run_metrics(run_rows)
    write_csv(OUT / 'run_metrics.csv', run_metric_rows)
    condition_rows = aggregate_condition_metrics(run_metric_rows)
    write_csv(OUT / 'condition_metrics.csv', condition_rows)
    lookup = build_lookup(condition_rows)
    outfiles = []
    for metric in METRICS:
        outfiles.append(save_metric_figure(lookup, metric))
    for metric in WINDOW_METRICS:
        for group_spec in GROUP_SPECS:
            outfiles.append(save_group_window_figure(lookup, metric, group_spec))
    print(OUT / 'run_status.csv')
    print(OUT / 'run_metrics.csv')
    print(OUT / 'condition_metrics.csv')
    for path in outfiles:
        print(path)


if __name__ == '__main__':
    main()
