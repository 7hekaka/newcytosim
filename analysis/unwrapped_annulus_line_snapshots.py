#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / 'analysis') not in sys.path:
    sys.path.insert(0, str(ROOT / 'analysis'))

from analysis import plot_timecourse_panels as base
from analysis import run_speckle_timecourse_batch as speck

MAP = ROOT / 'clu_fixed_global_force_sweep' / 'xlink_regime_map.csv'
MEDOID = ROOT / 'analysis' / 'results' / 'motor_force_sweep' / 'medoid_selection_summary.csv'
CONTROL = ROOT / 'control' / '2'
OUT = ROOT / 'analysis' / 'results' / 'motor_force_sweep' / 'unwrapped_annulus_lines'

XLINK_ORDER = ['1:4', '1:8', '1:16']
XLINK_SLUG = {'1:4': '1to4', '1:8': '1to8', '1:16': '1to16'}
FORCE_ORDER = [('c0_m0', 0.0), ('uf0p5', 0.5), ('uf1p0', 1.0), ('uf2p5', 2.5), ('uf5p0', 5.0)]
GROUP_SPECS = [
    ('c40_m12', '40 clusters, 12 motors/cluster'),
    ('c40_m50', '40 clusters, 50 motors/cluster'),
    ('c40_m100', '40 clusters, 100 motors/cluster'),
]
SNAPSHOTS = [
    ('t00', 0.00, '0%'),
    ('t25', 0.25, '25%'),
    ('t50', 0.50, '50%'),
    ('t100', 1.00, '100%'),
]
SPECKLE_NAME = 'speckle_i2.txt'
SPECKLE_INTERVAL = 2.0
SEAM_BINS = 180
TARGET_FINAL_MIN = 800.0 / 60.0


@dataclass
class SnapshotPanel:
    group: str
    regime: str
    case: str
    force_value: float
    run_dir: str
    run_path: str
    snapshot_slug: str
    snapshot_label: str
    time_min: float
    smax: float
    zmin: float
    zmax: float
    seam: float
    frame: np.ndarray


def load_run_lookup() -> dict[tuple[str, str, str, str], dict]:
    rows = list(csv.DictReader(MAP.open(encoding='utf-8')))
    lookup = {}
    for row in rows:
        if row.get('has_outputs') != '1':
            continue
        lookup[(row['group'], row['case'], row['xlink_regime'], row['run_dir'])] = dict(row)
    return lookup


def resolve_run_path(group: str, case: str, regime: str, run_dir: str, lookup: dict) -> Path:
    if case == 'c0_m0':
        return CONTROL / run_dir
    row = lookup.get((group, case, regime, run_dir))
    if row is None:
        raise KeyError(f'missing run lookup for {(group, case, regime, run_dir)}')
    return Path(row['run_path'])


def load_medoid_rows() -> list[dict]:
    return [dict(row) for row in csv.DictReader(MEDOID.open(encoding='utf-8'))]


def load_phase_frames(run_path: Path, report_bin: Path) -> tuple[list[np.ndarray], np.ndarray, dict]:
    speckle_path = run_path / SPECKLE_NAME
    speck.ensure_speckle(run_path, speckle_path, report_bin, SPECKLE_INTERVAL)
    cfg = base.parse_config_info(run_path / 'config.cym')
    runs = cfg['runs']
    final_run = runs[-1]
    nb_frames = final_run['nb_frames']
    dt_min = final_run['frame_dt'] / 60.0
    frames = speck.parse_frames_speckle(speckle_path)
    total_expected = sum(r['nb_frames'] for r in runs)
    total_expected_with_initial = sum(r['nb_frames'] + 1 for r in runs)
    if len(frames) == total_expected:
        phase_frames = frames[-nb_frames:]
        time_min = dt_min * np.arange(1, nb_frames + 1, dtype=float)
    elif len(frames) == total_expected_with_initial:
        phase_frames = frames[-(nb_frames + 1):]
        time_min = dt_min * np.arange(0, nb_frames + 1, dtype=float)
    else:
        raise ValueError(f'expected {total_expected} or {total_expected_with_initial} frames, parsed {len(frames)} in {speckle_path}')
    keep = int(np.count_nonzero(time_min <= TARGET_FINAL_MIN + 1e-9))
    if keep > 0:
        phase_frames = phase_frames[:keep]
        time_min = time_min[:keep]
    return phase_frames, time_min, cfg['space']


def pick_snapshot_index(time_min: np.ndarray, fraction: float) -> int:
    if time_min.size == 0:
        return 0
    target = float(time_min[-1]) * fraction
    return int(np.argmin(np.abs(time_min - target)))


def choose_theta_seam(frame: np.ndarray, n_bins: int = SEAM_BINS) -> float:
    if frame.size == 0:
        return 0.0
    theta = np.mod(np.arctan2(frame[:, 2], frame[:, 1]), 2.0 * math.pi)
    hist, edges = np.histogram(theta, bins=n_bins, range=(0.0, 2.0 * math.pi))
    idx = int(np.argmin(hist))
    return float(0.5 * (edges[idx] + edges[idx + 1]))


def collect_panels() -> list[SnapshotPanel]:
    report_bin = speck.pick_report_bin()
    run_lookup = load_run_lookup()
    medoid_rows = load_medoid_rows()
    wanted = {(group, case, regime) for group, _ in GROUP_SPECS for case, _ in FORCE_ORDER for regime in XLINK_ORDER}
    medoid_lookup = {}
    for row in medoid_rows:
        key = (row['group'], row['case'], row['xlink_regime'])
        if key in wanted:
            medoid_lookup[key] = row

    panels: list[SnapshotPanel] = []
    for group, _ in GROUP_SPECS:
        for regime in XLINK_ORDER:
            for case, force_value in FORCE_ORDER:
                medoid = medoid_lookup.get((group, case, regime))
                if medoid is None:
                    continue
                run_dir = medoid['medoid_run_dir']
                run_path = resolve_run_path(group, case, regime, run_dir, run_lookup)
                phase_frames, time_min, space = load_phase_frames(run_path, report_bin)
                seam = choose_theta_seam(phase_frames[0]) if phase_frames else 0.0
                inner = float(space.get('inner', 0.0))
                outer = float(space.get('outer', inner + 1.0))
                bottom = float(space.get('bottom', -20.0))
                top = float(space.get('top', 20.0))
                r_mid = 0.5 * (inner + outer)
                smax = 2.0 * math.pi * r_mid
                for snapshot_slug, fraction, snapshot_label in SNAPSHOTS:
                    idx = pick_snapshot_index(time_min, fraction)
                    panels.append(
                        SnapshotPanel(
                            group=group,
                            regime=regime,
                            case=case,
                            force_value=force_value,
                            run_dir=run_dir,
                            run_path=str(run_path),
                            snapshot_slug=snapshot_slug,
                            snapshot_label=snapshot_label,
                            time_min=float(time_min[idx]) if time_min.size else float('nan'),
                            smax=smax,
                            zmin=bottom,
                            zmax=top,
                            seam=seam,
                            frame=phase_frames[idx],
                        )
                    )
    return panels


def frame_segments(frame: np.ndarray, seam: float, smax: float, r_mid: float) -> list[tuple[np.ndarray, np.ndarray]]:
    segments: list[tuple[np.ndarray, np.ndarray]] = []
    if frame.size == 0:
        return segments
    ids = np.unique(frame[:, 0].astype(int))
    for fid in ids:
        pts = frame[frame[:, 0] == fid]
        pts = pts[np.argsort(pts[:, 4])]
        theta = np.mod(np.arctan2(pts[:, 2], pts[:, 1]) - seam, 2.0 * math.pi)
        s = theta * r_mid
        z = pts[:, 3]
        if s.size < 2:
            continue
        jumps = np.where(np.abs(np.diff(s)) > 0.5 * smax)[0]
        start = 0
        if jumps.size:
            for j in jumps:
                seg = slice(start, j + 1)
                if (j + 1 - start) >= 2:
                    segments.append((s[seg], z[seg]))
                start = j + 1
        if (s.size - start) >= 2:
            segments.append((s[start:], z[start:]))
    return segments


def style(ax) -> None:
    ax.tick_params(direction='out', width=1.2, labelsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.1)
    ax.spines['bottom'].set_linewidth(1.1)


def draw_panel(ax, panel: SnapshotPanel | None) -> None:
    if panel is not None:
        r_mid = panel.smax / (2.0 * math.pi)
        for s, z in frame_segments(panel.frame, panel.seam, panel.smax, r_mid):
            ax.plot(s, z, color='black', lw=0.65, alpha=0.45)
        ax.set_xlim(0.0, panel.smax)
        ax.set_ylim(panel.zmin, panel.zmax)
    style(ax)


def draw_group_snapshot(group: str, group_label: str, snapshot_slug: str, snapshot_label: str, lookup: dict[tuple[str, str, str, str], SnapshotPanel]) -> plt.Figure:
    fig, axes = plt.subplots(len(XLINK_ORDER), len(FORCE_ORDER), figsize=(17.0, 9.8), sharex=True, sharey=True)
    for i, regime in enumerate(XLINK_ORDER):
        for j, (case, force_value) in enumerate(FORCE_ORDER):
            ax = axes[i, j]
            panel = lookup.get((group, snapshot_slug, regime, case))
            draw_panel(ax, panel)
            if i == 0:
                title = '0' if case == 'c0_m0' else f'{force_value:g}'
                ax.set_title(title, fontsize=15, pad=8)
            if j == 0:
                ax.text(0.02, 0.94, regime, transform=ax.transAxes, ha='left', va='top', fontsize=13, color='black', weight='bold')
                ax.set_ylabel('z (um)', fontsize=14)
            if i == len(XLINK_ORDER) - 1:
                ax.set_xlabel('Unwrapped circumference (um)', fontsize=13)
    fig.text(0.5, 0.985, f'{group_label}   snapshot {snapshot_label}', ha='center', va='top', fontsize=16)
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.08, top=0.91, wspace=0.08, hspace=0.12)
    return fig


def draw_group_regime_timeline(group: str, group_label: str, regime: str, lookup: dict[tuple[str, str, str, str], SnapshotPanel]) -> plt.Figure:
    fig, axes = plt.subplots(len(FORCE_ORDER), len(SNAPSHOTS), figsize=(14.8, 12.5), sharex=True, sharey=True)
    for i, (case, force_value) in enumerate(FORCE_ORDER):
        for j, (snapshot_slug, _fraction, snapshot_label) in enumerate(SNAPSHOTS):
            ax = axes[i, j]
            panel = lookup.get((group, snapshot_slug, regime, case))
            draw_panel(ax, panel)
            if i == 0:
                time_txt = snapshot_label
                if panel is not None:
                    time_txt = f'{snapshot_label}\n{panel.time_min:.1f} min'
                ax.set_title(time_txt, fontsize=12, pad=7)
            if j == 0:
                force_txt = '0' if case == 'c0_m0' else f'{force_value:g}'
                ax.text(0.02, 0.94, force_txt, transform=ax.transAxes, ha='left', va='top', fontsize=13, color='black', weight='bold')
                ax.set_ylabel('z (um)', fontsize=13)
            if i == len(FORCE_ORDER) - 1:
                ax.set_xlabel('Unwrapped circumference (um)', fontsize=12)
    fig.text(0.5, 0.988, f'{group_label}   xlink regime {regime}', ha='center', va='top', fontsize=16)
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.07, top=0.90, wspace=0.08, hspace=0.10)
    return fig


def write_manifest(panels: list[SnapshotPanel]) -> Path:
    path = OUT / 'panel_manifest.csv'
    rows = []
    for p in panels:
        rows.append(
            {
                'group': p.group,
                'xlink_regime': p.regime,
                'case': p.case,
                'force_value': p.force_value,
                'snapshot_slug': p.snapshot_slug,
                'snapshot_label': p.snapshot_label,
                'time_min': p.time_min,
                'run_dir': p.run_dir,
                'run_path': p.run_path,
            }
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    plt.rcParams.update({
        'font.size': 12,
        'axes.linewidth': 1.1,
        'savefig.facecolor': 'white',
        'figure.facecolor': 'white',
    })
    (OUT / 'plots').mkdir(parents=True, exist_ok=True)
    panels = collect_panels()
    lookup = {(p.group, p.snapshot_slug, p.regime, p.case): p for p in panels}
    written = [write_manifest(panels)]
    for group, group_label in GROUP_SPECS:
        for snapshot_slug, _fraction, snapshot_label in SNAPSHOTS:
            fig = draw_group_snapshot(group, group_label, snapshot_slug, snapshot_label, lookup)
            path = OUT / 'plots' / f'{group}_{snapshot_slug}_unwrapped_annulus_lines.png'
            fig.savefig(path, dpi=220, bbox_inches='tight')
            plt.close(fig)
            written.append(path)
        for regime in XLINK_ORDER:
            fig = draw_group_regime_timeline(group, group_label, regime, lookup)
            path = OUT / 'plots' / f'{group}_{XLINK_SLUG[regime]}_unwrapped_annulus_timeline_lines.png'
            fig.savefig(path, dpi=220, bbox_inches='tight')
            plt.close(fig)
            written.append(path)
    for path in written:
        print(path)


if __name__ == '__main__':
    main()
