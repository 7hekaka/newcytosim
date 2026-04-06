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

from analysis import plot_timecourse_panels as base
from analysis import run_speckle_timecourse_batch as speck

MAP = ROOT / 'clu_fixed_global_force_sweep' / 'xlink_regime_map.csv'
MEDOID = ROOT / 'analysis' / 'results' / 'motor_force_sweep' / 'medoid_selection_summary.csv'
CONTROL = ROOT / 'control' / '2'
OUT = ROOT / 'analysis' / 'results' / 'motor_force_sweep' / 'unwrapped_annulus'

XLINK_ORDER = ['1:4', '1:8', '1:16']
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
NBINS_S = 96
NBINS_Z = 80
SEAM_BINS = 180
COLOR_CMAP = 'viridis'


@dataclass
class PanelData:
    group: str
    regime: str
    case: str
    force_value: float
    run_dir: str
    run_path: str
    snapshot_slug: str
    snapshot_label: str
    target_fraction: float
    time_min: float
    density: np.ndarray
    smax: float
    zmin: float
    zmax: float


def pick_report_bin() -> Path:
    return speck.pick_report_bin()


def load_run_lookup() -> dict[tuple[str, str, str, str], dict]:
    rows = list(csv.DictReader(MAP.open(encoding='utf-8')))
    lookup = {}
    for row in rows:
        if row.get('has_outputs') != '1':
            continue
        lookup[(row['group'], row['case'], row['xlink_regime'], row['run_dir'])] = dict(row)
    return lookup


def load_medoid_rows() -> list[dict]:
    return [dict(row) for row in csv.DictReader(MEDOID.open(encoding='utf-8'))]


def resolve_run_path(group: str, case: str, regime: str, run_dir: str, lookup: dict) -> Path:
    if case == 'c0_m0':
        return CONTROL / run_dir
    row = lookup.get((group, case, regime, run_dir))
    if row is None:
        raise KeyError(f'missing run lookup for {(group, case, regime, run_dir)}')
    return Path(row['run_path'])


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


def frame_to_density(frame: np.ndarray, *, seam: float, inner: float, outer: float, bottom: float, top: float, nbins_s: int, nbins_z: int) -> tuple[np.ndarray, float]:
    if frame.size == 0:
        r_mid = 0.5 * (inner + outer)
        smax = 2.0 * math.pi * r_mid
        return np.zeros((nbins_z, nbins_s), dtype=float), smax
    x = frame[:, 1]
    y = frame[:, 2]
    z = frame[:, 3]
    theta = np.mod(np.arctan2(y, x) - seam, 2.0 * math.pi)
    r_mid = 0.5 * (inner + outer)
    smax = 2.0 * math.pi * r_mid
    s = theta * r_mid
    hist, z_edges, s_edges = np.histogram2d(z, s, bins=[nbins_z, nbins_s], range=[[bottom, top], [0.0, smax]])
    mean = float(hist.mean())
    if mean > 0:
        hist = hist / mean
    return hist, smax


def collect_panels() -> tuple[list[PanelData], float]:
    report_bin = pick_report_bin()
    run_lookup = load_run_lookup()
    medoid_rows = load_medoid_rows()
    wanted = {(group, case, regime) for group, _ in GROUP_SPECS for case, _ in FORCE_ORDER for regime in XLINK_ORDER}
    medoid_lookup = {}
    for row in medoid_rows:
        key = (row['group'], row['case'], row['xlink_regime'])
        if key in wanted:
            medoid_lookup[key] = row

    panels: list[PanelData] = []
    all_values = []
    for group, _group_label in GROUP_SPECS:
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
                for snapshot_slug, fraction, snapshot_label in SNAPSHOTS:
                    idx = pick_snapshot_index(time_min, fraction)
                    density, smax = frame_to_density(
                        phase_frames[idx],
                        seam=seam,
                        inner=inner,
                        outer=outer,
                        bottom=bottom,
                        top=top,
                        nbins_s=NBINS_S,
                        nbins_z=NBINS_Z,
                    )
                    all_values.extend(density[np.isfinite(density)].ravel().tolist())
                    panels.append(
                        PanelData(
                            group=group,
                            regime=regime,
                            case=case,
                            force_value=force_value,
                            run_dir=run_dir,
                            run_path=str(run_path),
                            snapshot_slug=snapshot_slug,
                            snapshot_label=snapshot_label,
                            target_fraction=fraction,
                            time_min=float(time_min[idx]) if time_min.size else float('nan'),
                            density=density,
                            smax=smax,
                            zmin=bottom,
                            zmax=top,
                        )
                    )
    vmax = float(np.quantile(np.asarray(all_values, dtype=float), 0.995)) if all_values else 2.0
    vmax = min(max(vmax, 1.5), 3.0)
    return panels, vmax


def panel_lookup(panels: list[PanelData]) -> dict[tuple[str, str, str, str], PanelData]:
    return {(p.group, p.snapshot_slug, p.regime, p.case): p for p in panels}


def style(ax) -> None:
    ax.tick_params(direction='out', width=1.2, labelsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.1)
    ax.spines['bottom'].set_linewidth(1.1)


def draw_group_snapshot(group: str, group_label: str, snapshot_slug: str, snapshot_label: str, lookup: dict, vmax: float) -> tuple[plt.Figure, np.ndarray]:
    fig, axes = plt.subplots(len(XLINK_ORDER), len(FORCE_ORDER), figsize=(17.0, 9.8), sharex=True, sharey=True)
    image = None
    for i, regime in enumerate(XLINK_ORDER):
        for j, (case, force_value) in enumerate(FORCE_ORDER):
            ax = axes[i, j]
            panel = lookup.get((group, snapshot_slug, regime, case))
            if panel is not None:
                image = ax.imshow(
                    panel.density,
                    origin='lower',
                    aspect='auto',
                    extent=[0.0, panel.smax, panel.zmin, panel.zmax],
                    cmap=COLOR_CMAP,
                    vmin=0.0,
                    vmax=vmax,
                    interpolation='nearest',
                )
            ax.set_xlim(left=0.0)
            if i == 0:
                title = '0' if case == 'c0_m0' else f'{force_value:g}'
                ax.set_title(title, fontsize=15, pad=8)
            if j == 0:
                ax.text(0.02, 0.94, regime, transform=ax.transAxes, ha='left', va='top', fontsize=13, color='white', weight='bold')
                ax.set_ylabel('z (um)', fontsize=14)
            if i == len(XLINK_ORDER) - 1:
                ax.set_xlabel('Unwrapped circumference (um)', fontsize=13)
            style(ax)
    fig.text(0.5, 0.985, f'{group_label}   snapshot {snapshot_label}', ha='center', va='top', fontsize=17)
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.08, top=0.92, wspace=0.08, hspace=0.12)
    return fig, image


def save_colorbar(vmax: float) -> Path:
    path = OUT / 'plots' / 'unwrapped_annulus_colorbar_only.png'
    fig = plt.figure(figsize=(1.5, 5.2))
    ax = fig.add_axes([0.35, 0.05, 0.32, 0.9])
    norm = plt.Normalize(vmin=0.0, vmax=vmax)
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=COLOR_CMAP), cax=ax)
    cbar.set_label('Normalized areal density', fontsize=16)
    cbar.ax.tick_params(labelsize=12, width=1.1)
    cbar.outline.set_linewidth(1.2)
    fig.savefig(path, dpi=220, bbox_inches='tight')
    plt.close(fig)
    return path


def write_manifest(panels: list[PanelData]) -> Path:
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
                'target_fraction': p.target_fraction,
                'time_min': p.time_min,
                'run_dir': p.run_dir,
                'run_path': p.run_path,
            }
        )
    fieldnames = list(rows[0].keys()) if rows else []
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
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
    panels, vmax = collect_panels()
    lookup = panel_lookup(panels)
    manifest = write_manifest(panels)
    written = [manifest]
    for group, group_label in GROUP_SPECS:
        for snapshot_slug, _fraction, snapshot_label in SNAPSHOTS:
            fig, image = draw_group_snapshot(group, group_label, snapshot_slug, snapshot_label, lookup, vmax)
            png = OUT / 'plots' / f'{group}_{snapshot_slug}_unwrapped_annulus.png'
            fig.savefig(png, dpi=220, bbox_inches='tight')
            written.append(png)
            if image is not None:
                cfig = fig
                cax = cfig.add_axes([0.988, 0.17, 0.012, 0.68])
                cb = cfig.colorbar(image, cax=cax)
                cb.set_label('Normalized areal density', fontsize=13)
                cb.ax.tick_params(labelsize=10, width=1.0)
                with_cb = OUT / 'plots' / f'{group}_{snapshot_slug}_unwrapped_annulus_with_colorbar.png'
                cfig.savefig(with_cb, dpi=220, bbox_inches='tight')
                written.append(with_cb)
            plt.close(fig)
    written.append(save_colorbar(vmax))
    for path in written:
        print(path)


if __name__ == '__main__':
    main()
