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

ROT_MAP = ROOT / 'clu' / 'xlink_regime_map.csv'
FIX_MAP = ROOT / 'clu_fixed_global' / 'xlink_regime_map.csv'
MEDOID = ROOT / 'analysis' / 'results' / 'fixed_global_story' / 'transport_regime' / 'medoid_selection_summary.csv'
CONTROL = ROOT / 'control' / '2'
OUT = ROOT / 'analysis' / 'results' / 'fixed_global_story' / 'unwrapped_annulus'

FAMILIES = {'total480', 'mpc12'}
XLINK_ORDER = ['1:4', '1:8', '1:16']
XLINK_SLUG = {'1:4': '1to4', '1:8': '1to8', '1:16': '1to16'}
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
NBINS_S = 120
NBINS_Z = 84
SAMPLE_STEP = 0.35
HEATMAP_CMAP = 'viridis'


@dataclass
class Panel:
    family: str
    family_label: str
    regime: str
    row_model: str
    row_case: str
    row_label: str
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
    density: np.ndarray


def read_csv_rows(path: Path) -> list[dict]:
    return [dict(row) for row in csv.DictReader(path.open(encoding='utf-8'))]


def build_run_lookup(path: Path) -> dict[tuple[str, str, str, str], dict]:
    rows = read_csv_rows(path)
    lookup = {}
    for row in rows:
        if row.get('has_properties') == '0':
            continue
        lookup[(row['group'], row['case'], row['xlink_regime'], row['run_dir'])] = row
    return lookup


def load_medoid_lookup() -> dict[tuple[str, str, str, str], dict]:
    rows = read_csv_rows(MEDOID)
    return {(row['model'], row['group'], row['case'], row['xlink_regime']): row for row in rows}


def row_specs_for_family(family_spec: dict) -> list[dict]:
    rows = [{'model': 'control', 'group': 'control', 'case': 'c0_m0', 'label': '0 no motors'}]
    for case, label in family_spec['order']:
        rows.append({'model': 'rotatable', 'group': family_spec['group'], 'case': case, 'label': f'{label} rotatable'})
        rows.append({'model': 'fixed_global', 'group': family_spec['group'], 'case': case, 'label': f'{label} fixed-global'})
    return rows


def resolve_run_path(row_spec: dict, regime: str, medoid_lookup: dict, rot_lookup: dict, fix_lookup: dict) -> tuple[Path, str]:
    medoid = medoid_lookup.get((row_spec['model'], row_spec['group'], row_spec['case'], regime))
    if medoid is None:
        raise KeyError(f'missing medoid selection for {(row_spec["model"], row_spec["group"], row_spec["case"], regime)}')
    run_dir = medoid['medoid_run_dir']
    if row_spec['model'] == 'control':
        return CONTROL / run_dir, run_dir
    if row_spec['model'] == 'rotatable':
        row = rot_lookup.get((row_spec['group'], row_spec['case'], regime, run_dir))
    else:
        row = fix_lookup.get((row_spec['group'], row_spec['case'], regime, run_dir))
    if row is None:
        raise KeyError(f'missing run path for {(row_spec["group"], row_spec["case"], regime, run_dir)}')
    return Path(row['run_path']), run_dir


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


def sample_density_from_segments(segments: list[tuple[np.ndarray, np.ndarray]], *, smax: float, zmin: float, zmax: float) -> np.ndarray:
    if not segments:
        return np.zeros((NBINS_Z, NBINS_S), dtype=float)
    s_samples = []
    z_samples = []
    for s, z in segments:
        for i in range(len(s) - 1):
            ds = float(s[i + 1] - s[i])
            dz = float(z[i + 1] - z[i])
            seg_len = math.hypot(ds, dz)
            n = max(2, int(math.ceil(seg_len / SAMPLE_STEP)) + 1)
            t = np.linspace(0.0, 1.0, n)
            s_interp = s[i] + ds * t
            z_interp = z[i] + dz * t
            s_samples.append(s_interp)
            z_samples.append(z_interp)
    if not s_samples:
        return np.zeros((NBINS_Z, NBINS_S), dtype=float)
    s_all = np.concatenate(s_samples)
    z_all = np.concatenate(z_samples)
    hist, _, _ = np.histogram2d(z_all, s_all, bins=[NBINS_Z, NBINS_S], range=[[zmin, zmax], [0.0, smax]])
    hist = smooth2d(hist)
    mean = float(hist.mean())
    if mean > 0:
        hist = hist / mean
    return hist


def smooth2d(arr: np.ndarray, passes: int = 2) -> np.ndarray:
    kernel = np.array([1.0, 2.0, 1.0], dtype=float)
    kernel /= kernel.sum()
    out = np.asarray(arr, dtype=float)
    for _ in range(passes):
        out = np.apply_along_axis(lambda m: np.convolve(np.pad(m, 1, mode='edge'), kernel, mode='valid'), 1, out)
        out = np.apply_along_axis(lambda m: np.convolve(np.pad(m, 1, mode='edge'), kernel, mode='valid'), 0, out)
    return out


def collect_panels() -> tuple[dict[str, list[Panel]], dict[str, float]]:
    report_bin = speck.pick_report_bin()
    medoid_lookup = load_medoid_lookup()
    rot_lookup = build_run_lookup(ROT_MAP)
    fix_lookup = build_run_lookup(FIX_MAP)
    family_panels: dict[str, list[Panel]] = {}
    family_values: dict[str, list[float]] = {family: [] for family in FAMILIES}

    for family_spec in base.FAMILY_SPECS:
        family = family_spec['group']
        if family not in FAMILIES:
            continue
        panels: list[Panel] = []
        for regime in XLINK_ORDER:
            for row_spec in row_specs_for_family(family_spec):
                run_path, run_dir = resolve_run_path(row_spec, regime, medoid_lookup, rot_lookup, fix_lookup)
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
                    frame = phase_frames[idx]
                    segments = frame_segments(frame, seam, smax, r_mid)
                    density = sample_density_from_segments(segments, smax=smax, zmin=bottom, zmax=top)
                    family_values[family].extend(density[np.isfinite(density)].ravel().tolist())
                    panels.append(
                        Panel(
                            family=family,
                            family_label=family_spec['slug'],
                            regime=regime,
                            row_model=row_spec['model'],
                            row_case=row_spec['case'],
                            row_label=row_spec['label'],
                            run_dir=run_dir,
                            run_path=str(run_path),
                            snapshot_slug=snapshot_slug,
                            snapshot_label=snapshot_label,
                            time_min=float(time_min[idx]) if time_min.size else float('nan'),
                            smax=smax,
                            zmin=bottom,
                            zmax=top,
                            seam=seam,
                            frame=frame,
                            density=density,
                        )
                    )
        family_panels[family] = panels

    family_vmax = {}
    for family, vals in family_values.items():
        arr = np.asarray(vals, dtype=float)
        vmax = float(np.quantile(arr, 0.995)) if arr.size else 2.0
        family_vmax[family] = min(max(vmax, 1.5), 3.5)
    return family_panels, family_vmax


def panel_lookup(panels: list[Panel]) -> dict[tuple[str, str, str, str], Panel]:
    return {(p.regime, p.row_model, p.row_case, p.snapshot_slug): p for p in panels}


def style(ax) -> None:
    ax.tick_params(direction='out', width=1.0, labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.9)
    ax.spines['bottom'].set_linewidth(0.9)


def row_specs_for_lookup(family_spec: dict) -> list[dict]:
    return row_specs_for_family(family_spec)


def draw_line_panel(ax, panel: Panel | None) -> None:
    if panel is not None:
        r_mid = panel.smax / (2.0 * math.pi)
        for s, z in frame_segments(panel.frame, panel.seam, panel.smax, r_mid):
            ax.plot(s, z, color='black', lw=0.50, alpha=0.30)
        ax.set_xlim(0.0, panel.smax)
        ax.set_ylim(panel.zmin, panel.zmax)
    style(ax)


def draw_density_panel(ax, panel: Panel | None, *, vmax: float):
    im = None
    if panel is not None:
        im = ax.imshow(
            panel.density,
            origin='lower',
            aspect='auto',
            extent=[0.0, panel.smax, panel.zmin, panel.zmax],
            cmap=HEATMAP_CMAP,
            vmin=0.0,
            vmax=vmax,
            interpolation='nearest',
        )
    style(ax)
    return im


def make_timeline_figure(family_spec: dict, regime: str, panels: list[Panel], *, mode: str, vmax: float | None = None) -> tuple[plt.Figure, object | None]:
    rows = row_specs_for_lookup(family_spec)
    lookup = panel_lookup(panels)
    nrows = len(rows)
    ncols = len(SNAPSHOTS)
    fig, axes = plt.subplots(nrows, ncols, figsize=(12.5, 2.0 * nrows), sharex=True, sharey=True)
    if nrows == 1:
        axes = np.array([axes])
    im = None
    for i, row_spec in enumerate(rows):
        for j, (snapshot_slug, _fraction, snapshot_label) in enumerate(SNAPSHOTS):
            ax = axes[i, j]
            panel = lookup.get((regime, row_spec['model'], row_spec['case'], snapshot_slug))
            if mode == 'line':
                draw_line_panel(ax, panel)
            else:
                im = draw_density_panel(ax, panel, vmax=vmax)
            if i == 0:
                time_txt = snapshot_label
                if panel is not None:
                    time_txt = f'{snapshot_label}\n{panel.time_min:.1f} min'
                ax.set_title(time_txt, fontsize=11, pad=6)
            if j == 0:
                ax.set_ylabel('z (um)', fontsize=10)
                ax.text(0.01, 0.94, row_spec['label'], transform=ax.transAxes, ha='left', va='top', fontsize=8.5, color='white' if mode == 'density' and panel is not None else 'black', weight='bold')
            else:
                ax.set_yticklabels([])
            if i == nrows - 1:
                ax.set_xlabel('Unwrapped circumference (um)', fontsize=9.5)
            else:
                ax.set_xticklabels([])
    fig.text(0.5, 0.992, f'{family_spec["group"]}   xlink regime {regime}', ha='center', va='top', fontsize=15)
    right = 0.965 if mode == 'line' else 0.93
    fig.subplots_adjust(left=0.10, right=right, bottom=0.04, top=0.965, wspace=0.06, hspace=0.10)
    return fig, im


def write_manifest(family_panels: dict[str, list[Panel]]) -> Path:
    path = OUT / 'panel_manifest.csv'
    rows = []
    for family, panels in family_panels.items():
        for p in panels:
            rows.append(
                {
                    'family': family,
                    'xlink_regime': p.regime,
                    'row_model': p.row_model,
                    'row_case': p.row_case,
                    'row_label': p.row_label,
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


def save_colorbar(path: Path, vmax: float) -> None:
    fig = plt.figure(figsize=(1.6, 5.0), facecolor='white')
    cax = fig.add_axes([0.35, 0.05, 0.30, 0.9])
    norm = plt.Normalize(vmin=0.0, vmax=vmax)
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=HEATMAP_CMAP), cax=cax)
    cbar.set_label('Normalized areal density', fontsize=13)
    cbar.ax.tick_params(labelsize=10, width=0.9)
    cbar.outline.set_linewidth(1.0)
    fig.savefig(path, dpi=220, bbox_inches='tight')
    plt.close(fig)


def main() -> None:
    plt.rcParams.update({
        'font.size': 10,
        'axes.linewidth': 0.9,
        'savefig.facecolor': 'white',
        'figure.facecolor': 'white',
    })
    (OUT / 'line_plots').mkdir(parents=True, exist_ok=True)
    (OUT / 'heatmap_plots').mkdir(parents=True, exist_ok=True)

    family_panels, family_vmax = collect_panels()
    written = [write_manifest(family_panels)]

    family_specs = [f for f in base.FAMILY_SPECS if f['group'] in FAMILIES]
    for family_spec in family_specs:
        family = family_spec['group']
        panels = family_panels[family]
        vmax = family_vmax[family]
        save_colorbar(OUT / 'heatmap_plots' / f'{family}_colorbar_only.png', vmax)
        written.append(OUT / 'heatmap_plots' / f'{family}_colorbar_only.png')
        for regime in XLINK_ORDER:
            fig, _ = make_timeline_figure(family_spec, regime, panels, mode='line')
            path = OUT / 'line_plots' / f'{family}_{XLINK_SLUG[regime]}_timeline_lines.png'
            fig.savefig(path, dpi=220, bbox_inches='tight')
            plt.close(fig)
            written.append(path)

            fig, im = make_timeline_figure(family_spec, regime, panels, mode='density', vmax=vmax)
            path = OUT / 'heatmap_plots' / f'{family}_{XLINK_SLUG[regime]}_timeline_heatmap.png'
            fig.savefig(path, dpi=220, bbox_inches='tight')
            written.append(path)
            if im is not None:
                cax = fig.add_axes([0.94, 0.14, 0.014, 0.72])
                cbar = fig.colorbar(im, cax=cax)
                cbar.set_label('Normalized areal density', fontsize=10)
                cbar.ax.tick_params(labelsize=8, width=0.8)
                with_cb = OUT / 'heatmap_plots' / f'{family}_{XLINK_SLUG[regime]}_timeline_heatmap_with_colorbar.png'
                fig.savefig(with_cb, dpi=220, bbox_inches='tight')
                written.append(with_cb)
            plt.close(fig)
    for path in written:
        print(path)


if __name__ == '__main__':
    main()
