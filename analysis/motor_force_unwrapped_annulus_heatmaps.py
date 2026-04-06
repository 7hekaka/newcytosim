#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / 'analysis') not in sys.path:
    sys.path.insert(0, str(ROOT / 'analysis'))

from analysis import unwrapped_annulus_line_snapshots as linebase

OUT = ROOT / 'analysis' / 'results' / 'motor_force_sweep' / 'unwrapped_annulus'
NBINS_S = 120
NBINS_Z = 84
SAMPLE_STEP = 0.35
HEATMAP_CMAP = 'viridis'


def smooth2d(arr: np.ndarray, passes: int = 2) -> np.ndarray:
    kernel = np.array([1.0, 2.0, 1.0], dtype=float)
    kernel /= kernel.sum()
    out = np.asarray(arr, dtype=float)
    for _ in range(passes):
        out = np.apply_along_axis(lambda m: np.convolve(np.pad(m, 1, mode='edge'), kernel, mode='valid'), 1, out)
        out = np.apply_along_axis(lambda m: np.convolve(np.pad(m, 1, mode='edge'), kernel, mode='valid'), 0, out)
    return out


def sample_density_from_panel(panel) -> np.ndarray:
    r_mid = panel.smax / (2.0 * math.pi)
    segments = linebase.frame_segments(panel.frame, panel.seam, panel.smax, r_mid)
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
            s_samples.append(s[i] + ds * t)
            z_samples.append(z[i] + dz * t)
    if not s_samples:
        return np.zeros((NBINS_Z, NBINS_S), dtype=float)
    s_all = np.concatenate(s_samples)
    z_all = np.concatenate(z_samples)
    hist, _, _ = np.histogram2d(z_all, s_all, bins=[NBINS_Z, NBINS_S], range=[[panel.zmin, panel.zmax], [0.0, panel.smax]])
    hist = smooth2d(hist)
    mean = float(hist.mean())
    if mean > 0:
        hist = hist / mean
    return hist


def panel_lookup(panels):
    return {(p.group, p.snapshot_slug, p.regime, p.case): p for p in panels}


def density_lookup(panels):
    data = {}
    vals = []
    for p in panels:
        dens = sample_density_from_panel(p)
        data[(p.group, p.snapshot_slug, p.regime, p.case)] = dens
        vals.extend(dens[np.isfinite(dens)].ravel().tolist())
    vmax = float(np.quantile(np.asarray(vals, dtype=float), 0.995)) if vals else 2.0
    vmax = min(max(vmax, 1.5), 3.5)
    return data, vmax


def style(ax):
    ax.tick_params(direction='out', width=1.2, labelsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.1)
    ax.spines['bottom'].set_linewidth(1.1)


def draw_group_snapshot(group, group_label, snapshot_slug, snapshot_label, panel_lu, dens_lu, vmax):
    fig, axes = plt.subplots(len(linebase.XLINK_ORDER), len(linebase.FORCE_ORDER), figsize=(17.0, 9.8), sharex=True, sharey=True)
    im = None
    for i, regime in enumerate(linebase.XLINK_ORDER):
        for j, (case, force_value) in enumerate(linebase.FORCE_ORDER):
            ax = axes[i, j]
            panel = panel_lu.get((group, snapshot_slug, regime, case))
            density = dens_lu.get((group, snapshot_slug, regime, case))
            if panel is not None and density is not None:
                im = ax.imshow(density, origin='lower', aspect='auto', extent=[0.0, panel.smax, panel.zmin, panel.zmax], cmap=HEATMAP_CMAP, vmin=0.0, vmax=vmax, interpolation='nearest')
                ax.set_xlim(0.0, panel.smax)
                ax.set_ylim(panel.zmin, panel.zmax)
            if i == 0:
                title = '0' if case == 'c0_m0' else f'{force_value:g}'
                ax.set_title(title, fontsize=15, pad=8)
            if j == 0:
                ax.text(0.02, 0.94, regime, transform=ax.transAxes, ha='left', va='top', fontsize=13, color='white', weight='bold')
                ax.set_ylabel('z (um)', fontsize=14)
            if i == len(linebase.XLINK_ORDER) - 1:
                ax.set_xlabel('Unwrapped circumference (um)', fontsize=13)
            style(ax)
    fig.text(0.5, 0.985, f'{group_label}   snapshot {snapshot_label}', ha='center', va='top', fontsize=16)
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.08, top=0.91, wspace=0.08, hspace=0.12)
    return fig, im


def draw_group_regime_timeline(group, group_label, regime, panel_lu, dens_lu, vmax):
    fig, axes = plt.subplots(len(linebase.FORCE_ORDER), len(linebase.SNAPSHOTS), figsize=(14.8, 12.5), sharex=True, sharey=True)
    im = None
    for i, (case, force_value) in enumerate(linebase.FORCE_ORDER):
        for j, (snapshot_slug, _fraction, snapshot_label) in enumerate(linebase.SNAPSHOTS):
            ax = axes[i, j]
            panel = panel_lu.get((group, snapshot_slug, regime, case))
            density = dens_lu.get((group, snapshot_slug, regime, case))
            if panel is not None and density is not None:
                im = ax.imshow(density, origin='lower', aspect='auto', extent=[0.0, panel.smax, panel.zmin, panel.zmax], cmap=HEATMAP_CMAP, vmin=0.0, vmax=vmax, interpolation='nearest')
                ax.set_xlim(0.0, panel.smax)
                ax.set_ylim(panel.zmin, panel.zmax)
            if i == 0:
                title = snapshot_label
                if panel is not None:
                    title = f'{snapshot_label}\n{panel.time_min:.1f} min'
                ax.set_title(title, fontsize=12, pad=7)
            if j == 0:
                force_txt = '0' if case == 'c0_m0' else f'{force_value:g}'
                ax.text(0.02, 0.94, force_txt, transform=ax.transAxes, ha='left', va='top', fontsize=13, color='white', weight='bold')
                ax.set_ylabel('z (um)', fontsize=13)
            if i == len(linebase.FORCE_ORDER) - 1:
                ax.set_xlabel('Unwrapped circumference (um)', fontsize=12)
            style(ax)
    fig.text(0.5, 0.988, f'{group_label}   xlink regime {regime}', ha='center', va='top', fontsize=16)
    fig.subplots_adjust(left=0.07, right=0.92, bottom=0.07, top=0.90, wspace=0.08, hspace=0.10)
    return fig, im


def save_colorbar(path: Path, vmax: float):
    fig = plt.figure(figsize=(1.5, 5.2), facecolor='white')
    ax = fig.add_axes([0.35, 0.05, 0.32, 0.9])
    norm = plt.Normalize(vmin=0.0, vmax=vmax)
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=HEATMAP_CMAP), cax=ax)
    cbar.set_label('Normalized areal density', fontsize=16)
    cbar.ax.tick_params(labelsize=12, width=1.1)
    cbar.outline.set_linewidth(1.2)
    fig.savefig(path, dpi=220, bbox_inches='tight')
    plt.close(fig)


def write_manifest(path: Path, panels):
    rows = []
    for p in panels:
        rows.append({'group': p.group, 'xlink_regime': p.regime, 'case': p.case, 'force_value': p.force_value, 'snapshot_slug': p.snapshot_slug, 'snapshot_label': p.snapshot_label, 'time_min': p.time_min, 'run_dir': p.run_dir, 'run_path': p.run_path})
    with path.open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main():
    plt.rcParams.update({'font.size': 12, 'axes.linewidth': 1.1, 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
    (OUT / 'plots').mkdir(parents=True, exist_ok=True)
    panels = linebase.collect_panels()
    panel_lu = panel_lookup(panels)
    dens_lu, vmax = density_lookup(panels)
    write_manifest(OUT / 'panel_manifest.csv', panels)
    save_colorbar(OUT / 'plots' / 'unwrapped_annulus_colorbar_only.png', vmax)
    for group, group_label in linebase.GROUP_SPECS:
        for snapshot_slug, _fraction, snapshot_label in linebase.SNAPSHOTS:
            fig, im = draw_group_snapshot(group, group_label, snapshot_slug, snapshot_label, panel_lu, dens_lu, vmax)
            png = OUT / 'plots' / f'{group}_{snapshot_slug}_unwrapped_annulus.png'
            fig.savefig(png, dpi=220, bbox_inches='tight')
            if im is not None:
                cax = fig.add_axes([0.988, 0.17, 0.012, 0.68])
                cb = fig.colorbar(im, cax=cax)
                cb.set_label('Normalized areal density', fontsize=13)
                cb.ax.tick_params(labelsize=10, width=1.0)
                with_cb = OUT / 'plots' / f'{group}_{snapshot_slug}_unwrapped_annulus_with_colorbar.png'
                fig.savefig(with_cb, dpi=220, bbox_inches='tight')
            plt.close(fig)
        for regime in linebase.XLINK_ORDER:
            fig, im = draw_group_regime_timeline(group, group_label, regime, panel_lu, dens_lu, vmax)
            png = OUT / 'plots' / f'{group}_{linebase.XLINK_SLUG[regime]}_unwrapped_annulus_timeline_heatmap.png'
            fig.savefig(png, dpi=220, bbox_inches='tight')
            if im is not None:
                cax = fig.add_axes([0.94, 0.14, 0.014, 0.72])
                cb = fig.colorbar(im, cax=cax)
                cb.set_label('Normalized areal density', fontsize=10)
                cb.ax.tick_params(labelsize=8, width=0.8)
                with_cb = OUT / 'plots' / f'{group}_{linebase.XLINK_SLUG[regime]}_unwrapped_annulus_timeline_heatmap_with_colorbar.png'
                fig.savefig(with_cb, dpi=220, bbox_inches='tight')
            plt.close(fig)


if __name__ == '__main__':
    main()
