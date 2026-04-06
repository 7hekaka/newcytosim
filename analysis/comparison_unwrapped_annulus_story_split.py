#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / 'analysis') not in sys.path:
    sys.path.insert(0, str(ROOT / 'analysis'))

from analysis import comparison_unwrapped_annulus_story as story

OUT = story.OUT
PARTS = [
    ('low_counts', 3),
    ('high_counts', 2),
]


def split_row_specs(family_spec: dict, part_name: str) -> list[dict]:
    control = {'model': 'control', 'group': 'control', 'case': 'c0_m0', 'label': '0 no motors'}
    order = list(family_spec['order'])
    if part_name == 'low_counts':
        chosen = order[:3]
    else:
        chosen = order[3:]
    rows = [control]
    for case, label in chosen:
        rows.append({'model': 'rotatable', 'group': family_spec['group'], 'case': case, 'label': f'{label} rotatable'})
        rows.append({'model': 'fixed_global', 'group': family_spec['group'], 'case': case, 'label': f'{label} fixed-global'})
    return rows


def make_timeline_figure(family_spec: dict, regime: str, panels: list, *, mode: str, vmax: float | None, part_name: str):
    rows = split_row_specs(family_spec, part_name)
    lookup = story.panel_lookup(panels)
    nrows = len(rows)
    ncols = len(story.SNAPSHOTS)
    fig, axes = plt.subplots(nrows, ncols, figsize=(12.5, 2.15 * nrows), sharex=True, sharey=True)
    if nrows == 1:
        axes = np.array([axes])
    im = None
    for i, row_spec in enumerate(rows):
        for j, (snapshot_slug, _fraction, snapshot_label) in enumerate(story.SNAPSHOTS):
            ax = axes[i, j]
            panel = lookup.get((regime, row_spec['model'], row_spec['case'], snapshot_slug))
            if mode == 'line':
                story.draw_line_panel(ax, panel)
            else:
                im = story.draw_density_panel(ax, panel, vmax=vmax)
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
    pretty = 'low counts' if part_name == 'low_counts' else 'high counts'
    fig.text(0.5, 0.992, f'{family_spec["group"]}   xlink regime {regime}   {pretty}', ha='center', va='top', fontsize=15)
    right = 0.965 if mode == 'line' else 0.93
    fig.subplots_adjust(left=0.10, right=right, bottom=0.045, top=0.965, wspace=0.06, hspace=0.10)
    return fig, im


def main():
    plt.rcParams.update({'font.size': 10, 'axes.linewidth': 0.9, 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
    (OUT / 'line_plots_split').mkdir(parents=True, exist_ok=True)
    (OUT / 'heatmap_plots_split').mkdir(parents=True, exist_ok=True)
    family_panels, family_vmax = story.collect_panels()
    family_specs = [f for f in story.base.FAMILY_SPECS if f['group'] in story.FAMILIES]
    for family_spec in family_specs:
        family = family_spec['group']
        panels = family_panels[family]
        vmax = family_vmax[family]
        for regime in story.XLINK_ORDER:
            for part_name, _n in PARTS:
                fig, _ = make_timeline_figure(family_spec, regime, panels, mode='line', vmax=None, part_name=part_name)
                path = OUT / 'line_plots_split' / f'{family}_{story.XLINK_SLUG[regime]}_{part_name}_timeline_lines.png'
                fig.savefig(path, dpi=220, bbox_inches='tight')
                plt.close(fig)
                fig, im = make_timeline_figure(family_spec, regime, panels, mode='density', vmax=vmax, part_name=part_name)
                path = OUT / 'heatmap_plots_split' / f'{family}_{story.XLINK_SLUG[regime]}_{part_name}_timeline_heatmap.png'
                fig.savefig(path, dpi=220, bbox_inches='tight')
                if im is not None:
                    cax = fig.add_axes([0.94, 0.14, 0.014, 0.72])
                    cb = fig.colorbar(im, cax=cax)
                    cb.set_label('Normalized areal density', fontsize=10)
                    cb.ax.tick_params(labelsize=8, width=0.8)
                    with_cb = OUT / 'heatmap_plots_split' / f'{family}_{story.XLINK_SLUG[regime]}_{part_name}_timeline_heatmap_with_colorbar.png'
                    fig.savefig(with_cb, dpi=220, bbox_inches='tight')
                plt.close(fig)


if __name__ == '__main__':
    main()
