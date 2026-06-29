#!/usr/bin/env python3
from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path('/home/thekaka/project/cytosim')
BUNDLE = ROOT / 'analysis' / 'results' / 'full_story_2026-04-12'


def copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def copy_many(src_paths: list[Path], dst_dir: Path) -> None:
    dst_dir.mkdir(parents=True, exist_ok=True)
    for src in src_paths:
        if src.exists():
            copy(src, dst_dir / src.name)


def main() -> None:
    if BUNDLE.exists():
        shutil.rmtree(BUNDLE)

    random_transport_dir = BUNDLE / 'random_orientation' / 'transport_metrics'
    random_unwrapped_heatmaps = BUNDLE / 'random_orientation' / 'unwrapped_heatmaps'
    random_unwrapped_lines = BUNDLE / 'random_orientation' / 'unwrapped_lines'
    random_tables = BUNDLE / 'random_orientation' / 'tables'

    aligned_transport_dir = BUNDLE / 'aligned_minusz' / 'transport_metrics'
    aligned_unwrapped_heatmaps = BUNDLE / 'aligned_minusz' / 'unwrapped_heatmaps'
    aligned_unwrapped_lines = BUNDLE / 'aligned_minusz' / 'unwrapped_lines'
    aligned_tables = BUNDLE / 'aligned_minusz' / 'tables'

    notes_dir = BUNDLE / 'notes'

    random_transport = [
        ROOT / 'analysis/results/fixed_global_story/presentation_plots/vz_abs_mean_with_control.png',
        ROOT / 'analysis/results/fixed_global_story/transport_regime/plots/auc_vz_abs.png',
        ROOT / 'analysis/results/fixed_global_story/transport_regime/plots/peak_vz_abs.png',
        ROOT / 'analysis/results/fixed_global_story/transport_regime/plots/active_transport_fraction.png',
        ROOT / 'analysis/results/fixed_global_story/transport_regime/plots/abs_net_z_displacement.png',
        ROOT / 'analysis/results/fixed_global_story/transport_regime/plots/time_to_peak_min.png',
        ROOT / 'analysis/results/fixed_global_story/transport_regime/plots/total480_window_mean_vz_abs.png',
        ROOT / 'analysis/results/fixed_global_story/transport_regime/plots/mpc12_window_mean_vz_abs.png',
    ]
    copy_many(random_transport, random_transport_dir)

    aligned_transport = [
        ROOT / 'analysis/results/init_minusz_comparison/transport_regime/plots/mean_vz_abs.png',
        ROOT / 'analysis/results/init_minusz_comparison/transport_regime/plots/auc_vz_abs.png',
        ROOT / 'analysis/results/init_minusz_comparison/transport_regime/plots/peak_vz_abs.png',
        ROOT / 'analysis/results/init_minusz_comparison/transport_regime/plots/active_transport_fraction.png',
        ROOT / 'analysis/results/init_minusz_comparison/transport_regime/plots/abs_net_z_displacement.png',
        ROOT / 'analysis/results/init_minusz_comparison/transport_regime/plots/time_to_peak_min.png',
        ROOT / 'analysis/results/init_minusz_comparison/transport_regime/plots/total480_window_mean_vz_abs.png',
        ROOT / 'analysis/results/init_minusz_comparison/transport_regime/plots/mpc12_window_mean_vz_abs.png',
    ]
    copy_many(aligned_transport, aligned_transport_dir)

    for family in ['total480', 'mpc12']:
        for regime in ['1to4', '1to8', '1to16']:
            copy(
                ROOT / f'analysis/results/fixed_global_story/unwrapped_annulus/heatmap_plots/{family}_{regime}_timeline_heatmap_with_colorbar.png',
                random_unwrapped_heatmaps / f'{family}_{regime}_timeline_heatmap_with_colorbar.png',
            )
            copy(
                ROOT / f'analysis/results/fixed_global_story/unwrapped_annulus/line_plots/{family}_{regime}_timeline_lines.png',
                random_unwrapped_lines / f'{family}_{regime}_timeline_lines.png',
            )
            copy(
                ROOT / f'analysis/results/init_minusz_comparison/unwrapped_annulus/heatmap_plots/{family}_{regime}_timeline_heatmap_with_colorbar.png',
                aligned_unwrapped_heatmaps / f'{family}_{regime}_timeline_heatmap_with_colorbar.png',
            )
            copy(
                ROOT / f'analysis/results/init_minusz_comparison/unwrapped_annulus/line_plots/{family}_{regime}_timeline_lines.png',
                aligned_unwrapped_lines / f'{family}_{regime}_timeline_lines.png',
            )
        copy(
            ROOT / f'analysis/results/fixed_global_story/unwrapped_annulus/heatmap_plots/{family}_colorbar_only.png',
            random_unwrapped_heatmaps / f'{family}_colorbar_only.png',
        )
        copy(
            ROOT / f'analysis/results/init_minusz_comparison/unwrapped_annulus/heatmap_plots/{family}_colorbar_only.png',
            aligned_unwrapped_heatmaps / f'{family}_colorbar_only.png',
        )

    copy_many([
        ROOT / 'analysis/results/fixed_global_story/transport_regime/condition_metrics.csv',
        ROOT / 'analysis/results/fixed_global_story/transport_regime/run_metrics.csv',
    ], random_tables)

    copy_many([
        ROOT / 'analysis/results/init_minusz_comparison/transport_regime/condition_metrics.csv',
        ROOT / 'analysis/results/init_minusz_comparison/transport_regime/run_metrics.csv',
        ROOT / 'analysis/results/init_minusz_comparison/unwrapped_annulus/panel_manifest.csv',
    ], aligned_tables)

    notes_dir.mkdir(parents=True, exist_ok=True)
    (notes_dir / 'story_summary.md').write_text(
        """# Full Transport Story Bundle

## Random orientation: fixed vs rotatable
- Rotatable clusters consistently outperform fixed-global clusters across mean transport, peak transport, AUC, active-transport fraction, and windowed transport metrics.
- The random-orientation unwrapped annulus panels show stronger bright-track emergence in the rotatable cases, especially as motor loading increases.
- This bundle preserves the random-initialization story in the same plot language used for the aligned campaign.

## Aligned minus-z orientation: fixed vs rotatable
- Once the network is pre-aligned, both models transport more strongly and more cleanly than in the random-initialization case.
- In the total480 family, repartitioning a fixed total motor budget across 10 to 80 clusters has only a modest effect on mean transport and AUC.
- In the mpc12 family, transport rises strongly with cluster count because total motor number increases across the series.
- Rotatable still exceeds fixed-global under aligned conditions, but the aligned campaign shows that motor-budget and polarity are the dominant control knobs.

## Next comparison layer
- The natural next step is a 4-way comparison: random fixed, random rotatable, aligned fixed, aligned rotatable.
- The tables copied into this bundle are separated by campaign so we can build that comparison directly without another data hunt.
""",
        encoding='utf-8',
    )
    print(BUNDLE)


if __name__ == '__main__':
    main()
