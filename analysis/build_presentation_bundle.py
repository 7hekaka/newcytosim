#!/usr/bin/env python3
from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUNDLE = ROOT / 'analysis' / 'results' / 'presentation_ready_2026-04-06'


def copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def main() -> None:
    if BUNDLE.exists():
        shutil.rmtree(BUNDLE)
    (BUNDLE / 'fixed_vs_rotatable' / 'transport_metrics').mkdir(parents=True, exist_ok=True)
    (BUNDLE / 'fixed_vs_rotatable' / 'unwrapped_heatmaps').mkdir(parents=True, exist_ok=True)
    (BUNDLE / 'motor_force_sweep' / 'transport_metrics').mkdir(parents=True, exist_ok=True)
    (BUNDLE / 'motor_force_sweep' / 'unwrapped_heatmaps').mkdir(parents=True, exist_ok=True)
    (BUNDLE / 'total_motor_sweep').mkdir(parents=True, exist_ok=True)
    (BUNDLE / 'notes').mkdir(parents=True, exist_ok=True)

    fixed_transport = [
        'analysis/results/fixed_global_story/presentation_plots/vz_abs_mean_with_control.png',
        'analysis/results/fixed_global_story/presentation_plots/vz_std_presentation.png',
        'analysis/results/fixed_global_story/transport_regime/plots/peak_vz_abs.png',
        'analysis/results/fixed_global_story/transport_regime/plots/auc_vz_abs.png',
        'analysis/results/fixed_global_story/transport_regime/plots/active_transport_fraction.png',
        'analysis/results/fixed_global_story/transport_regime/plots/abs_net_z_displacement.png',
        'analysis/results/fixed_global_story/transport_regime/plots/time_to_peak_min.png',
    ]
    for rel in fixed_transport:
        src = ROOT / rel
        copy(src, BUNDLE / 'fixed_vs_rotatable' / 'transport_metrics' / src.name)

    for family in ['total480', 'mpc12']:
        copy(ROOT / f'analysis/results/fixed_global_story/unwrapped_annulus/heatmap_plots/{family}_colorbar_only.png', BUNDLE / 'fixed_vs_rotatable' / 'unwrapped_heatmaps' / f'{family}_colorbar_only.png')
        for regime in ['1to4', '1to8', '1to16']:
            for part in ['low_counts', 'high_counts']:
                src = ROOT / f'analysis/results/fixed_global_story/unwrapped_annulus/heatmap_plots_split/{family}_{regime}_{part}_timeline_heatmap_with_colorbar.png'
                copy(src, BUNDLE / 'fixed_vs_rotatable' / 'unwrapped_heatmaps' / src.name)

    motor_transport = [
        'analysis/results/motor_force_sweep/plots/mean_vz_abs.png',
        'analysis/results/motor_force_sweep/plots/std_vz_abs.png',
        'analysis/results/motor_force_sweep/plots/peak_vz_abs.png',
        'analysis/results/motor_force_sweep/plots/auc_vz_abs.png',
        'analysis/results/motor_force_sweep/plots/active_transport_fraction.png',
        'analysis/results/motor_force_sweep/plots/time_to_peak_min.png',
        'analysis/results/motor_force_sweep/plots/abs_net_z_displacement.png',
        'analysis/results/motor_force_sweep/plots/window_mean_vz_abs.png',
    ]
    for rel in motor_transport:
        src = ROOT / rel
        copy(src, BUNDLE / 'motor_force_sweep' / 'transport_metrics' / src.name)

    copy(ROOT / 'analysis/results/motor_force_sweep/unwrapped_annulus/plots/unwrapped_annulus_colorbar_only.png', BUNDLE / 'motor_force_sweep' / 'unwrapped_heatmaps' / 'motor_force_colorbar_only.png')
    for group in ['c40_m12', 'c40_m50', 'c40_m100']:
        for regime in ['1to4', '1to8', '1to16']:
            src = ROOT / f'analysis/results/motor_force_sweep/unwrapped_annulus/plots/{group}_{regime}_unwrapped_annulus_timeline_heatmap_with_colorbar.png'
            copy(src, BUNDLE / 'motor_force_sweep' / 'unwrapped_heatmaps' / src.name)

    for name in ['mean_vz_abs.png', 'peak_vz_abs.png', 'auc_vz_abs.png', 'active_transport_fraction.png', 'legend_only.png']:
        src = ROOT / 'analysis/results/fixed_global_story/total_motor_sweep/plots' / name
        copy(src, BUNDLE / 'total_motor_sweep' / name)
    copy(ROOT / 'analysis/results/fixed_global_story/total_motor_sweep/selected_conditions.csv', BUNDLE / 'total_motor_sweep' / 'selected_conditions.csv')

    summary = BUNDLE / 'notes' / 'story_summary.md'
    summary.write_text("""# Cytosim Fixed-Global Story Summary

## Main fixed-global vs rotatable comparison
- Fixed-global remains consistently slower than rotatable in the main transport metrics.
- The cleanest headline plot is `vz_abs_mean_with_control.png`.
- `vz_std_presentation.png` supports the same story: fixed-global is not only slower, but less dynamically bursty.
- The unwrapped-annulus heatmaps show the filament-material organization in an opened annulus view using medoid `fiber:speckle` trajectories over the same 800 s post-onset window used in the transport analysis.
- The split `low_counts` / `high_counts` heatmaps are the presentation-friendly versions of the large 11-row comparison panels.

## Motor unbinding-force sweep
- This sweep turned out to be scientifically useful even though it was not the originally intended crosslinker test.
- In code terms, motor `unbinding_force` controls how strongly a loaded motor resists force-dependent detachment.
- Lower motor `unbinding_force` means motors let go more easily under load; higher values keep them engaged longer under load.
- The data show that low-force cases, especially below the default `2.5`, can transport less than the no-motor control.
- Transport increases strongly with motor `unbinding_force`, especially when motors per cluster increase from `12` to `50` and `100`.
- The motor-force unwrapped-annulus heatmaps now use the same continuous filament-painted density style as the main story figures.

## Total motor sweep
- Using the baseline fixed-global `mpc12` series plus the new `40 x 50` and `40 x 100` default-force cases, transport rises strongly with total motor number.
- The compact total-motor figures show a clear crossover from below-control transport at low totals to strong transport at `2000` and `4000` motors.

## In progress
- The corrected crosslinker-force sweep is still in progress.
- Partial xlink-force results exist locally, but the full sweep has not yet been folded into the final comparison figures.
- Once the remaining crosslinker-force runs are transferred, the next key comparison will be motor persistence under load vs network rearrangeability under load.
""", encoding='utf-8')
    print(BUNDLE)


if __name__ == '__main__':
    main()
