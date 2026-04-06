# Cytosim Fixed-Global Story Summary

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
