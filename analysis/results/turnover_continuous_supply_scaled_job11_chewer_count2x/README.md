# Continuous top-supply scaled turnover analysis

## Dataset
- Runs parsed successfully: `20` / `20`.
- Cases: `small_long` and `large_long`, each with no-motor and rotatable-motor conditions.
- Movement and retention metrics are computed from `1.00 min` onward.
- Top/source-zone mass fraction uses the upper 30% of the annulus.
- Chewer-band mass fraction uses the bottom 30% chewer zone from the campaign manifest.

## Lowest final chewer-band mass fraction
- `small_long` + `nomotor_xlink`: final chewer-band fraction `0.018`, bottom-half fraction `0.216`.
- `large_long` + `nomotor_xlink`: final chewer-band fraction `0.022`, bottom-half fraction `0.180`.
- `large_long` + `rotatable_xlink`: final chewer-band fraction `0.124`, bottom-half fraction `0.327`.
- `small_long` + `rotatable_xlink`: final chewer-band fraction `0.297`, bottom-half fraction `0.495`.

## Lowest mass retention
- `small_long` + `rotatable_xlink`: mass retention `1.425`, AUC(|v_z|) `9.138 um`. final fiber length `203.7 um`, off length `32.3 um`.
- `small_long` + `nomotor_xlink`: mass retention `1.429`, AUC(|v_z|) `3.353 um`. final fiber length `207.7 um`, off length `29.3 um`.
- `large_long` + `rotatable_xlink`: mass retention `1.708`, AUC(|v_z|) `4.617 um`. final fiber length `1999.1 um`, off length `236.0 um`.
- `large_long` + `nomotor_xlink`: mass retention `1.715`, AUC(|v_z|) `1.893 um`. final fiber length `2045.7 um`, off length `190.9 um`.

## Files
- `run_metrics.csv`: one row per replicate.
- `condition_metrics.csv`: mean/SD/SEM across five replicates.
- `condition_timecourses.csv`: mean time traces.
- `fiber_length_frames.csv` and `fiber_length_condition_summary.csv`: exact `report fiber:length` summaries.
- `movement_matrices/`: metric maps.
- `timecourses/by_size/`: condition mean +/- SEM traces.
- `kymographs/full_timeline_kymograph_grid_medoids.png`: representative axial-density trajectories.
- `unwrapped_annulus/`: representative unwrapped heatmaps and line snapshots.
