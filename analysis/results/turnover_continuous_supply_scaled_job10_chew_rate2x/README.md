# Continuous top-supply scaled turnover analysis

## Dataset
- Runs parsed successfully: `20` / `20`.
- Cases: `small_long` and `large_long`, each with no-motor and rotatable-motor conditions.
- Movement and retention metrics are computed from `1.00 min` onward.
- Top/source-zone mass fraction uses the upper 30% of the annulus.
- Chewer-band mass fraction uses the bottom 30% chewer zone from the campaign manifest.

## Lowest final chewer-band mass fraction
- `large_long` + `nomotor_xlink`: final chewer-band fraction `0.018`, bottom-half fraction `0.166`.
- `small_long` + `nomotor_xlink`: final chewer-band fraction `0.029`, bottom-half fraction `0.195`.
- `large_long` + `rotatable_xlink`: final chewer-band fraction `0.156`, bottom-half fraction `0.366`.
- `small_long` + `rotatable_xlink`: final chewer-band fraction `0.322`, bottom-half fraction `0.481`.

## Lowest mass retention
- `small_long` + `rotatable_xlink`: mass retention `1.424`, AUC(|v_z|) `7.432 um`. final fiber length `200.4 um`, off length `36.9 um`.
- `small_long` + `nomotor_xlink`: mass retention `1.524`, AUC(|v_z|) `3.583 um`. final fiber length `206.0 um`, off length `31.1 um`.
- `large_long` + `rotatable_xlink`: mass retention `1.665`, AUC(|v_z|) `4.766 um`. final fiber length `1986.0 um`, off length `250.5 um`.
- `large_long` + `nomotor_xlink`: mass retention `1.754`, AUC(|v_z|) `1.805 um`. final fiber length `2048.2 um`, off length `188.1 um`.

## Files
- `run_metrics.csv`: one row per replicate.
- `condition_metrics.csv`: mean/SD/SEM across five replicates.
- `condition_timecourses.csv`: mean time traces.
- `fiber_length_frames.csv` and `fiber_length_condition_summary.csv`: exact `report fiber:length` summaries.
- `movement_matrices/`: metric maps.
- `timecourses/by_size/`: condition mean +/- SEM traces.
- `kymographs/full_timeline_kymograph_grid_medoids.png`: representative axial-density trajectories.
- `unwrapped_annulus/`: representative unwrapped heatmaps and line snapshots.
