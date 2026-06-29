# Continuous top-supply scaled turnover analysis

## Dataset
- Runs parsed successfully: `20` / `20`.
- Cases: `small_long` and `large_long`, each with no-motor and rotatable-motor conditions.
- Movement and retention metrics are computed from `1.00 min` onward.
- Top/source-zone mass fraction uses the upper 30% of the annulus.
- Chewer-band mass fraction uses the bottom 30% chewer zone from the campaign manifest.

## Lowest final chewer-band mass fraction
- `large_long` + `nomotor_xlink`: final chewer-band fraction `0.000`, bottom-half fraction `0.153`.
- `small_long` + `nomotor_xlink`: final chewer-band fraction `0.004`, bottom-half fraction `0.152`.
- `large_long` + `rotatable_xlink`: final chewer-band fraction `0.036`, bottom-half fraction `0.281`.
- `small_long` + `rotatable_xlink`: final chewer-band fraction `0.096`, bottom-half fraction `0.287`.

## Lowest mass retention
- `small_long` + `rotatable_xlink`: mass retention `1.503`, AUC(|v_z|) `6.149 um`. final fiber length `203.1 um`, off length `44.2 um`.
- `small_long` + `nomotor_xlink`: mass retention `1.540`, AUC(|v_z|) `3.704 um`. final fiber length `210.4 um`, off length `33.3 um`.
- `large_long` + `rotatable_xlink`: mass retention `1.622`, AUC(|v_z|) `4.183 um`. final fiber length `1882.9 um`, off length `381.1 um`.
- `large_long` + `nomotor_xlink`: mass retention `1.731`, AUC(|v_z|) `2.055 um`. final fiber length `2030.7 um`, off length `236.4 um`.

## Files
- `run_metrics.csv`: one row per replicate.
- `condition_metrics.csv`: mean/SD/SEM across five replicates.
- `condition_timecourses.csv`: mean time traces.
- `fiber_length_frames.csv` and `fiber_length_condition_summary.csv`: exact `report fiber:length` summaries.
- `movement_matrices/`: metric maps.
- `timecourses/by_size/`: condition mean +/- SEM traces.
- `kymographs/full_timeline_kymograph_grid_medoids.png`: representative axial-density trajectories.
- `unwrapped_annulus/`: representative unwrapped heatmaps and line snapshots.
