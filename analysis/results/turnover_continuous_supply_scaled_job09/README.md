# Continuous top-supply scaled turnover analysis

## Dataset
- Runs parsed successfully: `20` / `20`.
- Cases: `small_long` and `large_long`, each with no-motor and rotatable-motor conditions.
- Movement and retention metrics are computed from `1.00 min` onward.
- Top/source-zone mass fraction uses the upper 30% of the annulus.
- Chewer-band mass fraction uses the bottom 30% chewer zone from the campaign manifest.

## Lowest final chewer-band mass fraction
- `large_long` + `nomotor_xlink`: final chewer-band fraction `0.055`, bottom-half fraction `0.192`.
- `small_long` + `nomotor_xlink`: final chewer-band fraction `0.080`, bottom-half fraction `0.219`.
- `large_long` + `rotatable_xlink`: final chewer-band fraction `0.175`, bottom-half fraction `0.348`.
- `small_long` + `rotatable_xlink`: final chewer-band fraction `0.461`, bottom-half fraction `0.630`.

## Lowest mass retention
- `small_long` + `rotatable_xlink`: mass retention `1.558`, AUC(|v_z|) `9.614 um`. final fiber length `215.2 um`, off length `17.2 um`.
- `small_long` + `nomotor_xlink`: mass retention `1.625`, AUC(|v_z|) `3.274 um`. final fiber length `218.2 um`, off length `15.0 um`.
- `large_long` + `rotatable_xlink`: mass retention `1.752`, AUC(|v_z|) `4.826 um`. final fiber length `2093.6 um`, off length `125.4 um`.
- `large_long` + `nomotor_xlink`: mass retention `1.793`, AUC(|v_z|) `1.657 um`. final fiber length `2109.9 um`, off length `109.4 um`.

## Files
- `run_metrics.csv`: one row per replicate.
- `condition_metrics.csv`: mean/SD/SEM across five replicates.
- `condition_timecourses.csv`: mean time traces.
- `fiber_length_frames.csv` and `fiber_length_condition_summary.csv`: exact `report fiber:length` summaries.
- `movement_matrices/`: metric maps.
- `timecourses/by_size/`: condition mean +/- SEM traces.
- `kymographs/full_timeline_kymograph_grid_medoids.png`: representative axial-density trajectories.
- `unwrapped_annulus/`: representative unwrapped heatmaps and line snapshots.
