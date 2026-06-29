# Mini minus-end-chew turnover analysis

## Dataset
- Runs parsed successfully: `45` / `45`.
- Replicates: `5` per scenario/condition.
- Growth-only phase ends at `2.0 min`; crosslink conditioning ends and motor readout starts at `3.0 min`.
- Movement metrics are computed only from the post-motor window.
- Matrix plots use replicate means. Kymographs and unwrapped line snapshots use the replicate whose AUC is closest to that condition mean.

## Top conditions by mean AUC(|v_z|)
- `top_two_thirds_aligned` + `rotatable_xlink`: AUC `9.386 +/- 0.516 um SEM, mean |v_z| `0.01563 um/s`, medoid `r0030`.
- `top_two_thirds_mixed_polarity` + `rotatable_xlink`: AUC `9.115 +/- 0.400 um SEM, mean |v_z| `0.01518 um/s`, medoid `r0042`.
- `full_height_aligned` + `rotatable_xlink`: AUC `9.008 +/- 0.483 um SEM, mean |v_z| `0.01501 um/s`, medoid `r0014`.
- `full_height_aligned` + `fixed_global_xlink`: AUC `5.690 +/- 0.277 um SEM, mean |v_z| `0.00948 um/s`, medoid `r0005`.
- `top_two_thirds_aligned` + `fixed_global_xlink`: AUC `5.205 +/- 0.372 um SEM, mean |v_z| `0.00868 um/s`, medoid `r0018`.
- `top_two_thirds_mixed_polarity` + `fixed_global_xlink`: AUC `4.301 +/- 0.107 um SEM, mean |v_z| `0.00717 um/s`, medoid `r0034`.

## Top conditions by excess AUC above matched no-motor control
- `top_two_thirds_mixed_polarity` + `rotatable_xlink`: excess AUC `4.229 +/- 0.313 um SEM, matched-active fraction `0.453`.
- `top_two_thirds_aligned` + `rotatable_xlink`: excess AUC `4.029 +/- 0.408 um SEM, matched-active fraction `0.376`.
- `full_height_aligned` + `rotatable_xlink`: excess AUC `3.633 +/- 0.296 um SEM, matched-active fraction `0.361`.
- `full_height_aligned` + `fixed_global_xlink`: excess AUC `1.453 +/- 0.262 um SEM, matched-active fraction `0.186`.
- `top_two_thirds_aligned` + `fixed_global_xlink`: excess AUC `1.141 +/- 0.294 um SEM, matched-active fraction `0.167`.
- `top_two_thirds_mixed_polarity` + `fixed_global_xlink`: excess AUC `1.034 +/- 0.118 um SEM, matched-active fraction `0.151`.

## Lowest final chewer-band mass fraction
- `top_two_thirds_aligned` + `nomotor_xlink`: final chewer-band fraction `0.048`, bottom-half fraction `0.243`.
- `full_height_aligned` + `nomotor_xlink`: final chewer-band fraction `0.105`, bottom-half fraction `0.301`.
- `top_two_thirds_aligned` + `fixed_global_xlink`: final chewer-band fraction `0.139`, bottom-half fraction `0.364`.
- `top_two_thirds_mixed_polarity` + `rotatable_xlink`: final chewer-band fraction `0.147`, bottom-half fraction `0.367`.
- `full_height_aligned` + `fixed_global_xlink`: final chewer-band fraction `0.193`, bottom-half fraction `0.417`.
- `top_two_thirds_mixed_polarity` + `nomotor_xlink`: final chewer-band fraction `0.266`, bottom-half fraction `0.447`.

## Exact fiber:length at final frame
- `top_two_thirds_aligned` + `rotatable_xlink`: total length `1011.8 um`, off/chewed length `148.4 um`.
- `top_two_thirds_mixed_polarity` + `rotatable_xlink`: total length `1020.0 um`, off/chewed length `41.0 um`.
- `full_height_aligned` + `rotatable_xlink`: total length `1012.2 um`, off/chewed length `159.2 um`.
- `full_height_aligned` + `fixed_global_xlink`: total length `1011.8 um`, off/chewed length `161.5 um`.
- `top_two_thirds_aligned` + `fixed_global_xlink`: total length `1011.5 um`, off/chewed length `145.5 um`.
- `top_two_thirds_mixed_polarity` + `fixed_global_xlink`: total length `1020.0 um`, off/chewed length `43.9 um`.

## Files
- `run_metrics.csv`: one row per replicate.
- `condition_metrics.csv`: mean/SD/SEM across five replicates.
- `condition_timecourses.csv`: replicate-mean time traces.
- `fiber_length_frames.csv` and `fiber_length_condition_summary.csv`: exact `report fiber:length` summaries.
- `movement_matrices/`: condition-level metric maps.
- `timecourses/by_scenario/`: mean +/- SEM traces.
- `kymographs/full_timeline_kymograph_grid_medoids.png`: representative axial-density trajectories.
- `unwrapped_annulus/lines/`: representative line snapshots.
