# Large Continuous-Supply Turnover Presentation Plots

Inputs:
- Baseline: `analysis/results/turnover_continuous_supply_scaled_job09`
- Ultra aggressive chewing: `analysis/results/turnover_continuous_supply_scaled_job12_chewer_count2x_rate2x`

Included conditions:
- `large_long/nomotor_xlink`
- `large_long/rotatable_xlink`

Excluded:
- Small-system tests are intentionally omitted.
- `job13_sever_chew` is not plotted because all 20 runs were partial/incomplete in `run_status.csv`; representative run `r0001` stopped with Cytosim status 11 and `Segmentation fault`.

Generated figures:
- `large_story_summary.svg/.png`: compact comparison for the main slide.
- `large_velocity_endpoint_metrics.svg/.png`: endpoint velocity and directionality metrics.
- `large_velocity_timecourses.svg/.png`: speed and cumulative speed timecourses.
- `large_turnover_endpoint_metrics.svg/.png`: final bottom-clearing and mass readouts.
- `large_turnover_timecourses.svg/.png`: bottom-region mass timecourses.

Key large-system readout:
- Baseline motors increase AUC |v_z| by 2.91x over no motors.
- Ultra-aggressive chewing motors increase AUC |v_z| by 2.04x over no motors.
- Ultra-aggressive chewing drops the motor-condition final chewer-band mass fraction by 0.139 versus baseline.
- Motor AUC remains comparable after aggressive chewing: baseline 4.83 um, ultra 4.18 um.
