# One-Sided Nuclear-Aster Turnover Pilot

Purpose: quick Cytosim test for whether a short, nucleus-proximal actin aster with one-sided turnover can bias nuclear motion.

This is a follow-up to `ast/nuclear_aster_turnover_propulsion_pilot/`. The previous "polar" cases used `radius = 0.5, 0.1` without `type=angular`, which still distributed aster fibers across the full nuclear surface. This pilot explicitly uses `type=angular` for one-sided cap asters.

Key assumptions:

- Actin filaments are short-capped: `max_length = 2.5 um`, initial length `2.0 um`.
- Fast actin dynamics use `growing_speed = 0.75 um/s` and minus-end shrinkage `-0.35 um/s`.
- One-sided turnover factors are fixed `single` objects with `activity = fixed`, placed in a slab on the minus-x side of the starting nucleus.
- The fixed turnover slab is static in the lab frame. It is a first-pass proxy, not a fully nucleus-following ACTIN10 cap.
- The one-sided aster uses `type=angular`, `direction = -1 0 0`, and `aster_angle = 0.85 rad`.
- Delayed-turnover cases introduce the fixed one-sided turnover factors at `t = 4 s`, after the 2.0 um initial filaments have reached the 2.5 um length cap.

Cases:

- `01_symmetric_short_no_chewer`: symmetric short aster, no localized turnover.
- `02_symmetric_short_one_sided_chewer`: symmetric short aster, one-sided fixed turnover factors.
- `03_cap_short_no_chewer`: one-sided angular-cap aster, no localized turnover.
- `04_cap_short_one_sided_chewer`: one-sided angular-cap aster, same-side fixed turnover factors.
- `05_no_nucleator_fixed_chewer_control`: no actin nucleator, fixed turnover factors only.
- `06_symmetric_short_delayed_one_sided_chewer`: symmetric short aster, one-sided fixed turnover factors introduced at 4 s.
- `07_cap_short_delayed_one_sided_chewer`: one-sided angular-cap aster, same-side fixed turnover factors introduced at 4 s.

Run locally:

```bash
python3 ast/nuclear_aster_one_sided_turnover_pilot/generate_one_sided_turnover_pilot.py
bash ast/nuclear_aster_one_sided_turnover_pilot/run_all.sh
python3 ast/nuclear_aster_one_sided_turnover_pilot/analyze_one_sided_turnover.py
```

Outputs:

- `analysis/run_metrics.csv`
- `analysis/condition_summary.csv`
- `analysis/timecourse_run_metrics.csv`
- `analysis/timecourse_condition_summary.csv`
- `analysis/dx_timecourse.png`
