# Nuclear Propulsion Mechanism Mockups

Purpose: fast Cytosim mechanism screen for directed nuclear motion in the endosperm project.

The prior treadmilling-only pilot produced only weak signed displacement. These mockups add the missing physical ingredients one at a time:

- Geometry-only nucleus-bound asters.
- One-sided fixed non-motor binders, to test external anchoring/friction.
- One-sided fixed motors, to test actin-flow or cortical-traction-like force generation.
- Lab-frame growing actin filaments aimed at the nucleus, to test polymerization push with and without anchoring.

## Conditions

The production set contains 10 conditions x 3 replicates = 30 trajectories.

- `01_symmetric_aster_no_external`: symmetric nucleus-bound aster, no external traction.
- `02_plus_x_cap_no_external`: plus-x one-sided aster, geometry-only control.
- `03_plus_x_cap_fixed_binders`: one-sided aster plus fixed non-motor binders on the plus-x side.
- `04_plus_x_cap_fixed_plus_motors`: one-sided aster plus fixed plus-end motors on the plus-x side.
- `05_symmetric_aster_fixed_plus_motors`: symmetric aster plus one-sided fixed plus-end motors.
- `06_plus_x_cap_fixed_minus_motors`: one-sided aster plus fixed minus-end motors on the plus-x side.
- `07_symmetric_aster_fixed_minus_motors`: symmetric aster plus one-sided fixed minus-end motors.
- `08_push_filaments_unanchored`: lab-frame filaments grow toward the nucleus, no anchors.
- `09_push_filaments_fixed_binders`: lab-frame growing filaments plus fixed non-motor anchors.
- `10_push_filaments_fixed_plus_motors`: lab-frame growing filaments plus fixed motors.

## Generate

```bash
python3 ast/nuclear_propulsion_mechanism_mockups/generate_mechanism_mockups.py
```

## Smoke Test

```bash
bash ast/nuclear_propulsion_mechanism_mockups/run_all.sh smoke
```

## Run Locally

Run a 20 s quick-look set:

```bash
bash ast/nuclear_propulsion_mechanism_mockups/run_all.sh quick
```

Run one replicate per condition:

```bash
bash ast/nuclear_propulsion_mechanism_mockups/run_all.sh first-reps
```

Run all production trajectories locally:

```bash
bash ast/nuclear_propulsion_mechanism_mockups/run_all.sh all
```

## Submit To Cluster

```bash
cd ast/nuclear_propulsion_mechanism_mockups
DRY_RUN=1 ./submit_cluster.sh
./submit_cluster.sh
```

The submit helper expects `build_cluster/bin/sim` unless `SIM_EXE` is set.

## Analyze

For local production outputs:

```bash
python3 ast/nuclear_propulsion_mechanism_mockups/analyze_mechanism_mockups.py
```

For returned cluster outputs in `jobXX/save`:

```bash
python3 ast/nuclear_propulsion_mechanism_mockups/analyze_mechanism_mockups.py --input-root jobXX/save --output-dir analysis/jobXX
```

Main outputs:

- `analysis/run_metrics.csv`
- `analysis/condition_summary.csv`
- `analysis/timecourse_run_metrics.csv`
- `analysis/timecourse_condition_summary.csv`
- `analysis/mechanism_final_displacement.png`
- `analysis/mechanism_dx_timecourse.png`

## Interpretation Logic

- If only motor cases move strongly, the relevant ingredient is external motor traction or actin flow, not pure treadmilling.
- If plus- and minus-directed motors move in opposite directions, the sign is set by motor polarity and filament polarity.
- If `09` moves but `08` does not, polymerization push requires anchoring.
- If `03` moves nearly as much as or more than motor cases, fixed anchoring/friction is doing more than motor stepping.
- If `02` moves similarly to all active conditions, the screen is dominated by asymmetric geometry and needs a stricter force-transmission setup.

## Current Quick-Look Result

The 20 s one-replicate quick-look output is in `analysis/quick_20s/`.

Main observations:

- Geometry-only is weak: `02_plus_x_cap_no_external` gives only small +x motion.
- Polymerization-push mockups are near zero over 20 s, even with fixed binders/motors near the seeded filaments.
- Fixed non-motor binders can move the nucleus strongly, so external anchoring/friction alone is a major force-transmission route.
- Fixed plus-end motors and minus-end motors move the nucleus in opposite directions, consistent with gliding-assay polarity physics.
- The strongest +x motion is produced by minus-end-directed fixed motors on a plus-x actin cap, followed by the same motor field acting on a symmetric aster.

Working interpretation:

The useful mechanism is likely not treadmilling-only propulsion. It is more likely external actin anchoring/motor traction plus filament polarity, with turnover acting to remodel or localize the actin structure.
