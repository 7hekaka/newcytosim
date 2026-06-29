# Long Treadmilling-Only Nuclear-Aster Pilot

Purpose: test whether a dense, short, nucleus-bound actin aster can bias nuclear motion when turnover is driven only by treadmilling.

Major changes relative to `ast/nuclear_aster_one_sided_turnover_pilot/`:

- No chewers.
- No crosslinkers.
- Actin turnover only through `activity = treadmill`.
- Larger system: `rectangle` length `32, 16, 2`.
- Longer runs: 5 min, `60000` steps at `dt = 0.005 s`.
- Dense aster: `120` nucleus-bound filaments.
- Filaments start at `0.1 um`.
- Filaments are capped at `max_length = 1.0 um`.
- Polymer pool is non-limiting: `total_polymer = 12000`.
- Treadmilling rates are net-growing but turnover-active: plus-end growth `0.75 um/s`, minus-end shrinkage `-0.35 um/s`. This lets filaments grow from `0.1 um` toward the `1.0 um` cap while still losing material at the minus end.

Production cases:

- `01_symmetric_dense_treadmill`: symmetric dense aster.
- `02_minus_x_cap_dense_treadmill`: one-sided minus-x angular cap.
- `03_plus_x_cap_dense_treadmill`: one-sided plus-x angular cap, direction-flip control.
- `04_no_nucleator_control`: passive body baseline.

Generate configs:

```bash
python3 ast/nuclear_aster_treadmilling_long_pilot/generate_treadmilling_long_pilot.py
```

Smoke test:

```bash
bash ast/nuclear_aster_treadmilling_long_pilot/run_all.sh smoke
```

Run all locally:

```bash
bash ast/nuclear_aster_treadmilling_long_pilot/run_all.sh
```

Submit to cluster:

```bash
cd ast/nuclear_aster_treadmilling_long_pilot
DRY_RUN=1 ./submit_cluster.sh
./submit_cluster.sh
```

Analyze returned/local outputs:

```bash
python3 ast/nuclear_aster_treadmilling_long_pilot/analyze_treadmilling_long.py
```

Analyze the returned cluster outputs in `job03/save`:

```bash
python3 ast/nuclear_aster_treadmilling_long_pilot/analyze_treadmilling_long.py --input-root job03/save --output-dir analysis/job03
```

Analysis outputs:

- `analysis/run_metrics.csv`
- `analysis/condition_summary.csv`
- `analysis/timecourse_run_metrics.csv`
- `analysis/timecourse_condition_summary.csv`
- `analysis/dx_timecourse.png`

Returned `job03` result:

- `analysis/job03/run_metrics.csv`
- `analysis/job03/condition_summary.csv`
- `analysis/job03/timecourse_run_metrics.csv`
- `analysis/job03/timecourse_condition_summary.csv`
- `analysis/job03/dx_timecourse.png`
- `analysis/job03/representative_runs.txt`
