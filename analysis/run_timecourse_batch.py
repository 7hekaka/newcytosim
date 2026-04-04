#!/usr/bin/env python3
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt

import plot_timecourse_panels as base


def main() -> None:
    ap = argparse.ArgumentParser(description="Run batch-sized time-course comparisons for selected campaign groups.")
    ap.add_argument("--rotatable", type=Path, default=base.DEFAULT_ROTATABLE)
    ap.add_argument("--fixed", type=Path, default=base.DEFAULT_FIXED)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--group", action="append", default=[])
    ap.add_argument("--metric", action="append", default=[])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--nr", type=int, default=8)
    ap.add_argument("--nth", type=int, default=72)
    ap.add_argument("--nz", type=int, default=80)
    args = ap.parse_args()

    group_filter = set(args.group)
    metric_filter = set(args.metric)

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.linewidth": 1.1,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )

    tasks = []
    for row in base.load_metadata(args.rotatable, "rotatable"):
        if group_filter and row["group"] not in group_filter:
            continue
        tasks.append((row, args.nr, args.nth, args.nz))
    for row in base.load_metadata(args.fixed, "fixed_global"):
        if group_filter and row["group"] not in group_filter:
            continue
        tasks.append((row, args.nr, args.nth, args.nz))

    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            run_rows = list(ex.map(base.analyze_run, tasks))
    else:
        run_rows = [base.analyze_run(task) for task in tasks]

    args.outdir.mkdir(parents=True, exist_ok=True)
    base.write_csv(args.outdir / "timecourse_run_status.csv", run_rows)

    condition_rows = base.aggregate_runs(run_rows)
    base.write_csv(args.outdir / "condition_timecourses.csv", condition_rows)

    lookup = base.build_lookup(condition_rows)
    family_specs = [f for f in base.FAMILY_SPECS if not group_filter or f["group"] in group_filter]
    metrics = [m for m in base.METRICS if not metric_filter or m["field"] in metric_filter or m["slug"] in metric_filter]

    outfiles = []
    for family_spec in family_specs:
        for metric in metrics:
            outfiles.extend(base.make_timecourse_figure(lookup, family_spec, metric, args.outdir))

    print(args.outdir / "timecourse_run_status.csv")
    print(args.outdir / "condition_timecourses.csv")
    for path in outfiles:
        print(path)


if __name__ == "__main__":
    main()
