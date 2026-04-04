#!/usr/bin/env python3
"""
Multi-panel plotting helper.

Example:
python -m analysis.panel --root ./experiments/cyl --metric radial_flux --dt 1.0 --nfil 1024
"""
import argparse, pathlib, re
import numpy as np
import matplotlib.pyplot as plt

from analysis.compare import discover_runs
from .analyze_run import analyze_file

RUN_RE = re.compile(r"^r\d{4}$", re.IGNORECASE)

def _cli():
    ap = argparse.ArgumentParser(description="Multi-panel plot of multiple conditions side-by-side.")
    ap.add_argument("--root", required=True,
                    help="Directory with multiple condition subfolders (each containing r???? runs)")
    ap.add_argument("--metric", default="radial_flux",
                    choices=["top_bias", "azimuthal_throughput", "radial_flux", "nematic_order",
                             "annulus_envelope", "voxelized_union"],
                    help="Metric to plot")
    ap.add_argument("--dt", type=float, required=True, help="Frame interval (seconds)")
    ap.add_argument("--nfil", type=int, required=True, help="Number of filaments per run")
    ap.add_argument("--drop-first", action="store_true", help="Exclude frame 0 from each run")
    ap.add_argument("--minutes", action="store_true", help="Use minutes on x-axis")
    ap.add_argument("--out", default="panel.png", help="Output PNG filename")
    args = ap.parse_args()

    root = pathlib.Path(args.root).resolve()
    subconds = [d for d in sorted(root.iterdir()) if d.is_dir() and not RUN_RE.match(d.name)]
    if not subconds:
        ap.error(f"No condition subfolders found under {root}")

    n = len(subconds)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 3*nrows), sharey=True, constrained_layout=True)
    axes = np.atleast_2d(axes)

    for idx, cond in enumerate(subconds):
        ax = axes[idx // ncols, idx % ncols]
        ax.set_title(cond.name)

        runs = discover_runs(cond)
        for pt in runs:
            res = analyze_file(pt, nfil=args.nfil, dt=args.dt)
            t = res["t"][1:] if args.drop_first else res["t"]
            y = res[args.metric][1:] if args.drop_first else res[args.metric]
            x = t / 60.0 if args.minutes else t
            ax.plot(x, y, lw=1)

        ax.grid(True, alpha=0.3)
        if idx % ncols == 0:
            ax.set_ylabel(args.metric.replace("_", " ").title())

    for j in range(n, nrows*ncols):
        fig.delaxes(axes[j // ncols, j % ncols])

    axes[-1, min(ncols-1, n-1)].set_xlabel("Time (min)" if args.minutes else "Time (s)")
    fig.suptitle(f"Panel: {args.metric.replace('_',' ').title()}")
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"[saved] {args.out}")

if __name__ == "__main__":
    _cli()
