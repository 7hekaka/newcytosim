#!/usr/bin/env python3
import argparse, pathlib, re, math
from typing import List, Dict, Tuple
import numpy as np
import matplotlib.pyplot as plt

from analysis.analyze_run import analyze_file

RUN_RE = re.compile(r"^r\d{4}$", re.IGNORECASE)

METRICS = [
    ("top_bias", "Top Bias Index"),
    ("top_fraction_band", "Top Fraction (banded)"),       # NEW
    ("azimuthal_throughput", "Azimuthal Throughput (rad/s)"),
    ("radial_flux", "Radial Flux"),
    ("nematic_order", "Nematic Order (S₂)"),
    ("annulus_envelope", "Annulus Envelope (avg z-span)"),
    ("voxelized_union", "Voxelized Occupancy"),
    ("voxelized_union_per_fil", "Voxelized Occupancy / Filament"),  # NEW
]


# Accept legacy/alternate metric keys
METRIC_ALIASES = {
    "top_bias": [
        "top_bias", "topbias", "bias_top"
        ],
    "top_fraction_band": [
        "top_fraction_band", "top_fraction"
        ],
    "azimuthal_throughput": [
        "azimuthal_throughput", "azimuthal_throughput_net", "azimuthal_net",
        "theta_flux_net", "azimuthal"  # last resort names if you used them
    ],
    "radial_flux": [
        "radial_flux", "radial_flux_net", "radial_net", "radial"
    ],
    "nematic_order": [
        "nematic_order", "nematic_xy", "S2", "nematic"
    ],
    "annulus_envelope": [
        "annulus_envelope", "annulus_envelope_zspan", "annulus_zspan"
    ],
    "voxelized_union": [
        "voxelized_union", "voxel_union", "voxelized_occupancy"
    ],
    "voxelized_union_per_fil": [
        "voxelized_union_per_fil", "voxel_union_per_fil", "voxelized_occupancy_per_fil"
    ],

}

def get_series(res: dict, key: str):
    """Return (array, found_key) for metric 'key' using aliases. None if missing."""
    for k in METRIC_ALIASES.get(key, [key]):
        if k in res:
            return np.asarray(res[k], dtype=float), k
    return None, None


def _normalize_result(out, dt: float):
    """
    Accepts analyze_file(...) outputs in several legacy formats and returns a dict with 't'.
    Supported:
      - dict with 't'
      - dict with 'time'
      - (t, dict)
    """
    if isinstance(out, tuple) and len(out) == 2:
        t, d = out
        if not isinstance(d, dict):
            raise TypeError("analyze_file returned (t, <not dict>)")
        d = dict(d)
        d["t"] = np.asarray(t, float)
        return d

    if isinstance(out, dict):
        d = dict(out)
        if "t" in d:
            d["t"] = np.asarray(d["t"], float)
            return d
        if "time" in d:
            d["t"] = np.asarray(d["time"], float)
            return d
        # Last-resort fallback: build t from any metric length if present
        for k, v in d.items():
            if hasattr(v, "__len__"):
                d["t"] = np.arange(len(v), dtype=float) * float(dt)
                return d
        raise KeyError("No 't' or 'time' in analyze_file output")
    raise TypeError("analyze_file returned unsupported type")


def discover_runs(root: pathlib.Path) -> List[pathlib.Path]:
    """Find all r????/points.txt under a condition root."""
    runs = []
    for child in sorted(root.iterdir()):
        if child.is_dir() and RUN_RE.match(child.name):
            pt = child / "points.txt"
            if pt.exists():
                runs.append(pt)
    if not runs:
        raise FileNotFoundError(f"No r????/points.txt found under {root}")
    return runs

def interpolate_to_grid(t_src: np.ndarray, y_src: np.ndarray, t_grid: np.ndarray) -> np.ndarray:
    """Linear interpolate onto t_grid; returns NaNs outside range."""
    if len(t_src) < 2:
        return np.full_like(t_grid, np.nan, dtype=float)
    y = np.full_like(t_grid, np.nan, dtype=float)
    mask = (t_grid >= t_src[0]) & (t_grid <= t_src[-1])
    if mask.any():
        y[mask] = np.interp(t_grid[mask], t_src, y_src)
    return y

def aggregate_series(series_list: List[np.ndarray]) -> Dict[str, np.ndarray]:
    """Stack (n_runs, n_time) with NaNs allowed, return mean and 95% CI."""
    M = np.vstack(series_list) if series_list else np.zeros((0, len(series_list[0])))
    mean = np.nanmean(M, axis=0)
    n = np.sum(~np.isnan(M), axis=0).astype(float)
    std = np.nanstd(M, axis=0, ddof=1, where=~np.isnan(M))
    se = np.divide(std, np.sqrt(np.maximum(n, 1.0)), out=np.zeros_like(std), where=n > 0)
    ci = 1.96 * se
    lo = mean - ci
    hi = mean + ci
    return dict(mean=mean, lo=lo, hi=hi, n=n)

def load_condition(root: pathlib.Path, nfil: int, dt: float, drop_first: bool) -> Tuple[np.ndarray, Dict[str, List[np.ndarray]]]:
    """Analyze all runs under a condition root. Returns t_grid and metric stacks."""
    run_files = discover_runs(root)
    results = []
    for pt in run_files:
        res = _normalize_result(analyze_file(pt, nfil=nfil, dt=dt), dt=dt)
        t = np.asarray(res["t"], dtype=float)
        if drop_first and len(t) >= 2:
            t = t[1:] - t[1]
            res = {k: (np.asarray(v)[1:] if hasattr(v, "__len__") and len(v) == len(t)+1 else v)
                   for k, v in res.items()}
            res["t"] = t
        results.append(res)

    t0s = [r["t"][0] for r in results if len(r["t"]) > 0]
    tEnds = [r["t"][-1] for r in results if len(r["t"]) > 0]
    if not t0s:
        raise RuntimeError(f"No time data found under {root}")

    t0 = max(min(t0s), 0.0)
    t1 = min(tEnds)
    if t1 <= t0:
        t0 = min(t0s); t1 = max(tEnds)

    nsteps = int(math.floor((t1 - t0) / dt)) + 1
    t_grid = t0 + np.arange(nsteps) * dt

    stacked: Dict[str, List[np.ndarray]] = {k: [] for k, _ in METRICS}
    for r in results:
        for k, _label in METRICS:
            y, found = get_series(r, k)
            t_run = np.asarray(r["t"], dtype=float)

            if y is None:
                # Metric missing for this run; keep NaNs so aggregation still works.
                stacked[k].append(np.full_like(t_grid, np.nan, dtype=float))
                continue

            # length-guard in case of off-by-one
            m = min(len(y), len(t_run))
            y = y[:m]
            t_run = t_run[:m]

            stacked[k].append(interpolate_to_grid(t_run, y, t_grid))

    return t_grid, stacked

def plot_overlay(t: np.ndarray, aggregates: Dict[str, Dict[str, np.ndarray]], label_map: Dict[str, str],
                 minutes: bool, out_png: pathlib.Path):
    x = t / 60.0 if minutes else t
    xlab = "Time (min)" if minutes else "Time (s)"

    nrows = len(METRICS)
    fig, axes = plt.subplots(nrows, 1, figsize=(9, 1.9*nrows), sharex=True, constrained_layout=True)
    if nrows == 1:
        axes = [axes]

    for ax, (key, title) in zip(axes, METRICS):
        for label, agg in aggregates.items():
            ax.plot(x, agg[key]["mean"], lw=2, label=f"{label_map[label]} (mean)")
            ax.fill_between(x, agg[key]["lo"], agg[key]["hi"], alpha=0.2)
        ax.set_ylabel(title)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel(xlab)
    axes[0].legend(ncol=2, fontsize=9)

    fig.suptitle("Condition Comparison — transport metrics", y=1.02, fontsize=12)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    print(f"[saved] {out_png}")

def _cli():
    ap = argparse.ArgumentParser(description="Compare conditions across r???? runs.")
    ap.add_argument("--cond", action="append", required=True,
                    help="Label=PATH[,nfil=NNN]  (can be given multiple times)")
    ap.add_argument("--nfil", type=int, default=None,
                    help="Fallback number of filaments if a condition doesn't specify nfil")
    ap.add_argument("--dt", type=float, required=True, help="Frame interval (seconds)")
    ap.add_argument("--minutes", action="store_true", help="Plot time in minutes on x-axis")
    ap.add_argument("--drop-first", action="store_true", help="Exclude frame 0 and re-zero time")
    ap.add_argument("--out", default="compare_metrics.png", help="Output PNG file")
    args = ap.parse_args()

    def parse_label_path(s: str):
        if "=" not in s:
            raise ValueError("Use Label=PATH[,nfil=NNN] for --cond")
        label, rest = s.split("=", 1)
        label = label.strip()
        nfil_opt = None
        parts = [p.strip() for p in rest.split(",") if p.strip()]
        path_str = parts[0]
        for opt in parts[1:]:
            if opt.lower().startswith("nfil="):
                nfil_opt = int(opt.split("=", 1)[1])
            else:
                raise ValueError(f"Unknown option '{opt}' (supported: nfil=)")
        return label, pathlib.Path(path_str).resolve(), nfil_opt

    parsed = [parse_label_path(s) for s in args.cond]
    label_map = {lab: lab for (lab, _, _) in parsed}

    # Load each condition
    t_grids = []
    stacks = {}
    for lab, path, nfil_opt in parsed:
        nfil_use = nfil_opt if nfil_opt is not None else args.nfil
        if nfil_use is None:
            raise SystemExit(f"{lab} missing nfil and no global --nfil provided.")
        print(f"[scan] {lab}: {path}  (nfil={nfil_use})")
        t, stack = load_condition(path, nfil=nfil_use, dt=args.dt, drop_first=args.drop_first)
        t_grids.append(t)
        stacks[lab] = stack

    # Build common grid across all conditions (overlap)
    t0 = max(t[0] for t in t_grids)
    t1 = min(t[-1] for t in t_grids)
    dt = args.dt
    nsteps = int(max(1, (t1 - t0) // dt) + 1)
    t_grid = t0 + np.arange(nsteps) * dt

    def agg_from_stack(t_src, stack):
        aggs = {}
        for key, _ in METRICS:
            series_interp = [interpolate_to_grid(t_src, y, t_grid) for y in stack[key]]
            aggs[key] = aggregate_series(series_interp)
        return aggs

    aggregates = {lab: agg_from_stack(t_grids[i], stacks[lab]) for i, (lab, *_ ) in enumerate(parsed)}

    plot_overlay(t_grid, aggregates, label_map, args.minutes, pathlib.Path(args.out))

if __name__ == "__main__":
    _cli()
