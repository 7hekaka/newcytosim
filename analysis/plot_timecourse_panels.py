#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.metrics import compute_top_bottom_bands, nematic_order_xy
from trust.analyze_axial import central_diff, parse_frames_point, reduce_frame_point, unwrap_angles


DEFAULT_ROTATABLE = ROOT / "clu" / "xlink_regime_map.csv"
DEFAULT_FIXED = ROOT / "clu_fixed_global" / "xlink_regime_map.csv"
DEFAULT_OUTDIR = ROOT / "comparison_timecourses"

TIME_STEP_RE = re.compile(r"^\s*time_step\s*=\s*([0-9.eE+-]+)", re.M)
RUN_BLOCK_RE = re.compile(r"run\s+(\d+)\s+system\s*\{(.*?)\}", re.S)
NB_FRAMES_RE = re.compile(r"nb_frames\s*=\s*(\d+)")
SPACE_PARAM_RE = re.compile(r"^\s*(inner|outer|top|bottom)\s*=\s*([0-9.eE+-]+)", re.M)

XLINK_ORDER = ["1:4", "1:8", "1:16"]

MODEL_SPECS = [
    {
        "key": "rotatable",
        "label": "Rotatable",
        "linestyle": "--",
    },
    {
        "key": "fixed_global",
        "label": "Fixed global",
        "linestyle": "-",
    },
]

COUNT_COLORS = {
    "c10_m48": "#1b9e77",
    "c20_m24": "#d95f02",
    "c40_m12": "#7570b3",
    "c60_m8": "#e7298a",
    "c80_m6": "#66a61e",
    "c10_m12": "#1b9e77",
    "c20_m12": "#d95f02",
    "c40_m12": "#7570b3",
    "c60_m12": "#e7298a",
    "c80_m12": "#66a61e",
}

FOCUS_COLORS = {
    "parallel_only": "#1b9e77",
    "parallel_antiparallel_50_50": "#d95f02",
    "offrate_0p5": "#7570b3",
    "offrate_1p0": "#e7298a",
}

FAMILY_SPECS = [
    {
        "group": "total480",
        "slug": "total480",
        "order": [
            ("c10_m48", "10 clusters"),
            ("c20_m24", "20 clusters"),
            ("c40_m12", "40 clusters"),
            ("c60_m8", "60 clusters"),
            ("c80_m6", "80 clusters"),
        ],
        "colors": COUNT_COLORS,
    },
    {
        "group": "mpc12",
        "slug": "mpc12",
        "order": [
            ("c10_m12", "10 clusters"),
            ("c20_m12", "20 clusters"),
            ("c40_m12", "40 clusters"),
            ("c60_m12", "60 clusters"),
            ("c80_m12", "80 clusters"),
        ],
        "colors": COUNT_COLORS,
    },
    {
        "group": "focus40x12",
        "slug": "focus40x12",
        "order": [
            ("parallel_only", "parallel"),
            ("parallel_antiparallel_50_50", "50:50"),
            ("offrate_0p5", "off-rate 0.5"),
            ("offrate_1p0", "off-rate 1.0"),
        ],
        "colors": FOCUS_COLORS,
    },
]

METRICS = [
    {"field": "vz_abs", "ylabel": "|v_z| (um/s)", "slug": "vz_abs", "clip_zero": True},
    {"field": "Dz", "ylabel": "Depletion index", "slug": "depletion", "clip_zero": True},
    {"field": "shell_occ", "ylabel": "Occupied shell fraction", "slug": "shell_occupancy", "clip_zero": True},
    {"field": "nematic_xy", "ylabel": "XY nematic order", "slug": "nematic_xy", "clip_zero": True},
    {"field": "top_bias", "ylabel": "Top-bottom bias", "slug": "top_bias", "clip_zero": False},
    {"field": "swirl_rate_abs", "ylabel": "|swirl rate| (rad/s)", "slug": "swirl_rate_abs", "clip_zero": True},
]


def parse_config_info(config_path: Path) -> dict:
    text = config_path.read_text(encoding="utf-8", errors="replace")

    m = TIME_STEP_RE.search(text)
    if not m:
        raise ValueError(f"could not find time_step in {config_path}")
    time_step = float(m.group(1))

    runs = []
    for steps_s, body in RUN_BLOCK_RE.findall(text):
        nb = NB_FRAMES_RE.search(body)
        if not nb:
            continue
        steps = int(steps_s)
        nb_frames = int(nb.group(1))
        duration = steps * time_step
        frame_dt = duration / nb_frames if nb_frames else float("nan")
        runs.append(
            {
                "steps": steps,
                "nb_frames": nb_frames,
                "duration": duration,
                "frame_dt": frame_dt,
            }
        )

    if not runs:
        raise ValueError(f"could not find run blocks in {config_path}")

    space = {}
    for key, value in SPACE_PARAM_RE.findall(text):
        space[key] = float(value)

    return {"time_step": time_step, "runs": runs, "space": space}


def shell_occupancy_fraction(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    *,
    inner: float,
    outer: float,
    bottom: float,
    top: float,
    nr: int,
    nth: int,
    nz: int,
) -> float:
    if x.size == 0:
        return float("nan")

    r = np.hypot(x, y)
    mask = (r >= inner) & (r <= outer) & (z >= bottom) & (z <= top)
    if not np.any(mask):
        return 0.0

    r = r[mask]
    theta = np.arctan2(y[mask], x[mask])
    z = z[mask]

    rbins = np.linspace(inner, outer + 1e-9, nr + 1)
    tbins = np.linspace(-np.pi, np.pi + 1e-9, nth + 1)
    zbins = np.linspace(bottom, top + 1e-9, nz + 1)

    ir = np.clip(np.digitize(r, rbins) - 1, 0, nr - 1)
    it = np.clip(np.digitize(theta, tbins) - 1, 0, nth - 1)
    iz = np.clip(np.digitize(z, zbins) - 1, 0, nz - 1)

    occ = np.zeros((nr, nth, nz), dtype=bool)
    occ[ir, it, iz] = True
    return float(occ.mean())


def load_metadata(path: Path, model: str) -> list[dict]:
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    keep = []
    for row in rows:
        row = dict(row)
        row["model"] = model
        if row.get("has_point") != "1":
            continue
        keep.append(row)
    return keep


def analyze_run(task: tuple[dict, int, int, int]) -> dict:
    row, nr, nth, nz = task
    run_path = Path(row["run_path"])
    point_path = run_path / "point.txt"
    config_path = run_path / "config.cym"

    out = {
        "model": row["model"],
        "group": row["group"],
        "case": row["case"],
        "xlink_regime": row["xlink_regime"],
        "run_dir": row["run_dir"],
        "run_path": row["run_path"],
        "status": "skipped_missing_point",
    }
    if not point_path.exists():
        return out

    try:
        cfg = parse_config_info(config_path)
        runs = cfg["runs"]
        final_run = runs[-1]
        nb_frames = final_run["nb_frames"]
        dt = final_run["frame_dt"]
        frames = parse_frames_point(point_path)

        total_expected = sum(r["nb_frames"] for r in runs)
        total_expected_with_initial = sum(r["nb_frames"] + 1 for r in runs)
        if len(frames) == total_expected:
            phase_frames = frames[-nb_frames:]
        elif len(frames) == total_expected_with_initial:
            phase_frames = frames[-(nb_frames + 1):][1:]
        else:
            out["status"] = f"error:expected {total_expected} or {total_expected_with_initial} frames but parsed {len(frames)}"
            return out

        space = cfg["space"]
        zmin = float(space.get("bottom", -20.0))
        zmax = float(space.get("top", 20.0))
        inner = float(space.get("inner", 0.0))
        outer = float(space.get("outer", 1.0))

        z_com = []
        theta_com = []
        Dz = []
        top_bias = []
        nematic_xy = []
        shell_occ = []

        for frame in phase_frames:
            reduced = reduce_frame_point(frame, zmin=zmin, zmax=zmax, nbins_z=80, zmid=0.0, Lz=(zmax - zmin))
            x = frame[:, 1]
            y = frame[:, 2]
            z = frame[:, 3]
            bands = compute_top_bottom_bands(z, top_band_frac=1 / 3)

            z_com.append(reduced["z_com"])
            theta_com.append(reduced["theta_com"])
            Dz.append(reduced["Dz"])
            top_bias.append(bands.top_bias)
            nematic_xy.append(nematic_order_xy(frame[:, :3]))
            shell_occ.append(
                shell_occupancy_fraction(
                    x,
                    y,
                    z,
                    inner=inner,
                    outer=outer,
                    bottom=zmin,
                    top=zmax,
                    nr=nr,
                    nth=nth,
                    nz=nz,
                )
            )

        z_com = np.asarray(z_com, dtype=float)
        theta_com = np.asarray(theta_com, dtype=float)
        vz_abs = np.abs(central_diff(z_com, dt))
        swirl_rate_abs = np.abs(central_diff(unwrap_angles(theta_com), dt))
        time_min = (dt * np.arange(1, nb_frames + 1, dtype=float)) / 60.0

        out.update(
            {
                "status": "ok",
                "n_frames": int(nb_frames),
                "dt_s": float(dt),
                "time_min": time_min.tolist(),
                "vz_abs": vz_abs.tolist(),
                "Dz": np.asarray(Dz, dtype=float).tolist(),
                "shell_occ": np.asarray(shell_occ, dtype=float).tolist(),
                "nematic_xy": np.asarray(nematic_xy, dtype=float).tolist(),
                "top_bias": np.asarray(top_bias, dtype=float).tolist(),
                "swirl_rate_abs": swirl_rate_abs.tolist(),
            }
        )
        return out
    except Exception as exc:
        out["status"] = f"error:{exc}"
        return out


def aggregate_runs(run_rows: list[dict]) -> list[dict]:
    grouped: dict[tuple[str, str, str, str], list[dict]] = defaultdict(list)
    for row in run_rows:
        if row.get("status") == "ok":
            grouped[(row["model"], row["group"], row["case"], row["xlink_regime"])].append(row)

    out_rows = []
    for (model, group, case, regime), items in sorted(grouped.items()):
        n_runs = len(items)
        n_frames = min(len(item["time_min"]) for item in items)
        time_min = np.asarray(items[0]["time_min"][:n_frames], dtype=float)

        metrics = {}
        for metric in [m["field"] for m in METRICS]:
            stack = np.vstack([np.asarray(item[metric][:n_frames], dtype=float) for item in items])
            mean = np.nanmean(stack, axis=0)
            if stack.shape[0] > 1:
                sem = np.nanstd(stack, axis=0, ddof=1) / math.sqrt(stack.shape[0])
            else:
                sem = np.zeros_like(mean)
            metrics[metric] = (mean, sem)

        for idx, t in enumerate(time_min, start=1):
            row = {
                "model": model,
                "group": group,
                "case": case,
                "xlink_regime": regime,
                "frame_index": idx,
                "time_min": float(t),
                "n_runs_used": n_runs,
            }
            for metric in [m["field"] for m in METRICS]:
                mean, sem = metrics[metric]
                row[f"{metric}_mean"] = float(mean[idx - 1])
                row[f"{metric}_sem"] = float(sem[idx - 1])
            out_rows.append(row)
    return out_rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.1)
    ax.spines["bottom"].set_linewidth(1.1)
    ax.tick_params(direction="out", width=1.1, labelsize=10)
    ax.grid(axis="y", color="#d7d7d7", linewidth=0.6, alpha=0.8)


def build_lookup(rows: list[dict]) -> dict[tuple[str, str, str, str], list[dict]]:
    grouped: dict[tuple[str, str, str, str], list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(row["model"], row["group"], row["case"], row["xlink_regime"])].append(row)
    for key in grouped:
        grouped[key].sort(key=lambda r: int(r["frame_index"]))
    return grouped


def compute_ylim(lookup: dict, family_spec: dict, metric: dict) -> tuple[float, float]:
    vals = []
    for regime in XLINK_ORDER:
        for case, _ in family_spec["order"]:
            for model in [m["key"] for m in MODEL_SPECS]:
                rows = lookup.get((model, family_spec["group"], case, regime), [])
                for row in rows:
                    mean = float(row[f'{metric["field"]}_mean'])
                    if math.isfinite(mean):
                        vals.append(mean)
    if not vals:
        return 0.0, 1.0
    lo = min(vals)
    hi = max(vals)
    pad = 0.1 * max(abs(hi), 1.0) if hi <= lo else 0.08 * (hi - lo)
    lo -= pad
    hi += pad
    if metric["clip_zero"]:
        lo = max(0.0, lo)
    return lo, hi


def make_timecourse_figure(lookup: dict, family_spec: dict, metric: dict, outdir: Path) -> list[Path]:
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.3), sharey=True)
    ylo, yhi = compute_ylim(lookup, family_spec, metric)

    for ax, regime in zip(axes, XLINK_ORDER):
        style_axes(ax)
        for case, label in family_spec["order"]:
            color = family_spec["colors"][case]
            for model_spec in MODEL_SPECS:
                rows = lookup.get((model_spec["key"], family_spec["group"], case, regime), [])
                if not rows:
                    continue
                t = np.array([float(r["time_min"]) for r in rows], dtype=float)
                y = np.array([float(r[f'{metric["field"]}_mean']) for r in rows], dtype=float)
                ax.plot(
                    t,
                    y,
                    color=color,
                    linestyle=model_spec["linestyle"],
                    linewidth=1.8,
                    alpha=0.95,
                )
        ax.set_ylim(ylo, yhi)
        ax.set_xlabel("Time after motor onset (min)", fontsize=11)
        ax.text(0.03, 0.94, regime, transform=ax.transAxes, ha="left", va="top", fontsize=10, color="#444444")

    axes[0].set_ylabel(metric["ylabel"], fontsize=11)

    case_handles = []
    case_labels = []
    for case, label in family_spec["order"]:
        case_handles.append(plt.Line2D([0], [0], color=family_spec["colors"][case], lw=2.0))
        case_labels.append(label)

    model_handles = []
    model_labels = []
    for model_spec in MODEL_SPECS:
        model_handles.append(plt.Line2D([0], [0], color="#444444", lw=2.0, linestyle=model_spec["linestyle"]))
        model_labels.append(model_spec["label"])

    fig.legend(case_handles, case_labels, loc="lower center", ncol=min(len(case_labels), 5), frameon=False, bbox_to_anchor=(0.5, -0.03), fontsize=10)
    fig.legend(model_handles, model_labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02), fontsize=10)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.83, bottom=0.18, wspace=0.12)

    outbase = outdir / f'{family_spec["slug"]}_{metric["slug"]}_time'
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return [outbase.with_suffix(".png"), outbase.with_suffix(".pdf")]


def main() -> None:
    ap = argparse.ArgumentParser(description="Aggregate and plot final-phase time-course comparisons for rotatable vs fixed-global campaigns.")
    ap.add_argument("--rotatable", type=Path, default=DEFAULT_ROTATABLE)
    ap.add_argument("--fixed", type=Path, default=DEFAULT_FIXED)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--nr", type=int, default=8)
    ap.add_argument("--nth", type=int, default=72)
    ap.add_argument("--nz", type=int, default=80)
    args = ap.parse_args()

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.linewidth": 1.1,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )

    tasks = []
    for row in load_metadata(args.rotatable, "rotatable"):
        tasks.append((row, args.nr, args.nth, args.nz))
    for row in load_metadata(args.fixed, "fixed_global"):
        tasks.append((row, args.nr, args.nth, args.nz))

    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            run_rows = list(ex.map(analyze_run, tasks))
    else:
        run_rows = [analyze_run(task) for task in tasks]

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "timecourse_run_status.csv", run_rows)

    condition_rows = aggregate_runs(run_rows)
    write_csv(args.outdir / "condition_timecourses.csv", condition_rows)

    lookup = build_lookup(condition_rows)
    outfiles = []
    for family_spec in FAMILY_SPECS:
        for metric in METRICS:
            outfiles.extend(make_timecourse_figure(lookup, family_spec, metric, args.outdir))

    print(args.outdir / "timecourse_run_status.csv")
    print(args.outdir / "condition_timecourses.csv")
    for path in outfiles:
        print(path)


if __name__ == "__main__":
    main()
