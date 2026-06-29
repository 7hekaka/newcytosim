#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS_DIR = ROOT / "analysis"
for path in (str(ROOT), str(ANALYSIS_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)

from analysis.metrics import compute_top_bottom_bands
from analysis import plot_timecourse_panels as base
from analysis import run_speckle_timecourse_batch as speck
from trust.analyze_axial import central_diff, reduce_frame_point, unwrap_angles


NBINS_Z = 80
NBINS_S = 96
SEAM_BINS = 180
SNAPSHOT_TIMES_MIN = [
    ("t00", 0.0, "0.0 min"),
    ("t100", 100.0 / 60.0, "1.7 min"),
    ("t120", 120.0 / 60.0, "2.0 min"),
    ("t360", 360.0 / 60.0, "6.0 min"),
    ("t720", 720.0 / 60.0, "12.0 min"),
]
SUMMARY_METRICS = [
    {"field": "mean_vz_abs", "label": "Mean |v_z| (um/s)", "slug": "mean_vz_abs", "clip_zero": True, "fmt": "{:.4f}", "cmap": "magma"},
    {"field": "auc_vz_abs", "label": "AUC(|v_z|) (um)", "slug": "auc_vz_abs", "clip_zero": True, "fmt": "{:.2f}", "cmap": "viridis"},
    {"field": "excess_auc_vz_abs", "label": "Excess AUC above matched control (um)", "slug": "excess_auc_vz_abs", "clip_zero": True, "fmt": "{:.2f}", "cmap": "viridis"},
    {"field": "abs_net_z_displacement", "label": "|Net z shift| (um)", "slug": "abs_net_z_displacement", "clip_zero": True, "fmt": "{:.2f}", "cmap": "plasma"},
    {"field": "directionality_index", "label": "Directionality index", "slug": "directionality_index", "clip_zero": True, "fmt": "{:.2f}", "cmap": "cividis"},
    {"field": "active_transport_fraction", "label": "Active transport fraction", "slug": "active_transport_fraction", "clip_zero": True, "fmt": "{:.2f}", "cmap": "magma"},
    {"field": "matched_active_transport_fraction", "label": "Matched active transport fraction", "slug": "matched_active_transport_fraction", "clip_zero": True, "fmt": "{:.2f}", "cmap": "magma"},
    {"field": "longest_active_streak_min", "label": "Longest matched-active streak (min)", "slug": "longest_active_streak_min", "clip_zero": True, "fmt": "{:.2f}", "cmap": "magma"},
    {"field": "mean_vz_signed", "label": "Mean v_z (um/s)", "slug": "mean_vz_signed", "clip_zero": False, "fmt": "{:.4f}", "cmap": "coolwarm"},
    {"field": "mass_normalized_auc", "label": "Mass-normalized AUC", "slug": "mass_normalized_auc", "clip_zero": True, "fmt": "{:.4f}", "cmap": "viridis"},
    {"field": "mass_retention", "label": "Mass retention (final/onset)", "slug": "mass_retention", "clip_zero": True, "fmt": "{:.2f}", "cmap": "cividis"},
    {"field": "swirl_penalty", "label": "Swirl penalty", "slug": "swirl_penalty", "clip_zero": True, "fmt": "{:.2f}", "cmap": "magma"},
    {"field": "chewer_band_mass_fraction_final", "label": "Final chewer-band mass fraction", "slug": "chewer_band_mass_fraction_final", "clip_zero": True, "fmt": "{:.2f}", "cmap": "YlGnBu"},
    {"field": "bottom_half_mass_fraction_final", "label": "Final bottom-half mass fraction", "slug": "bottom_half_mass_fraction_final", "clip_zero": True, "fmt": "{:.2f}", "cmap": "YlGnBu"},
    {"field": "seed_zone_retention", "label": "Seed-zone retention (final/onset)", "slug": "seed_zone_retention", "clip_zero": True, "fmt": "{:.2f}", "cmap": "cividis"},
    {"field": "axial_spread_entropy_mean", "label": "Mean axial spread entropy", "slug": "axial_spread_entropy_mean", "clip_zero": True, "fmt": "{:.2f}", "cmap": "viridis"},
    {"field": "speckle_length_proxy_final_um", "label": "Final filament mass proxy (um)", "slug": "speckle_length_proxy_final_um", "clip_zero": True, "fmt": "{:.0f}", "cmap": "viridis"},
    {"field": "speckle_length_proxy_motor_onset_um", "label": "Motor-onset mass proxy (um)", "slug": "speckle_length_proxy_motor_onset_um", "clip_zero": True, "fmt": "{:.0f}", "cmap": "viridis"},
]
OVERVIEW_METRICS = [
    "auc_vz_abs",
    "excess_auc_vz_abs",
    "matched_active_transport_fraction",
    "mass_retention",
    "chewer_band_mass_fraction_final",
    "axial_spread_entropy_mean",
]
TIMECOURSE_FIELDS = [
    ("vz_abs", "Mean |v_z| (um/s)", True),
    ("speckle_length_proxy_um", "Filament mass proxy (um)", True),
    ("chewer_band_mass_fraction", "Chewer-band mass fraction", True),
    ("axial_spread_entropy", "Axial spread entropy", True),
]
PHASE_NAMES = ["growth", "xlink", "motor"]
CONDITION_COLORS = {
    "c00_nomotor_noxlink": "#7f7f7f",
    "c01_nomotor_xlink": "#4d4d4d",
    "c02_rotatable_noxlink": "#5ab4ac",
    "c03_rotatable_xlink": "#1b9e77",
    "c04_fixed_global_noxlink": "#fc8d62",
    "c05_fixed_global_xlink": "#d95f02",
}

SPACE_START_RE = re.compile(r"^\s*set\s+space\s+([A-Za-z0-9_]+)")
SPACE_INSTANCE_RE = re.compile(r"^\s*new\s+([A-Za-z0-9_]+)")
SPACE_VALUE_RE = re.compile(r"^\s*(inner|outer|top|bottom)\s*=\s*([0-9.eE+-]+)")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Analyze returned turnover campaign runs with full-timeline speckle metrics, kymographs, and unwrapped snapshots.")
    ap.add_argument("--campaign-root", required=True, help="Campaign directory containing manifest.csv and job/save outputs")
    ap.add_argument("--save-dir", required=True, help="Directory containing r0001, r0002, ... run folders")
    ap.add_argument("--outdir", required=True, help="Output directory for plots and CSVs")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--interval", type=float, default=2.0)
    ap.add_argument("--speckle-name", default="speckle_i2.txt")
    ap.add_argument("--report-bin", default=None)
    return ap.parse_args()


def pick_report_bin(explicit: str | None) -> Path:
    if explicit:
        path = Path(explicit).expanduser().resolve()
        if path.exists():
            return path
        raise FileNotFoundError(f"report binary not found: {path}")
    for cand in [ROOT / "vesicles" / "report", ROOT / "build" / "bin" / "report", ROOT / "bin" / "report"]:
        if cand.exists() and os.access(cand, os.X_OK):
            return cand
    raise FileNotFoundError("Could not find report executable in vesicles/report, build/bin/report, or bin/report")


def load_manifest(campaign_root: Path, save_dir: Path) -> list[dict]:
    manifest_path = campaign_root / "manifest.csv"
    rows = []
    for row in csv.DictReader(manifest_path.open(encoding="utf-8")):
        item = dict(row)
        index = int(item["index"])
        run_dir = f"r{index + 1:04d}"
        run_path = save_dir / run_dir
        item["index"] = index
        item["run_dir"] = run_dir
        item["run_path"] = str(run_path)
        item["include_xlinks"] = int(item["include_xlinks"])
        item["n_clusters"] = int(item["n_clusters"])
        item["motors_per_cluster"] = int(item["motors_per_cluster"])
        item["total_motors"] = int(item["total_motors"])
        item["scenario_short"] = shorten_scenario_label(item["scenario_label"])
        item["condition_short"] = shorten_condition_label(item["condition_label"])
        rows.append(item)
    return rows


def shorten_scenario_label(label: str) -> str:
    mapping = {
        "Throughout system, mixed polarity": "Full height\nmixed",
        "Throughout system, aligned upward": "Full height\naligned",
        "Top two-thirds, aligned upward": "Top 2/3\naligned",
        "Top one-third, aligned upward": "Top 1/3\naligned",
        "Top cap, aligned upward": "Top cap\naligned",
    }
    return mapping.get(label, label)


def shorten_condition_label(label: str) -> str:
    mapping = {
        "No motors, no crosslinkers": "No motors\nNo xlinks",
        "No motors, with crosslinkers": "No motors\n+ xlinks",
        "Rotatable motors, no crosslinkers": "Rotatable\nNo xlinks",
        "Rotatable motors, with crosslinkers": "Rotatable\n+ xlinks",
        "Fixed-global motors, no crosslinkers": "Fixed global\nNo xlinks",
        "Fixed-global motors, with crosslinkers": "Fixed global\n+ xlinks",
    }
    return mapping.get(label, label)


def phase_segments(frames: list[np.ndarray], runs: list[dict]) -> tuple[list[dict], str]:
    expected_no = sum(run["nb_frames"] for run in runs)
    expected_with = sum(run["nb_frames"] + 1 for run in runs)
    segments = []
    start_min = 0.0
    cursor = 0

    def append_segment(idx: int, seg_frames: list[np.ndarray], dt_s: float, duration_min: float, mode_label: str) -> None:
        nonlocal start_min, segments
        dt_min = dt_s / 60.0
        if mode_label in {"with_initial", "truncated_with_initial"} and idx == 0:
            time_min = start_min + dt_min * np.arange(0, len(seg_frames), dtype=float)
        else:
            offset = 1 if mode_label in {"with_initial", "truncated_with_initial"} else 1
            time_min = start_min + dt_min * np.arange(offset, offset + len(seg_frames), dtype=float)
        actual_duration_min = dt_min * len(seg_frames)
        segments.append(
            {
                "phase": PHASE_NAMES[idx] if idx < len(PHASE_NAMES) else f"phase_{idx}",
                "frames": seg_frames,
                "time_min": time_min,
                "dt_s": float(dt_s),
                "start_min": start_min,
                "end_min": start_min + (duration_min if len(seg_frames) and actual_duration_min >= duration_min - 1e-12 else actual_duration_min),
            }
        )
        start_min += duration_min if len(seg_frames) and actual_duration_min >= duration_min - 1e-12 else actual_duration_min

    if len(frames) == expected_no:
        mode = "no_initial"
        for idx, run in enumerate(runs):
            nb = int(run["nb_frames"])
            seg_frames = frames[cursor:cursor + nb]
            cursor += nb
            append_segment(idx, seg_frames, float(run["frame_dt"]), float(run["duration"]) / 60.0, mode)
        return segments, mode

    if len(frames) == expected_with:
        mode = "with_initial"
        for idx, run in enumerate(runs):
            nb = int(run["nb_frames"])
            raw_frames = frames[cursor:cursor + nb + 1]
            cursor += nb + 1
            seg_frames = raw_frames if idx == 0 else raw_frames[1:]
            append_segment(idx, seg_frames, float(run["frame_dt"]), float(run["duration"]) / 60.0, mode)
        return segments, mode

    full_no_before_last = sum(int(run["nb_frames"]) for run in runs[:-1])
    if full_no_before_last < len(frames) < expected_no:
        mode = "truncated_no_initial"
        for idx, run in enumerate(runs[:-1]):
            nb = int(run["nb_frames"])
            seg_frames = frames[cursor:cursor + nb]
            cursor += nb
            append_segment(idx, seg_frames, float(run["frame_dt"]), float(run["duration"]) / 60.0, "no_initial")
        remaining = frames[cursor:]
        if remaining:
            run = runs[-1]
            append_segment(len(runs) - 1, remaining, float(run["frame_dt"]), float(run["duration"]) / 60.0, mode)
            return segments, mode

    full_with_before_last = sum((int(run["nb_frames"]) + 1) for run in runs[:-1])
    if full_with_before_last + 1 <= len(frames) < expected_with:
        mode = "truncated_with_initial"
        for idx, run in enumerate(runs[:-1]):
            nb = int(run["nb_frames"])
            raw_frames = frames[cursor:cursor + nb + 1]
            cursor += nb + 1
            seg_frames = raw_frames if idx == 0 else raw_frames[1:]
            append_segment(idx, seg_frames, float(run["frame_dt"]), float(run["duration"]) / 60.0, "with_initial")
        remaining_raw = frames[cursor:]
        if len(remaining_raw) >= 1:
            run = runs[-1]
            seg_frames = remaining_raw[1:]
            append_segment(len(runs) - 1, seg_frames, float(run["frame_dt"]), float(run["duration"]) / 60.0, mode)
            return segments, mode

    raise ValueError(f"expected {expected_no} or {expected_with} frames, parsed {len(frames)}")

def stable_diff(values: np.ndarray, dt_s: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return values
    if values.size == 1:
        return np.zeros_like(values)
    return np.asarray(central_diff(values, dt_s), dtype=float)


def optional_float(row: dict, key: str, default: float = float("nan")) -> float:
    value = row.get(key, "")
    if value in {"", None}:
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_named_spaces(config_path: Path) -> dict[str, dict[str, float]]:
    spaces: dict[str, dict[str, float]] = {}
    current: str | None = None
    for line in config_path.read_text(encoding="utf-8", errors="replace").splitlines():
        start = SPACE_START_RE.match(line)
        if start:
            spaces.setdefault(start.group(1), {})
            current = None
            continue
        instance = SPACE_INSTANCE_RE.match(line)
        if instance and instance.group(1) in spaces:
            current = instance.group(1)
            spaces.setdefault(current, {})
            continue
        if current is None:
            continue
        if line.strip().startswith("}"):
            current = None
            continue
        value = SPACE_VALUE_RE.match(line)
        if value:
            spaces[current][value.group(1)] = float(value.group(2))
    return spaces


def pick_analysis_space(config_path: Path, fallback: dict) -> dict[str, float]:
    spaces = parse_named_spaces(config_path)
    if "stripbox" in spaces:
        return spaces["stripbox"]
    for name, space in spaces.items():
        if "chewer" not in name and {"bottom", "top"}.issubset(space):
            return space
    return dict(fallback)


def nanmean_or_nan(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return float("nan")
    return float(np.mean(finite))


def ratio_or_nan(numerator: float, denominator: float) -> float:
    if not np.isfinite(numerator) or not np.isfinite(denominator) or denominator <= 0:
        return float("nan")
    return float(numerator / denominator)


def finite_value_at(values: np.ndarray, idx: int) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0 or idx < 0 or idx >= arr.size:
        return float("nan")
    value = float(arr[idx])
    return value if np.isfinite(value) else float("nan")


def band_fraction_from_rho(rho_z: np.ndarray, z_edges: np.ndarray, z_low: float, z_high: float) -> np.ndarray:
    rho = np.asarray(rho_z, dtype=float)
    if rho.size == 0 or not np.isfinite(z_low) or not np.isfinite(z_high):
        rows = rho.shape[0] if rho.ndim >= 2 else 1
        return np.full(rows, np.nan, dtype=float)
    if z_high < z_low:
        z_low, z_high = z_high, z_low
    centers = 0.5 * (np.asarray(z_edges[:-1], dtype=float) + np.asarray(z_edges[1:], dtype=float))
    band_mask = (centers >= z_low) & (centers <= z_high)
    if not np.any(band_mask):
        rows = rho.shape[0] if rho.ndim >= 2 else 1
        return np.full(rows, np.nan, dtype=float)
    rho2 = np.atleast_2d(np.maximum(rho, 0.0))
    total = np.nansum(rho2, axis=1)
    band = np.nansum(rho2[:, band_mask], axis=1)
    out = np.full(rho2.shape[0], np.nan, dtype=float)
    valid = total > 0
    out[valid] = band[valid] / total[valid]
    return out


def axial_spread_entropy_from_rho(rho_z: np.ndarray) -> np.ndarray:
    rho = np.atleast_2d(np.maximum(np.asarray(rho_z, dtype=float), 0.0))
    total = np.nansum(rho, axis=1)
    entropy = np.full(rho.shape[0], np.nan, dtype=float)
    valid = total > 0
    if not np.any(valid):
        return entropy
    probs = rho[valid] / total[valid, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(probs > 0, probs * np.log(probs), 0.0)
    denom = math.log(rho.shape[1]) if rho.shape[1] > 1 else 1.0
    entropy[valid] = -np.nansum(terms, axis=1) / denom
    return entropy


def longest_true_streak_min(time_min: np.ndarray, mask: np.ndarray) -> float:
    flags = np.asarray(mask, dtype=bool)
    if flags.size == 0:
        return float("nan")
    longest = 0
    current = 0
    for flag in flags:
        if flag:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    if longest == 0:
        return 0.0
    times = np.asarray(time_min, dtype=float)
    diffs = np.diff(times[np.isfinite(times)])
    finite_diffs = diffs[diffs > 0]
    dt_min = float(np.median(finite_diffs)) if finite_diffs.size else 0.0
    return float(longest * dt_min)


def analyze_frame(frame: np.ndarray, *, zmin: float, zmax: float, inner: float, outer: float, interval: float, nr: int = 8, nth: int = 72, nz: int = 80) -> dict:
    if frame.size == 0:
        return {
            "z_com": float("nan"),
            "theta_com": float("nan"),
            "Dz": float("nan"),
            "top_bias": float("nan"),
            "nematic_xy": float("nan"),
            "shell_occ": float("nan"),
            "rho_z": np.zeros(NBINS_Z, dtype=float),
            "speckle_length_proxy_um": 0.0,
        }
    reduced = reduce_frame_point(frame[:, :4], zmin=zmin, zmax=zmax, nbins_z=NBINS_Z, zmid=0.0, Lz=(zmax - zmin))
    x = frame[:, 1]
    y = frame[:, 2]
    z = frame[:, 3]
    bands = compute_top_bottom_bands(z, top_band_frac=1 / 3)
    return {
        "z_com": float(reduced["z_com"]),
        "theta_com": float(reduced["theta_com"]),
        "Dz": float(reduced["Dz"]),
        "top_bias": float(bands.top_bias),
        "nematic_xy": float(speck.nematic_order_xy_from_speckles(frame)),
        "shell_occ": float(base.shell_occupancy_fraction(x, y, z, inner=inner, outer=outer, bottom=zmin, top=zmax, nr=nr, nth=nth, nz=nz)),
        "rho_z": np.asarray(reduced["rho_z"], dtype=float),
        "speckle_length_proxy_um": float(interval * frame.shape[0]),
    }


def choose_theta_seam(frame: np.ndarray, bins: int = SEAM_BINS) -> float:
    if frame.size == 0:
        return 0.0
    theta = np.mod(np.arctan2(frame[:, 2], frame[:, 1]), 2.0 * math.pi)
    hist, edges = np.histogram(theta, bins=bins, range=(0.0, 2.0 * math.pi))
    idx = int(np.argmin(hist))
    return float(0.5 * (edges[idx] + edges[idx + 1]))


def frame_to_density(frame: np.ndarray, *, seam: float, inner: float, outer: float, bottom: float, top: float, nbins_s: int, nbins_z: int) -> tuple[np.ndarray, float]:
    if frame.size == 0:
        r_mid = 0.5 * (inner + outer)
        return np.zeros((nbins_z, nbins_s), dtype=float), 2.0 * math.pi * r_mid
    x = frame[:, 1]
    y = frame[:, 2]
    z = frame[:, 3]
    theta = np.mod(np.arctan2(y, x) - seam, 2.0 * math.pi)
    r_mid = 0.5 * (inner + outer)
    smax = 2.0 * math.pi * r_mid
    s = theta * r_mid
    hist, _z, _s = np.histogram2d(z, s, bins=[nbins_z, nbins_s], range=[[bottom, top], [0.0, smax]])
    mean = float(hist.mean())
    if mean > 0:
        hist = hist / mean
    return hist, smax


def frame_segments(frame: np.ndarray, seam: float, smax: float, r_mid: float) -> list[tuple[np.ndarray, np.ndarray]]:
    segments = []
    if frame.size == 0:
        return segments
    fiber_ids = np.unique(frame[:, 0].astype(int))
    for fid in fiber_ids:
        pts = frame[frame[:, 0] == fid]
        pts = pts[np.argsort(pts[:, 4])]
        theta = np.mod(np.arctan2(pts[:, 2], pts[:, 1]) - seam, 2.0 * math.pi)
        s = theta * r_mid
        z = pts[:, 3]
        if s.size < 2:
            continue
        jumps = np.where(np.abs(np.diff(s)) > 0.5 * smax)[0]
        start = 0
        if jumps.size:
            for jump in jumps:
                seg = slice(start, jump + 1)
                if (jump + 1 - start) >= 2:
                    segments.append((s[seg], z[seg]))
                start = jump + 1
        if (s.size - start) >= 2:
            segments.append((s[start:], z[start:]))
    return segments


def nearest_index(time_min: np.ndarray, target_min: float) -> int:
    if time_min.size == 0:
        return 0
    return int(np.argmin(np.abs(time_min - target_min)))


def analyze_run(task: tuple[dict, str, float, str]) -> dict:
    row, report_bin_s, interval, speckle_name = task
    run_path = Path(row["run_path"])
    out = {
        "index": int(row["index"]),
        "scenario": row["scenario"],
        "scenario_label": row["scenario_label"],
        "scenario_short": row["scenario_short"],
        "condition": row["condition"],
        "condition_label": row["condition_label"],
        "condition_short": row["condition_short"],
        "motor_mode": row["motor_mode"],
        "include_xlinks": int(row["include_xlinks"]),
        "n_clusters": int(row["n_clusters"]),
        "motors_per_cluster": int(row["motors_per_cluster"]),
        "total_motors": int(row["total_motors"]),
        "run_dir": row["run_dir"],
        "run_path": str(run_path),
        "status": "pending",
    }
    try:
        speckle_path = run_path / speckle_name
        speck.ensure_speckle(run_path, speckle_path, Path(report_bin_s), interval)
        cfg = base.parse_config_info(run_path / "config.cym")
        frames = speck.parse_frames_speckle(speckle_path)
        segments, frame_mode = phase_segments(frames, cfg["runs"])

        space = pick_analysis_space(run_path / "config.cym", cfg["space"])
        zmin = float(space.get("bottom", -20.0))
        zmax = float(space.get("top", 20.0))
        inner = float(space.get("inner", 0.0))
        outer = float(space.get("outer", inner + 1.0))

        time_all = []
        phase_all = []
        z_com_all = []
        theta_all = []
        vz_signed_all = []
        vz_abs_all = []
        Dz_all = []
        top_bias_all = []
        nematic_xy_all = []
        shell_occ_all = []
        swirl_all = []
        rho_all = []
        length_proxy_all = []
        frames_all = []

        for segment in segments:
            seg_features = [analyze_frame(frame, zmin=zmin, zmax=zmax, inner=inner, outer=outer, interval=interval) for frame in segment["frames"]]
            z_com = np.asarray([feat["z_com"] for feat in seg_features], dtype=float)
            theta = np.asarray([feat["theta_com"] for feat in seg_features], dtype=float)
            vz_signed = stable_diff(z_com, float(segment["dt_s"]))
            swirl_abs = np.abs(stable_diff(unwrap_angles(theta), float(segment["dt_s"])))

            time_all.extend(np.asarray(segment["time_min"], dtype=float).tolist())
            phase_all.extend([segment["phase"]] * len(seg_features))
            z_com_all.extend(z_com.tolist())
            theta_all.extend(theta.tolist())
            vz_signed_all.extend(vz_signed.tolist())
            vz_abs_all.extend(np.abs(vz_signed).tolist())
            Dz_all.extend([feat["Dz"] for feat in seg_features])
            top_bias_all.extend([feat["top_bias"] for feat in seg_features])
            nematic_xy_all.extend([feat["nematic_xy"] for feat in seg_features])
            shell_occ_all.extend([feat["shell_occ"] for feat in seg_features])
            swirl_all.extend(swirl_abs.tolist())
            rho_all.extend([feat["rho_z"] for feat in seg_features])
            length_proxy_all.extend([feat["speckle_length_proxy_um"] for feat in seg_features])
            frames_all.extend(segment["frames"])

        time_min = np.asarray(time_all, dtype=float)
        z_com = np.asarray(z_com_all, dtype=float)
        vz_abs = np.asarray(vz_abs_all, dtype=float)
        vz_signed = np.asarray(vz_signed_all, dtype=float)
        length_proxy = np.asarray(length_proxy_all, dtype=float)

        growth_end_min = float(segments[0]["end_min"])
        motor_onset_min = float(segments[1]["end_min"]) if len(segments) > 1 else growth_end_min

        onset_idx = nearest_index(time_min, motor_onset_min)
        final_idx = len(time_min) - 1 if time_min.size else 0
        seam = choose_theta_seam(frames_all[0]) if frames_all else 0.0
        rho_arr = np.stack(rho_all, axis=0)
        z_edges = np.linspace(zmin, zmax, NBINS_Z + 1, dtype=float)
        scenario_z_low = optional_float(row, "scenario_z_low", zmin)
        scenario_z_high = optional_float(row, "scenario_z_high", zmax)
        if not np.isfinite(scenario_z_low):
            scenario_z_low = zmin
        if not np.isfinite(scenario_z_high):
            scenario_z_high = zmax
        chewer_zone_bottom = optional_float(row, "chewer_zone_bottom")
        chewer_zone_top = optional_float(row, "chewer_zone_top")
        bottom_half_top = zmin + 0.5 * (zmax - zmin)
        seed_zone_mass_fraction = band_fraction_from_rho(rho_arr, z_edges, scenario_z_low, scenario_z_high)
        chewer_band_mass_fraction = band_fraction_from_rho(rho_arr, z_edges, chewer_zone_bottom, chewer_zone_top)
        bottom_half_mass_fraction = band_fraction_from_rho(rho_arr, z_edges, zmin, bottom_half_top)
        axial_spread_entropy = axial_spread_entropy_from_rho(rho_arr)

        out.update(
            {
                "status": "ok",
                "frame_mode": frame_mode,
                "time_min": time_min.tolist(),
                "phase_name": phase_all,
                "z_com": z_com.tolist(),
                "theta_com": np.asarray(theta_all, dtype=float).tolist(),
                "vz_abs": vz_abs.tolist(),
                "vz_signed": vz_signed.tolist(),
                "Dz": np.asarray(Dz_all, dtype=float).tolist(),
                "top_bias": np.asarray(top_bias_all, dtype=float).tolist(),
                "nematic_xy": np.asarray(nematic_xy_all, dtype=float).tolist(),
                "shell_occ": np.asarray(shell_occ_all, dtype=float).tolist(),
                "swirl_rate_abs": np.asarray(swirl_all, dtype=float).tolist(),
                "speckle_length_proxy_um": length_proxy.tolist(),
                "seed_zone_mass_fraction": seed_zone_mass_fraction.tolist(),
                "chewer_band_mass_fraction": chewer_band_mass_fraction.tolist(),
                "bottom_half_mass_fraction": bottom_half_mass_fraction.tolist(),
                "axial_spread_entropy": axial_spread_entropy.tolist(),
                "rho_z": rho_arr.tolist(),
                "z_edges": z_edges.tolist(),
                "growth_end_min": growth_end_min,
                "motor_onset_min": motor_onset_min,
                "speckle_length_proxy_motor_onset_um": float(length_proxy[onset_idx]) if length_proxy.size else float("nan"),
                "speckle_length_proxy_final_um": float(length_proxy[final_idx]) if length_proxy.size else float("nan"),
                "seed_zone_mass_fraction_motor_onset": finite_value_at(seed_zone_mass_fraction, onset_idx),
                "seed_zone_mass_fraction_final": finite_value_at(seed_zone_mass_fraction, final_idx),
                "chewer_band_mass_fraction_motor_onset": finite_value_at(chewer_band_mass_fraction, onset_idx),
                "chewer_band_mass_fraction_final": finite_value_at(chewer_band_mass_fraction, final_idx),
                "bottom_half_mass_fraction_motor_onset": finite_value_at(bottom_half_mass_fraction, onset_idx),
                "bottom_half_mass_fraction_final": finite_value_at(bottom_half_mass_fraction, final_idx),
                "axial_spread_entropy_motor_onset": finite_value_at(axial_spread_entropy, onset_idx),
                "axial_spread_entropy_final": finite_value_at(axial_spread_entropy, final_idx),
                "top_bias_final": float(np.asarray(top_bias_all, dtype=float)[final_idx]) if top_bias_all else float("nan"),
                "nematic_xy_final": float(np.asarray(nematic_xy_all, dtype=float)[final_idx]) if nematic_xy_all else float("nan"),
                "inner_radius": inner,
                "outer_radius": outer,
                "zmin": zmin,
                "zmax": zmax,
                "scenario_z_low": scenario_z_low,
                "scenario_z_high": scenario_z_high,
                "chewer_zone_bottom": chewer_zone_bottom,
                "chewer_zone_top": chewer_zone_top,
                "theta_seam": seam,
                "snapshot_indices": {slug: nearest_index(time_min, target_min) for slug, target_min, _ in SNAPSHOT_TIMES_MIN},
                "snapshot_times": {slug: float(time_min[nearest_index(time_min, target_min)]) for slug, target_min, _ in SNAPSHOT_TIMES_MIN},
                "snapshot_densities": {},
                "snapshot_frames": {},
            }
        )

        for slug, _target_min, _label in SNAPSHOT_TIMES_MIN:
            idx = out["snapshot_indices"][slug]
            density, smax = frame_to_density(frames_all[idx], seam=seam, inner=inner, outer=outer, bottom=zmin, top=zmax, nbins_s=NBINS_S, nbins_z=NBINS_Z)
            out["snapshot_densities"][slug] = density.tolist()
            out["snapshot_frames"][slug] = np.asarray(frames_all[idx], dtype=float).tolist()
            out[f"{slug}_smax"] = smax
        return out
    except Exception as exc:
        out["status"] = f"error:{exc}"
        return out


def compute_active_thresholds(run_rows: list[dict]) -> dict[int, float]:
    grouped: dict[int, list[float]] = defaultdict(list)
    for row in run_rows:
        if row.get("status") != "ok" or row.get("motor_mode") != "none":
            continue
        time = np.asarray(row["time_min"], dtype=float)
        vz_abs = np.asarray(row["vz_abs"], dtype=float)
        mask = time >= float(row["motor_onset_min"]) - 1e-9
        grouped[int(row["include_xlinks"])].extend(vz_abs[mask].tolist())
    thresholds = {}
    for include_xlinks, values in grouped.items():
        arr = np.asarray(values, dtype=float)
        thresholds[include_xlinks] = float(np.nanpercentile(arr, 95.0)) if arr.size else float("nan")
    return thresholds


def compute_scenario_matched_thresholds(run_rows: list[dict]) -> dict[tuple[str, int], float]:
    grouped: dict[tuple[str, int], list[float]] = defaultdict(list)
    for row in run_rows:
        if row.get("status") != "ok" or row.get("motor_mode") != "none":
            continue
        time = np.asarray(row["time_min"], dtype=float)
        vz_abs = np.asarray(row["vz_abs"], dtype=float)
        mask = time >= float(row["motor_onset_min"]) - 1e-9
        grouped[(row["scenario"], int(row["include_xlinks"]))].extend(vz_abs[mask].tolist())
    thresholds = {}
    for key, values in grouped.items():
        arr = np.asarray(values, dtype=float)
        thresholds[key] = float(np.nanpercentile(arr, 95.0)) if arr.size else float("nan")
    return thresholds


def first_sustained_time(time_min: np.ndarray, signal: np.ndarray, threshold: float, consecutive: int = 3) -> float:
    if not np.isfinite(threshold):
        return float("nan")
    mask = signal > threshold
    count = 0
    for idx, flag in enumerate(mask):
        count = count + 1 if flag else 0
        if count >= consecutive:
            return float(time_min[idx - consecutive + 1])
    return float("nan")


def compute_run_metrics(run_rows: list[dict], thresholds: dict[int, float], matched_thresholds: dict[tuple[str, int], float]) -> list[dict]:
    rows = []
    for row in run_rows:
        if row.get("status") != "ok":
            continue
        time_min = np.asarray(row["time_min"], dtype=float)
        time_s = time_min * 60.0
        vz_abs = np.asarray(row["vz_abs"], dtype=float)
        vz_signed = np.asarray(row["vz_signed"], dtype=float)
        z_com = np.asarray(row["z_com"], dtype=float)
        swirl_rate_abs = np.asarray(row["swirl_rate_abs"], dtype=float)
        mass_proxy = np.asarray(row["speckle_length_proxy_um"], dtype=float)
        seed_zone = np.asarray(row["seed_zone_mass_fraction"], dtype=float)
        chewer_band = np.asarray(row["chewer_band_mass_fraction"], dtype=float)
        bottom_half = np.asarray(row["bottom_half_mass_fraction"], dtype=float)
        axial_entropy = np.asarray(row["axial_spread_entropy"], dtype=float)
        mask = time_min >= float(row["motor_onset_min"]) - 1e-9
        if not np.any(mask):
            continue
        time_min_post = time_min[mask]
        time_s_post = time_s[mask]
        vz_abs_post = vz_abs[mask]
        vz_signed_post = vz_signed[mask]
        z_post = z_com[mask]
        swirl_post = swirl_rate_abs[mask]
        seed_zone_post = seed_zone[mask]
        chewer_band_post = chewer_band[mask]
        bottom_half_post = bottom_half[mask]
        axial_entropy_post = axial_entropy[mask]
        peak_idx = int(np.nanargmax(vz_abs_post)) if vz_abs_post.size else 0
        net_z_displacement = float(z_post[-1] - z_post[0]) if z_post.size else float("nan")
        abs_net_z_displacement = abs(net_z_displacement) if np.isfinite(net_z_displacement) else float("nan")
        auc_vz_abs = float(np.trapezoid(vz_abs_post, time_s_post)) if time_s_post.size else float("nan")
        threshold = thresholds.get(int(row["include_xlinks"]), float("nan"))
        matched_threshold = matched_thresholds.get((row["scenario"], int(row["include_xlinks"])), float("nan"))
        matched_active_mask = vz_abs_post > matched_threshold if np.isfinite(matched_threshold) else np.full(vz_abs_post.shape, False, dtype=bool)
        excess_auc = float(np.trapezoid(np.maximum(vz_abs_post - matched_threshold, 0.0), time_s_post)) if np.isfinite(matched_threshold) and time_s_post.size else float("nan")
        mean_vz_abs = float(np.nanmean(vz_abs_post))
        mean_swirl = nanmean_or_nan(swirl_post)
        mass_onset = float(row["speckle_length_proxy_motor_onset_um"])
        mass_final = float(row["speckle_length_proxy_final_um"])
        seed_onset = float(row["seed_zone_mass_fraction_motor_onset"])
        seed_final = float(row["seed_zone_mass_fraction_final"])
        chewer_onset = float(row["chewer_band_mass_fraction_motor_onset"])
        chewer_final = float(row["chewer_band_mass_fraction_final"])
        bottom_half_onset = float(row["bottom_half_mass_fraction_motor_onset"])
        bottom_half_final = float(row["bottom_half_mass_fraction_final"])
        entropy_onset = float(row["axial_spread_entropy_motor_onset"])
        entropy_final = float(row["axial_spread_entropy_final"])
        rows.append(
            {
                "scenario": row["scenario"],
                "scenario_label": row["scenario_label"],
                "scenario_short": row["scenario_short"],
                "condition": row["condition"],
                "condition_label": row["condition_label"],
                "condition_short": row["condition_short"],
                "motor_mode": row["motor_mode"],
                "include_xlinks": int(row["include_xlinks"]),
                "n_clusters": int(row["n_clusters"]),
                "motors_per_cluster": int(row["motors_per_cluster"]),
                "total_motors": int(row["total_motors"]),
                "run_dir": row["run_dir"],
                "run_path": row["run_path"],
                "growth_end_min": float(row["growth_end_min"]),
                "motor_onset_min": float(row["motor_onset_min"]),
                "threshold_vz_abs": threshold,
                "threshold_vz_abs_matched": matched_threshold,
                "mean_vz_abs": mean_vz_abs,
                "mean_vz_signed": float(np.nanmean(vz_signed_post)),
                "peak_vz_abs": float(np.nanmax(vz_abs_post)),
                "time_to_peak_min": float(time_min_post[peak_idx] - row["motor_onset_min"]) if time_min_post.size else float("nan"),
                "auc_vz_abs": auc_vz_abs,
                "excess_auc_vz_abs": excess_auc,
                "directionality_index": (abs_net_z_displacement / auc_vz_abs) if np.isfinite(abs_net_z_displacement) and np.isfinite(auc_vz_abs) and auc_vz_abs > 0 else float("nan"),
                "active_transport_fraction": float(np.mean(vz_abs_post > threshold)) if np.isfinite(threshold) else float("nan"),
                "matched_active_transport_fraction": float(np.mean(matched_active_mask)) if np.isfinite(matched_threshold) else float("nan"),
                "longest_active_streak_min": longest_true_streak_min(time_min_post, matched_active_mask) if np.isfinite(matched_threshold) else float("nan"),
                "onset_time_after_motor_min": first_sustained_time(time_min_post - float(row["motor_onset_min"]), vz_abs_post, threshold),
                "net_z_displacement": net_z_displacement,
                "abs_net_z_displacement": abs_net_z_displacement,
                "mass_normalized_auc": ratio_or_nan(auc_vz_abs, mass_onset),
                "mass_retention": ratio_or_nan(mass_final, mass_onset),
                "mean_swirl_rate_abs": mean_swirl,
                "swirl_penalty": ratio_or_nan(mean_swirl, mean_vz_abs),
                "seed_zone_mass_fraction_motor_onset": seed_onset,
                "seed_zone_mass_fraction_final": seed_final,
                "seed_zone_retention": ratio_or_nan(seed_final, seed_onset),
                "seed_zone_mass_fraction_mean": nanmean_or_nan(seed_zone_post),
                "chewer_band_mass_fraction_motor_onset": chewer_onset,
                "chewer_band_mass_fraction_final": chewer_final,
                "chewer_band_mass_fraction_gain": chewer_final - chewer_onset if np.isfinite(chewer_final) and np.isfinite(chewer_onset) else float("nan"),
                "chewer_band_mass_fraction_mean": nanmean_or_nan(chewer_band_post),
                "bottom_half_mass_fraction_motor_onset": bottom_half_onset,
                "bottom_half_mass_fraction_final": bottom_half_final,
                "bottom_half_mass_fraction_gain": bottom_half_final - bottom_half_onset if np.isfinite(bottom_half_final) and np.isfinite(bottom_half_onset) else float("nan"),
                "bottom_half_mass_fraction_mean": nanmean_or_nan(bottom_half_post),
                "axial_spread_entropy_motor_onset": entropy_onset,
                "axial_spread_entropy_final": entropy_final,
                "axial_spread_entropy_gain": entropy_final - entropy_onset if np.isfinite(entropy_final) and np.isfinite(entropy_onset) else float("nan"),
                "axial_spread_entropy_mean": nanmean_or_nan(axial_entropy_post),
                "speckle_length_proxy_motor_onset_um": mass_onset,
                "speckle_length_proxy_final_um": mass_final,
                "top_bias_final": float(row["top_bias_final"]),
                "nematic_xy_final": float(row["nematic_xy_final"]),
                "zmin": float(row["zmin"]),
                "zmax": float(row["zmax"]),
                "scenario_z_low": float(row["scenario_z_low"]),
                "scenario_z_high": float(row["scenario_z_high"]),
                "chewer_zone_bottom": float(row["chewer_zone_bottom"]),
                "chewer_zone_top": float(row["chewer_zone_top"]),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_timecourses_long(path: Path, run_rows: list[dict]) -> None:
    fieldnames = [
        "scenario",
        "condition",
        "run_dir",
        "time_min",
        "phase_name",
        "z_com",
        "theta_com",
        "vz_abs",
        "vz_signed",
        "Dz",
        "top_bias",
        "nematic_xy",
        "shell_occ",
        "swirl_rate_abs",
        "speckle_length_proxy_um",
        "seed_zone_mass_fraction",
        "chewer_band_mass_fraction",
        "bottom_half_mass_fraction",
        "axial_spread_entropy",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in run_rows:
            if row.get("status") != "ok":
                continue
            for idx, tval in enumerate(row["time_min"]):
                writer.writerow(
                    {
                        "scenario": row["scenario"],
                        "condition": row["condition"],
                        "run_dir": row["run_dir"],
                        "time_min": float(tval),
                        "phase_name": row["phase_name"][idx],
                        "z_com": float(row["z_com"][idx]),
                        "theta_com": float(row["theta_com"][idx]),
                        "vz_abs": float(row["vz_abs"][idx]),
                        "vz_signed": float(row["vz_signed"][idx]),
                        "Dz": float(row["Dz"][idx]),
                        "top_bias": float(row["top_bias"][idx]),
                        "nematic_xy": float(row["nematic_xy"][idx]),
                        "shell_occ": float(row["shell_occ"][idx]),
                        "swirl_rate_abs": float(row["swirl_rate_abs"][idx]),
                        "speckle_length_proxy_um": float(row["speckle_length_proxy_um"][idx]),
                        "seed_zone_mass_fraction": float(row["seed_zone_mass_fraction"][idx]),
                        "chewer_band_mass_fraction": float(row["chewer_band_mass_fraction"][idx]),
                        "bottom_half_mass_fraction": float(row["bottom_half_mass_fraction"][idx]),
                        "axial_spread_entropy": float(row["axial_spread_entropy"][idx]),
                    }
                )


def matrix_limits(values: np.ndarray, clip_zero: bool) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return (0.0, 1.0) if clip_zero else (-1.0, 1.0)
    lo = float(np.nanmin(finite))
    hi = float(np.nanmax(finite))
    if math.isclose(lo, hi):
        pad = 0.10 * max(abs(hi), 1.0)
    else:
        pad = 0.08 * (hi - lo)
    lo -= pad
    hi += pad
    if clip_zero:
        lo = max(0.0, lo)
    return lo, hi


def build_metric_matrix(metric_rows: list[dict], scenarios: list[str], conditions: list[str], field: str) -> np.ndarray:
    lookup = {(row["scenario"], row["condition"]): row for row in metric_rows}
    arr = np.full((len(scenarios), len(conditions)), np.nan, dtype=float)
    for i, scenario in enumerate(scenarios):
        for j, condition in enumerate(conditions):
            row = lookup.get((scenario, condition))
            if row is not None:
                arr[i, j] = float(row[field])
    return arr


def metric_spec_lookup() -> dict[str, dict]:
    return {spec["field"]: spec for spec in SUMMARY_METRICS}


def plot_metric_matrix(metric_rows: list[dict], scenarios: list[str], conditions: list[str], scenario_labels: dict[str, str], condition_labels: dict[str, str], metric_spec: dict, outpath: Path) -> Path:
    arr = build_metric_matrix(metric_rows, scenarios, conditions, metric_spec["field"])
    vmin, vmax = matrix_limits(arr, metric_spec["clip_zero"])
    fig, ax = plt.subplots(figsize=(1.8 * len(conditions) + 2.6, 0.95 * len(scenarios) + 2.0))
    image = ax.imshow(arr, cmap=metric_spec["cmap"], aspect="auto", vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(len(conditions)))
    ax.set_xticklabels([condition_labels[c] for c in conditions], rotation=35, ha="right", fontsize=11)
    ax.set_yticks(np.arange(len(scenarios)))
    ax.set_yticklabels([scenario_labels[s] for s in scenarios], fontsize=11)
    ax.set_title(metric_spec["label"], fontsize=16, pad=10)
    ax.set_xlabel("Condition", fontsize=13)
    ax.set_ylabel("Scenario", fontsize=13)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if not np.isfinite(arr[i, j]):
                text = "n/a"
            else:
                text = metric_spec["fmt"].format(arr[i, j])
            color = "white" if np.isfinite(arr[i, j]) and arr[i, j] > (vmin + vmax) / 2 else "black"
            ax.text(j, i, text, ha="center", va="center", fontsize=9.5, color=color)
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    cbar.ax.tick_params(labelsize=10)
    fig.subplots_adjust(left=0.24, right=0.96, top=0.90, bottom=0.20)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outpath


def plot_metric_overview(metric_rows: list[dict], scenarios: list[str], conditions: list[str], scenario_labels: dict[str, str], condition_labels: dict[str, str], outpath: Path) -> Path:
    spec_by_field = metric_spec_lookup()
    specs = [spec_by_field[field] for field in OVERVIEW_METRICS]
    fig, axes = plt.subplots(2, 3, figsize=(22, 12))
    axes = axes.ravel()
    for ax, spec in zip(axes, specs):
        arr = build_metric_matrix(metric_rows, scenarios, conditions, spec["field"])
        vmin, vmax = matrix_limits(arr, spec["clip_zero"])
        image = ax.imshow(arr, cmap=spec["cmap"], aspect="auto", vmin=vmin, vmax=vmax)
        ax.set_xticks(np.arange(len(conditions)))
        ax.set_xticklabels([condition_labels[c] for c in conditions], rotation=35, ha="right", fontsize=9)
        ax.set_yticks(np.arange(len(scenarios)))
        ax.set_yticklabels([scenario_labels[s] for s in scenarios], fontsize=10)
        ax.set_title(spec["label"], fontsize=14, pad=8)
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                if not np.isfinite(arr[i, j]):
                    text = "n/a"
                else:
                    text = spec["fmt"].format(arr[i, j])
                color = "white" if np.isfinite(arr[i, j]) and arr[i, j] > (vmin + vmax) / 2 else "black"
                ax.text(j, i, text, ha="center", va="center", fontsize=7.5, color=color)
        cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
        cbar.ax.tick_params(labelsize=8)
    fig.suptitle("Turnover campaign overview metrics (post-motor window)", fontsize=18, y=0.99)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.93, bottom=0.16, wspace=0.28, hspace=0.30)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outpath


def style(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.tick_params(direction="out", width=1.1, labelsize=11)
    ax.grid(axis="y", color="#dddddd", linewidth=0.7, alpha=0.9)


def add_phase_lines(ax, growth_end_min: float, motor_onset_min: float) -> None:
    ax.axvline(growth_end_min, color="#666666", linestyle="--", linewidth=1.1)
    ax.axvline(motor_onset_min, color="#111111", linestyle=":", linewidth=1.3)


def plot_timecourses_by_scenario(run_rows: list[dict], scenarios: list[str], scenario_labels: dict[str, str], conditions: list[str], condition_labels: dict[str, str], outdir: Path) -> list[Path]:
    outpaths = []
    lookup = {(row["scenario"], row["condition"]): row for row in run_rows if row.get("status") == "ok"}
    growth_end = float(next(row["growth_end_min"] for row in run_rows if row.get("status") == "ok"))
    motor_onset = float(next(row["motor_onset_min"] for row in run_rows if row.get("status") == "ok"))
    for scenario in scenarios:
        fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
        axes = axes.ravel()
        for ax, (field, ylabel, clip_zero) in zip(axes, TIMECOURSE_FIELDS):
            ymin = float("inf")
            ymax = float("-inf")
            for condition in conditions:
                row = lookup.get((scenario, condition))
                if row is None:
                    continue
                x = np.asarray(row["time_min"], dtype=float)
                y = np.asarray(row[field], dtype=float)
                ax.plot(x, y, color=CONDITION_COLORS.get(condition, "#333333"), lw=1.9, label=condition_labels[condition])
                finite = y[np.isfinite(y)]
                if finite.size:
                    ymin = min(ymin, float(np.nanmin(finite)))
                    ymax = max(ymax, float(np.nanmax(finite)))
            add_phase_lines(ax, growth_end, motor_onset)
            if np.isfinite(ymin) and np.isfinite(ymax):
                if clip_zero:
                    ymin = max(0.0, ymin - 0.08 * (ymax - ymin if ymax > ymin else max(abs(ymax), 1.0)))
                else:
                    ymin -= 0.08 * (ymax - ymin if ymax > ymin else max(abs(ymax), 1.0))
                ymax += 0.08 * (ymax - ymin if ymax > ymin else max(abs(ymax), 1.0))
                ax.set_ylim(ymin, ymax)
            ax.set_ylabel(ylabel, fontsize=12)
            style(ax)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=11)
        for ax in axes[-2:]:
            ax.set_xlabel("Time (min)", fontsize=12)
        fig.suptitle(f"{scenario_labels[scenario]}: full-timeline turnover traces", fontsize=17, y=0.98)
        fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.14, wspace=0.24, hspace=0.26)
        outpath = outdir / f"{scenario}_timecourses.png"
        outpath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outpath, dpi=300, bbox_inches="tight")
        plt.close(fig)
        outpaths.append(outpath)
    return outpaths


def plot_kymograph_grid(run_rows: list[dict], scenarios: list[str], scenario_labels: dict[str, str], conditions: list[str], condition_labels: dict[str, str], outpath: Path) -> Path:
    lookup = {(row["scenario"], row["condition"]): row for row in run_rows if row.get("status") == "ok"}
    values = []
    for row in lookup.values():
        values.extend(np.asarray(row["rho_z"], dtype=float).ravel().tolist())
    vmax = float(np.nanpercentile(np.asarray(values, dtype=float), 99.5)) if values else 2.0
    vmax = max(vmax, 1.0)
    fig, axes = plt.subplots(len(scenarios), len(conditions), figsize=(3.4 * len(conditions) + 1.3, 2.2 * len(scenarios) + 1.3), sharex=True, sharey=True)
    if len(scenarios) == 1:
        axes = np.array([axes])
    if len(conditions) == 1:
        axes = axes[:, None]
    image = None
    growth_end = float(next(row["growth_end_min"] for row in run_rows if row.get("status") == "ok"))
    motor_onset = float(next(row["motor_onset_min"] for row in run_rows if row.get("status") == "ok"))
    for i, scenario in enumerate(scenarios):
        for j, condition in enumerate(conditions):
            ax = axes[i, j]
            row = lookup.get((scenario, condition))
            if row is None:
                ax.set_facecolor("#f3f3f3")
                ax.set_xticks([])
                ax.set_yticks([])
                ax.text(0.5, 0.5, "n/a", ha="center", va="center", transform=ax.transAxes, fontsize=12, color="#777777")
                continue
            time = np.asarray(row["time_min"], dtype=float)
            z_edges = np.asarray(row["z_edges"], dtype=float)
            rho = np.asarray(row["rho_z"], dtype=float).T
            image = ax.imshow(rho, aspect="auto", origin="upper", extent=[time[0], time[-1], z_edges[-1], z_edges[0]], cmap="viridis", vmin=0.0, vmax=vmax)
            add_phase_lines(ax, growth_end, motor_onset)
            if i == 0:
                ax.set_title(condition_labels[condition], fontsize=12, pad=8)
            if j == 0:
                ax.set_ylabel(f"{scenario_labels[scenario]}\nz (um)", fontsize=11)
            if i == len(scenarios) - 1:
                ax.set_xlabel("Time (min)", fontsize=11)
            ax.tick_params(direction="out", width=1.0, labelsize=10)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_linewidth(1.0)
            ax.spines["bottom"].set_linewidth(1.0)
    cbar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.018, pad=0.02)
    cbar.set_label("Normalized axial density", fontsize=13)
    cbar.ax.tick_params(labelsize=10)
    fig.suptitle("Full-timeline axial-density kymographs", fontsize=18, y=0.995)
    fig.subplots_adjust(left=0.08, right=0.93, top=0.91, bottom=0.08, wspace=0.08, hspace=0.10)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outpath


def snapshot_limits(run_rows: list[dict]) -> float:
    values = []
    for row in run_rows:
        if row.get("status") != "ok":
            continue
        for slug, _target, _label in SNAPSHOT_TIMES_MIN:
            values.extend(np.asarray(row["snapshot_densities"][slug], dtype=float).ravel().tolist())
    vmax = float(np.nanpercentile(np.asarray(values, dtype=float), 99.5)) if values else 2.0
    return min(max(vmax, 1.2), 3.5)


def plot_snapshot_heatmaps(run_rows: list[dict], scenarios: list[str], scenario_labels: dict[str, str], conditions: list[str], condition_labels: dict[str, str], outdir: Path) -> list[Path]:
    lookup = {(row["scenario"], row["condition"]): row for row in run_rows if row.get("status") == "ok"}
    vmax = snapshot_limits(run_rows)
    outpaths = []
    for slug, _target, label in SNAPSHOT_TIMES_MIN:
        fig, axes = plt.subplots(len(scenarios), len(conditions), figsize=(3.4 * len(conditions) + 1.3, 2.2 * len(scenarios) + 1.3), sharex=True, sharey=True)
        if len(scenarios) == 1:
            axes = np.array([axes])
        if len(conditions) == 1:
            axes = axes[:, None]
        image = None
        for i, scenario in enumerate(scenarios):
            for j, condition in enumerate(conditions):
                ax = axes[i, j]
                row = lookup.get((scenario, condition))
                if row is None:
                    ax.set_facecolor("#f3f3f3")
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.text(0.5, 0.5, "n/a", ha="center", va="center", transform=ax.transAxes, fontsize=12, color="#777777")
                    continue
                density = np.asarray(row["snapshot_densities"][slug], dtype=float)
                smax = float(row[f"{slug}_smax"])
                image = ax.imshow(density, origin="lower", aspect="auto", extent=[0.0, smax, float(row["zmin"]), float(row["zmax"])], cmap="viridis", vmin=0.0, vmax=vmax, interpolation="nearest")
                ax.set_ylim(float(row["zmax"]), float(row["zmin"]))
                if i == 0:
                    ax.set_title(condition_labels[condition], fontsize=12, pad=8)
                if j == 0:
                    ax.set_ylabel(f"{scenario_labels[scenario]}\nz (um)", fontsize=11)
                if i == len(scenarios) - 1:
                    ax.set_xlabel("Unwrapped circumference (um)", fontsize=11)
                ax.tick_params(direction="out", width=1.0, labelsize=9)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
        cbar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.018, pad=0.02)
        cbar.set_label("Normalized areal density", fontsize=13)
        cbar.ax.tick_params(labelsize=10)
        fig.suptitle(f"Unwrapped annulus heatmaps at {label}", fontsize=18, y=0.995)
        fig.subplots_adjust(left=0.08, right=0.93, top=0.91, bottom=0.08, wspace=0.08, hspace=0.10)
        outpath = outdir / f"{slug}_heatmaps.png"
        outpath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outpath, dpi=300, bbox_inches="tight")
        plt.close(fig)
        outpaths.append(outpath)
    return outpaths


def plot_snapshot_lines(run_rows: list[dict], scenarios: list[str], scenario_labels: dict[str, str], conditions: list[str], condition_labels: dict[str, str], outdir: Path) -> list[Path]:
    lookup = {(row["scenario"], row["condition"]): row for row in run_rows if row.get("status") == "ok"}
    outpaths = []
    for slug, _target, label in SNAPSHOT_TIMES_MIN:
        fig, axes = plt.subplots(len(scenarios), len(conditions), figsize=(3.4 * len(conditions) + 1.3, 2.2 * len(scenarios) + 1.3), sharex=True, sharey=True)
        if len(scenarios) == 1:
            axes = np.array([axes])
        if len(conditions) == 1:
            axes = axes[:, None]
        for i, scenario in enumerate(scenarios):
            for j, condition in enumerate(conditions):
                ax = axes[i, j]
                row = lookup.get((scenario, condition))
                if row is None:
                    ax.set_facecolor("#f3f3f3")
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.text(0.5, 0.5, "n/a", ha="center", va="center", transform=ax.transAxes, fontsize=12, color="#777777")
                    continue
                frame = np.asarray(row["snapshot_frames"][slug], dtype=float)
                smax = float(row[f"{slug}_smax"])
                r_mid = smax / (2.0 * math.pi)
                seam = float(row["theta_seam"])
                for s, z in frame_segments(frame, seam, smax, r_mid):
                    ax.plot(s, z, color="black", lw=0.55, alpha=0.40)
                ax.set_xlim(0.0, smax)
                ax.set_ylim(float(row["zmax"]), float(row["zmin"]))
                if i == 0:
                    ax.set_title(condition_labels[condition], fontsize=12, pad=8)
                if j == 0:
                    ax.set_ylabel(f"{scenario_labels[scenario]}\nz (um)", fontsize=11)
                if i == len(scenarios) - 1:
                    ax.set_xlabel("Unwrapped circumference (um)", fontsize=11)
                ax.tick_params(direction="out", width=1.0, labelsize=9)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
        fig.suptitle(f"Unwrapped annulus line snapshots at {label}", fontsize=18, y=0.995)
        fig.subplots_adjust(left=0.08, right=0.98, top=0.91, bottom=0.08, wspace=0.08, hspace=0.10)
        outpath = outdir / f"{slug}_lines.png"
        outpath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outpath, dpi=300, bbox_inches="tight")
        plt.close(fig)
        outpaths.append(outpath)
    return outpaths


def write_readme(path: Path, campaign_name: str, run_rows: list[dict], metric_rows: list[dict], thresholds: dict[int, float]) -> None:
    growth_end = float(next(row["growth_end_min"] for row in run_rows if row.get("status") == "ok"))
    motor_onset = float(next(row["motor_onset_min"] for row in run_rows if row.get("status") == "ok"))
    ranked_speed = sorted(metric_rows, key=lambda row: row["mean_vz_abs"], reverse=True)[:5]
    ranked_auc = sorted(metric_rows, key=lambda row: row["auc_vz_abs"], reverse=True)[:5]
    ranked_excess_auc = sorted(
        [row for row in metric_rows if math.isfinite(row.get("excess_auc_vz_abs", float("nan")))],
        key=lambda row: row["excess_auc_vz_abs"],
        reverse=True,
    )[:5]
    ranked_chewer_band_low = sorted(
        [row for row in metric_rows if math.isfinite(row.get("chewer_band_mass_fraction_final", float("nan")))],
        key=lambda row: row["chewer_band_mass_fraction_final"],
    )[:5]
    ranked_bottom_half_low = sorted(
        [row for row in metric_rows if math.isfinite(row.get("bottom_half_mass_fraction_final", float("nan")))],
        key=lambda row: row["bottom_half_mass_fraction_final"],
    )[:5]
    lines = [
        f"# {campaign_name} analysis bundle",
        "",
        "## Phase markers",
        f"- Growth-only ends at `{growth_end:.3f} min`.",
        f"- Motor observation starts at `{motor_onset:.3f} min`.",
        "- Movement metrics are computed only on the post-motor window, so the growth and crosslink-conditioning phases do not dilute the transport readout.",
        "",
        "## Active-transport thresholds",
        f"- No-xlink control threshold: `{thresholds.get(0, float('nan')):.6f} um/s`",
        f"- Xlink control threshold: `{thresholds.get(1, float('nan')):.6f} um/s`",
        "- Scenario-matched thresholds are written to `scenario_matched_active_thresholds.csv`; `excess_auc_vz_abs` and `matched_active_transport_fraction` use those fair baselines.",
        "",
        "## Best conditions by mean |v_z|",
    ]
    for row in ranked_speed:
        lines.append(
            f"- `{row['scenario']}` + `{row['condition']}`: mean |v_z| = `{row['mean_vz_abs']:.6f} um/s`, AUC = `{row['auc_vz_abs']:.3f} um`, directionality = `{row['directionality_index']:.3f}`"
        )
    lines.extend(["", "## Best conditions by AUC(|v_z|)"])
    for row in ranked_auc:
        lines.append(
            f"- `{row['scenario']}` + `{row['condition']}`: AUC = `{row['auc_vz_abs']:.3f} um`, |Net z shift| = `{row['abs_net_z_displacement']:.3f} um`, final filament mass proxy = `{row['speckle_length_proxy_final_um']:.1f} um`"
        )
    lines.extend(["", "## Best conditions by excess AUC above matched no-motor control"])
    for row in ranked_excess_auc:
        lines.append(
            f"- `{row['scenario']}` + `{row['condition']}`: excess AUC = `{row['excess_auc_vz_abs']:.3f} um`, matched-active fraction = `{row['matched_active_transport_fraction']:.3f}`, longest active streak = `{row['longest_active_streak_min']:.3f} min`"
        )
    if ranked_chewer_band_low:
        lines.extend(["", "## Lowest final chewer-band mass fraction"])
        for row in ranked_chewer_band_low:
            lines.append(
                f"- `{row['scenario']}` + `{row['condition']}`: final chewer-band fraction = `{row['chewer_band_mass_fraction_final']:.3f}`, gain from motor onset = `{row['chewer_band_mass_fraction_gain']:.3f}`"
            )
    if ranked_bottom_half_low:
        lines.extend(["", "## Lowest final common bottom-half mass fraction"])
        for row in ranked_bottom_half_low:
            lines.append(
                f"- `{row['scenario']}` + `{row['condition']}`: final bottom-half fraction = `{row['bottom_half_mass_fraction_final']:.3f}`, gain from motor onset = `{row['bottom_half_mass_fraction_gain']:.3f}`"
            )
    lines.extend(
        [
            "",
            "## Files",
            "- `movement_matrices/`: 5x6 scenario-condition metric maps.",
            "- `scenario_matched_active_thresholds.csv`: scenario- and xlink-matched no-motor movement thresholds.",
            "- `timecourses/by_scenario/`: full-timeline metric traces with phase markers.",
            "- `kymographs/full_timeline_kymograph_grid.png`: axial-density evolution for all 30 conditions.",
            "- `unwrapped_annulus/heatmaps/` and `unwrapped_annulus/lines/`: snapshot grids at 0, 100, 120, 360, and 720 s.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
    plt.rcParams.update({"font.size": 11, "axes.linewidth": 1.1, "savefig.facecolor": "white", "figure.facecolor": "white"})

    campaign_root = Path(args.campaign_root).expanduser().resolve()
    save_dir = Path(args.save_dir).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    report_bin = pick_report_bin(args.report_bin)

    manifest_rows = load_manifest(campaign_root, save_dir)
    tasks = [(row, str(report_bin), args.interval, args.speckle_name) for row in manifest_rows]
    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            run_rows = list(executor.map(analyze_run, tasks))
    else:
        run_rows = [analyze_run(task) for task in tasks]

    thresholds = compute_active_thresholds(run_rows)
    matched_thresholds = compute_scenario_matched_thresholds(run_rows)
    metric_rows = compute_run_metrics(run_rows, thresholds, matched_thresholds)

    scenario_order = [row["scenario"] for row in manifest_rows if row["condition"] == manifest_rows[0]["condition"]]
    seen = set()
    scenario_order = [s for s in scenario_order if not (s in seen or seen.add(s))]
    condition_order = [row["condition"] for row in manifest_rows if row["scenario"] == manifest_rows[0]["scenario"]]
    seen = set()
    condition_order = [c for c in condition_order if not (c in seen or seen.add(c))]
    scenario_labels = {row["scenario"]: row["scenario_short"] for row in manifest_rows}
    condition_labels = {row["condition"]: row["condition_short"] for row in manifest_rows}

    write_csv(outdir / "run_metrics.csv", metric_rows)
    write_csv(outdir / "active_thresholds.csv", [{"include_xlinks": key, "threshold_vz_abs": value} for key, value in sorted(thresholds.items())])
    write_csv(
        outdir / "scenario_matched_active_thresholds.csv",
        [
            {"scenario": scenario, "include_xlinks": include_xlinks, "threshold_vz_abs_matched": value}
            for (scenario, include_xlinks), value in sorted(matched_thresholds.items())
        ],
    )
    write_timecourses_long(outdir / "run_timecourses_long.csv", run_rows)

    metric_paths = []
    for spec in SUMMARY_METRICS:
        metric_paths.append(
            plot_metric_matrix(
                metric_rows,
                scenario_order,
                condition_order,
                scenario_labels,
                condition_labels,
                spec,
                outdir / "movement_matrices" / f"{spec['slug']}.png",
            )
        )
    overview_path = plot_metric_overview(metric_rows, scenario_order, condition_order, scenario_labels, condition_labels, outdir / "movement_matrices" / "movement_metric_overview.png")
    timecourse_paths = plot_timecourses_by_scenario(run_rows, scenario_order, scenario_labels, condition_order, condition_labels, outdir / "timecourses" / "by_scenario")
    kymograph_path = plot_kymograph_grid(run_rows, scenario_order, scenario_labels, condition_order, condition_labels, outdir / "kymographs" / "full_timeline_kymograph_grid.png")
    heatmap_paths = plot_snapshot_heatmaps(run_rows, scenario_order, scenario_labels, condition_order, condition_labels, outdir / "unwrapped_annulus" / "heatmaps")
    line_paths = plot_snapshot_lines(run_rows, scenario_order, scenario_labels, condition_order, condition_labels, outdir / "unwrapped_annulus" / "lines")
    write_readme(outdir / "README.md", campaign_root.name, run_rows, metric_rows, thresholds)

    print(outdir / "README.md")
    print(outdir / "run_metrics.csv")
    print(outdir / "active_thresholds.csv")
    print(outdir / "scenario_matched_active_thresholds.csv")
    print(outdir / "run_timecourses_long.csv")
    print(overview_path)
    print(kymograph_path)
    for path in metric_paths + timecourse_paths + heatmap_paths + line_paths:
        print(path)


if __name__ == "__main__":
    main()
