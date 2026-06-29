#!/usr/bin/env python3
"""Prototype movie-level actin motion estimates for egg-cell IMS files.

This is a first-pass QC tool, not a final quantitative pipeline. It reads
Imaris `.ims` files, max-projects selected z slices, estimates local
frame-to-frame shifts with block cross-correlation, and writes overlays plus
CSV summaries. Timing remains provisional until the collaborator frame-interval
conflict is resolved, so `um_per_frame` is the primary unit.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import correlate2d


DEFAULT_FILES = [
    "Kwaku and Berry/Wild type/1_5ADVMLE.ims",
    "Kwaku and Berry/CA ROP9/Image0004_5ADVMLE.oib.ims",
]


@dataclass(frozen=True)
class MovieMeta:
    relative_path: str
    line: str
    voxel_um_x: float
    voxel_um_y: float
    median_dt_s: float | None
    frame_count: int
    height: int
    width: int
    z_count: int


def decode_value(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    dtype = getattr(value, "dtype", None)
    if dtype is not None and getattr(dtype, "kind", None) == "S":
        return value.tobytes().decode("utf-8", errors="replace").rstrip("\x00")
    if hasattr(value, "shape") and getattr(value, "shape", None) == ():
        return decode_value(value.item())
    if hasattr(value, "tolist"):
        return decode_value(value.tolist())
    return value


def to_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_timepoint(value: Any) -> datetime | None:
    value = decode_value(value)
    if not isinstance(value, str):
        return None
    for fmt in ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S"):
        try:
            return datetime.strptime(value, fmt)
        except ValueError:
            pass
    return None


def frame_intervals_s(handle: h5py.File) -> list[float]:
    time_info = handle.get("DataSetInfo/TimeInfo")
    if time_info is None:
        return []
    points: list[tuple[int, datetime]] = []
    for key, value in time_info.attrs.items():
        if not key.startswith("TimePoint"):
            continue
        suffix = key.replace("TimePoint", "")
        if not suffix.isdigit():
            continue
        parsed = parse_timepoint(value)
        if parsed is not None:
            points.append((int(suffix), parsed))
    points.sort(key=lambda item: item[0])
    return [
        (points[idx + 1][1] - points[idx][1]).total_seconds()
        for idx in range(len(points) - 1)
    ]


def image_voxels_um(handle: h5py.File, shape_zyx: tuple[int, int, int]) -> tuple[float, float]:
    image = handle.get("DataSetInfo/Image")
    if image is None:
        return (1.0, 1.0)
    attrs = {key: decode_value(value) for key, value in image.attrs.items()}
    sx = int(attrs.get("ImageSizeX") or shape_zyx[2])
    sy = int(attrs.get("ImageSizeY") or shape_zyx[1])
    xmin = to_float(attrs.get("ExtMin0"))
    xmax = to_float(attrs.get("ExtMax0"))
    ymin = to_float(attrs.get("ExtMin1"))
    ymax = to_float(attrs.get("ExtMax1"))
    vx = ((xmax - xmin) / sx) if xmin is not None and xmax is not None and sx else 1.0
    vy = ((ymax - ymin) / sy) if ymin is not None and ymax is not None and sy else 1.0
    return (float(vx), float(vy))


def timepoint_names(handle: h5py.File) -> list[str]:
    root = handle["DataSet/ResolutionLevel 0"]
    names = [name for name in root.keys() if name.startswith("TimePoint ")]
    return sorted(names, key=lambda name: int(name.split()[-1]))


def projection_stack(path: Path, z_start: int, z_end: int) -> tuple[np.ndarray, MovieMeta]:
    with h5py.File(path, "r") as handle:
        names = timepoint_names(handle)
        if not names:
            raise ValueError(f"No time points found in {path}")
        first = handle[f"DataSet/ResolutionLevel 0/{names[0]}/Channel 0/Data"]
        z_count, height, width = first.shape
        lo = max(0, z_start - 1)
        hi = min(z_count, z_end)
        if lo >= hi:
            raise ValueError(f"Invalid z slice range {z_start}-{z_end} for {z_count} slices")
        frames: list[np.ndarray] = []
        for name in names:
            data = handle[f"DataSet/ResolutionLevel 0/{name}/Channel 0/Data"]
            frames.append(np.asarray(data[lo:hi], dtype=np.float32).max(axis=0))
        voxel_x, voxel_y = image_voxels_um(handle, first.shape)
        intervals = frame_intervals_s(handle)
        median_dt = float(np.median(intervals)) if intervals else None
    meta = MovieMeta(
        relative_path=path.as_posix(),
        line=path.parent.name,
        voxel_um_x=voxel_x,
        voxel_um_y=voxel_y,
        median_dt_s=median_dt,
        frame_count=len(frames),
        height=height,
        width=width,
        z_count=z_count,
    )
    return np.stack(frames), meta


def normalize_frame(frame: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(frame, [1, 99.5])
    if hi <= lo:
        return np.zeros_like(frame, dtype=np.float32)
    out = (frame - lo) / (hi - lo)
    return np.clip(out, 0, 1).astype(np.float32)


def estimate_pair_flow(
    frame_a: np.ndarray,
    frame_b: np.ndarray,
    window: int,
    step: int,
    search: int,
    min_std: float,
    min_mean: float,
) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    half = window // 2
    for y in range(search + half, frame_a.shape[0] - search - half + 1, step):
        for x in range(search + half, frame_a.shape[1] - search - half + 1, step):
            template = frame_a[y - half : y + half, x - half : x + half]
            if template.mean() < min_mean or template.std() < min_std:
                continue
            target = frame_b[
                y - half - search : y + half + search + 1,
                x - half - search : x + half + search + 1,
            ]
            if target.shape[0] < template.shape[0] or target.shape[1] < template.shape[1]:
                continue
            template = template - template.mean()
            target = target - target.mean()
            corr = correlate2d(target, template, mode="valid")
            peak_y, peak_x = np.unravel_index(int(np.argmax(corr)), corr.shape)
            dy = float(peak_y - search)
            dx = float(peak_x - search)
            peak = float(corr[peak_y, peak_x])
            rows.append({"x_pix": float(x), "y_pix": float(y), "dx_pix": dx, "dy_pix": dy, "peak": peak})
    return rows


def aggregate_flow(
    frames: np.ndarray,
    window: int,
    step: int,
    search: int,
    min_std: float,
    min_mean: float,
) -> list[dict[str, float]]:
    norm = np.stack([normalize_frame(frame) for frame in frames])
    keyed: dict[tuple[float, float], list[dict[str, float]]] = {}
    for idx in range(norm.shape[0] - 1):
        for row in estimate_pair_flow(norm[idx], norm[idx + 1], window, step, search, min_std, min_mean):
            keyed.setdefault((row["x_pix"], row["y_pix"]), []).append(row)
    out: list[dict[str, float]] = []
    for (x, y), rows in keyed.items():
        dx = np.array([row["dx_pix"] for row in rows], dtype=float)
        dy = np.array([row["dy_pix"] for row in rows], dtype=float)
        out.append(
            {
                "x_pix": x,
                "y_pix": y,
                "dx_pix_per_frame": float(np.median(dx)),
                "dy_pix_per_frame": float(np.median(dy)),
                "n_pairs": float(len(rows)),
                "valid_fraction": float(len(rows) / max(1, frames.shape[0] - 1)),
            }
        )
    return out


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def save_projection_preview(path: Path, frames: np.ndarray, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    indices = sorted({0, frames.shape[0] // 2, frames.shape[0] - 1})
    fig, axes = plt.subplots(1, len(indices), figsize=(4 * len(indices), 4), constrained_layout=True)
    if len(indices) == 1:
        axes = [axes]
    for ax, idx in zip(axes, indices):
        ax.imshow(normalize_frame(frames[idx]), cmap="gray")
        ax.set_title(f"frame {idx}")
        ax.axis("off")
    fig.suptitle(title)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_flow_overlay(path: Path, frames: np.ndarray, rows: list[dict[str, Any]], title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    base = normalize_frame(np.median(frames, axis=0))
    fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
    ax.imshow(base, cmap="gray")
    if rows:
        x = np.array([row["x_pix"] for row in rows])
        y = np.array([row["y_pix"] for row in rows])
        dx = np.array([row["dx_pix_per_frame"] for row in rows])
        dy = np.array([row["dy_pix_per_frame"] for row in rows])
        speed = np.hypot(dx, dy)
        keep = speed > 0
        ax.quiver(
            x[keep],
            y[keep],
            dx[keep],
            dy[keep],
            speed[keep],
            angles="xy",
            scale_units="xy",
            scale=0.35,
            cmap="magma",
            width=0.004,
        )
    ax.set_title(title)
    ax.axis("off")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def summarize(meta: MovieMeta, rows: list[dict[str, Any]], source: Path, z_start: int, z_end: int) -> dict[str, Any]:
    dx = np.array([row["dx_pix_per_frame"] for row in rows], dtype=float)
    dy = np.array([row["dy_pix_per_frame"] for row in rows], dtype=float)
    speed_um_frame = np.hypot(dx * meta.voxel_um_x, dy * meta.voxel_um_y) if len(rows) else np.array([])
    med_speed = float(np.median(speed_um_frame)) if speed_um_frame.size else None
    return {
        "source_file": source.as_posix(),
        "line": meta.line,
        "frames": meta.frame_count,
        "z_slices_1based": f"{z_start}-{z_end}",
        "xy_voxel_um": f"{meta.voxel_um_x:.6g},{meta.voxel_um_y:.6g}",
        "median_metadata_dt_s": meta.median_dt_s,
        "vectors": len(rows),
        "median_speed_um_per_frame": med_speed,
        "median_speed_um_per_s_metadata_dt": med_speed / meta.median_dt_s if med_speed is not None and meta.median_dt_s else None,
        "median_speed_um_per_s_assuming_30s": med_speed / 30.0 if med_speed is not None else None,
        "median_dx_um_per_frame": float(np.median(dx * meta.voxel_um_x)) if len(rows) else None,
        "median_dy_um_per_frame": float(np.median(dy * meta.voxel_um_y)) if len(rows) else None,
    }


def slug(path: Path) -> str:
    return "_".join(path.with_suffix("").parts[-2:]).replace(" ", "_").replace("=", "")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="analysis/results/egg_cell_flow_prototype_2026-06-15")
    parser.add_argument("--z-start", type=int, default=5, help="First z slice, 1-based inclusive.")
    parser.add_argument("--z-end", type=int, default=7, help="Last z slice, 1-based inclusive.")
    parser.add_argument("--window", type=int, default=48)
    parser.add_argument("--step", type=int, default=32)
    parser.add_argument("--search", type=int, default=8)
    parser.add_argument("--min-std", type=float, default=0.025)
    parser.add_argument("--min-mean", type=float, default=0.06)
    parser.add_argument("files", nargs="*", default=DEFAULT_FILES)
    args = parser.parse_args()

    outdir = Path(args.output)
    overlays = outdir / "overlays"
    tables = outdir / "tables"
    summaries: list[dict[str, Any]] = []

    for file_arg in args.files:
        source = Path(file_arg)
        frames, meta = projection_stack(source, args.z_start, args.z_end)
        rows = aggregate_flow(frames, args.window, args.step, args.search, args.min_std, args.min_mean)
        name = slug(source)

        for row in rows:
            row["dx_um_per_frame"] = row["dx_pix_per_frame"] * meta.voxel_um_x
            row["dy_um_per_frame"] = row["dy_pix_per_frame"] * meta.voxel_um_y
            row["speed_um_per_frame"] = float(np.hypot(row["dx_um_per_frame"], row["dy_um_per_frame"]))

        write_csv(
            tables / f"{name}_vectors.csv",
            rows,
            [
                "x_pix",
                "y_pix",
                "dx_pix_per_frame",
                "dy_pix_per_frame",
                "dx_um_per_frame",
                "dy_um_per_frame",
                "speed_um_per_frame",
                "n_pairs",
                "valid_fraction",
            ],
        )
        save_projection_preview(overlays / f"{name}_projection_preview.png", frames, source.as_posix())
        save_flow_overlay(overlays / f"{name}_flow_overlay.png", frames, rows, source.as_posix())
        summaries.append(summarize(meta, rows, source, args.z_start, args.z_end))

    write_csv(
        tables / "prototype_summary.csv",
        summaries,
        [
            "source_file",
            "line",
            "frames",
            "z_slices_1based",
            "xy_voxel_um",
            "median_metadata_dt_s",
            "vectors",
            "median_speed_um_per_frame",
            "median_speed_um_per_s_metadata_dt",
            "median_speed_um_per_s_assuming_30s",
            "median_dx_um_per_frame",
            "median_dy_um_per_frame",
        ],
    )
    print(tables / "prototype_summary.csv")
    print(overlays)


if __name__ == "__main__":
    main()
