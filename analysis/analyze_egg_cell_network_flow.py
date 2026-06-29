#!/usr/bin/env python3
"""Consistent movie-level actin flow analysis for egg-cell movies.

This script is the more reproducible alternative to hand-picked kymograph
slopes. It estimates local frame-to-frame actin motion across the movie using
masked block cross-correlation, with optional global drift correction.

Primary output units are um/frame. Use um/s only when the frame interval has
been confirmed or explicitly supplied with --assumed-dt-s.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
from scipy.signal import correlate2d

try:
    from PIL import Image, ImageSequence
except ImportError:  # pragma: no cover
    Image = None
    ImageSequence = None

try:
    import h5py
except ImportError:  # pragma: no cover
    h5py = None


@dataclass(frozen=True)
class MovieMeta:
    source_file: str
    line: str
    reader: str
    frames: int
    height: int
    width: int
    z_count: int | None
    z_slices_1based: str
    voxel_um_x: float
    voxel_um_y: float
    metadata_dt_s: float | None
    calibration_note: str


def decode_value(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    dtype = getattr(value, "dtype", None)
    if dtype is not None and getattr(dtype, "kind", None) == "S":
        return value.tobytes().decode("utf-8", errors="replace").rstrip("\x00")
    if dtype is not None and getattr(dtype, "kind", None) == "U":
        return "".join(value.tolist())
    if hasattr(value, "shape") and getattr(value, "shape", None) == ():
        return decode_value(value.item())
    if hasattr(value, "tolist"):
        return decode_value(value.tolist())
    return value


def to_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(decode_value(value))
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


def timepoint_names(handle: h5py.File) -> list[str]:
    root = handle["DataSet/ResolutionLevel 0"]
    names = [name for name in root.keys() if name.startswith("TimePoint ")]
    return sorted(names, key=lambda name: int(name.split()[-1]))


def ims_voxels_um(handle: h5py.File, shape_zyx: tuple[int, int, int]) -> tuple[float, float, str]:
    image = handle.get("DataSetInfo/Image")
    if image is None:
        return (1.0, 1.0, "No DataSetInfo/Image metadata; using 1 pixel = 1 unit.")

    attrs = {key: decode_value(value) for key, value in image.attrs.items()}
    xmin = to_float(attrs.get("ExtMin0"))
    xmax = to_float(attrs.get("ExtMax0"))
    ymin = to_float(attrs.get("ExtMin1"))
    ymax = to_float(attrs.get("ExtMax1"))
    width = shape_zyx[2]
    height = shape_zyx[1]
    if xmin is None or xmax is None or ymin is None or ymax is None:
        return (1.0, 1.0, "Missing ExtMin/ExtMax metadata; using 1 pixel = 1 unit.")
    note = (
        "Direct Imaris HDF5 calibration from ExtMin/ExtMax and ResolutionLevel 0 "
        "array shape. Verify against Bio-Formats/Fiji for final absolute units."
    )
    return ((xmax - xmin) / width, (ymax - ymin) / height, note)


def read_ims_projection(path: Path, z_start: int, z_end: int) -> tuple[np.ndarray, MovieMeta]:
    if h5py is None:
        raise RuntimeError(
            "Reading .ims files requires h5py. For collaborator use, export "
            "projected TIFF stacks from Fiji/Bio-Formats and analyze those instead."
        )
    with h5py.File(path, "r") as handle:
        names = timepoint_names(handle)
        if not names:
            raise ValueError(f"No time points found in {path}")
        first = handle[f"DataSet/ResolutionLevel 0/{names[0]}/Channel 0/Data"]
        z_count, height, width = first.shape
        requested_slices = max(1, z_end - z_start + 1)
        lo = max(0, z_start - 1)
        hi = min(z_count, z_end)
        if lo >= z_count:
            lo = max(0, z_count - requested_slices)
            hi = z_count
        if lo >= hi:
            raise ValueError(f"Invalid z slice range {z_start}-{z_end} for {z_count} slices")

        frames = []
        for name in names:
            data = handle[f"DataSet/ResolutionLevel 0/{name}/Channel 0/Data"]
            frames.append(np.asarray(data[lo:hi], dtype=np.float32).max(axis=0))

        voxel_x, voxel_y, note = ims_voxels_um(handle, first.shape)
        intervals = frame_intervals_s(handle)
        metadata_dt = float(np.median(intervals)) if intervals else None

    meta = MovieMeta(
        source_file=path.as_posix(),
        line=path.parent.name,
        reader="direct_ims_hdf5",
        frames=len(frames),
        height=height,
        width=width,
        z_count=z_count,
        z_slices_1based=f"{lo + 1}-{hi}",
        voxel_um_x=float(voxel_x),
        voxel_um_y=float(voxel_y),
        metadata_dt_s=metadata_dt,
        calibration_note=note,
    )
    return np.stack(frames), meta


def read_tiff_stack(
    path: Path,
    pixel_size_um: float,
    line: str | None,
) -> tuple[np.ndarray, MovieMeta]:
    if Image is None or ImageSequence is None:
        raise RuntimeError("Pillow is required to read TIFF stacks.")
    frames = []
    with Image.open(path) as image:
        for frame in ImageSequence.Iterator(image):
            frames.append(np.asarray(frame, dtype=np.float32))
    if not frames:
        raise ValueError(f"No frames found in {path}")
    stack = np.stack(frames)
    if stack.ndim != 3:
        raise ValueError(f"Expected a 2D time stack from {path}, got shape {stack.shape}")
    meta = MovieMeta(
        source_file=path.as_posix(),
        line=line or path.parent.name,
        reader="tiff_stack",
        frames=stack.shape[0],
        height=stack.shape[1],
        width=stack.shape[2],
        z_count=None,
        z_slices_1based="already_projected",
        voxel_um_x=pixel_size_um,
        voxel_um_y=pixel_size_um,
        metadata_dt_s=None,
        calibration_note="Pixel size supplied by --pixel-size-um for exported TIFF stack.",
    )
    return stack, meta


def read_movie(
    path: Path,
    z_start: int,
    z_end: int,
    pixel_size_um: float | None,
    line: str | None,
) -> tuple[np.ndarray, MovieMeta]:
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".ims"):
        return read_ims_projection(path, z_start, z_end)
    if path.suffix.lower() in {".tif", ".tiff"}:
        if pixel_size_um is None:
            raise ValueError("--pixel-size-um is required for TIFF input")
        return read_tiff_stack(path, pixel_size_um, line)
    raise ValueError(f"Unsupported input format: {path}")


def normalize_frame(frame: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(frame, [1, 99.7])
    if hi <= lo:
        return np.zeros_like(frame, dtype=np.float32)
    out = (frame - lo) / (hi - lo)
    return np.clip(out, 0, 1).astype(np.float32)


def otsu_threshold(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0
    hist, edges = np.histogram(values, bins=256, range=(0, 1))
    centers = (edges[:-1] + edges[1:]) / 2
    weight1 = np.cumsum(hist)
    weight2 = np.cumsum(hist[::-1])[::-1]
    mean1 = np.cumsum(hist * centers) / np.maximum(weight1, 1)
    mean2 = (np.cumsum((hist * centers)[::-1]) / np.maximum(weight2[::-1], 1))[::-1]
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2
    if variance12.size == 0:
        return float(np.percentile(values, 75))
    return float(centers[:-1][np.argmax(variance12)])


def foreground_mask(norm_stack: np.ndarray, percentile_floor: float) -> np.ndarray:
    base = np.median(norm_stack, axis=0)
    threshold = max(otsu_threshold(base), float(np.percentile(base, percentile_floor)))
    mask = base > threshold
    mask = ndimage.binary_opening(mask, iterations=1)
    mask = ndimage.binary_closing(mask, iterations=2)
    mask = ndimage.binary_dilation(mask, iterations=3)
    labels, count = ndimage.label(mask)
    if count == 0:
        return np.ones_like(base, dtype=bool)
    sizes = ndimage.sum(mask, labels, index=np.arange(1, count + 1))
    keep_labels = np.where(sizes >= max(50, 0.001 * mask.size))[0] + 1
    cleaned = np.isin(labels, keep_labels)
    return cleaned.astype(bool)


def phase_correlation_shift(frame_a: np.ndarray, frame_b: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    a = (frame_a * mask).astype(np.float32)
    b = (frame_b * mask).astype(np.float32)
    a -= float(a.mean())
    b -= float(b.mean())
    fft_product = np.fft.fft2(b) * np.conj(np.fft.fft2(a))
    denom = np.abs(fft_product)
    cross_power = fft_product / np.maximum(denom, 1e-12)
    corr = np.fft.ifft2(cross_power).real
    peak_y, peak_x = np.unravel_index(int(np.argmax(corr)), corr.shape)
    if peak_y > corr.shape[0] // 2:
        peak_y -= corr.shape[0]
    if peak_x > corr.shape[1] // 2:
        peak_x -= corr.shape[1]
    return float(peak_y), float(peak_x)


def drift_correct(norm_stack: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, list[dict[str, float]]]:
    corrected = [norm_stack[0]]
    rows = [{"frame": 0, "raw_shift_dy": 0.0, "raw_shift_dx": 0.0, "cumulative_dy": 0.0, "cumulative_dx": 0.0}]
    cumulative_y = 0.0
    cumulative_x = 0.0
    previous_corrected = norm_stack[0]
    for idx in range(1, norm_stack.shape[0]):
        dy, dx = phase_correlation_shift(previous_corrected, norm_stack[idx], mask)
        cumulative_y += dy
        cumulative_x += dx
        aligned = ndimage.shift(norm_stack[idx], shift=(-cumulative_y, -cumulative_x), order=1, mode="nearest")
        corrected.append(aligned.astype(np.float32))
        previous_corrected = aligned
        rows.append(
            {
                "frame": idx,
                "raw_shift_dy": dy,
                "raw_shift_dx": dx,
                "cumulative_dy": cumulative_y,
                "cumulative_dx": cumulative_x,
            }
        )
    return np.stack(corrected), rows


def estimate_pair_vectors(
    frame_a: np.ndarray,
    frame_b: np.ndarray,
    mask: np.ndarray,
    pair_index: int,
    window: int,
    step: int,
    search: int,
    min_std: float,
    min_mean: float,
    min_mask_fraction: float,
) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    half = window // 2
    for y in range(search + half, frame_a.shape[0] - search - half + 1, step):
        for x in range(search + half, frame_a.shape[1] - search - half + 1, step):
            block_mask = mask[y - half : y + half, x - half : x + half]
            mask_fraction = float(block_mask.mean()) if block_mask.size else 0.0
            if mask_fraction < min_mask_fraction:
                continue

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
            corr_median = float(np.median(corr))
            corr_std = float(np.std(corr))
            peak_z = (peak - corr_median) / (corr_std + 1e-12)

            rows.append(
                {
                    "pair_index": float(pair_index),
                    "x_pix": float(x),
                    "y_pix": float(y),
                    "dx_pix_per_frame": dx,
                    "dy_pix_per_frame": dy,
                    "speed_pix_per_frame": float(math.hypot(dx, dy)),
                    "peak": peak,
                    "peak_z": float(peak_z),
                    "mask_fraction": mask_fraction,
                    "template_mean": float(template.mean()),
                    "template_std": float(template.std()),
                }
            )
    return rows


def estimate_vectors(
    stack: np.ndarray,
    mask: np.ndarray,
    window: int,
    step: int,
    search: int,
    min_std: float,
    min_mean: float,
    min_mask_fraction: float,
) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for idx in range(stack.shape[0] - 1):
        rows.extend(
            estimate_pair_vectors(
                stack[idx],
                stack[idx + 1],
                mask,
                idx,
                window,
                step,
                search,
                min_std,
                min_mean,
                min_mask_fraction,
            )
        )
    return rows


def add_physical_units(rows: list[dict[str, float]], meta: MovieMeta) -> None:
    for row in rows:
        row["dx_um_per_frame"] = row["dx_pix_per_frame"] * meta.voxel_um_x
        row["dy_um_per_frame"] = row["dy_pix_per_frame"] * meta.voxel_um_y
        row["speed_um_per_frame"] = float(
            math.hypot(row["dx_um_per_frame"], row["dy_um_per_frame"])
        )


def iqr(values: np.ndarray) -> float | None:
    if values.size == 0:
        return None
    q25, q75 = np.percentile(values, [25, 75])
    return float(q75 - q25)


def summarize_movie(
    meta: MovieMeta,
    rows: list[dict[str, float]],
    drift_rows: list[dict[str, float]],
    mask: np.ndarray,
    assumed_dt_s: float | None,
) -> dict[str, Any]:
    speeds = np.array([row["speed_um_per_frame"] for row in rows], dtype=float)
    dx = np.array([row["dx_um_per_frame"] for row in rows], dtype=float)
    dy = np.array([row["dy_um_per_frame"] for row in rows], dtype=float)
    peak_z = np.array([row["peak_z"] for row in rows], dtype=float)
    cumulative_dx = drift_rows[-1]["cumulative_dx"] if drift_rows else 0.0
    cumulative_dy = drift_rows[-1]["cumulative_dy"] if drift_rows else 0.0
    median_speed = float(np.median(speeds)) if speeds.size else None
    return {
        "source_file": meta.source_file,
        "line": meta.line,
        "reader": meta.reader,
        "frames": meta.frames,
        "height": meta.height,
        "width": meta.width,
        "z_count": meta.z_count,
        "z_slices_1based": meta.z_slices_1based,
        "voxel_um_x": meta.voxel_um_x,
        "voxel_um_y": meta.voxel_um_y,
        "metadata_dt_s": meta.metadata_dt_s,
        "assumed_dt_s": assumed_dt_s,
        "mask_fraction": float(mask.mean()),
        "vectors": len(rows),
        "vectors_per_pair": len(rows) / max(1, meta.frames - 1),
        "median_speed_um_per_frame": median_speed,
        "iqr_speed_um_per_frame": iqr(speeds),
        "mean_speed_um_per_frame": float(np.mean(speeds)) if speeds.size else None,
        "median_dx_um_per_frame": float(np.median(dx)) if dx.size else None,
        "median_dy_um_per_frame": float(np.median(dy)) if dy.size else None,
        "median_peak_z": float(np.median(peak_z)) if peak_z.size else None,
        "median_speed_um_per_s_assumed_dt": median_speed / assumed_dt_s
        if median_speed is not None and assumed_dt_s
        else None,
        "median_speed_um_per_s_metadata_dt": median_speed / meta.metadata_dt_s
        if median_speed is not None and meta.metadata_dt_s
        else None,
        "total_drift_dx_pix": cumulative_dx,
        "total_drift_dy_pix": cumulative_dy,
        "calibration_note": meta.calibration_note,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def slug(path: Path) -> str:
    return "_".join(path.with_suffix("").parts[-2:]).replace(" ", "_").replace("=", "")


def save_preview(path: Path, raw_stack: np.ndarray, corrected_stack: np.ndarray, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    indices = sorted({0, raw_stack.shape[0] // 2, raw_stack.shape[0] - 1})
    fig, axes = plt.subplots(2, len(indices), figsize=(4 * len(indices), 7), constrained_layout=True)
    if len(indices) == 1:
        axes = np.array(axes).reshape(2, 1)
    for col, idx in enumerate(indices):
        axes[0, col].imshow(raw_stack[idx], cmap="gray", vmin=0, vmax=1)
        axes[0, col].set_title(f"raw frame {idx}")
        axes[0, col].axis("off")
        axes[1, col].imshow(corrected_stack[idx], cmap="gray", vmin=0, vmax=1)
        axes[1, col].set_title(f"drift-corrected frame {idx}")
        axes[1, col].axis("off")
    fig.suptitle(title)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_mask_overlay(path: Path, stack: np.ndarray, mask: np.ndarray, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    base = np.median(stack, axis=0)
    fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
    ax.imshow(base, cmap="gray", vmin=0, vmax=1)
    ax.contour(mask.astype(float), levels=[0.5], colors=["cyan"], linewidths=0.8)
    ax.set_title(title)
    ax.axis("off")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_flow_overlay(path: Path, stack: np.ndarray, rows: list[dict[str, float]], title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    base = np.median(stack, axis=0)
    fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
    ax.imshow(base, cmap="gray", vmin=0, vmax=1)
    if rows:
        x = np.array([row["x_pix"] for row in rows])
        y = np.array([row["y_pix"] for row in rows])
        dx = np.array([row["dx_pix_per_frame"] for row in rows])
        dy = np.array([row["dy_pix_per_frame"] for row in rows])
        speed = np.array([row["speed_um_per_frame"] for row in rows])
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


def save_drift_plot(path: Path, drift_rows: list[dict[str, float]], title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frames = [row["frame"] for row in drift_rows]
    dx = [row["cumulative_dx"] for row in drift_rows]
    dy = [row["cumulative_dy"] for row in drift_rows]
    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax.plot(frames, dx, marker="o", label="dx")
    ax.plot(frames, dy, marker="o", label="dy")
    ax.set_xlabel("Frame")
    ax.set_ylabel("Estimated cumulative drift (pixels)")
    ax.set_title(title)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def discover_ims(root: Path) -> list[Path]:
    return sorted(root.glob("*/*.ims"))


def optional_float(value: str | None, default: float | None = None) -> float | None:
    if value is None or value == "":
        return default
    return float(value)


def optional_int(value: str | None, default: int) -> int:
    if value is None or value == "":
        return default
    return int(value)


def manifest_path(value: str, manifest_file: Path) -> Path:
    path = Path(value)
    if path.exists() or path.is_absolute():
        return path
    candidate = manifest_file.parent / path
    if candidate.exists():
        return candidate
    return path


def process_manifest(manifest_file: Path, args: argparse.Namespace, outdir: Path) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    with manifest_file.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            source = row.get("file") or row.get("path") or row.get("projected_tiff")
            if not source:
                raise ValueError("Manifest rows require `file`, `path`, or `projected_tiff`.")

            row_args = argparse.Namespace(**vars(args))
            row_args.line = row.get("line") or row.get("condition") or args.line
            row_args.pixel_size_um = optional_float(
                row.get("pixel_size_um") or row.get("pixel_width_um"),
                args.pixel_size_um,
            )
            row_args.assumed_dt_s = optional_float(
                row.get("frame_interval_s") or row.get("assumed_dt_s"),
                args.assumed_dt_s,
            )
            row_args.z_start = optional_int(row.get("z_start"), args.z_start)
            row_args.z_end = optional_int(row.get("z_stop") or row.get("z_end"), args.z_end)

            path = manifest_path(source, manifest_file)
            print(f"Processing {path}")
            summary = process_one(path, row_args, outdir)
            summary["manifest_movie_id"] = row.get("movie_id") or Path(source).stem
            summary["manifest_raw_file"] = row.get("raw_file") or ""
            summary["manifest_frame_interval_source"] = row.get("frame_interval_source") or ""
            summary["manifest_notes"] = row.get("notes") or ""
            summaries.append(summary)
    return summaries


def process_one(path: Path, args: argparse.Namespace, outdir: Path) -> dict[str, Any]:
    stack_raw, meta = read_movie(path, args.z_start, args.z_end, args.pixel_size_um, args.line)
    norm_stack = np.stack([normalize_frame(frame) for frame in stack_raw])
    mask = foreground_mask(norm_stack, args.mask_percentile_floor)
    if args.no_drift_correction:
        corrected = norm_stack
        drift_rows = [
            {
                "frame": idx,
                "raw_shift_dy": 0.0,
                "raw_shift_dx": 0.0,
                "cumulative_dy": 0.0,
                "cumulative_dx": 0.0,
            }
            for idx in range(norm_stack.shape[0])
        ]
    else:
        corrected, drift_rows = drift_correct(norm_stack, mask)

    vectors = estimate_vectors(
        corrected,
        mask,
        args.window,
        args.step,
        args.search,
        args.min_std,
        args.min_mean,
        args.min_mask_fraction,
    )
    add_physical_units(vectors, meta)
    summary = summarize_movie(meta, vectors, drift_rows, mask, args.assumed_dt_s)

    name = slug(path)
    vector_fields = [
        "pair_index",
        "x_pix",
        "y_pix",
        "dx_pix_per_frame",
        "dy_pix_per_frame",
        "speed_pix_per_frame",
        "dx_um_per_frame",
        "dy_um_per_frame",
        "speed_um_per_frame",
        "peak",
        "peak_z",
        "mask_fraction",
        "template_mean",
        "template_std",
    ]
    write_csv(outdir / "vectors" / f"{name}_vectors.csv", vectors, vector_fields)
    write_csv(outdir / "drift" / f"{name}_drift.csv", drift_rows)
    save_preview(outdir / "qc" / f"{name}_preview.png", norm_stack, corrected, meta.source_file)
    save_mask_overlay(outdir / "qc" / f"{name}_mask.png", corrected, mask, meta.source_file)
    save_flow_overlay(outdir / "qc" / f"{name}_flow.png", corrected, vectors, meta.source_file)
    save_drift_plot(outdir / "qc" / f"{name}_drift.png", drift_rows, meta.source_file)
    return summary


def write_readme(outdir: Path, args: argparse.Namespace, summaries: list[dict[str, Any]]) -> None:
    lines = [
        "# Egg-Cell Network Flow Analysis",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "This is the movie-level analysis branch intended to replace hand-picked kymograph slopes.",
        "Primary units are `um/frame`; `um/s` is only reported when an assumed or metadata frame interval is present.",
        "",
        "## Parameters",
        "",
        f"- z projection: slices {args.z_start}-{args.z_end}, max intensity for direct `.ims` input",
        f"- drift correction: {'off' if args.no_drift_correction else 'phase-correlation global drift correction'}",
        f"- PIV window/step/search: {args.window}/{args.step}/{args.search} pixels",
        f"- mask percentile floor: {args.mask_percentile_floor}",
        f"- minimum mask fraction per window: {args.min_mask_fraction}",
        f"- assumed dt: {args.assumed_dt_s}",
        "",
        "## Outputs",
        "",
        "- `movie_summary.csv`: one row per movie.",
        "- `vectors/*_vectors.csv`: local displacement vectors.",
        "- `drift/*_drift.csv`: estimated cumulative drift correction.",
        "- `qc/*_preview.png`: raw vs drift-corrected frame previews.",
        "- `qc/*_mask.png`: foreground mask overlay.",
        "- `qc/*_flow.png`: vector overlay.",
        "- `qc/*_drift.png`: cumulative drift plot.",
        "",
        "## Calibration Warning",
        "",
        "Direct `.ims` reading uses Imaris HDF5 metadata and may not exactly match Bio-Formats/Fiji series calibration.",
        "For final reporting, prefer Fiji/Bio-Formats exported projected TIFF/OME-TIFF stacks with explicit pixel size.",
        "",
        "## Movies Processed",
        "",
    ]
    for row in summaries:
        lines.append(
            f"- {row['source_file']}: {row['vectors']} vectors, "
            f"median {row['median_speed_um_per_frame']} um/frame"
        )
    (outdir / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="*", help="Movie files. Defaults to all .ims under --root.")
    parser.add_argument("--manifest", default=None, help="CSV with per-movie file, line, pixel size, and frame interval.")
    parser.add_argument("--root", default="Kwaku and Berry")
    parser.add_argument("--output", default="analysis/results/egg_cell_network_flow_2026-06-15")
    parser.add_argument("--z-start", type=int, default=5, help="First z slice, 1-based inclusive.")
    parser.add_argument("--z-end", type=int, default=7, help="Last z slice, 1-based inclusive.")
    parser.add_argument("--pixel-size-um", type=float, default=None, help="Required for TIFF input.")
    parser.add_argument("--line", default=None, help="Optional line/genotype label for TIFF input.")
    parser.add_argument("--assumed-dt-s", type=float, default=None)
    parser.add_argument("--window", type=int, default=48)
    parser.add_argument("--step", type=int, default=32)
    parser.add_argument("--search", type=int, default=8)
    parser.add_argument("--min-std", type=float, default=0.025)
    parser.add_argument("--min-mean", type=float, default=0.06)
    parser.add_argument("--mask-percentile-floor", type=float, default=65)
    parser.add_argument("--min-mask-fraction", type=float, default=0.12)
    parser.add_argument("--no-drift-correction", action="store_true")
    args = parser.parse_args()

    root = Path(args.root)
    paths = [Path(item) for item in args.files] if args.files else discover_ims(root)
    if not args.manifest and not paths:
        raise SystemExit("No input movies found.")

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.manifest:
        summaries = process_manifest(Path(args.manifest), args, outdir)
    else:
        summaries = []
        for path in paths:
            print(f"Processing {path}")
            summaries.append(process_one(path, args, outdir))

    summary_fields = list(summaries[0].keys())
    write_csv(outdir / "movie_summary.csv", summaries, summary_fields)
    write_readme(outdir, args, summaries)
    print(outdir / "movie_summary.csv")


if __name__ == "__main__":
    main()
