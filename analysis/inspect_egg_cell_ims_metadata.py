#!/usr/bin/env python3
"""Inspect Imaris .ims metadata for the egg-cell actin velocity dataset.

The collaborator folder contains Imaris .ims files and Olympus .oib files.
The .ims files are HDF5 containers, so this script extracts the metadata we
need before building a velocity pipeline: image dimensions, physical extents,
time points, channels, and the primary data-array shape.

Dependency:
    python -m pip install h5py
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import json
from statistics import median
from pathlib import Path
from typing import Any

try:
    import h5py
except ImportError as exc:  # pragma: no cover - exercised when dependency absent.
    raise SystemExit(
        "Missing dependency: h5py. Install with `python -m pip install h5py`."
    ) from exc


INTERESTING_ATTRS = {
    "ImageSizeX",
    "ImageSizeY",
    "ImageSizeZ",
    "ExtMin0",
    "ExtMin1",
    "ExtMin2",
    "ExtMax0",
    "ExtMax1",
    "ExtMax2",
    "Unit",
    "RecordingDate",
    "LensPower",
    "NumericalAperture",
    "MicroscopeModality",
    "MicroscopeType",
    "Description",
    "Name",
}


def decode_value(value: Any) -> Any:
    """Convert HDF5 scalar/byte values into JSON-friendly Python values."""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    dtype = getattr(value, "dtype", None)
    if dtype is not None and getattr(dtype, "kind", None) == "S":
        return value.tobytes().decode("utf-8", errors="replace").rstrip("\x00")
    if dtype is not None and getattr(dtype, "kind", None) == "U":
        return "".join(value.tolist())
    if isinstance(value, dict):
        return {str(key): decode_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [decode_value(item) for item in value]
    if hasattr(value, "shape") and getattr(value, "shape", None) == ():
        return decode_value(value.item())
    if hasattr(value, "tolist"):
        return decode_value(value.tolist())
    return value


def flatten_attrs(group: h5py.Group) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in group.attrs.items():
        if key in INTERESTING_ATTRS or key.startswith("TimePoint"):
            out[key] = decode_value(value)
    return out


def first_data_shape(handle: h5py.File) -> list[int] | None:
    shape: list[int] | None = None

    def visitor(_name: str, obj: h5py.Dataset | h5py.Group) -> None:
        nonlocal shape
        if shape is not None:
            return
        if isinstance(obj, h5py.Dataset) and obj.name.endswith("/Data"):
            shape = list(obj.shape)

    handle.visititems(visitor)
    return shape


def collect_channel_names(handle: h5py.File) -> list[str]:
    info = handle.get("DataSetInfo")
    if info is None:
        return []
    return sorted(name for name in info.keys() if name.lower().startswith("channel"))


def inspect_ims(path: Path, root: Path) -> dict[str, Any]:
    with h5py.File(path, "r") as handle:
        info = handle.get("DataSetInfo")
        image_attrs = {}
        time_attrs = {}
        channel_attrs: dict[str, dict[str, Any]] = {}

        if info is not None:
            image = info.get("Image")
            if image is not None:
                image_attrs = flatten_attrs(image)

            time_info = info.get("TimeInfo")
            if time_info is not None:
                time_attrs = flatten_attrs(time_info)

            for channel in collect_channel_names(handle):
                channel_group = info.get(channel)
                if channel_group is not None:
                    channel_attrs[channel] = flatten_attrs(channel_group)

        data_shape = first_data_shape(handle)
        ext = physical_extent_um(image_attrs)
        size = image_size(image_attrs, data_shape)
        intervals = frame_intervals_s(time_attrs)

        return {
            "path": path.as_posix(),
            "relative_path": path.relative_to(root).as_posix(),
            "line": path.parent.name,
            "data_shape": data_shape,
            "image_size_xyz": size,
            "physical_extent_um_xyz": ext,
            "voxel_size_um_xyz": voxel_size_um(size, ext),
            "channels": sorted(channel_attrs.keys()),
            "channel_attrs": channel_attrs,
            "image_attrs": image_attrs,
            "time_attrs": time_attrs,
            "timepoint_count": len([key for key in time_attrs if key.startswith("TimePoint")]),
            "frame_intervals_s": intervals,
            "median_frame_interval_s": median(intervals) if intervals else None,
            "duration_s": sum(intervals) if intervals else None,
        }


def image_size(
    attrs: dict[str, Any], data_shape: list[int] | None = None
) -> list[int | None]:
    size = [to_int(attrs.get(f"ImageSize{axis}")) for axis in ("X", "Y", "Z")]
    if all(value is not None for value in size):
        return size
    if data_shape is not None and len(data_shape) >= 3:
        # Imaris stores the image data as Z, Y, X.
        return [data_shape[-1], data_shape[-2], data_shape[-3]]
    return size


def physical_extent_um(attrs: dict[str, Any]) -> list[float | None]:
    extents: list[float | None] = []
    for idx in range(3):
        min_value = to_float(attrs.get(f"ExtMin{idx}"))
        max_value = to_float(attrs.get(f"ExtMax{idx}"))
        if min_value is None or max_value is None:
            extents.append(None)
        else:
            extents.append(max_value - min_value)
    return extents


def voxel_size_um(
    size_xyz: list[int | None], extent_xyz: list[float | None]
) -> list[float | None]:
    voxels: list[float | None] = []
    for size, extent in zip(size_xyz, extent_xyz):
        if size in (None, 0) or extent is None:
            voxels.append(None)
        else:
            voxels.append(extent / float(size))
    return voxels


def to_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def to_int(value: Any) -> int | None:
    if value is None:
        return None


def frame_intervals_s(attrs: dict[str, Any]) -> list[float]:
    points: list[tuple[int, datetime]] = []
    for key, value in attrs.items():
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


def parse_timepoint(value: Any) -> datetime | None:
    if not isinstance(value, str):
        return None
    for fmt in ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S"):
        try:
            return datetime.strptime(value, fmt)
        except ValueError:
            pass
    return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "relative_path",
        "line",
        "data_shape",
        "image_size_xyz",
        "physical_extent_um_xyz",
        "voxel_size_um_xyz",
        "timepoint_count",
        "median_frame_interval_s",
        "duration_s",
        "channels",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: json.dumps(row[field]) for field in fields})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "root",
        nargs="?",
        default="Kwaku and Berry",
        help="Root folder containing collaborator files.",
    )
    parser.add_argument(
        "--output",
        default="analysis/results/egg_cell_velocity_audit_2026-06-15/ims_metadata",
        help="Output directory for JSON/CSV summaries.",
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)

    rows = [inspect_ims(path, root) for path in sorted(root.rglob("*.ims"))]
    (output / "ims_metadata.json").write_text(json.dumps(rows, indent=2) + "\n")
    write_csv(rows, output / "ims_metadata_summary.csv")

    print(f"Inspected {len(rows)} .ims files")
    print(f"Wrote {output / 'ims_metadata_summary.csv'}")
    print(f"Wrote {output / 'ims_metadata.json'}")


if __name__ == "__main__":
    main()
