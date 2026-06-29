#!/usr/bin/env python3
"""Render timecourse Cytosim snapshots for 80-cluster, 12-mpc minus-z annulus runs.

The output is designed to replace overly wide unwrapped-annulus panels with
actual Cytosim play snapshots.  Two views are rendered:

- zoomed oblique view: shows the annulus/cylindrical shell context.
- upright side view: renders the side view, then rotates the image upright to
  mimic a cell-like vertical projection.
"""

from __future__ import annotations

import math
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "analysis" / "results" / "clu_minusz_slides_2026-06-02" / "mpc12_c80_timecourse_play_snapshots"
RAW = OUT / "raw"
PROCESSED = OUT / "processed"
PANELS = OUT / "panels"

PLAY = ROOT / "clu_init_minusz" / "mpc12" / "c80_m12" / "play"

RUNS = {
    "rotatable": {
        "label": "Rotatable",
        "run": ROOT / "clu_init_minusz" / "mpc12" / "c80_m12" / "r0024",
    },
    "fixed_orientation": {
        "label": "Fixed",
        "run": ROOT / "clu_fixed_global_init_minusz" / "mpc12" / "c80_m12" / "r0024",
    },
}

# Configuration for these runs:
# run 10000, nb_frames 40; run 25000, nb_frames 100; motorized run 200000,
# nb_frames 400; timestep 0.004 s.  Therefore the motorized phase is 800 s
# with 2 s between saved frames.  Global frame 139 is the last pre-motor frame,
# and global frame 539 is the final frame.
PRE_MOTOR_FRAME = 139
FINAL_FRAME = 539
MOTORIZED_DURATION_MIN = 800.0 / 60.0
TIMEPOINTS = [
    ("t00", 0.0, PRE_MOTOR_FRAME, "0 min"),
    ("t04", 4.0, PRE_MOTOR_FRAME + round(4.0 * 60.0 / 2.0), "4 min"),
    ("t08", 8.0, PRE_MOTOR_FRAME + round(8.0 * 60.0 / 2.0), "8 min"),
    ("t12", 12.0, PRE_MOTOR_FRAME + round(12.0 * 60.0 / 2.0), "12 min"),
    ("tfinal", MOTORIZED_DURATION_MIN, FINAL_FRAME, "13.3 min"),
]

RENDER_SIZE = 2400
SAMPLES = 12


@dataclass(frozen=True)
class Camera:
    key: str
    label: str
    rotation: str
    zoom: float
    view_scale: float
    postprocess_rot90: int = 0


def quat(axis: str, degrees: float) -> tuple[float, float, float, float]:
    angle = math.radians(degrees) / 2.0
    c = math.cos(angle)
    s = math.sin(angle)
    if axis == "x":
        return c, s, 0.0, 0.0
    if axis == "y":
        return c, 0.0, s, 0.0
    if axis == "z":
        return c, 0.0, 0.0, s
    raise KeyError(axis)


def qmul(q2: tuple[float, float, float, float], q1: tuple[float, float, float, float]) -> tuple[float, float, float, float]:
    w2, x2, y2, z2 = q2
    w1, x1, y1, z1 = q1
    return (
        w2 * w1 - x2 * x1 - y2 * y1 - z2 * z1,
        w2 * x1 + x2 * w1 + y2 * z1 - z2 * y1,
        w2 * y1 - x2 * z1 + y2 * w1 + z2 * x1,
        w2 * z1 + x2 * y1 - y2 * x1 + z2 * w1,
    )


def qstr(q: tuple[float, float, float, float]) -> str:
    norm = math.sqrt(sum(value * value for value in q))
    return " ".join(f"{value / norm:.8f}" for value in q)


CAMERAS = [
    Camera(
        "oblique_full_annulus",
        "Oblique view",
        qstr(qmul(quat("z", 35), quat("y", 70))),
        zoom=0.95,
        view_scale=62.0,
        postprocess_rot90=2,
    ),
    Camera(
        "upright_side",
        "Side view, upright",
        qstr(quat("y", 90)),
        zoom=0.88,
        view_scale=68.0,
        postprocess_rot90=-1,
    ),
]


def raw_path(condition: str, camera: Camera, time_slug: str, frame: int) -> Path:
    return RAW / f"{condition}_{camera.key}_{time_slug}_frame{frame:04d}_raw_{RENDER_SIZE}px.png"


def processed_path(condition: str, camera: Camera, time_slug: str, frame: int, style: str) -> Path:
    return PROCESSED / f"{condition}_{camera.key}_{time_slug}_frame{frame:04d}_{style}_{RENDER_SIZE}px.png"


def render_raw(condition: str, run_dir: Path, camera: Camera, time_slug: str, frame: int) -> tuple[Path, str, int]:
    RAW.mkdir(parents=True, exist_ok=True)
    output = raw_path(condition, camera, time_slug, frame)
    if output.exists():
        return output, "skipped existing output", 0

    before = {path.name for path in RAW.glob("image*.png")}
    cmd = [
        str(PLAY),
        "on",
        "image",
        f"frame={frame}",
        f"size={RENDER_SIZE},{RENDER_SIZE}",
        f"samples={SAMPLES}",
        "auto_scale=0",
        f"view_scale={camera.view_scale}",
        f"zoom={camera.zoom}",
        f"rotation={camera.rotation}",
        "focus=0 0 0",
        "draw_memo=0",
        "label=none",
        "couple_select=0",
        "single_select=0",
        "draw_links=0",
        "back_color=white",
        "front_color=black",
        f"image_dir={RAW}",
        f"image_name={condition}_{camera.key}_{time_slug}_frame{frame:04d}_raw_{RENDER_SIZE}px",
        "image_format=png",
    ]
    proc = subprocess.run(cmd, cwd=run_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if not output.exists():
        new_images = [path for path in RAW.glob("image*.png") if path.name not in before]
        if len(new_images) == 1:
            shutil.move(new_images[0], output)
    if not output.exists():
        raise RuntimeError(f"play did not produce {output}\nreturn={proc.returncode}\n{proc.stdout}")
    return output, " ".join(cmd), proc.returncode


def signal_from_raw(path: Path) -> np.ndarray:
    arr = np.asarray(Image.open(path).convert("RGB"), dtype=float) / 255.0
    lum = 0.299 * arr[..., 0] + 0.587 * arr[..., 1] + 0.114 * arr[..., 2]
    signal = np.clip(1.0 - lum, 0.0, 1.0)
    edge = max(60, signal.shape[0] // 32)
    signal[:edge, :] = 0.0
    signal[-edge:, :] = 0.0
    return signal


def union_crop(signals: list[np.ndarray], *, margin_fraction: float = 0.08) -> tuple[slice, slice]:
    mask = np.zeros_like(signals[0], dtype=bool)
    for signal in signals:
        mask |= signal > 0.018
    ys, xs = np.where(mask)
    if ys.size == 0 or xs.size == 0:
        return slice(0, signals[0].shape[0]), slice(0, signals[0].shape[1])
    height, width = signals[0].shape
    y0, y1 = int(ys.min()), int(ys.max()) + 1
    x0, x1 = int(xs.min()), int(xs.max()) + 1
    pad_y = int(round((y1 - y0) * margin_fraction))
    pad_x = int(round((x1 - x0) * margin_fraction))
    y0 = max(0, y0 - pad_y)
    y1 = min(height, y1 + pad_y)
    x0 = max(0, x0 - pad_x)
    x1 = min(width, x1 + pad_x)
    return slice(y0, y1), slice(x0, x1)


def fluorescence_transform(signal: np.ndarray, style: str) -> np.ndarray:
    fine = gaussian_filter(signal, sigma=1.2)
    glow = gaussian_filter(signal, sigma=7.0)
    image = fine + 0.43 * glow
    positive = image[image > 0]
    scale = float(np.quantile(positive, 0.996)) if positive.size else 1.0
    image = np.power(np.clip(image / max(scale, 1e-12), 0.0, 1.0), 0.70)

    if style == "cyan_fluorescence":
        return np.dstack([0.035 * image, 0.92 * image, image])
    if style == "inverted_fluorescence":
        gray = 1.0 - 0.90 * image
        return np.dstack([gray, gray, gray])
    if style == "line_render":
        gray = 1.0 - 0.98 * np.power(signal, 0.82)
        return np.dstack([gray, gray, gray])
    raise KeyError(style)


def orient_image(rgb: np.ndarray, camera: Camera) -> np.ndarray:
    return np.rot90(rgb, k=camera.postprocess_rot90)


def process_images(raw_records: list[dict]) -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    for camera in CAMERAS:
        camera_records = [record for record in raw_records if record["camera"].key == camera.key]
        signals = [signal_from_raw(record["raw_path"]) for record in camera_records]
        crop_y, crop_x = union_crop(signals)
        for record, signal in zip(camera_records, signals):
            cropped = signal[crop_y, crop_x]
            for style in ("cyan_fluorescence", "inverted_fluorescence", "line_render"):
                rgb = orient_image(fluorescence_transform(cropped, style), camera)
                output = processed_path(record["condition"], camera, record["time_slug"], record["frame"], style)
                Image.fromarray(np.uint8(np.clip(rgb, 0.0, 1.0) * 255)).save(output)
                record.setdefault("processed", {})[style] = output


def add_panel_image(ax: plt.Axes, path: Path, *, background: str) -> None:
    ax.imshow(mpimg.imread(path))
    ax.set_axis_off()
    ax.set_facecolor(background)


def make_timecourse_panel(camera: Camera, style: str) -> Path:
    PANELS.mkdir(parents=True, exist_ok=True)
    background = "black" if style == "cyan_fluorescence" else "white"
    text_color = "white" if background == "black" else "black"
    ncols = len(TIMEPOINTS)
    fig_height = 7.8 if camera.key == "upright_side" else 5.4
    fig, axes = plt.subplots(2, ncols, figsize=(15.6, fig_height), constrained_layout=False)
    fig.patch.set_facecolor(background)

    for row, condition in enumerate(("rotatable", "fixed_orientation")):
        row_label = RUNS[condition]["label"]
        for col, (time_slug, _time_min, frame, label) in enumerate(TIMEPOINTS):
            ax = axes[row, col]
            path = processed_path(condition, camera, time_slug, frame, style)
            add_panel_image(ax, path, background=background)
            if row == 0:
                ax.text(
                    0.5,
                    1.035,
                    label,
                    transform=ax.transAxes,
                    ha="center",
                    va="bottom",
                    color=text_color,
                    fontsize=12.5,
                    weight="bold",
                )
            if col == 0:
                ax.text(
                    0.03,
                    0.96,
                    row_label,
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    color=text_color,
                    fontsize=11.5,
                    weight="bold",
                    bbox={"boxstyle": "round,pad=0.20", "facecolor": background, "edgecolor": "none", "alpha": 0.70},
                )

    fig.text(
        0.5,
        0.992,
        "80 clusters, 12 motors/cluster, 1:16 crosslinkers",
        ha="center",
        va="top",
        color=text_color,
        fontsize=14.0,
    )
    fig.subplots_adjust(left=0.012, right=0.988, top=0.885, bottom=0.025, wspace=0.025, hspace=0.07)
    output = PANELS / f"mpc12_c80_1to16_timecourse_{camera.key}_{style}.png"
    fig.savefig(output, dpi=600, bbox_inches="tight", facecolor=background)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight", facecolor=background)
    plt.close(fig)
    return output


def write_readme(raw_records: list[dict], panels: list[Path]) -> None:
    lines = [
        "# 80-Cluster 12-mpc Timecourse Play Snapshots",
        "",
        "Actual Cytosim `play` snapshots for the initially minus-z-aligned annulus simulations.",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "Condition:",
        "- Family: `mpc12`",
        "- Case: `c80_m12`",
        "- Crosslinker regime: `1:16`",
        "- Comparison: rotatable versus fixed motor clusters",
        "",
        "Runs:",
    ]
    for condition, spec in RUNS.items():
        lines.append(f"- `{condition}`: `{spec['run'].relative_to(ROOT)}`")
    lines.extend(
        [
            "",
            "Timepoints:",
            "- Global frame 139 is the last pre-motor frame and is labeled `0 min`.",
            "- The motorized phase is 800 s with 2 s between saved frames.",
        ]
    )
    for time_slug, time_min, frame, label in TIMEPOINTS:
        lines.append(f"- `{time_slug}`: {label}; frame `{frame}`; motorized time `{time_min:.3f} min`")
    lines.extend(["", "Cameras:"])
    for camera in CAMERAS:
        lines.append(
            f"- `{camera.key}`: {camera.label}; rotation=`{camera.rotation}`, zoom={camera.zoom}, "
            f"view_scale={camera.view_scale}, postprocess rot90={camera.postprocess_rot90}"
        )
    lines.extend(["", "Panels:"])
    for panel in panels:
        lines.append(f"- `{panel.relative_to(ROOT)}`")
    lines.extend(["", "Render commands:"])
    for record in raw_records:
        lines.append(
            f"- `{record['condition']}` `{record['camera'].key}` `{record['time_slug']}` "
            f"frame={record['frame']} return={record['returncode']} output=`{record['raw_path'].relative_to(ROOT)}`"
        )
        lines.append(f"  `{record['command']}`")
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    raw_records: list[dict] = []
    for condition, spec in RUNS.items():
        run_dir = spec["run"]
        if not run_dir.exists():
            raise FileNotFoundError(run_dir)
        for camera in CAMERAS:
            for time_slug, _time_min, frame, _label in TIMEPOINTS:
                raw, command, returncode = render_raw(condition, run_dir, camera, time_slug, frame)
                raw_records.append(
                    {
                        "condition": condition,
                        "camera": camera,
                        "time_slug": time_slug,
                        "frame": frame,
                        "raw_path": raw,
                        "command": command,
                        "returncode": returncode,
                    }
                )

    process_images(raw_records)
    panels = []
    for camera in CAMERAS:
        for style in ("cyan_fluorescence", "inverted_fluorescence", "line_render"):
            panels.append(make_timecourse_panel(camera, style))
    write_readme(raw_records, panels)
    print(OUT)


if __name__ == "__main__":
    main()
