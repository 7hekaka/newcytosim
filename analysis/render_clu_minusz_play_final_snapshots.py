#!/usr/bin/env python3
"""Render final microscopy-style Cytosim snapshots for clu_minusz simulations.

This script is intentionally narrower than render_clu_minusz_play_camera_tests.py:
it renders only the two camera views that were useful in the camera test pass,
then builds compact cyan and inverted-fluorescence panels for slides/manuscript
discussion.
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
OUT = ROOT / "analysis" / "results" / "clu_minusz_slides_2026-06-02" / "play_final_side_snapshots"
RAW = OUT / "raw"
PROCESSED = OUT / "processed"

RUNS = {
    "rotatable": ROOT / "clu_init_minusz" / "total480" / "c80_m6" / "r0028",
    "fixed_orientation": ROOT / "clu_fixed_global_init_minusz" / "total480" / "c80_m6" / "r0027",
}

PLAY = ROOT / "clu_init_minusz" / "total480" / "c80_m6" / "play"
FRAME = 539
RENDER_SIZE = 3000
SAMPLES = 16


@dataclass(frozen=True)
class Camera:
    key: str
    label: str
    rotation: str
    zoom: float
    view_scale: float


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
    Camera("side_y_90", "Side view", qstr(quat("y", 90)), 1.25, 47.0),
    Camera("oblique_y70_z35", "Oblique side view", qstr(qmul(quat("z", 35), quat("y", 70))), 1.25, 50.0),
]


def raw_path(model: str, camera: Camera) -> Path:
    return RAW / f"{model}_{camera.key}_raw_{RENDER_SIZE}px.png"


def processed_path(model: str, camera: Camera, style: str) -> Path:
    return PROCESSED / f"{model}_{camera.key}_{style}_{RENDER_SIZE}px.png"


def render_raw(run_dir: Path, model: str, camera: Camera, *, overwrite: bool = False) -> tuple[Path, str, int]:
    RAW.mkdir(parents=True, exist_ok=True)
    output = raw_path(model, camera)
    if output.exists() and not overwrite:
        return output, "skipped existing output", 0

    before = {path.name for path in RAW.glob("image*.png")}
    cmd = [
        str(PLAY),
        "on",
        "image",
        f"frame={FRAME}",
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
        f"image_name={model}_{camera.key}_raw_{RENDER_SIZE}px",
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
    # Cytosim can write small overlay labels near image edges; suppress those.
    edge = max(80, signal.shape[0] // 30)
    signal[:edge, :] = 0.0
    signal[-edge:, :] = 0.0
    return signal


def fluorescence_transform(signal: np.ndarray, mode: str) -> np.ndarray:
    fine = gaussian_filter(signal, sigma=1.2)
    glow = gaussian_filter(signal, sigma=7.5)
    image = fine + 0.45 * glow
    positive = image[image > 0]
    scale = float(np.quantile(positive, 0.996)) if positive.size else 1.0
    image = np.clip(image / max(scale, 1e-12), 0.0, 1.0)
    image = np.power(image, 0.70)

    if mode == "cyan_fluorescence":
        return np.dstack([0.035 * image, 0.90 * image, image])
    if mode == "inverted_fluorescence":
        gray = 1.0 - 0.90 * image
        return np.dstack([gray, gray, gray])
    if mode == "line_render":
        gray = 1.0 - 0.98 * np.power(signal, 0.85)
        return np.dstack([gray, gray, gray])
    raise KeyError(mode)


def process_raw(raw: Path, model: str, camera: Camera) -> list[Path]:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    signal = signal_from_raw(raw)
    outputs: list[Path] = []
    for style in ("cyan_fluorescence", "inverted_fluorescence", "line_render"):
        output = processed_path(model, camera, style)
        rgb = fluorescence_transform(signal, style)
        Image.fromarray(np.uint8(np.clip(rgb, 0.0, 1.0) * 255)).save(output)
        outputs.append(output)
    return outputs


def add_image(ax: plt.Axes, path: Path, label: str, *, background: str, text_color: str) -> None:
    ax.imshow(mpimg.imread(path))
    ax.set_axis_off()
    ax.set_facecolor(background)
    ax.text(
        0.03,
        0.96,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        color=text_color,
        fontsize=11,
        weight="bold",
        bbox={"boxstyle": "round,pad=0.20", "facecolor": background, "edgecolor": "none", "alpha": 0.72},
    )


def make_panel(style: str) -> Path:
    background = "black" if style == "cyan_fluorescence" else "white"
    text_color = "white" if background == "black" else "black"
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 10.2), constrained_layout=False)
    fig.patch.set_facecolor(background)

    for row, (model, model_label) in enumerate((("rotatable", "Rotatable"), ("fixed_orientation", "Fixed orientation"))):
        for col, camera in enumerate(CAMERAS):
            path = processed_path(model, camera, style)
            add_image(axes[row, col], path, f"{model_label}: {camera.label}", background=background, text_color=text_color)

    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.02, hspace=0.02)
    output = OUT / f"clu_minusz_final_side_snapshots_{style}.png"
    fig.savefig(output, dpi=600, bbox_inches="tight", facecolor=background)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight", facecolor=background)
    plt.close(fig)
    return output


def write_readme(commands: list[tuple[str, str, Path, str, int]]) -> None:
    lines = [
        "# Final Cytosim Side-View Snapshots",
        "",
        "High-resolution microscopy-style renders for the initially minus-z-aligned annulus simulations.",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"Frame: `{FRAME}`",
        f"Render size: `{RENDER_SIZE} x {RENDER_SIZE}`",
        f"Samples: `{SAMPLES}`",
        "",
        "Runs:",
    ]
    for model, run_dir in RUNS.items():
        lines.append(f"- `{model}`: `{run_dir.relative_to(ROOT)}`")
    lines.extend(["", "Cameras:"])
    for camera in CAMERAS:
        lines.append(
            f"- `{camera.key}`: {camera.label}; rotation=`{camera.rotation}`, zoom={camera.zoom}, view_scale={camera.view_scale}"
        )
    lines.extend(
        [
            "",
            "Styles:",
            "- `cyan_fluorescence`: black background with cyan glow, closest to confocal-style presentation.",
            "- `inverted_fluorescence`: white background with gray/black signal, closest to inverted fluorescence figures.",
            "- `line_render`: white background line-style render for structural comparison.",
            "",
            "Notes:",
            "- `play on image` writes the PNG before exiting with a segmentation fault on this WSL/OpenGL setup; output-file existence is treated as success.",
            "- The camera-test diagnostics remain in `analysis/results/clu_minusz_slides_2026-06-02/play_snapshot_camera_tests/`.",
            "",
            "Render commands:",
        ]
    )
    for model, camera_key, output, command, returncode in commands:
        lines.append(f"- `{model}` `{camera_key}` return={returncode} output=`{output.relative_to(ROOT)}`")
        lines.append(f"  `{command}`")
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    commands: list[tuple[str, str, Path, str, int]] = []
    for model, run_dir in RUNS.items():
        if not run_dir.exists():
            raise FileNotFoundError(run_dir)
        for camera in CAMERAS:
            raw, command, returncode = render_raw(run_dir, model, camera)
            process_raw(raw, model, camera)
            commands.append((model, camera.key, raw, command, returncode))

    for style in ("cyan_fluorescence", "inverted_fluorescence", "line_render"):
        make_panel(style)

    write_readme(commands)
    print(OUT)


if __name__ == "__main__":
    main()
