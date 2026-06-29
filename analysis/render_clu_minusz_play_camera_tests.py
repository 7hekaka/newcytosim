#!/usr/bin/env python3
from __future__ import annotations

import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "analysis" / "results" / "clu_minusz_slides_2026-06-02" / "play_snapshot_camera_tests"
RAW = OUT / "raw"
PROCESSED = OUT / "processed"

RUNS = {
    "rotatable": ROOT / "clu_init_minusz" / "total480" / "c80_m6" / "r0028",
    "fixed_orientation": ROOT / "clu_fixed_global_init_minusz" / "total480" / "c80_m6" / "r0027",
}
PLAY = ROOT / "clu_init_minusz" / "total480" / "c80_m6" / "play"
FRAME = 539


@dataclass(frozen=True)
class Camera:
    key: str
    label: str
    rotation: str | None
    zoom: float
    view_scale: float


def quat(axis: str, degrees: float) -> tuple[float, float, float, float]:
    a = math.radians(degrees) / 2.0
    c = math.cos(a)
    s = math.sin(a)
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
    n = math.sqrt(sum(v * v for v in q))
    return " ".join(f"{v / n:.8f}" for v in q)


CAMERAS = [
    Camera("top_default", "Top/default", None, 1.45, 42.0),
    Camera("side_x_90", "Side view, X 90", qstr(quat("x", 90)), 1.35, 44.0),
    Camera("side_y_90", "Side view, Y 90", qstr(quat("y", 90)), 1.35, 44.0),
    Camera("tilt_x_70", "Tilted side, X 70", qstr(quat("x", 70)), 1.35, 44.0),
    Camera("oblique_x70_z35", "Oblique side, X 70 + Z 35", qstr(qmul(quat("z", 35), quat("x", 70))), 1.35, 44.0),
    Camera("oblique_y70_z35", "Oblique side, Y 70 + Z 35", qstr(qmul(quat("z", 35), quat("y", 70))), 1.35, 44.0),
]


def image_path(model: str, camera: Camera, style: str) -> Path:
    root = RAW if style == "raw" else PROCESSED
    return root / f"{model}_{camera.key}_{style}.png"


def run_play(run_dir: Path, model: str, camera: Camera, *, size: int = 1800) -> Path:
    RAW.mkdir(parents=True, exist_ok=True)
    out = image_path(model, camera, "raw")
    if out.exists():
        return out
    before = {p.name for p in RAW.glob("image*.png")}
    cmd = [
        str(PLAY),
        "on",
        "image",
        f"frame={FRAME}",
        f"size={size},{size}",
        "samples=8",
        "auto_scale=0",
        f"view_scale={camera.view_scale}",
        f"zoom={camera.zoom}",
        "focus=0 0 0",
        "draw_memo=0",
        "label=none",
        "couple_select=0",
        "single_select=0",
        "draw_links=0",
        "back_color=white",
        "front_color=black",
        f"image_dir={RAW}",
        f"image_name={model}_{camera.key}_raw",
        "image_format=png",
    ]
    if camera.rotation:
        cmd.insert(9, f"rotation={camera.rotation}")
    proc = subprocess.run(cmd, cwd=run_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    produced = RAW / f"{model}_{camera.key}_raw.png"
    if not produced.exists():
        after = [p for p in RAW.glob("image*.png") if p.name not in before]
        if len(after) == 1:
            shutil.move(after[0], produced)
    if not produced.exists():
        raise RuntimeError(f"play did not produce {produced}\nreturn={proc.returncode}\n{proc.stdout}")
    return produced


def signal_from_raw(path: Path) -> np.ndarray:
    arr = np.asarray(Image.open(path).convert("RGB"), dtype=float) / 255.0
    lum = 0.299 * arr[..., 0] + 0.587 * arr[..., 1] + 0.114 * arr[..., 2]
    signal = np.clip(1.0 - lum, 0.0, 1.0)
    signal[:80, :] = 0.0
    signal[-100:, :] = 0.0
    return signal


def fluorescence_transform(signal: np.ndarray, mode: str) -> np.ndarray:
    fine = gaussian_filter(signal, sigma=1.0)
    glow = gaussian_filter(signal, sigma=6.0)
    img = fine + 0.40 * glow
    positive = img[img > 0]
    scale = float(np.quantile(positive, 0.995)) if positive.size else 1.0
    img = np.power(np.clip(img / max(scale, 1e-12), 0.0, 1.0), 0.72)
    if mode == "cyan":
        return np.dstack([0.05 * img, 0.95 * img, img])
    if mode == "inverted":
        gray = 1.0 - 0.92 * img
        return np.dstack([gray, gray, gray])
    raise KeyError(mode)


def process_image(raw: Path, model: str, camera: Camera) -> list[Path]:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    signal = signal_from_raw(raw)
    paths = []
    for mode in ("cyan", "inverted"):
        out = image_path(model, camera, mode)
        rgb = fluorescence_transform(signal, mode)
        Image.fromarray(np.uint8(np.clip(rgb, 0.0, 1.0) * 255)).save(out)
        paths.append(out)
    return paths


def draw_contact(paths_by_style: dict[str, list[Path]], *, model: str) -> list[Path]:
    outputs = []
    for style, paths in paths_by_style.items():
        fig, axes = plt.subplots(2, 3, figsize=(12.0, 8.1), constrained_layout=False)
        axes = axes.ravel()
        background = "black" if style == "cyan" else "white"
        text_color = "white" if background == "black" else "black"
        fig.patch.set_facecolor(background)
        for ax, camera, path in zip(axes, CAMERAS, paths):
            ax.set_axis_off()
            ax.set_facecolor(background)
            ax.imshow(mpimg.imread(path))
            ax.text(
                0.02,
                0.98,
                camera.label,
                transform=ax.transAxes,
                ha="left",
                va="top",
                color=text_color,
                fontsize=10.5,
                weight="bold",
            )
        fig.text(
            0.5,
            0.985,
            f"{model.replace('_', ' ')} camera presets, {style}",
            ha="center",
            va="top",
            color=text_color,
            fontsize=14,
        )
        fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.01, wspace=0.02, hspace=0.05)
        out = OUT / f"{model}_camera_contact_{style}.png"
        fig.savefig(out, dpi=360, bbox_inches="tight", facecolor=background)
        plt.close(fig)
        outputs.append(out)
    return outputs


def write_readme(outputs: list[Path]) -> None:
    lines = [
        "# Cytosim Play Camera Tests",
        "",
        "Camera/style contact sheets for initially minus-z-aligned annulus simulations.",
        "",
        f"Frame rendered: `{FRAME}`.",
        "Raw `play` renders were generated using on-screen image mode because off-screen read-pixels failed on this WSL/OpenGL setup.",
        "`play` exits with a segmentation fault after writing images on this setup; the script treats the image file as the success criterion.",
        "",
        "Camera presets:",
    ]
    for cam in CAMERAS:
        lines.append(f"- `{cam.key}`: {cam.label}; rotation=`{cam.rotation or 'default'}`, zoom={cam.zoom}, view_scale={cam.view_scale}")
    lines.extend(["", "Generated files:"])
    for path in outputs:
        lines.append(f"- `{path.relative_to(ROOT)}`")
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    outputs: list[Path] = []
    for model, run_dir in RUNS.items():
        style_paths = {"raw": [], "cyan": [], "inverted": []}
        for camera in CAMERAS:
            raw = run_play(run_dir, model, camera)
            style_paths["raw"].append(raw)
            cyan, inverted = process_image(raw, model, camera)
            style_paths["cyan"].append(cyan)
            style_paths["inverted"].append(inverted)
        outputs.extend(draw_contact(style_paths, model=model))
    write_readme(outputs)
    print(OUT)


if __name__ == "__main__":
    main()
