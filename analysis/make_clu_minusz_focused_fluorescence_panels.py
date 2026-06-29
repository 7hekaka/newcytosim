#!/usr/bin/env python3
from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import gaussian_filter


ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "analysis") not in sys.path:
    sys.path.insert(0, str(ROOT / "analysis"))

from analysis import build_clu_minusz_slide_graphs as slides
from analysis import comparison_unwrapped_annulus_story as story


OUT = ROOT / "analysis" / "results" / "clu_minusz_slides_2026-06-02" / "fluorescence_focused_1to16"
REGIME = "1:16"
CASES = {
    "total480": [
        ("c20_m24", "20 clusters"),
        ("c40_m12", "40 clusters"),
        ("c80_m6", "80 clusters"),
    ],
    "mpc12": [
        ("c20_m12", "20 clusters"),
        ("c40_m12", "40 clusters"),
        ("c80_m12", "80 clusters"),
    ],
}
GROUP_LABELS = {
    "total480": "Total motors = 480",
    "mpc12": "12 motors per cluster",
}
MODEL_LABELS = {
    "rotatable": "rotatable",
    "fixed_global": "fixed orientation",
}
STYLE_LABELS = {
    "inverted_fluorescence": "inverted fluorescence",
    "cyan_fluorescence": "cyan fluorescence",
    "line_render": "line render",
}


def line_intensity(panel: story.Panel, *, width: int = 1500, height: int = 450, supersample: int = 4) -> np.ndarray:
    w = width * supersample
    h = height * supersample
    canvas = Image.new("L", (w, h), 0)
    draw = ImageDraw.Draw(canvas)
    r_mid = panel.smax / (2.0 * math.pi)

    for s, z in story.frame_segments(panel.frame, panel.seam, panel.smax, r_mid):
        if len(s) < 2:
            continue
        x = np.clip(s / panel.smax, 0.0, 1.0) * (w - 1)
        y = np.clip((z - panel.zmin) / (panel.zmax - panel.zmin), 0.0, 1.0) * (h - 1)
        draw.line(
            list(zip(x.astype(float), y.astype(float))),
            fill=255,
            width=max(1, int(round(1.15 * supersample))),
        )

    arr = np.asarray(canvas, dtype=float) / 255.0
    fine = gaussian_filter(arr, sigma=0.55 * supersample)
    glow = gaussian_filter(arr, sigma=3.6 * supersample)
    img = fine + 0.40 * glow
    img = img / max(float(img.max()), 1e-12)
    resized = Image.fromarray(np.uint8(np.clip(img, 0.0, 1.0) * 255)).resize((width, height), Image.Resampling.LANCZOS)
    return np.asarray(resized, dtype=float) / 255.0


def clean_line_rgb(panel: story.Panel, *, width: int = 1500, height: int = 450, supersample: int = 4) -> np.ndarray:
    w = width * supersample
    h = height * supersample
    canvas = Image.new("L", (w, h), 255)
    draw = ImageDraw.Draw(canvas)
    r_mid = panel.smax / (2.0 * math.pi)

    for s, z in story.frame_segments(panel.frame, panel.seam, panel.smax, r_mid):
        if len(s) < 2:
            continue
        x = np.clip(s / panel.smax, 0.0, 1.0) * (w - 1)
        y = np.clip((z - panel.zmin) / (panel.zmax - panel.zmin), 0.0, 1.0) * (h - 1)
        draw.line(
            list(zip(x.astype(float), y.astype(float))),
            fill=20,
            width=max(1, int(round(0.85 * supersample))),
        )

    img = Image.fromarray(np.asarray(canvas, dtype=np.uint8)).resize((width, height), Image.Resampling.LANCZOS)
    arr = np.asarray(img, dtype=float) / 255.0
    return np.dstack([arr, arr, arr])


def normalize_intensity(images: list[np.ndarray], *, gamma: float = 0.72, clip_quantile: float = 0.995) -> list[np.ndarray]:
    positive = np.concatenate([img[img > 0].ravel() for img in images if np.any(img > 0)])
    scale = float(np.quantile(positive, clip_quantile)) if positive.size else 1.0
    scale = max(scale, 1e-12)
    out = []
    for img in images:
        val = np.clip(img / scale, 0.0, 1.0)
        out.append(np.power(val, gamma))
    return out


def add_noise(img: np.ndarray, *, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return np.clip(img + rng.normal(0.0, 0.010, img.shape), 0.0, 1.0)


def fluorescence_rgb(img: np.ndarray, style: str, *, seed: int) -> np.ndarray:
    img = add_noise(img, seed=seed)
    if style == "inverted_fluorescence":
        gray = 1.0 - 0.92 * img
        return np.dstack([gray, gray, gray])
    if style == "cyan_fluorescence":
        return np.dstack([0.05 * img, 0.95 * img, 1.00 * img])
    raise KeyError(style)


def build_image_grid(
    lookup: dict[tuple[str, str, str, str], story.Panel],
    *,
    group: str,
    model: str,
    style: str,
) -> tuple[list[list[np.ndarray]], list[list[story.Panel]]]:
    panels: list[list[story.Panel]] = []
    raw: list[np.ndarray] = []
    for case, _case_label in CASES[group]:
        panel_row: list[story.Panel] = []
        for snapshot_slug, _fraction, _snapshot_label in story.SNAPSHOTS:
            panel = lookup[(REGIME, model, case, snapshot_slug)]
            panel_row.append(panel)
            if style != "line_render":
                raw.append(line_intensity(panel))
        panels.append(panel_row)

    if style == "line_render":
        return (
            [[clean_line_rgb(panel) for panel in row] for row in panels],
            panels,
        )

    normalized = iter(normalize_intensity(raw))
    seed_base = abs(hash((group, model, style))) % 100000
    images: list[list[np.ndarray]] = []
    for ridx, row in enumerate(panels):
        image_row: list[np.ndarray] = []
        for cidx, _panel in enumerate(row):
            image_row.append(fluorescence_rgb(next(normalized), style, seed=seed_base + 100 * ridx + cidx))
        images.append(image_row)
    return images, panels


def draw_figure(
    images: list[list[np.ndarray]],
    panels: list[list[story.Panel]],
    *,
    group: str,
    model: str,
    style: str,
) -> Path:
    nrows = len(CASES[group])
    ncols = len(story.SNAPSHOTS)
    background = "black" if style == "cyan_fluorescence" else "white"
    text_color = "white" if background == "black" else "black"

    fig, axes = plt.subplots(nrows, ncols, figsize=(14.8, 5.1), constrained_layout=False)
    axes = np.asarray(axes)
    fig.patch.set_facecolor(background)

    for ridx, (_case, case_label) in enumerate(CASES[group]):
        for cidx, (_snapshot_slug, _fraction, _snapshot_label) in enumerate(story.SNAPSHOTS):
            ax = axes[ridx, cidx]
            ax.set_axis_off()
            ax.set_facecolor(background)
            ax.imshow(images[ridx][cidx], interpolation="nearest")
            if ridx == 0:
                ax.text(
                    0.5,
                    1.035,
                    f"{panels[ridx][cidx].time_min:.1f} min",
                    transform=ax.transAxes,
                    ha="center",
                    va="bottom",
                    color=text_color,
                    fontsize=12.5,
                )
            if cidx == 0:
                ax.text(
                    0.015,
                    0.93,
                    case_label,
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    color=text_color,
                    fontsize=10.5,
                    weight="bold",
                )

    title = f"{GROUP_LABELS[group]}, {MODEL_LABELS[model]}, {REGIME} | {STYLE_LABELS[style]}"
    fig.text(0.5, 0.987, title, ha="center", va="top", color=text_color, fontsize=15.0)
    fig.subplots_adjust(left=0.010, right=0.990, top=0.900, bottom=0.030, wspace=0.018, hspace=0.085)

    OUT.mkdir(parents=True, exist_ok=True)
    stem = f"{group}_{model}_1to16_20_40_80_{style}"
    png = OUT / f"{stem}.png"
    pdf = OUT / f"{stem}.pdf"
    fig.savefig(png, dpi=600, bbox_inches="tight", facecolor=background)
    fig.savefig(pdf, bbox_inches="tight", facecolor=background)
    plt.close(fig)
    return png


def write_readme(paths: list[Path]) -> None:
    lines = [
        "# Focused 1:16 Fluorescence-Style Unwrapped Annulus Panels",
        "",
        "Focused visualization panels for initially minus-z-aligned annulus simulations.",
        "",
        "Selection:",
        "- Crosslinker ratio: `1:16` only.",
        "- Cluster counts: `20`, `40`, and `80` only.",
        "- Rotatable and fixed-orientation conditions are plotted separately.",
        "- Families: `total480` and `mpc12`.",
        "",
        "Styles:",
        "- `inverted_fluorescence`: line-derived black/gray intensity on white background.",
        "- `cyan_fluorescence`: line-derived cyan intensity on black background.",
        "- `line_render`: clean black filament lines on white background.",
        "",
        "Generated files:",
    ]
    for path in paths:
        lines.append(f"- `{path.relative_to(ROOT)}`")
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    outputs: list[Path] = []
    slides.configure_story_module()
    family_panels, _ = story.collect_panels()
    for group in ("total480", "mpc12"):
        lookup = story.panel_lookup(family_panels[group])
        for model in ("rotatable", "fixed_global"):
            for style in ("inverted_fluorescence", "cyan_fluorescence", "line_render"):
                images, panels = build_image_grid(lookup, group=group, model=model, style=style)
                outputs.append(draw_figure(images, panels, group=group, model=model, style=style))
    write_readme(outputs)
    print(OUT)


if __name__ == "__main__":
    main()
