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


OUT = ROOT / "analysis" / "results" / "clu_minusz_slides_2026-06-02" / "fluorescence_style_prototypes"

STYLE_SPECS = {
    "cyan_confocal": {
        "label": "cyan confocal-style intensity",
        "background": "black",
        "text_color": "white",
        "mode": "line",
    },
    "inverted_fluorescence": {
        "label": "inverted fluorescence-style intensity",
        "background": "white",
        "text_color": "black",
        "mode": "line",
    },
    "density_cyan": {
        "label": "density-derived cyan intensity",
        "background": "black",
        "text_color": "white",
        "mode": "density",
    },
}

PROTOTYPES = [
    ("total480", "1:16", "high_counts", "total480_1to16_high_counts"),
    ("mpc12", "1:16", "high_counts", "mpc12_1to16_high_counts"),
    ("total480", "1:4", "low_counts", "total480_1to4_low_counts"),
]


def family_spec(group: str) -> dict:
    for spec in slides.family_specs():
        if spec["group"] == group:
            return spec
    raise KeyError(group)


def row_specs(spec: dict, part: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for case, label in slides.split_cases(spec, part):
        rows.append({"model": "rotatable", "case": case, "label": f"{label}\nrotatable"})
        rows.append({"model": "fixed_global", "case": case, "label": f"{label}\nfixed orientation"})
    return rows


def panel_key(row: dict[str, str], regime: str, snapshot_slug: str) -> tuple[str, str, str, str]:
    return (regime, row["model"], row["case"], snapshot_slug)


def raw_line_intensity(panel: story.Panel, *, width: int = 900, height: int = 260, supersample: int = 3) -> np.ndarray:
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
        points = list(zip(x.astype(float), y.astype(float)))
        draw.line(points, fill=255, width=max(1, int(round(1.2 * supersample))))

    arr = np.asarray(canvas, dtype=float) / 255.0
    fine = gaussian_filter(arr, sigma=0.65 * supersample)
    glow = gaussian_filter(arr, sigma=4.0 * supersample)
    img = fine + 0.48 * glow
    img = np.asarray(Image.fromarray(np.uint8(np.clip(img / max(img.max(), 1e-12), 0, 1) * 255)).resize((width, height), Image.Resampling.LANCZOS), dtype=float) / 255.0
    return img


def raw_density_intensity(panel: story.Panel, *, width: int = 900, height: int = 260) -> np.ndarray:
    density = np.asarray(panel.density, dtype=float)
    density = np.clip(density, 0.0, None)
    if density.max() > 0:
        density = density / density.max()
    img = Image.fromarray(np.uint8(density * 255)).resize((width, height), Image.Resampling.BICUBIC)
    arr = np.asarray(img, dtype=float) / 255.0
    return gaussian_filter(arr, sigma=1.1)


def normalize_images(raw_images: list[np.ndarray], *, gamma: float = 0.72, clip_quantile: float = 0.995) -> list[np.ndarray]:
    positive = np.concatenate([img[img > 0].ravel() for img in raw_images if np.any(img > 0)])
    scale = float(np.quantile(positive, clip_quantile)) if positive.size else 1.0
    scale = max(scale, 1e-12)
    normalized = []
    for img in raw_images:
        val = np.clip(img / scale, 0.0, 1.0)
        normalized.append(np.power(val, gamma))
    return normalized


def add_microscope_noise(img: np.ndarray, *, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    noisy = img + rng.normal(0.0, 0.012, img.shape)
    background = rng.normal(0.0, 0.004, img.shape)
    return np.clip(noisy + background, 0.0, 1.0)


def to_rgb(img: np.ndarray, style_key: str, *, seed: int) -> np.ndarray:
    img = add_microscope_noise(img, seed=seed)
    if style_key == "inverted_fluorescence":
        gray = 1.0 - 0.92 * img
        return np.dstack([gray, gray, gray])
    if style_key in {"cyan_confocal", "density_cyan"}:
        return np.dstack([0.05 * img, 0.95 * img, 1.00 * img])
    raise KeyError(style_key)


def collect_grid_images(
    panels: dict[tuple[str, str, str, str], story.Panel],
    rows: list[dict[str, str]],
    regime: str,
    style_key: str,
) -> tuple[list[list[np.ndarray | None]], list[list[story.Panel | None]]]:
    raw_grid: list[list[np.ndarray | None]] = []
    panel_grid: list[list[story.Panel | None]] = []
    raw_images: list[np.ndarray] = []
    mode = STYLE_SPECS[style_key]["mode"]

    for row in rows:
        raw_row: list[np.ndarray | None] = []
        panel_row: list[story.Panel | None] = []
        for snapshot_slug, _fraction, _snapshot_label in story.SNAPSHOTS:
            panel = panels.get(panel_key(row, regime, snapshot_slug))
            panel_row.append(panel)
            if panel is None:
                raw_row.append(None)
                continue
            raw = raw_density_intensity(panel) if mode == "density" else raw_line_intensity(panel)
            raw_row.append(raw)
            raw_images.append(raw)
        raw_grid.append(raw_row)
        panel_grid.append(panel_row)

    normalized = normalize_images(raw_images)
    it = iter(normalized)
    image_grid: list[list[np.ndarray | None]] = []
    seed_base = abs(hash((regime, style_key))) % 100000
    for ridx, raw_row in enumerate(raw_grid):
        image_row: list[np.ndarray | None] = []
        for cidx, raw in enumerate(raw_row):
            if raw is None:
                image_row.append(None)
            else:
                image_row.append(to_rgb(next(it), style_key, seed=seed_base + 100 * ridx + cidx))
        image_grid.append(image_row)
    return image_grid, panel_grid


def draw_grid(
    image_grid: list[list[np.ndarray | None]],
    panel_grid: list[list[story.Panel | None]],
    rows: list[dict[str, str]],
    *,
    group_label: str,
    regime: str,
    part: str,
    style_key: str,
    outstem: str,
) -> Path:
    nrows = len(rows)
    ncols = len(story.SNAPSHOTS)
    fig, axes = plt.subplots(nrows, ncols, figsize=(14.2, 1.75 * nrows), constrained_layout=False)
    axes = np.asarray(axes)
    if nrows == 1:
        axes = axes[np.newaxis, :]

    spec = STYLE_SPECS[style_key]
    fig.patch.set_facecolor(spec["background"])
    for ridx, row in enumerate(rows):
        for cidx, (snapshot_slug, _fraction, snapshot_label) in enumerate(story.SNAPSHOTS):
            ax = axes[ridx, cidx]
            ax.set_facecolor(spec["background"])
            ax.set_axis_off()
            img = image_grid[ridx][cidx]
            panel = panel_grid[ridx][cidx]
            if img is not None:
                ax.imshow(img, interpolation="nearest")
            if ridx == 0:
                time_label = snapshot_label
                if panel is not None and np.isfinite(panel.time_min):
                    time_label = f"{panel.time_min:.1f} min"
                ax.text(
                    0.5,
                    1.035,
                    time_label,
                    transform=ax.transAxes,
                    ha="center",
                    va="bottom",
                    color=spec["text_color"],
                    fontsize=11.5,
                )
            if cidx == 0:
                ax.text(
                    0.012,
                    0.94,
                    row["label"],
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    color=spec["text_color"],
                    fontsize=9.8,
                    weight="bold",
                )

    fig.text(
        0.5,
        0.986,
        f"{group_label}, {regime}, {part.replace('_', ' ')} | {spec['label']}",
        ha="center",
        va="top",
        color=spec["text_color"],
        fontsize=14,
    )
    fig.subplots_adjust(left=0.012, right=0.988, top=0.92, bottom=0.02, wspace=0.018, hspace=0.10)
    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / f"{outstem}_{style_key}.png"
    pdf = OUT / f"{outstem}_{style_key}.pdf"
    fig.savefig(png, dpi=420, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(pdf, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return png


def draw_style_comparison(
    panels: dict[tuple[str, str, str, str], story.Panel],
    *,
    regime: str = "1:16",
    row: dict[str, str],
    outstem: str = "single_condition_style_comparison",
) -> Path:
    styles = ["inverted_fluorescence", "cyan_confocal", "density_cyan"]
    short_labels = {
        "inverted_fluorescence": "Inverted\nline render",
        "cyan_confocal": "Cyan\nline render",
        "density_cyan": "Cyan\ndensity render",
    }
    fig, axes = plt.subplots(len(styles), len(story.SNAPSHOTS), figsize=(13.8, 5.5), constrained_layout=False)
    fig.patch.set_facecolor("white")
    for sidx, style_key in enumerate(styles):
        image_grid, panel_grid = collect_grid_images(panels, [row], regime, style_key)
        for cidx, (_snapshot_slug, _fraction, snapshot_label) in enumerate(story.SNAPSHOTS):
            ax = axes[sidx, cidx]
            ax.set_axis_off()
            ax.imshow(image_grid[0][cidx])
            panel = panel_grid[0][cidx]
            if sidx == 0:
                time_label = snapshot_label
                if panel is not None and np.isfinite(panel.time_min):
                    time_label = f"{panel.time_min:.1f} min"
                ax.text(0.5, 1.04, time_label, transform=ax.transAxes, ha="center", va="bottom", fontsize=11.5)
            if cidx == 0:
                ax.text(
                    -0.055,
                    0.5,
                    short_labels[style_key],
                    transform=ax.transAxes,
                    ha="right",
                    va="center",
                    rotation=0,
                    fontsize=10.8,
                )
    fig.text(0.5, 0.98, f"Style comparison: {row['label'].replace(chr(10), ', ')} at {regime}", ha="center", va="top", fontsize=13.5)
    fig.subplots_adjust(left=0.12, right=0.99, top=0.88, bottom=0.04, wspace=0.025, hspace=0.11)
    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / f"{outstem}.png"
    pdf = OUT / f"{outstem}.pdf"
    fig.savefig(png, dpi=420, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png


def write_readme(paths: list[Path]) -> None:
    lines = [
        "# Fluorescence-Style CLU Minus-Z Prototypes",
        "",
        "These images are visualization prototypes generated from the unwrapped Cytosim filament panels.",
        "They are not additional analysis metrics. A fixed microscopy-style transfer function is applied to all panels within each figure.",
        "",
        "Rendering modes:",
        "- `inverted_fluorescence`: filament geometry rasterized, blurred, gamma-compressed, lightly noised, and inverted to black/gray signal on white background.",
        "- `cyan_confocal`: filament geometry rasterized, blurred, gamma-compressed, lightly noised, and displayed as cyan signal on black background.",
        "- `density_cyan`: precomputed unwrapped filament density converted to cyan fluorescence-style intensity.",
        "",
        "Source modules:",
        f"- `{Path(__file__).relative_to(ROOT)}`",
        "- `analysis/comparison_unwrapped_annulus_story.py`",
        "- `analysis/build_clu_minusz_slide_graphs.py`",
        "",
        "Generated files:",
    ]
    for path in paths:
        lines.append(f"- `{path.relative_to(ROOT)}`")
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    slides.configure_story_module()
    family_panels, _family_vmax = story.collect_panels()
    outputs: list[Path] = []

    for group, regime, part, outstem in PROTOTYPES:
        spec = family_spec(group)
        rows = row_specs(spec, part)
        lookup = story.panel_lookup(family_panels[group])
        group_label = "Total motors = 480" if group == "total480" else "12 motors per cluster"
        for style_key in ("inverted_fluorescence", "cyan_confocal", "density_cyan"):
            image_grid, panel_grid = collect_grid_images(lookup, rows, regime, style_key)
            outputs.append(
                draw_grid(
                    image_grid,
                    panel_grid,
                    rows,
                    group_label=group_label,
                    regime=regime,
                    part=part,
                    style_key=style_key,
                    outstem=outstem,
                )
            )

    comparison_lookup = story.panel_lookup(family_panels["total480"])
    outputs.append(
        draw_style_comparison(
            comparison_lookup,
            regime="1:16",
            row={"model": "fixed_global", "case": "c80_m6", "label": "80 clusters\nfixed orientation"},
        )
    )
    write_readme(outputs)
    print(OUT)


if __name__ == "__main__":
    main()
