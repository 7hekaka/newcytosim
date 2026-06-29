#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "analysis") not in sys.path:
    sys.path.insert(0, str(ROOT / "analysis"))

from analysis import comparison_unwrapped_annulus_story as story

OUT = ROOT / "analysis" / "results" / "presentation_ready_2026-04-06" / "slide_assets"
HEATMAP_OUT = OUT / "fixed_vs_rotatable_unwrapped_slides"
NOTES_PATH = OUT / "slide_notes.md"

FAMILY_ORDER = ("total480", "mpc12")
REGIMES = story.XLINK_ORDER


def family_specs() -> dict[str, dict]:
    return {spec["group"]: spec for spec in story.base.FAMILY_SPECS if spec["group"] in FAMILY_ORDER}


def lookup_panel_map(panels: list[story.Panel]) -> dict[tuple[str, str, str, str], story.Panel]:
    return story.panel_lookup(panels)


def draw_mini_timeline(
    fig: plt.Figure,
    subplot_spec,
    lookup: dict,
    regime: str,
    row_spec: dict,
    vmax: float,
    *,
    row_label: bool = False,
    show_titles: bool = False,
) -> None:
    sub = subplot_spec.subgridspec(1, len(story.SNAPSHOTS), wspace=0.02)
    for j, (snapshot_slug, _fraction, snapshot_label) in enumerate(story.SNAPSHOTS):
        ax = fig.add_subplot(sub[0, j])
        panel = lookup.get((regime, row_spec["model"], row_spec["case"], snapshot_slug))
        story.draw_density_panel(ax, panel, vmax=vmax)
        if show_titles:
            title = snapshot_label
            if panel is not None:
                title = f"{snapshot_label}\n{panel.time_min:.1f} min"
            ax.set_title(title, fontsize=8.5, pad=3)
        else:
            ax.set_title("")
        if j == 0 and row_label:
            ax.set_ylabel("z (um)", fontsize=8.5)
            ax.text(
                -0.42,
                0.50,
                row_spec["label"],
                transform=ax.transAxes,
                rotation=90,
                ha="center",
                va="center",
                fontsize=9,
                weight="bold",
            )
        else:
            ax.set_ylabel("")
            ax.set_yticklabels([])
        if j == 1:
            ax.set_xlabel("Time after motor onset", fontsize=8.2, labelpad=2)
        else:
            ax.set_xlabel("")
        if j != len(story.SNAPSHOTS) - 1:
            ax.set_xticklabels([])
        ax.tick_params(labelsize=7, width=0.7, length=2.5)


def make_slide_heatmap_figure(family: str, regime: str, panels: list[story.Panel], vmax: float) -> Path:
    specs = family_specs()
    family_spec = specs[family]
    lookup = lookup_panel_map(panels)
    top_counts = family_spec["order"][:3]
    bottom_counts = family_spec["order"][3:]

    fig = plt.figure(figsize=(13.333, 7.5), facecolor="white")
    outer = fig.add_gridspec(
        4,
        2,
        width_ratios=[5.2, 1.1],
        height_ratios=[1.0, 2.0, 2.0, 0.25],
        left=0.04,
        right=0.98,
        top=0.94,
        bottom=0.06,
        wspace=0.10,
        hspace=0.16,
    )

    top = outer[0, 0].subgridspec(1, 1)
    mid = outer[1, 0].subgridspec(2, 3, wspace=0.10, hspace=0.18)
    bottom = outer[2, 0].subgridspec(2, 2, wspace=0.10, hspace=0.18)
    notes_ax = fig.add_subplot(outer[:3, 1])
    notes_ax.axis("off")

    control_spec = {"model": "control", "group": "control", "case": "c0_m0", "label": "No motors"}
    draw_mini_timeline(fig, top[0, 0], lookup, regime, control_spec, vmax, row_label=True, show_titles=True)

    title_prefix = "Total motors = 480" if family == "total480" else "12 motors per cluster"
    fig.text(0.06, 0.965, f"{title_prefix}   |   xlink ratio {regime}", ha="left", va="top", fontsize=18, weight="bold")

    for col, (_case, label) in enumerate(top_counts):
        fig.text(0.15 + 0.18 * col, 0.71, label, ha="center", va="bottom", fontsize=11.5, weight="bold")
    for col, (_case, label) in enumerate(bottom_counts):
        fig.text(0.22 + 0.28 * col, 0.36, label, ha="center", va="bottom", fontsize=11.5, weight="bold")

    for col, (case, _label) in enumerate(top_counts):
        draw_mini_timeline(
            fig,
            mid[0, col],
            lookup,
            regime,
            {"model": "rotatable", "group": family, "case": case, "label": "Rotatable"},
            vmax,
            row_label=(col == 0),
        )
        draw_mini_timeline(
            fig,
            mid[1, col],
            lookup,
            regime,
            {"model": "fixed_global", "group": family, "case": case, "label": "Fixed-global"},
            vmax,
            row_label=(col == 0),
        )

    for col, (case, _label) in enumerate(bottom_counts):
        draw_mini_timeline(
            fig,
            bottom[0, col],
            lookup,
            regime,
            {"model": "rotatable", "group": family, "case": case, "label": "Rotatable"},
            vmax,
            row_label=(col == 0),
        )
        draw_mini_timeline(
            fig,
            bottom[1, col],
            lookup,
            regime,
            {"model": "fixed_global", "group": family, "case": case, "label": "Fixed-global"},
            vmax,
            row_label=(col == 0),
        )

    notes_ax.text(0.02, 0.98, "Comments", ha="left", va="top", fontsize=14, weight="bold")
    notes_ax.text(
        0.02,
        0.90,
        "- Compare each fixed row directly against the rotatable row above it.\n"
        "- Use the control row as the no-motor baseline texture.\n"
        "- Look for sharpening, drift, or band breakup as time progresses.",
        ha="left",
        va="top",
        fontsize=11,
        linespacing=1.7,
    )

    cax = fig.add_axes([0.91, 0.18, 0.015, 0.62])
    dummy = np.linspace(0, vmax, 64).reshape(-1, 1)
    im = cax.imshow(dummy, aspect="auto", cmap=story.HEATMAP_CMAP, origin="lower", vmin=0.0, vmax=vmax)
    cax.set_visible(False)
    cb_ax = fig.add_axes([0.92, 0.18, 0.015, 0.62])
    cb = fig.colorbar(im, cax=cb_ax)
    cb.set_label("Normalized areal density", fontsize=11)
    cb.ax.tick_params(labelsize=9)

    HEATMAP_OUT.mkdir(parents=True, exist_ok=True)
    path = HEATMAP_OUT / f"{family}_{story.XLINK_SLUG[regime]}_slide_layout.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def write_notes(paths: list[Path]) -> None:
    lines = [
        "# Slide Notes",
        "",
        "## Fixed vs rotatable unwrapped annulus heatmaps",
        "- Use one slide per xlink ratio and family.",
        "- Control sits at the top as the no-motor baseline texture.",
        "- Rotatable should be read against fixed-global row-by-row at the same cluster count.",
        "- The 10/20/40 block is placed above the 60/80 block to keep each panel large enough for slides.",
        "",
        "## Transport metric slides",
        "- Give each metric its own slide when the comparison is obvious.",
        "- Recommended bullets:",
        "- Mean |v_z|: headline transport magnitude.",
        "- Peak |v_z|: strongest transport state reached.",
        "- AUC(|v_z|): total transport over the motor-on window.",
        "- Active transport fraction: how often the condition is truly above control-like noise.",
        "",
        "## Motor-force sweep slides",
        "- Use one slide per xlink ratio.",
        "- The 5x4 summary plots already work well as a compact family view.",
        "- Pair the summary plot with one sentence about low-force cases underperforming the no-motor control and one sentence about the rescue at higher unbinding force / higher motors per cluster.",
        "",
        "## Generated slide assets",
    ]
    lines.extend(f"- {path.relative_to(ROOT)}" for path in paths)
    NOTES_PATH.parent.mkdir(parents=True, exist_ok=True)
    NOTES_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.linewidth": 0.9,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )
    family_panels, family_vmax = story.collect_panels()
    made: list[Path] = []
    for family in FAMILY_ORDER:
        panels = family_panels[family]
        vmax = family_vmax[family]
        for regime in REGIMES:
            made.append(make_slide_heatmap_figure(family, regime, panels, vmax))
    write_notes(made)


if __name__ == "__main__":
    main()
