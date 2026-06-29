#!/usr/bin/env python3
"""Summarize movie-level F-actin flow results by line/genotype.

Input is the `movie_summary.csv` file written by
`analyze_egg_cell_network_flow.py`.

The replicate unit is one movie. Pairwise p-values are exact permutation tests
when the number of label permutations is modest; otherwise the script uses a
fixed-seed Monte Carlo permutation test.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
from pathlib import Path
import random
import statistics
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_movies(path: Path, metric: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if not row.get("line"):
                raise ValueError("movie_summary.csv must contain a `line` column.")
            if metric not in row or row[metric] in ("", "None"):
                raise ValueError(f"Metric `{metric}` is missing for at least one movie.")
            row["metric_value"] = float(row[metric])
            rows.append(row)
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return rows


def line_order(rows: list[dict[str, Any]], preferred: list[str] | None = None) -> list[str]:
    seen = []
    for row in rows:
        line = row["line"]
        if line not in seen:
            seen.append(line)
    if not preferred:
        return seen
    ordered = [line for line in preferred if line in seen]
    ordered.extend(line for line in seen if line not in ordered)
    return ordered


def exact_permutation_p(a: list[float], b: list[float]) -> tuple[float, str]:
    obs = abs(statistics.mean(a) - statistics.mean(b))
    combined = a + b
    n = len(a)
    total_combos = math.comb(len(combined), n)
    if total_combos > 200000:
        return monte_carlo_permutation_p(a, b), "monte_carlo_100000_seed_1"

    ge = 0
    total = 0
    for idxs in itertools.combinations(range(len(combined)), n):
        idx = set(idxs)
        aa = [combined[i] for i in range(len(combined)) if i in idx]
        bb = [combined[i] for i in range(len(combined)) if i not in idx]
        diff = abs(statistics.mean(aa) - statistics.mean(bb))
        ge += diff >= obs - 1e-12
        total += 1
    return ge / total, "exact"


def monte_carlo_permutation_p(a: list[float], b: list[float], samples: int = 100000) -> float:
    rng = random.Random(1)
    obs = abs(statistics.mean(a) - statistics.mean(b))
    combined = a + b
    n = len(a)
    ge = 0
    for _ in range(samples):
        shuffled = combined[:]
        rng.shuffle(shuffled)
        aa = shuffled[:n]
        bb = shuffled[n:]
        ge += abs(statistics.mean(aa) - statistics.mean(bb)) >= obs - 1e-12
    return (ge + 1) / (samples + 1)


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def summarize_lines(rows: list[dict[str, Any]], order: list[str]) -> list[dict[str, Any]]:
    out = []
    for line in order:
        vals = [row["metric_value"] for row in rows if row["line"] == line]
        out.append(
            {
                "line": line,
                "n_movies": len(vals),
                "mean": statistics.mean(vals),
                "sd": statistics.stdev(vals) if len(vals) > 1 else 0.0,
                "sem": statistics.stdev(vals) / math.sqrt(len(vals)) if len(vals) > 1 else 0.0,
                "median": statistics.median(vals),
                "min": min(vals),
                "max": max(vals),
            }
        )
    return out


def pairwise_tests(rows: list[dict[str, Any]], order: list[str]) -> list[dict[str, Any]]:
    out = []
    for idx, a_line in enumerate(order):
        for b_line in order[idx + 1 :]:
            a = [row["metric_value"] for row in rows if row["line"] == a_line]
            b = [row["metric_value"] for row in rows if row["line"] == b_line]
            p_value, method = exact_permutation_p(a, b)
            out.append(
                {
                    "comparison": f"{a_line} vs {b_line}",
                    "n_a": len(a),
                    "n_b": len(b),
                    "mean_a": statistics.mean(a),
                    "mean_b": statistics.mean(b),
                    "mean_diff_b_minus_a": statistics.mean(b) - statistics.mean(a),
                    "two_sided_permutation_p_mean_diff": p_value,
                    "test_method": method,
                }
            )
    return out


def plot_lines(rows: list[dict[str, Any]], order: list[str], output: Path, ylabel: str) -> None:
    palette = [
        "#e76f51",
        "#2a9d8f",
        "#8ab17d",
        "#4cc9f0",
        "#b56576",
        "#6d597a",
        "#f4a261",
        "#457b9d",
    ]
    fig, ax = plt.subplots(figsize=(max(6.5, 1.25 * len(order)), 4.8), constrained_layout=True)
    for i, line in enumerate(order):
        vals = np.array([row["metric_value"] for row in rows if row["line"] == line], dtype=float)
        color = palette[i % len(palette)]
        ax.boxplot(
            [vals],
            positions=[i],
            widths=0.5,
            showfliers=False,
            patch_artist=True,
            boxprops={"facecolor": color, "alpha": 0.35, "edgecolor": "black"},
            medianprops={"color": "black", "linewidth": 2},
            whiskerprops={"color": "black"},
            capprops={"color": "black"},
        )
        jitter = np.linspace(-0.16, 0.16, len(vals)) if len(vals) > 1 else np.array([0.0])
        ax.scatter(
            np.full(len(vals), i) + jitter,
            vals,
            s=42,
            color=color,
            edgecolor="black",
            linewidth=0.5,
            zorder=3,
        )
        ax.scatter(i, vals.mean(), marker="D", s=64, color="white", edgecolor="black", zorder=4)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=25, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title("Network-scale F-actin flow, movie-level replicates")
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(output, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("movie_summary", help="Path to movie_summary.csv from analyze_egg_cell_network_flow.py")
    parser.add_argument("--output", default=None, help="Output folder. Defaults beside movie_summary.csv.")
    parser.add_argument(
        "--metric",
        default="median_speed_um_per_s_assumed_dt",
        help="Metric column to summarize.",
    )
    parser.add_argument(
        "--line-order",
        default="",
        help="Optional comma-separated line order, e.g. 'Wild type,CA ROP9,DN ROP9,scar4-1'.",
    )
    parser.add_argument("--ylabel", default="Movie median F-actin speed (um/s)")
    args = parser.parse_args()

    movie_summary = Path(args.movie_summary)
    outdir = Path(args.output) if args.output else movie_summary.parent / "summary_by_line"
    rows = read_movies(movie_summary, args.metric)
    preferred = [item.strip() for item in args.line_order.split(",") if item.strip()]
    order = line_order(rows, preferred)

    line_rows = summarize_lines(rows, order)
    pair_rows = pairwise_tests(rows, order)
    write_csv(outdir / "line_summary.csv", line_rows)
    write_csv(
        outdir / "pairwise_permutation_tests.csv",
        pair_rows,
        [
            "comparison",
            "n_a",
            "n_b",
            "mean_a",
            "mean_b",
            "mean_diff_b_minus_a",
            "two_sided_permutation_p_mean_diff",
            "test_method",
        ],
    )
    plot_lines(rows, order, outdir / "line_summary_plot.png", args.ylabel)

    print(outdir / "line_summary.csv")
    print(outdir / "pairwise_permutation_tests.csv")
    print(outdir / "line_summary_plot.png")


if __name__ == "__main__":
    main()
