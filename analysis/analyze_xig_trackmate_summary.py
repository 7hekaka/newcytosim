#!/usr/bin/env python3
"""Analyze collaborator TrackMate XIG track-summary workbooks.

The attached workbook is an xlsx file but this script intentionally avoids
openpyxl/pandas so it can run in the current lightweight environment.
It parses the sheet XML directly, normalizes TrackMate summary rows, and
compares WT versus abcc8 at both track and movie/file levels.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import random
import statistics
import sys
import xml.etree.ElementTree as ET
import zipfile


NS = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}


def cell_col(ref: str) -> int:
    out = 0
    for ch in "".join(c for c in ref if c.isalpha()):
        out = out * 26 + ord(ch.upper()) - 64
    return out


def cell_row(ref: str) -> int:
    digits = "".join(c for c in ref if c.isdigit())
    return int(digits or 0)


def xml_text(node: ET.Element | None) -> str:
    if node is None:
        return ""
    return "".join(t.text or "" for t in node.findall(".//m:t", NS))


def parse_cell_value(cell: ET.Element, shared_strings: list[str]) -> object | None:
    typ = cell.attrib.get("t")
    value = cell.find("m:v", NS)
    if typ == "s" and value is not None:
        return shared_strings[int(value.text)]
    if typ == "inlineStr":
        return xml_text(cell.find("m:is", NS))
    if value is None:
        return None
    text = value.text
    try:
        number = float(text)
        return int(number) if number.is_integer() else number
    except (TypeError, ValueError):
        return text


def parse_workbook(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with zipfile.ZipFile(path) as zf:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            shared_strings = [xml_text(si) for si in root.findall("m:si", NS)]

        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        relmap = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}

        for sheet in workbook.findall(".//m:sheet", NS):
            sheet_name = sheet.attrib["name"]
            rel_id = sheet.attrib[f"{{{NS['r']}}}id"]
            sheet_path = "xl/" + relmap[rel_id].lstrip("/")
            root = ET.fromstring(zf.read(sheet_path))

            cells: dict[tuple[int, int], object] = {}
            max_row = 0
            for cell in root.findall(".//m:c", NS):
                ref = cell.attrib["r"]
                row = cell_row(ref)
                col = cell_col(ref)
                max_row = max(max_row, row)
                cells[(row, col)] = parse_cell_value(cell, shared_strings)

            # TrackMate headers occupy row 1. The collaborator inserted
            # Sample/File columns before the data, so TrackMate data begin in col 3.
            headers = [cells.get((1, col)) for col in range(1, 29)]
            headers = [str(h) for h in headers if h is not None]

            current_sample = None
            current_file = None
            genotype = "WT" if "WT" in sheet_name.upper() else "abcc8"

            for row_idx in range(6, max_row + 1):
                if cells.get((row_idx, 1)) is not None:
                    current_sample = cells.get((row_idx, 1))
                if cells.get((row_idx, 2)) is not None:
                    current_file = cells.get((row_idx, 2))
                if cells.get((row_idx, 3)) is None:
                    continue

                record: dict[str, object] = {
                    "sheet": sheet_name,
                    "genotype": genotype,
                    "sample": current_sample,
                    "file": current_file,
                    "source_row": row_idx,
                }
                for idx, header in enumerate(headers, start=1):
                    record[header] = cells.get((row_idx, idx + 2))
                rows.append(record)
    return rows


def as_float(value: object) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return None


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    idx = (len(ordered) - 1) * probability
    lo = math.floor(idx)
    hi = math.ceil(idx)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - idx) + ordered[hi] * (idx - lo)


def summarize(values: list[float]) -> dict[str, float]:
    return {
        "n": len(values),
        "mean": statistics.mean(values),
        "median": statistics.median(values),
        "q1": quantile(values, 0.25),
        "q3": quantile(values, 0.75),
        "min": min(values),
        "max": max(values),
    }


def exact_or_sampled_permutation(
    wt_values: list[float],
    abcc8_values: list[float],
    *,
    statistic: str,
    max_exact: int,
    samples: int,
    seed: int,
) -> dict[str, object]:
    def stat(vals: list[float]) -> float:
        if statistic == "mean":
            return statistics.mean(vals)
        if statistic == "median":
            return statistics.median(vals)
        raise ValueError(statistic)

    values = wt_values + abcc8_values
    n_wt = len(wt_values)
    observed = stat(abcc8_values) - stat(wt_values)
    combinations = math.comb(len(values), n_wt)

    if combinations <= max_exact:
        iterator = itertools.combinations(range(len(values)), n_wt)
        total_expected = combinations
        mode = "exact"
    else:
        rng = random.Random(seed)
        iterator = (rng.sample(range(len(values)), n_wt) for _ in range(samples))
        total_expected = samples
        mode = "sampled"

    total = 0
    extreme = 0
    for wt_idx_iter in iterator:
        wt_idx = set(wt_idx_iter)
        perm_wt = [values[idx] for idx in wt_idx]
        perm_abcc8 = [values[idx] for idx in range(len(values)) if idx not in wt_idx]
        diff = stat(perm_abcc8) - stat(perm_wt)
        total += 1
        if abs(diff) >= abs(observed) - 1e-15:
            extreme += 1

    return {
        "statistic": statistic,
        "observed_abcc8_minus_wt": observed,
        "two_sided_p": extreme / total,
        "mode": mode,
        "n_permutations": total,
        "n_possible_permutations": total_expected,
    }


def build_file_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, object], list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault((str(row["genotype"]), row["file"]), []).append(row)

    out: list[dict[str, object]] = []
    for (genotype, file_id), group_rows in sorted(grouped.items()):
        speeds = [as_float(row.get("TRACK_MEAN_SPEED")) for row in group_rows]
        speeds = [v for v in speeds if v is not None]
        durations = [as_float(row.get("TRACK_DURATION")) for row in group_rows]
        durations = [v for v in durations if v is not None]
        out.append(
            {
                "genotype": genotype,
                "file": file_id,
                "n_tracks": len(speeds),
                "mean_track_mean_speed": statistics.mean(speeds),
                "median_track_mean_speed": statistics.median(speeds),
                "min_track_mean_speed": min(speeds),
                "max_track_mean_speed": max(speeds),
                "mean_track_duration_s": statistics.mean(durations) if durations else "",
            }
        )
    return out


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(
    path: Path,
    workbook: Path,
    track_summary: dict[str, dict[str, float]],
    file_summary_stats: dict[str, dict[str, float]],
    tests: dict[str, list[dict[str, object]]],
    rows: list[dict[str, object]],
    file_rows: list[dict[str, object]],
) -> None:
    lines: list[str] = []
    lines.append("# XIG TrackMate Summary Analysis")
    lines.append("")
    lines.append(f"Workbook: `{workbook}`")
    lines.append("")
    lines.append("## Data Type")
    lines.append("")
    lines.append(
        "This workbook contains TrackMate track-summary metrics, not full per-frame spot trajectories. "
        "It supports mean-speed, total-distance, displacement, and confinement summaries, but not direct "
        "velocity autocorrelation, reversal frequency, MSD, or directional persistence calculations."
    )
    lines.append("")
    lines.append("## Track-Level Mean Speed")
    lines.append("")
    for genotype, stats in track_summary.items():
        lines.append(
            f"- {genotype}: n_tracks={int(stats['n'])}, mean={stats['mean']:.4g}, "
            f"median={stats['median']:.4g}, IQR={stats['q1']:.4g}-{stats['q3']:.4g}, "
            f"range={stats['min']:.4g}-{stats['max']:.4g}"
        )
    track_diff = track_summary["abcc8"]["mean"] - track_summary["WT"]["mean"]
    track_fold = track_summary["abcc8"]["mean"] / track_summary["WT"]["mean"]
    lines.append(f"- Mean difference, abcc8 minus WT: {track_diff:.4g}; fold-change: {track_fold:.3g}x.")
    lines.append("")
    lines.append("## Movie/File-Level Mean Speed")
    lines.append("")
    for genotype, stats in file_summary_stats.items():
        lines.append(
            f"- {genotype}: n_files={int(stats['n'])}, mean={stats['mean']:.4g}, "
            f"median={stats['median']:.4g}, IQR={stats['q1']:.4g}-{stats['q3']:.4g}, "
            f"range={stats['min']:.4g}-{stats['max']:.4g}"
        )
    file_diff = file_summary_stats["abcc8"]["mean"] - file_summary_stats["WT"]["mean"]
    file_fold = file_summary_stats["abcc8"]["mean"] / file_summary_stats["WT"]["mean"]
    lines.append(f"- Mean difference, abcc8 minus WT: {file_diff:.4g}; fold-change: {file_fold:.3g}x.")
    lines.append("")
    lines.append("## Permutation Tests")
    lines.append("")
    for level, level_tests in tests.items():
        lines.append(f"- {level}:")
        for result in level_tests:
            lines.append(
                f"  - {result['statistic']} difference={result['observed_abcc8_minus_wt']:.4g}, "
                f"two-sided p={result['two_sided_p']:.4g} "
                f"({result['mode']}, n={result['n_permutations']})"
            )
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append(
        "The collaborator's qualitative conclusion is supported directionally: abcc8 tracks have higher "
        "TrackMate mean speed than WT in this workbook. However, treating each track as an independent "
        "replicate gives stronger statistics than treating each movie/file as the replicate. The file-level "
        "comparison still trends faster in abcc8, but with only five movies per condition it should be "
        "reported cautiously unless more movies or full trajectory exports are available."
    )
    lines.append("")
    lines.append("## Needed For Deeper Motion Behavior")
    lines.append("")
    lines.append(
        "To compute persistence, back-and-forth/reversal behavior, confined diffusion, velocity "
        "autocorrelation, and MSD, request the TrackMate spots/edges/tracks export with at least: "
        "track ID, frame/time, x, y, z if available, and quality/intensity. Also request pixel size, "
        "frame interval, whether the movie is z-projected or single-plane/3D, and whether drift correction "
        "was applied."
    )
    lines.append("")
    lines.append("## Files Parsed")
    lines.append("")
    lines.append(f"- Track rows: {len(rows)}")
    lines.append(f"- File-level rows: {len(file_rows)}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("workbook", type=Path)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("analysis/results/xig_cluster_speed_2026-06-15"),
    )
    parser.add_argument("--max-exact-permutations", type=int, default=3_000_000)
    parser.add_argument("--sampled-permutations", type=int, default=200_000)
    args = parser.parse_args()

    rows = parse_workbook(args.workbook)
    if not rows:
        raise SystemExit("No TrackMate rows found")

    file_rows = build_file_summary(rows)
    args.outdir.mkdir(parents=True, exist_ok=True)

    track_fields = [
        "sheet",
        "genotype",
        "sample",
        "file",
        "source_row",
        "LABEL",
        "TRACK_INDEX",
        "TRACK_ID",
        "NUMBER_SPOTS",
        "TRACK_DURATION",
        "TRACK_START",
        "TRACK_STOP",
        "TRACK_DISPLACEMENT",
        "TRACK_MEAN_SPEED",
        "TRACK_MAX_SPEED",
        "TRACK_MEDIAN_SPEED",
        "TRACK_STD_SPEED",
        "TOTAL_DISTANCE_TRAVELED",
        "CONFINEMENT_RATIO",
        "MEAN_STRAIGHT_LINE_SPEED",
        "LINEARITY_OF_FORWARD_PROGRESSION",
        "MEAN_DIRECTIONAL_CHANGE_RATE",
    ]
    write_csv(args.outdir / "track_summary_normalized.csv", rows, track_fields)
    write_csv(
        args.outdir / "file_summary.csv",
        file_rows,
        [
            "genotype",
            "file",
            "n_tracks",
            "mean_track_mean_speed",
            "median_track_mean_speed",
            "min_track_mean_speed",
            "max_track_mean_speed",
            "mean_track_duration_s",
        ],
    )

    track_by_genotype: dict[str, list[float]] = {}
    for row in rows:
        speed = as_float(row.get("TRACK_MEAN_SPEED"))
        if speed is not None:
            track_by_genotype.setdefault(str(row["genotype"]), []).append(speed)

    file_by_genotype: dict[str, list[float]] = {}
    for row in file_rows:
        file_by_genotype.setdefault(str(row["genotype"]), []).append(float(row["mean_track_mean_speed"]))

    track_summary = {key: summarize(vals) for key, vals in sorted(track_by_genotype.items())}
    file_summary_stats = {key: summarize(vals) for key, vals in sorted(file_by_genotype.items())}

    tests = {
        "track-level": [
            exact_or_sampled_permutation(
                track_by_genotype["WT"],
                track_by_genotype["abcc8"],
                statistic="mean",
                max_exact=args.max_exact_permutations,
                samples=args.sampled_permutations,
                seed=1,
            ),
            exact_or_sampled_permutation(
                track_by_genotype["WT"],
                track_by_genotype["abcc8"],
                statistic="median",
                max_exact=args.max_exact_permutations,
                samples=args.sampled_permutations,
                seed=2,
            ),
        ],
        "file-level": [
            exact_or_sampled_permutation(
                file_by_genotype["WT"],
                file_by_genotype["abcc8"],
                statistic="mean",
                max_exact=args.max_exact_permutations,
                samples=args.sampled_permutations,
                seed=3,
            ),
            exact_or_sampled_permutation(
                file_by_genotype["WT"],
                file_by_genotype["abcc8"],
                statistic="median",
                max_exact=args.max_exact_permutations,
                samples=args.sampled_permutations,
                seed=4,
            ),
        ],
    }

    summary_payload = {
        "workbook": str(args.workbook),
        "track_summary": track_summary,
        "file_summary": file_summary_stats,
        "permutation_tests": tests,
    }
    (args.outdir / "summary.json").write_text(
        json.dumps(summary_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_summary(
        args.outdir / "README.md",
        args.workbook,
        track_summary,
        file_summary_stats,
        tests,
        rows,
        file_rows,
    )

    print(f"Wrote {args.outdir}")
    print((args.outdir / "README.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
