#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SRC_ROOT = ROOT / "clu_init_minusz"
OUT_ROOT = ROOT / "clu_init_minusz_diffuse"

SRC_MODEL = "rotatable"
OUT_MODEL = "diffuse_init_minusz"

MYOSIN_BLOCK_RE = re.compile(r"(set single myosin1\s*\{\n.*?\n\})", re.S)
ACTIVITY_RE = re.compile(r"(^\s*activity\s*=\s*)fixed\s*$", re.M)
CLUSTER_COMMENT = "% stationary surface clusters on INNER wall"
DIFFUSE_CLUSTER_COMMENT = "% initial motor clusters near INNER wall; myosin1 activity=diffuse"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def make_myosin_diffuse(text: str) -> str:
    def replace_block(match: re.Match[str]) -> str:
        block = match.group(1)
        block, count = ACTIVITY_RE.subn(r"\1diffuse", block, count=1)
        if count != 1:
            raise RuntimeError("set single myosin1 block did not contain `activity = fixed`")
        if "confine" not in block:
            block = block[:-2].rstrip() + "\n    confine = inside\n}"
        return block

    text, count = MYOSIN_BLOCK_RE.subn(replace_block, text, count=1)
    if count != 1:
        raise RuntimeError("could not find `set single myosin1` block")
    return text.replace(CLUSTER_COMMENT, DIFFUSE_CLUSTER_COMMENT)


def write_configs(manifest_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out_rows: list[dict[str, str]] = []
    for row in manifest_rows:
        src_config = ROOT / row["path"]
        if not src_config.is_file():
            raise RuntimeError(f"missing source config: {src_config}")

        rel = Path(row["group"]) / row["case"] / row["run"] / "config.cym"
        out_config = OUT_ROOT / rel
        out_config.parent.mkdir(parents=True, exist_ok=True)
        out_config.write_text(make_myosin_diffuse(src_config.read_text(encoding="utf-8")), encoding="utf-8")

        out_row = dict(row)
        out_row["model"] = OUT_MODEL
        out_row["notes"] = row["notes"] + "; initial clusters, myosin1 diffuse after initialization"
        out_row["path"] = str(out_config.relative_to(ROOT))
        out_rows.append(out_row)
    return out_rows


def write_xlink_map(src_rows: list[dict[str, str]]) -> None:
    out_rows: list[dict[str, str]] = []
    for row in src_rows:
        out_row = dict(row)
        out_row["model"] = OUT_MODEL
        out_row["run_path"] = str(OUT_ROOT / row["group"] / row["case"] / row["run_dir"])
        for field in ("has_point", "has_outputs"):
            if field in out_row:
                out_row[field] = "0"
        out_rows.append(out_row)
    write_csv(OUT_ROOT / "xlink_regime_map.csv", out_rows, list(src_rows[0].keys()))


def verify_configs(config_rows: list[dict[str, str]]) -> None:
    fixed_hits = 0
    diffuse_hits = 0
    for row in config_rows:
        text = (ROOT / row["path"]).read_text(encoding="utf-8")
        block_match = MYOSIN_BLOCK_RE.search(text)
        if block_match is None:
            raise RuntimeError(f"missing myosin1 block in {row['path']}")
        block = block_match.group(1)
        fixed_hits += int("activity = fixed" in block)
        diffuse_hits += int("activity = diffuse" in block)
    if fixed_hits:
        raise RuntimeError(f"{fixed_hits} generated configs still have myosin1 activity=fixed")
    if diffuse_hits != len(config_rows):
        raise RuntimeError(f"expected {len(config_rows)} diffuse myosin1 blocks, saw {diffuse_hits}")


def main() -> None:
    manifest_path = SRC_ROOT / "manifest.csv"
    xlink_path = SRC_ROOT / "xlink_regime_map.csv"
    manifest_rows = read_csv(manifest_path)
    xlink_rows = read_csv(xlink_path)

    if not manifest_rows:
        raise RuntimeError(f"empty manifest: {manifest_path}")
    if not xlink_rows:
        raise RuntimeError(f"empty xlink map: {xlink_path}")
    if {row["model"] for row in manifest_rows} != {SRC_MODEL}:
        raise RuntimeError("source manifest is not the expected rotatable init-minus-z campaign")

    out_manifest_rows = write_configs(manifest_rows)
    write_csv(OUT_ROOT / "manifest.csv", out_manifest_rows, list(manifest_rows[0].keys()))
    write_xlink_map(xlink_rows)
    verify_configs(out_manifest_rows)

    print(f"Wrote {len(out_manifest_rows)} diffuse init-minus-z configs under {OUT_ROOT}")
    print(f"Manifest: {OUT_ROOT / 'manifest.csv'}")
    print(f"Map: {OUT_ROOT / 'xlink_regime_map.csv'}")


if __name__ == "__main__":
    main()
