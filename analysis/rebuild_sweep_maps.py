#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def regime_from_run(run_dir: str) -> str:
    idx = int(run_dir[1:])
    if idx <= 10:
        return '1:4'
    if idx <= 20:
        return '1:8'
    return '1:16'


def main() -> None:
    ap = argparse.ArgumentParser(description='Rebuild xlink_regime_map.csv from a sweep manifest.')
    ap.add_argument('campaign_root', type=Path)
    ap.add_argument('manifest', type=Path)
    ap.add_argument('model_name')
    ap.add_argument('--out', type=Path, default=None)
    args = ap.parse_args()

    campaign_root = args.campaign_root.resolve()
    manifest = args.manifest.resolve()
    out = (args.out.resolve() if args.out else campaign_root / 'xlink_regime_map.csv')

    rows = []
    with manifest.open(encoding='utf-8', newline='') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rel = Path(row['path'])
            run_dir = rel.parent.name
            case = row['case']
            abs_run = (campaign_root / rel.parent.relative_to(campaign_root.name if rel.parts and rel.parts[0] == campaign_root.name else rel.parent)).resolve() if False else (campaign_root / row['group'] / case / run_dir).resolve()
            if not abs_run.exists():
                # Fall back to the path encoded in the manifest if campaign layout differs.
                abs_run = (campaign_root / row['group'] / case / run_dir).resolve()
            has_outputs = int((abs_run / 'log.txt').exists() and (abs_run / 'objects.cmo').exists())
            rows.append({
                'model': args.model_name,
                'group': row['group'],
                'case': case,
                'xlink_regime': regime_from_run(run_dir),
                'run_dir': run_dir,
                'run_path': str(abs_run),
                'has_outputs': has_outputs,
            })

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=['model', 'group', 'case', 'xlink_regime', 'run_dir', 'run_path', 'has_outputs'])
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
