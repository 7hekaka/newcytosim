#!/usr/bin/env python3
from pathlib import Path
import csv
import shutil

BASE = Path('/home/thekaka/project/cytosim/clu_fixed_global_init_minusz')
MAPPING = {
    'job00': ['mpc12/c10_m12'],
    'job01': ['mpc12/c20_m12'],
    'job02': ['mpc12/c40_m12', 'total480/c40_m12'],
    'job03': ['mpc12/c60_m12'],
    'job04': ['mpc12/c80_m12'],
    'job05': ['total480/c10_m48'],
    'job06': ['total480/c20_m24'],
    'job07': ['total480/c60_m8'],
    'job09': ['total480/c80_m6'],
}

rows = []
for job, targets in MAPPING.items():
    save_dir = BASE / job / 'save'
    run_dirs = sorted([p for p in save_dir.glob('r*') if p.is_dir()]) if save_dir.exists() else []
    for target in targets:
        target_dir = BASE / target
        imported = 0
        for src_run in run_dirs:
            dst_run = target_dir / src_run.name
            if not dst_run.exists():
                continue
            for item in src_run.iterdir():
                dst_item = dst_run / item.name
                if item.is_dir():
                    shutil.copytree(item, dst_item, dirs_exist_ok=True)
                else:
                    shutil.copy2(item, dst_item)
            imported += 1
        rows.append({
            'job': job,
            'target': str(target_dir),
            'runs_in_save': len(run_dirs),
            'runs_imported': imported,
        })

summary = BASE / 'job_import_summary.csv'
with summary.open('w', newline='', encoding='utf-8') as fh:
    writer = csv.DictWriter(fh, fieldnames=['job', 'target', 'runs_in_save', 'runs_imported'])
    writer.writeheader()
    writer.writerows(rows)

print(summary)
for row in rows:
    print(f"{row['job']} -> {row['target']} ({row['runs_imported']}/{row['runs_in_save']})")
