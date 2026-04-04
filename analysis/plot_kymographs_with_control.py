#!/usr/bin/env python3
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import plot_timecourse_panels as base
import run_speckle_kymograph_batch as kymo
import run_speckle_timecourse_batch as speck

ROOT = Path(__file__).resolve().parent.parent
CONTROL = ROOT / 'control' / '2'
ROT = base.DEFAULT_ROTATABLE
FIX = base.DEFAULT_FIXED
OUT = ROOT / 'comparison_kymographs_counts_n0'
MATCH_FINAL_MIN = 800.0 / 60.0

FAMILY_SPECS = [
    {
        'group': 'total480',
        'slug': 'total480',
        'order': [
            ('c0_m0', '0'),
            ('c10_m48', '10'),
            ('c20_m24', '20'),
            ('c40_m12', '40'),
            ('c60_m8', '60'),
            ('c80_m6', '80'),
        ],
    },
    {
        'group': 'mpc12',
        'slug': 'mpc12',
        'order': [
            ('c0_m0', '0'),
            ('c10_m12', '10'),
            ('c20_m12', '20'),
            ('c40_m12', '40'),
            ('c60_m12', '60'),
            ('c80_m12', '80'),
        ],
    },
]


def control_metadata() -> list[dict]:
    rows = []
    for idx in range(1, 31):
        if idx <= 10:
            regime = '1:4'
        elif idx <= 20:
            regime = '1:8'
        else:
            regime = '1:16'
        run_dir = f'r{idx:04d}'
        run_path = CONTROL / run_dir
        rows.append(
            {
                'model': 'control',
                'group': 'control',
                'case': 'c0_m0',
                'xlink_regime': regime,
                'run_dir': run_dir,
                'run_path': str(run_path),
            }
        )
    return rows


def duplicate_control_rows(rows: list[dict]) -> list[dict]:
    out = list(rows)
    control_rows = [row for row in rows if row['model'] == 'control']
    for row in control_rows:
        for model in ['rotatable', 'fixed_global']:
            for group in ['total480', 'mpc12']:
                clone = dict(row)
                clone['model'] = model
                clone['group'] = group
                out.append(clone)
    return out


def main() -> None:
    report_bin = speck.pick_report_bin()
    tasks = []
    for row in control_metadata():
        tasks.append((row, 2.0, 'speckle_i2.txt', str(report_bin), 80))
    for row in base.load_metadata(ROT, 'rotatable'):
        if row['group'] in {'total480', 'mpc12'}:
            tasks.append((row, 2.0, 'speckle_i2.txt', str(report_bin), 80))
    for row in base.load_metadata(FIX, 'fixed_global'):
        if row['group'] in {'total480', 'mpc12'}:
            tasks.append((row, 2.0, 'speckle_i2.txt', str(report_bin), 80))

    with ProcessPoolExecutor(max_workers=4) as ex:
        run_rows = list(ex.map(kymo.analyze_run, tasks))

    for row in run_rows:
        if row.get('model') != 'control' or row.get('status') != 'ok':
            continue
        time = row['time_min']
        keep = sum(1 for t in time if float(t) <= MATCH_FINAL_MIN + 1e-9)
        row['time_min'] = row['time_min'][:keep]
        row['rho_z'] = row['rho_z'][:keep]

    OUT.mkdir(parents=True, exist_ok=True)
    kymo.write_status_csv(OUT / 'kymograph_run_status.csv', run_rows)

    condition_rows = kymo.aggregate_runs(run_rows)
    condition_rows = duplicate_control_rows(condition_rows)
    kymo.write_condition_csv(OUT / 'kymograph_condition_summary.csv', condition_rows)
    lookup = kymo.build_lookup(condition_rows)

    vmin, vmax = 0.0, 2.0
    outputs = []
    for family_spec in FAMILY_SPECS:
        outputs.append(kymo.plot_family_model(lookup, family_spec, model='rotatable', outdir=OUT, vmin=vmin, vmax=vmax))
        outputs.append(kymo.plot_family_model(lookup, family_spec, model='fixed_global', outdir=OUT, vmin=vmin, vmax=vmax))
    cbar_path = kymo.save_standalone_colorbar(OUT, vmin=vmin, vmax=vmax)

    with (OUT / 'kymograph_scale.txt').open('w', encoding='utf-8') as fh:
        fh.write(f'vmin={vmin}\n')
        fh.write(f'vmax={vmax}\n')

    print(OUT / 'kymograph_run_status.csv')
    print(OUT / 'kymograph_condition_summary.csv')
    print(OUT / 'kymograph_scale.txt')
    for path in outputs:
        print(path)
    print(cbar_path)


if __name__ == '__main__':
    main()

