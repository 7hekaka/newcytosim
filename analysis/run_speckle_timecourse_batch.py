#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import plot_timecourse_panels as base
from analysis.metrics import compute_top_bottom_bands
from trust.analyze_axial import central_diff, reduce_frame_point, unwrap_angles

FRAME_RE = re.compile(r'^% frame\s+(\d+)')
FIBER_RE = re.compile(r'^%\s+fiber\s+\S+:(\d+)')
DATA_RE = re.compile(r'^\s*([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s*$')


def pick_report_bin() -> Path:
    for cand in [base.ROOT / 'build' / 'bin' / 'report', base.ROOT / 'bin' / 'report']:
        if cand.exists():
            return cand
    raise FileNotFoundError('Could not find report executable in build/bin or bin')


def parse_frames_speckle(path: Path) -> list[np.ndarray]:
    frames = []
    cur = []
    current_fid = None
    seen_frame = False
    with path.open('r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if FRAME_RE.match(line):
                if seen_frame and cur:
                    frames.append(np.array(cur, dtype=float))
                    cur = []
                seen_frame = True
                current_fid = None
                continue
            m_fib = FIBER_RE.match(line)
            if m_fib:
                current_fid = int(m_fib.group(1))
                continue
            m = DATA_RE.match(line)
            if m and current_fid is not None:
                x = float(m.group(1)); y = float(m.group(2)); z = float(m.group(3)); a = float(m.group(4))
                cur.append((current_fid, x, y, z, a))
    if cur:
        frames.append(np.array(cur, dtype=float))
    if not frames:
        raise ValueError(f'No speckle frames parsed from {path}')
    return frames


def nematic_order_xy_from_speckles(frame: np.ndarray) -> float:
    if frame.size == 0:
        return 0.0
    order = np.lexsort((frame[:, 4], frame[:, 0]))
    data = frame[order]
    c2_sum = 0.0
    s2_sum = 0.0
    n = 0
    for k in range(1, data.shape[0]):
        if data[k, 0] != data[k - 1, 0]:
            continue
        dx = data[k, 1] - data[k - 1, 1]
        dy = data[k, 2] - data[k - 1, 2]
        norm = math.hypot(dx, dy)
        if norm <= 0:
            continue
        th = math.atan2(dy, dx)
        c2_sum += math.cos(2 * th)
        s2_sum += math.sin(2 * th)
        n += 1
    if n == 0:
        return 0.0
    return float(math.hypot(c2_sum, s2_sum) / n)


def ensure_speckle(run_path: Path, speckle_path: Path, report_bin: Path, interval: float) -> None:
    if speckle_path.exists() and speckle_path.stat().st_size > 0:
        return
    with speckle_path.open('w', encoding='utf-8') as fh:
        subprocess.run(
            [str(report_bin), 'fiber:speckle', f'interval={interval}'],
            cwd=run_path,
            stdout=fh,
            stderr=subprocess.DEVNULL,
            check=True,
        )


def analyze_speckle_run(task: tuple[dict, float, str, str, int, int, int]) -> dict:
    row, interval, speckle_name, report_bin_s, nr, nth, nz = task
    run_path = Path(row['run_path'])
    config_path = run_path / 'config.cym'
    speckle_path = run_path / speckle_name
    out = {
        'model': row['model'],
        'group': row['group'],
        'case': row['case'],
        'xlink_regime': row['xlink_regime'],
        'run_dir': row['run_dir'],
        'run_path': row['run_path'],
        'status': 'pending',
    }
    try:
        ensure_speckle(run_path, speckle_path, Path(report_bin_s), interval)
        cfg = base.parse_config_info(config_path)
        runs = cfg['runs']
        final_run = runs[-1]
        nb_frames = final_run['nb_frames']
        dt = final_run['frame_dt']
        frames = parse_frames_speckle(speckle_path)
        total_expected = sum(r['nb_frames'] for r in runs)
        total_expected_with_initial = sum(r['nb_frames'] + 1 for r in runs)
        if len(frames) == total_expected:
            phase_frames = frames[-nb_frames:]
        elif len(frames) == total_expected_with_initial:
            phase_frames = frames[-(nb_frames + 1):][1:]
        else:
            out['status'] = f'error:expected {total_expected} or {total_expected_with_initial} frames but parsed {len(frames)}'
            return out
        space = cfg['space']
        zmin = float(space.get('bottom', -20.0))
        zmax = float(space.get('top', 20.0))
        inner = float(space.get('inner', 0.0))
        outer = float(space.get('outer', 1.0))
        z_com = []
        theta_com = []
        Dz = []
        top_bias = []
        nematic_xy = []
        shell_occ = []
        for frame in phase_frames:
            reduced = reduce_frame_point(frame[:, :4], zmin=zmin, zmax=zmax, nbins_z=80, zmid=0.0, Lz=(zmax - zmin))
            x = frame[:, 1]; y = frame[:, 2]; z = frame[:, 3]
            bands = compute_top_bottom_bands(z, top_band_frac=1 / 3)
            z_com.append(reduced['z_com'])
            theta_com.append(reduced['theta_com'])
            Dz.append(reduced['Dz'])
            top_bias.append(bands.top_bias)
            nematic_xy.append(nematic_order_xy_from_speckles(frame))
            shell_occ.append(base.shell_occupancy_fraction(x, y, z, inner=inner, outer=outer, bottom=zmin, top=zmax, nr=nr, nth=nth, nz=nz))
        z_com = np.asarray(z_com, dtype=float)
        theta_com = np.asarray(theta_com, dtype=float)
        out.update({
            'status': 'ok',
            'n_frames': int(nb_frames),
            'dt_s': float(dt),
            'time_min': ((dt * np.arange(1, nb_frames + 1, dtype=float)) / 60.0).tolist(),
            'vz_abs': np.abs(central_diff(z_com, dt)).tolist(),
            'Dz': np.asarray(Dz, dtype=float).tolist(),
            'shell_occ': np.asarray(shell_occ, dtype=float).tolist(),
            'nematic_xy': np.asarray(nematic_xy, dtype=float).tolist(),
            'top_bias': np.asarray(top_bias, dtype=float).tolist(),
            'swirl_rate_abs': np.abs(central_diff(unwrap_angles(theta_com), dt)).tolist(),
        })
        return out
    except Exception as exc:
        out['status'] = f'error:{exc}'
        return out


def main() -> None:
    ap = argparse.ArgumentParser(description='Run sparse-speckle time-course batches for selected campaign groups.')
    ap.add_argument('--rotatable', type=Path, default=base.DEFAULT_ROTATABLE)
    ap.add_argument('--fixed', type=Path, default=base.DEFAULT_FIXED)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--group', action='append', default=[])
    ap.add_argument('--metric', action='append', default=[])
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--interval', type=float, default=2.0)
    ap.add_argument('--speckle-name', default='speckle_i2.txt')
    ap.add_argument('--nr', type=int, default=8)
    ap.add_argument('--nth', type=int, default=72)
    ap.add_argument('--nz', type=int, default=80)
    args = ap.parse_args()

    plt.rcParams.update({'font.size': 10, 'axes.linewidth': 1.1, 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})

    group_filter = set(args.group)
    metric_filter = set(args.metric)
    report_bin = pick_report_bin()

    tasks = []
    for row in base.load_metadata(args.rotatable, 'rotatable'):
        if group_filter and row['group'] not in group_filter:
            continue
        tasks.append((row, args.interval, args.speckle_name, str(report_bin), args.nr, args.nth, args.nz))
    for row in base.load_metadata(args.fixed, 'fixed_global'):
        if group_filter and row['group'] not in group_filter:
            continue
        tasks.append((row, args.interval, args.speckle_name, str(report_bin), args.nr, args.nth, args.nz))

    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            run_rows = list(ex.map(analyze_speckle_run, tasks))
    else:
        run_rows = [analyze_speckle_run(task) for task in tasks]

    args.outdir.mkdir(parents=True, exist_ok=True)
    base.write_csv(args.outdir / 'timecourse_run_status.csv', run_rows)
    condition_rows = base.aggregate_runs(run_rows)
    base.write_csv(args.outdir / 'condition_timecourses.csv', condition_rows)
    lookup = base.build_lookup(condition_rows)
    family_specs = [f for f in base.FAMILY_SPECS if not group_filter or f['group'] in group_filter]
    metrics = [m for m in base.METRICS if not metric_filter or m['field'] in metric_filter or m['slug'] in metric_filter]
    outfiles = []
    for family_spec in family_specs:
        for metric in metrics:
            outfiles.extend(base.make_timecourse_figure(lookup, family_spec, metric, args.outdir))
    print(args.outdir / 'timecourse_run_status.csv')
    print(args.outdir / 'condition_timecourses.csv')
    for path in outfiles:
        print(path)


if __name__ == '__main__':
    main()