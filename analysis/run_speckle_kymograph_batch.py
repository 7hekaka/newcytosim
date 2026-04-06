#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

import plot_timecourse_panels as base
import run_speckle_timecourse_batch as speck
from trust.analyze_axial import reduce_frame_point


def analyze_run(task: tuple[dict, float, str, str, int]) -> dict:
    row, interval, speckle_name, report_bin_s, nbins_z = task
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
        speck.ensure_speckle(run_path, speckle_path, Path(report_bin_s), interval)
        cfg = base.parse_config_info(config_path)
        runs = cfg['runs']
        final_run = runs[-1]
        nb_frames = final_run['nb_frames']
        dt = final_run['frame_dt']
        frames = speck.parse_frames_speckle(speckle_path)

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
        rho = []
        for frame in phase_frames:
            reduced = reduce_frame_point(frame[:, :4], zmin=zmin, zmax=zmax, nbins_z=nbins_z, zmid=0.0, Lz=(zmax - zmin))
            rho.append(np.asarray(reduced['rho_z'], dtype=float))

        out.update(
            {
                'status': 'ok',
                'time_min': ((dt * np.arange(1, nb_frames + 1, dtype=float)) / 60.0).tolist(),
                'z_edges': np.linspace(zmin, zmax, nbins_z + 1, dtype=float).tolist(),
                'rho_z': np.stack(rho, axis=0).tolist(),
            }
        )
        return out
    except Exception as exc:
        out['status'] = f'error:{exc}'
        return out


def aggregate_runs(run_rows: list[dict]) -> list[dict]:
    grouped: dict[tuple[str, str, str, str], list[dict]] = defaultdict(list)
    for row in run_rows:
        if row.get('status') == 'ok':
            grouped[(row['model'], row['group'], row['case'], row['xlink_regime'])].append(row)

    out_rows = []
    for (model, group, case, regime), items in sorted(grouped.items()):
        n_frames = min(len(item['time_min']) for item in items)
        t = np.asarray(items[0]['time_min'][:n_frames], dtype=float)
        z_edges = np.asarray(items[0]['z_edges'], dtype=float)
        rho_stack = np.stack([np.asarray(item['rho_z'][:n_frames], dtype=float) for item in items], axis=0)
        rho_mean = np.nanmean(rho_stack, axis=0)
        out_rows.append(
            {
                'model': model,
                'group': group,
                'case': case,
                'xlink_regime': regime,
                'n_runs_used': len(items),
                'time_min': t,
                'z_edges': z_edges,
                'rho_z_mean': rho_mean,
            }
        )
    return out_rows


def build_lookup(rows: list[dict]) -> dict[tuple[str, str, str, str], dict]:
    return {(r['model'], r['group'], r['case'], r['xlink_regime']): r for r in rows}


def compute_vrange(rows: list[dict], percentile: float = 99.0) -> tuple[float, float]:
    vals = []
    for row in rows:
        arr = np.asarray(row['rho_z_mean'], dtype=float)
        vals.append(arr.ravel())
    if not vals:
        return 0.0, 1.0
    all_vals = np.concatenate(vals)
    vmax = float(np.nanpercentile(all_vals, percentile))
    return 0.0, vmax


def plot_family_model(
    lookup: dict,
    family_spec: dict,
    *,
    model: str,
    outdir: Path,
    vmin: float,
    vmax: float,
    with_colorbar: bool = False,
) -> Path:
    ncols = len(family_spec['order'])
    fig, axes = plt.subplots(3, ncols, figsize=(4.0 * ncols + 0.8, 9.8), constrained_layout=False)
    if ncols == 1:
        axes = np.array(axes).reshape(3, 1)

    im = None
    for ridx, regime in enumerate(base.XLINK_ORDER):
        for cidx, (case, label) in enumerate(family_spec['order']):
            ax = axes[ridx, cidx]
            row = lookup.get((model, family_spec['group'], case, regime))
            ax.tick_params(direction='out', width=1.1, labelsize=13)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(1.1)
            ax.spines['bottom'].set_linewidth(1.1)

            if ridx == 0:
                ax.set_title(label, fontsize=16, pad=10)

            if row is None:
                ax.set_facecolor('#f3f3f3')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.text(0.5, 0.5, 'n/a', ha='center', va='center', transform=ax.transAxes, fontsize=14, color='#777777')
            else:
                rho = np.asarray(row['rho_z_mean'], dtype=float)
                t = np.asarray(row['time_min'], dtype=float)
                z_edges = np.asarray(row['z_edges'], dtype=float)
                if len(t) > 1:
                    mids = 0.5 * (t[1:] + t[:-1])
                    t_edges = np.concatenate(([t[0] - (mids[0] - t[0])], mids, [t[-1] + (t[-1] - mids[-1])]))
                else:
                    t_edges = np.array([0.0, t[0] + 1e-6], dtype=float)
                im = ax.pcolormesh(t_edges, z_edges, rho.T, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax)
                ax.set_xlim(t_edges[0], t_edges[-1])

            if cidx == 0:
                ax.text(
                    0.03,
                    0.94,
                    regime,
                    transform=ax.transAxes,
                    ha='left',
                    va='top',
                    fontsize=15,
                    color='white' if row is not None else '#555555',
                )
            else:
                ax.set_yticklabels([])
            if ridx < 2:
                ax.set_xticklabels([])

    if with_colorbar and im is not None:
        cax = fig.add_axes([0.92, 0.14, 0.018, 0.72])
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label('Normalized axial density', fontsize=16)
        cbar.ax.tick_params(labelsize=12, width=1.0)

    fig.text(0.5, 0.04, 'Time after motor onset (min)', ha='center', va='center', fontsize=18)
    fig.text(0.025, 0.5, 'z (um)', ha='center', va='center', rotation='vertical', fontsize=18)
    right = 0.89 if with_colorbar else 0.98
    fig.subplots_adjust(left=0.09, right=right, top=0.92, bottom=0.10, wspace=0.10, hspace=0.10)
    suffix = '_kymograph_with_colorbar.png' if with_colorbar else '_kymograph.png'
    outpath = outdir / f"{family_spec['slug']}_{model}{suffix}"
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return outpath


def save_standalone_colorbar(outdir: Path, *, vmin: float, vmax: float) -> Path:
    fig = plt.figure(figsize=(1.8, 4.8), facecolor='white')
    cax = fig.add_axes([0.40, 0.08, 0.28, 0.88])
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(norm=norm, cmap='viridis')
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label('Normalized axial density', fontsize=15)
    cbar.ax.tick_params(labelsize=11, width=1.0)
    outpath = outdir / 'kymograph_colorbar_only.png'
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return outpath


def write_status_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fieldnames = []
    slim_rows = []
    for row in rows:
        slim = {k: v for k, v in row.items() if k not in {'time_min', 'z_edges', 'rho_z'}}
        slim_rows.append(slim)
        for key in slim.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open('w', newline='', encoding='utf-8') as fh:
        import csv
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(slim_rows)


def write_condition_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fieldnames = ['model', 'group', 'case', 'xlink_regime', 'n_runs_used']
    with path.open('w', newline='', encoding='utf-8') as fh:
        import csv
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in fieldnames})


def main() -> None:
    ap = argparse.ArgumentParser(description='Aggregate and plot sparse-speckle axial-density kymographs for one campaign family.')
    ap.add_argument('--rotatable', type=Path, default=base.DEFAULT_ROTATABLE)
    ap.add_argument('--fixed', type=Path, default=base.DEFAULT_FIXED)
    ap.add_argument('--group', required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--jobs', type=int, default=3)
    ap.add_argument('--interval', type=float, default=2.0)
    ap.add_argument('--speckle-name', default='speckle_i2.txt')
    ap.add_argument('--nbins-z', type=int, default=80)
    ap.add_argument('--scale-max', type=float, default=0.0)
    args = ap.parse_args()

    plt.rcParams.update({'font.size': 12, 'axes.linewidth': 1.1, 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})

    family_specs = [f for f in base.FAMILY_SPECS if f['group'] == args.group]
    if not family_specs:
        raise ValueError(f'Unknown group {args.group}')
    family_spec = family_specs[0]

    report_bin = speck.pick_report_bin()
    tasks = []
    for row in base.load_metadata(args.rotatable, 'rotatable'):
        if row['group'] == args.group:
            tasks.append((row, args.interval, args.speckle_name, str(report_bin), args.nbins_z))
    for row in base.load_metadata(args.fixed, 'fixed_global'):
        if row['group'] == args.group:
            tasks.append((row, args.interval, args.speckle_name, str(report_bin), args.nbins_z))

    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            run_rows = list(ex.map(analyze_run, tasks))
    else:
        run_rows = [analyze_run(task) for task in tasks]

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_status_csv(args.outdir / 'kymograph_run_status.csv', run_rows)

    condition_rows = aggregate_runs(run_rows)
    write_condition_csv(args.outdir / 'kymograph_condition_summary.csv', condition_rows)
    lookup = build_lookup(condition_rows)

    if args.scale_max > 0:
        vmin, vmax = 0.0, args.scale_max
    else:
        vmin, vmax = compute_vrange(condition_rows)

    fixed_path = plot_family_model(lookup, family_spec, model='fixed_global', outdir=args.outdir, vmin=vmin, vmax=vmax)
    rot_path = plot_family_model(lookup, family_spec, model='rotatable', outdir=args.outdir, vmin=vmin, vmax=vmax)
    fixed_cb_path = plot_family_model(lookup, family_spec, model='fixed_global', outdir=args.outdir, vmin=vmin, vmax=vmax, with_colorbar=True)
    rot_cb_path = plot_family_model(lookup, family_spec, model='rotatable', outdir=args.outdir, vmin=vmin, vmax=vmax, with_colorbar=True)
    cbar_path = save_standalone_colorbar(args.outdir, vmin=vmin, vmax=vmax)

    with (args.outdir / 'kymograph_scale.txt').open('w', encoding='utf-8') as fh:
        fh.write(f'vmin={vmin}\n')
        fh.write(f'vmax={vmax}\n')
    print(args.outdir / 'kymograph_run_status.csv')
    print(args.outdir / 'kymograph_condition_summary.csv')
    print(args.outdir / 'kymograph_scale.txt')
    print(rot_path)
    print(fixed_path)
    print(rot_cb_path)
    print(fixed_cb_path)
    print(cbar_path)


if __name__ == '__main__':
    main()
