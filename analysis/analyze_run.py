# analysis/analyze_run.py
import os
import argparse
from .reader import parse_frames
from .engine import analyze_from_frames
from .plotting import make_plots

def analyze_file(path, nfil, dt, drop_first=False):
    frames = list(parse_frames(path, nfil=nfil))
    print(f"Parsed {len(frames)} frames | Duration ≈ {(len(frames)-1)*dt:.2f} s")
    return analyze_from_frames(frames, dt=dt, drop_first=drop_first)

def _cli():
    p = argparse.ArgumentParser()
    p.add_argument("--file", required=True)
    p.add_argument("--nfil", type=int, required=True)
    p.add_argument("--dt", type=float, default=1.0)
    p.add_argument("--drop-first", action="store_true")
    p.add_argument("--minutes", action="store_true")
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    res = analyze_file(args.file, nfil=args.nfil, dt=args.dt, drop_first=args.drop_first)

    # base name for outputs
    base = os.path.splitext(os.path.basename(args.file))[0]
    os.makedirs(args.outdir, exist_ok=True)
    save_prefix = os.path.join(args.outdir, base)

    make_plots(res, minutes=args.minutes, save_prefix=save_prefix)

if __name__ == "__main__":
    _cli()
