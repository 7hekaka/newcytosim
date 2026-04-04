import re
from pathlib import Path
from typing import List
import numpy as np

ROW_RE = re.compile(
    r"^\s*(\d+)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s*$"
)

def parse_frames(path: Path, nfil: int) -> List[np.ndarray]:
    """
    Parse a NON-VERBOSE Cytosim point.txt into frames.
    Each frame is a numpy array with columns [fid, x, y, z].

    Frame boundary = filament id wraps from nfil -> 1.
    """
    if not isinstance(path, Path):
        path = Path(path)

    frames, cur = [], []
    last_fid = None

    with path.open("r") as fh:
        for line in fh:
            m = ROW_RE.match(line)
            if not m:
                # ignore malformed lines quietly (secure: no exec/eval)
                continue
            fid = int(m.group(1))
            x = float(m.group(2)); y = float(m.group(3)); z = float(m.group(4))
            # curv = float(m.group(5))  # parsed but unused

            # boundary: last filament id followed by 1 => new frame
            if last_fid is not None and last_fid == nfil and fid == 1 and len(cur) > 0:
                frames.append(np.array(cur, dtype=float))
                cur = []

            cur.append((fid, x, y, z))
            last_fid = fid

    if cur:
        frames.append(np.array(cur, dtype=float))

    if not frames:
        raise ValueError("No frames parsed. Check --nfil or file format.")

    return frames
