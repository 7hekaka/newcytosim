#!/usr/bin/env python3
import sys
from pathlib import Path

root = Path('/home/thekaka/project/cytosim')
if str(root) not in sys.path:
    sys.path.insert(0, str(root))
if str(root / 'analysis') not in sys.path:
    sys.path.insert(0, str(root / 'analysis'))

from analysis import comparison_unwrapped_annulus_story as s

s.ROT_MAP = root / 'clu_init_minusz' / 'xlink_regime_map.csv'
s.FIX_MAP = root / 'clu_fixed_global_init_minusz' / 'xlink_regime_map.csv'
s.MEDOID = root / 'analysis' / 'results' / 'init_minusz_comparison' / 'transport_regime' / 'medoid_selection_summary.csv'
s.OUT = root / 'analysis' / 'results' / 'init_minusz_comparison' / 'unwrapped_annulus'
s.FAMILIES = {'total480', 'mpc12'}


def build_run_lookup(path):
    rows = s.read_csv_rows(path)
    lookup = {}
    for row in rows:
        avail = row.get('has_properties')
        if avail is None:
            avail = row.get('has_point')
        if avail is None:
            avail = row.get('has_outputs')
        if str(avail).strip().lower() not in {'1', 'true', 'yes'}:
            continue
        lookup[(row['group'], row['case'], row['xlink_regime'], row['run_dir'])] = row
    return lookup


s.build_run_lookup = build_run_lookup

s.main()
