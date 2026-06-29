#!/usr/bin/env python3
import sys
from pathlib import Path

root = Path('/home/thekaka/project/cytosim')
if str(root) not in sys.path:
    sys.path.insert(0, str(root))

import analysis.transport_regime_batch as t

t.ROT = root / 'clu_init_minusz' / 'xlink_regime_map.csv'
t.FIX = root / 'clu_fixed_global_init_minusz' / 'xlink_regime_map.csv'
t.OUT = root / 'analysis' / 'results' / 'init_minusz_comparison' / 'transport_regime'
t.PLOTS = t.OUT / 'plots'
t.MODEL_SPECS = [
    ('rotatable', 'Rotatable', '#222222', 's', 'white'),
    ('fixed_global', 'Fixed global', '#d95f02', 'o', '#d95f02'),
]

t.main()
