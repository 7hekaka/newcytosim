import sys
from pathlib import Path

root = Path('/home/thekaka/project/cytosim')
if str(root) not in sys.path:
    sys.path.insert(0, str(root))

import analysis.transport_regime_batch as t

t.ROT = root / 'clu_init_minusz' / 'xlink_regime_map.csv'
t.FIX = root / 'clu_init_minusz' / 'empty_fixed_map.csv'
t.OUT = root / 'analysis' / 'results' / 'init_minusz_rotatable' / 'transport_regime'
t.PLOTS = t.OUT / 'plots'
t.MODEL_SPECS = [
    ('rotatable', 'Rotatable', '#222222', 's', 'white'),
]

t.main()
