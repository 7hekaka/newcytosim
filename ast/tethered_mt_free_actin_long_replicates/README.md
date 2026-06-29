# Replicated Long Tethered-MT / Free-Actin Campaign

This is the cluster-ready follow-up to `ast/tethered_mt_free_actin_long_pilot`.

Purpose:

- Run 5 replicates per condition.
- Compare fixed-length tethered microtubules against growing tethered microtubules.
- Keep actin free in the volume with no direct actin-MT binding.
- Use longer actin growth: 5 min simulations, actin total polymer at least doubled, and `ACTIN max_length = 16 um`.
- Add a larger-system version so the small system can finish first while the larger system continues.

Generated systems:

- `small_system`: `12 x 12 x 1.2 um`, nuclei at approximately the first pilot positions, fixed MT length `6 um`, actin total polymer `1200 um`.
- `large_system`: `24 x 24 x 1.2 um`, nuclei at the earlier large-test positions, fixed MT length `12 um`, material scale `4`, actin total polymer `4800 um`.

Conditions:

- `01_actin_only_perinuclear_long`
- `02_fixed_mt_uniform_actin_long`
- `03_growing_mt_uniform_actin_long`
- `04_fixed_mt_perinuclear_actin_long`
- `05_growing_mt_perinuclear_actin_long`
- `06_division_fixed_mt_perinuclear_actin_long`
- `07_division_growing_mt_perinuclear_actin_long`

Each condition has `r0001` to `r0005`.

Regenerate configs:

```sh
cd /home/thekaka/project/cytosim
python3 ast/tethered_mt_free_actin_long_replicates/generate_tethered_mt_free_actin_long_replicates.py
```

Smoke tests:

```sh
cd /home/thekaka/project/cytosim/ast/tethered_mt_free_actin_long_replicates/00_smoke_small_system
../../../sim config.cym

cd /home/thekaka/project/cytosim/ast/tethered_mt_free_actin_long_replicates/00_smoke_large_system
../../../sim config.cym
```

Submit all replicated runs:

```sh
cd /home/thekaka/project/cytosim/ast/tethered_mt_free_actin_long_replicates
./submit_cluster.sh
```

Submit only one system if needed:

```sh
SYSTEM_FILTER=small ./submit_cluster.sh
SYSTEM_FILTER=large ./submit_cluster.sh
```

Check the submission list without submitting:

```sh
DRY_RUN=1 ./submit_cluster.sh
DRY_RUN=1 SYSTEM_FILTER=small ./submit_cluster.sh
DRY_RUN=1 SYSTEM_FILTER=large ./submit_cluster.sh
```

The submit helper defaults to `build_cluster/bin/sim`, `condo-sabel1`, account `ACF-UTK0049`, QoS `condo`, Bergamo nodes `ber1528,ber1529`, 72 h, 8192 MB, and 4 CPUs. Override with environment variables, for example:

```sh
HOURS=48 MEM=12288 ./submit_cluster.sh
```

First movies to inspect:

- Small system `04` vs `05`: fixed vs growing MTs with perinuclear free actin.
- Small system `06` vs `07`: division-stage fixed vs growing MTs.
- Large system `04` vs `05`: checks whether the qualitative behavior survives the larger geometry.
- Large system `06` vs `07`: larger division-stage test.
