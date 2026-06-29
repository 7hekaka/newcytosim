# Continuous Top-Supply Full-Comparable Campaign

Purpose: rerun the aggressive bottom-turnover model at the older annulus scale so it can be compared directly to the earlier initialized-minus-z c80/m12 simulations.

## Geometry and density
- Annulus: inner radius `10.5`, outer radius `11.0`, z range `[-20.0, 20.0]`.
- Initial actin: `256` filaments, length `15.0`.
- Top supply: `64` new aligned filaments every 60 s for `13` pulses.
- Explicit filament length introduced by config: `16320.0 um`.
- Actin polymer pool: `19008.0 um`.
- Crosslinkers: `1024` initially, `256` added per pulse.

## Turnover and motor model
- Bottom chewer zone: lower 30% of the annulus, z `[-20.0, -8.0]`.
- Chewers: `2048` diffuse bottom chewers, equal to 2x the full-size base chewer count.
- `chewing_speed = 1.6` and `max_chewing_speed = 3.0`.
- Motor case: `80` inner-wall clusters x `12` motors/cluster = `960` motors.
- Conditions: no motors with crosslinkers; rotatable wall motors with crosslinkers.

## Runtime
- Initial block: 60 s.
- Post-onset supply blocks: `13` x 60 s = `780` s.
- Total simulated time: `840` s.
- Replicates: `10` per condition.
- Total configs: `20`.

## Submit on cluster
Compile first on a Bergamo node, then:
```bash
cd /lustre/isaac24/scratch/kacheamp/newcytosim/turnover_continuous_supply_full_comparable
./submit_cluster.sh
```

The submit script defaults to `NODELIST=ber1528,ber1529` to avoid the Milan illegal-instruction issue.
