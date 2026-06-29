# Cytosim Project Planner

Last updated: 2026-06-25

Use this file as the shared working TODO list for collaborator-facing simulation work. Daily rule: pick one item from `Active Next Actions`, do or advance it, then update this file.

Operational path note: the active WSL repo is `/home/thekaka/project/cytosim`, and the Windows UNC path is `\\wsl.localhost\Ubuntu-22.04\home\thekaka\project\cytosim`. Avoid `C:\home\thekaka\project\cytosim`; that was a stale bad Codex path.

## Active Next Actions

- [x] Re-anchor the planner after the 2026-06-11 PI/collaborator priority reset.
- [x] Identify the existing bottom-chewer full-height, top-2/3, and top-1/3 aligned runs with crosslinkers, with and without motors.
- [x] Use selected existing diffuse bottom-chewer results as the starting reference, but treat them as growth plus bottom-localized binders rather than true turnover because old fiber-end chewing was compiled out.
- [x] Design a miniature top-polymerization / bottom-depolymerization F-actin pilot with crosslinkers, with and without membrane-bound rotatable motors.
- [x] Generate and smoke-test local real-minus-end-chewer pilot configs before deciding what to submit on the cluster.
- [x] Rerun the cluster compile after the 2026-06-12 OpenMP/range-for fix and confirm `sim` and `report` build.
- [x] Submit/start the 45-config `turnover_mini_minus_chew/cluster_runs` production set on the cluster.
- [x] Wait for the 45 mini-turnover cluster jobs to finish, then copy/check outputs and run the first analysis pass.
- [x] Generate the scaled continuous top-supply follow-up campaign with small and large elongated systems, matched motor/no-motor cases, and density-scaled filament/crosslinker/chewer counts.
- [x] Move/force-add the scaled continuous-supply campaign to the cluster, submit the 20 configs, then analyze whether sustained top supply prevents depletion while bottom minus-end chewing clears transported material. Job09 returned and was analyzed on 2026-06-24 in `analysis/results/turnover_continuous_supply_scaled_job09/`.
- [x] Submit/analyze the aggressive continuous-supply depolymerization follow-up. Jobs 10-12 completed and were analyzed; job13 sever-chew crashed early with segmentation faults after partial output. Combined comparison is in `analysis/results/turnover_continuous_supply_scaled_aggressive_comparison/`.
- [ ] Decide whether to run a gentler severing branch, because the first true-cutter setup was too aggressive/unstable. Candidate change: reduce severer count and/or `cutting_rate`, then smoke-test before cluster submission.
- [x] Microscopy branch: inspect the initial collaborator TrackMate XIG summary workbook and separate what can be concluded from what still needs raw trajectories.
- [x] Microscopy branch: draft a short XIG movement meeting brief for the planned Tuesday 2026-06-16 postdoc discussion.
- [x] Microscopy branch: inspect collaborator confocal image/movie folder, identify file formats/metadata, and draft the actin-speed audit.
- [ ] Microscopy branch: reproduce the collaborator's Fiji/ImageJ kymograph measurements on representative files, saving ROIs/kymographs/results and testing the time-axis ambiguity.
- [x] Microscopy branch: first WT Fiji/ImageJ reproduction pass completed on `1_0001`, `1_0006`, and `1_0007`, with four manual traces per movie and comparison against the prototype flow metric.
- [x] Microscopy branch: first independent network-scale actin-flow workflow implemented and run on all readable WT/CA/DN `.ims` files, with movie-level summaries and QC overlays.
- [x] Microscopy branch: export and analyze scar4-1 `.oib` movies through Bio-Formats/Fiji Java, then run the same network-scale flow analysis on projected TIFF stacks.
- [x] Microscopy branch: standardize the collaborator-facing F-actin flow protocol with Fiji export steps, manifest format, Python command, QC checks, and reporting rules.
- [ ] Microscopy branch: validate the network-scale workflow against Bio-Formats/Fiji exported projected TIFF/OME-TIFF stacks so absolute calibration matches the series used manually.
- [x] Endosperm pivot: design a tethered-microtubule plus free-actin pilot and choose the first small parameter set.
- [x] Endosperm pivot: summarize returned first tethered-MT/free-actin pilot fiber lengths.
- [x] Endosperm pivot: generate and smoke-test the 5 min fixed-vs-growing MT/free-actin long pilot with doubled actin polymer pool.
- [x] Endosperm pivot: generate and smoke-test the 5-replicate small/large fixed-vs-growing MT/free-actin cluster campaign.
- [x] Endosperm pivot: full 5-replicate small/large fixed-vs-growing MT campaign returned in `ast/tethered_mt_free_actin_pilot/job02/save` and analyzed.
- [ ] Ask collaborators for egg-cell geometry inputs from images: approximate dome/oval dimensions, boundary shape, and whether there is a large vacuole.
- [ ] Convert the collaborator F-actin setup questions into a short data-request table.
- [ ] Submit canonical paper-geometry comparison pair when it becomes relevant again.
- [ ] When PI feedback arrives, update the manuscript plan and figure wording before making new manuscript plots.

## 2026-06-11 Meeting Priority Reset

New priority order:

- F-actin: the aligned-filament results convinced the group that polarity/alignment matters. Next question is whether top-biased polymerization plus bottom-biased depolymerization can generate or maintain the relevant aligned network, both with and without motors.
- Existing bottom-chewer campaigns should be mined before creating new runs. Keep crosslinkers; ignore no-crosslinker cases for this question. Include full-height aligned, top-2/3 aligned, and top-1/3 aligned cases.
- Start with a miniature system containing the relevant ingredients before scaling to the larger central-cell runs.
- Endosperm: pause the actin-only tethered branch as the main route. Test tethered microtubules instead, with actin placed freely in the system and not directly bound to microtubules.
- Collaborator support: help replace rough kymograph-based speed estimates with a reproducible image-analysis workflow for actin movement.

Selected existing turnover runs:

- Table: `analysis/results/meeting_2026-06-11_polymerization_priority/selected_existing_runs.csv`
- Quick comparison plot: `analysis/results/meeting_2026-06-11_polymerization_priority/selected_existing_runs_metrics.png`
- Main existing figure bundle: `analysis/results/pi_meeting_turnover_2026-04-23/`

## Collaborator TODO: F-actin Simulation Setup

Status: open, with central-cell geometry resolved

Their listed items:

- [x] 3D volume for central cell simulation.
- [ ] F-actin length and numbers.
- [ ] F-actin bundling orientations.
- [ ] Linker details.
- [ ] Anything else needed for setup.
- [ ] What information is needed from collaborators.

Current geometry decisions:

- Central cell: current annular geometry is acceptable and matches the central-cell setup.
- Egg cell: likely needs a dome-shaped oval geometry rather than the annular central-cell geometry.
- Egg cell probably does not have the large vacuole used in the central-cell model, but this should be confirmed from collaborator images.

Information needed from collaborators:

- Egg-cell geometry: image-derived dimensions, dome/oval shape constraints, boundary assumptions, and whether a large vacuole should be included.
- Initial F-actin network: filament count or density, length distribution, polarity/orientation information, bundle width, bundle persistence, and any known spatial enrichment.
- F-actin dynamics: turnover rate, severing/disassembly evidence, nucleation zones, growth/shrinkage rates, and whether filaments should treadmill.
- Linkers/crosslinkers: molecular identity if known, concentration or ratio to actin, binding/unbinding rates if available, stiffness estimate, and whether bundling is static or remodels.
- Boundaries and anchoring: whether actin is anchored at membrane, around nuclei, around organelles, or free in cytoplasm.
- Readouts they care about: nuclear spacing, actin density, organelle exclusion, bundle orientation, velocity, or qualitative movie match.

Our next useful action:

- [ ] Make a one-page collaborator request table with columns: parameter, why needed, acceptable estimate if unknown, current assumption. Prioritize egg-cell dimensions first.

## Collaborator TODO: Myosin

Status: open

Their listed items:

- [ ] Speed of myosin movement.
- [ ] What data is needed from collaborators.

Information needed from collaborators:

- Myosin class or best biological proxy.
- Measured or literature speed in this system.
- Directionality relative to F-actin polarity.
- Localization: membrane-bound, organelle-bound, cytosolic, or clustered.
- Cluster size/distribution if membrane-bound.
- Binding lifetime/processivity and force sensitivity if known.
- Whether myosin activity should be fixed, turn over, diffuse, or reorganize.

Current model reference:

- Fixed-global simulations used `preferred_direction = 0, 0, -1` and `angular_tolerance = 0.35` rad, about 20 degrees.
- Rotatable simulations allow motor orientation to follow local filament geometry.

Our next useful action:

- [ ] Draft a myosin-parameter table that separates known data, literature assumptions, and free sweep parameters.

## Collaborator TODO: Organelles

Status: open

Their listed items:

- [ ] Simulate F-actin dynamics when movement is not affected by peroxisomes.
- [ ] Simulate F-actin dynamics when movement is not affected by mitochondria.
- [ ] ER.

Interpretation to clarify:

- "Movement is not affected by organelle" could mean organelles are ignored, included as passive excluded volume, or included visually but without force coupling.
- Need to know whether peroxisomes, mitochondria, and ER should act as obstacles, cargo, boundaries, or just analysis masks.

Information needed from collaborators:

- Organelle size/shape distributions and approximate counts.
- Spatial maps or representative images for peroxisome, mitochondria, and ER localization.
- Whether actin contacts or avoids each organelle type.
- Whether organelle motion should be measured, constrained, or ignored.
- Which organelle condition is highest priority.

Our next useful action:

- [ ] Propose three model levels: no organelles, passive excluded-volume organelles, and force-coupled organelles.

## Endosperm Nuclear Positioning Project

Status: active, pivoting from free-actin/MT scaffold tests toward nucleus-bound actin asters and asymmetric turnover-driven nuclear motion

Known notes from `ast/endosperm_method_notes.md`:

- Local activation/assembly zones may matter more than uniformly random free filaments.
- Active turnover should be tested explicitly in division cases.
- Collaborator interpretation now favors actin asters that are nucleated or anchored at nuclei, with microtubules playing a weaker role in initial aster formation than assumed in the first pilot.
- New working hypothesis: rapid, spatially asymmetric actin turnover near one side of a nucleus may help move nuclei toward the chalazal/ACTIN10-rich region after the aster morphology breaks down.
- Current ACTIN8 has `growing_force = inf`, so force-dependent polymerization is simplified.
- Steric/crowding and mesh mechanics may matter, but steric-on should follow promising smoke tests.
- Static crosslinker counts may miss elastic gel-like force transmission.
- Nuclear division should be staged as two nuclei becoming four, with actin assembly continuing through division.

Immediate TODOs:

- [x] Inspect local Cytosim microtubule examples and choose a starting microtubule rigidity. Current local examples use microtubule-like rigidities around `10` to `30`; current actin runs use `0.075`. First MT/free-actin pilot uses MT rigidity `20`.
- [x] Build a miniature tethered-microtubule plus free-actin pilot: microtubules tethered/anchored, actin free in the volume, no direct actin-microtubule binding initially. Current pilot is `ast/tethered_mt_free_actin_pilot/`.
- [x] First returned pilot summarized in `ast/tethered_mt_free_actin_pilot/analysis/job01_final_fiber_lengths.csv`. This confirmed the first pilot was short and actin-limited: the division run ended with ACTIN total length `369.1` um, ACTIN max length `8.0` um, and MT max length `15.855` um.
- [x] Built the second-pass long pilot in `ast/tethered_mt_free_actin_long_pilot/`: 5 min runs, doubled actin total polymer pool `1200`, actin `max_length = 16`, fixed MT length `6.0` um controls, and matched growing-MT cases. `00_smoke_parse/config.cym` passed locally on 2026-06-16; `submit_cluster.sh` selects the seven non-smoke configs.
- [x] Built the replicated cluster campaign in `ast/tethered_mt_free_actin_long_replicates/`: 70 production configs = 2 systems x 7 conditions x 5 replicates. `small_system` is `12 x 12 x 1.2` um with fixed MT length `6` um and actin total polymer `1200` um. `large_system` is `24 x 24 x 1.2` um with fixed MT length `12` um, material scale `4`, and actin total polymer `4800` um. Both smoke configs passed locally on 2026-06-16.
- [x] First returned long-replicate batch analyzed in `ast/tethered_mt_free_actin_long_replicates/analysis/returned_save_first_pass/`. Returned outputs were copied into `ast/tethered_mt_free_actin_pilot/save`. Small system was complete; large system was partial.
- [x] Full `job02` long-replicate sweep analyzed in `ast/tethered_mt_free_actin_long_replicates/analysis/job02_full_sweep/`: 70/70 completed runs = 2 systems x 7 conditions x 5 replicates. Main result: fixed MT cases preserve MT length and compact nuclear/body spacing in both system sizes, while the current growing-MT parameterization collapses MTs to short lengths and produces larger body spacing. This is a real qualitative difference, but also indicates the growing-MT parameter set is too depolymerizing if the intended comparison is fixed long MTs versus growing long/spanning MTs.
- [x] Built and ran a 15-trajectory quick nuclear-aster turnover-propulsion pilot in `ast/nuclear_aster_turnover_propulsion_pilot/`: 5 conditions x 3 replicates, 60 s each. Initial result: asymmetric fast polymerization without strong chewing gives the clearest signed nucleus drift (`dx = +0.704 +/- 0.054 um` SEM), while the current rapid-turnover/chewer cases drift negative or are weaker (`dx = -0.377 +/- 0.158 um` for symmetric fast turnover; `dx = -0.161 +/- 0.091 um` for polar fast turnover). The no-nucleator control shows large passive wandering, so future propulsion tests need a better drag/anchoring baseline and more biological turnover geometry before claiming turnover-driven propulsion.
- [x] Built and ran a stricter one-sided short-filament turnover pilot in `ast/nuclear_aster_one_sided_turnover_pilot/`: 7 conditions x 3 replicates, 60 s each. This fixes the earlier geometry issue by using `type=angular`, caps actin at `max_length = 2.5 um`, and tests fixed one-sided chewers either from t=0 or delayed until t=4 s. Main result: the one-sided cap geometry itself dominates the displacement (`dx = -0.672 +/- 0.078 um` without chewers; `dx = -0.817 +/- 0.169 um` with chewers). After the 4 s turnover onset, delayed cap+chewer motion is not stronger than the cap-only relaxation (`dx_4to60 = +0.172 +/- 0.070 um` versus `+0.224 +/- 0.064 um`). Chewers bind locally and the filament length cap holds, but this still does not give a clean turnover-only propulsion signal.
- [x] Built a long chewer-free treadmilling-only pilot in `ast/nuclear_aster_treadmilling_long_pilot/`: 12 production configs = symmetric dense aster, minus-x cap, plus-x cap direction-flip control, and no-nucleator baseline, each with 3 replicates. Setup: larger `32 x 16 x 2` rectangle, 120 nucleus-bound actin filaments, initial length `0.1 um`, `max_length = 1.0 um`, `total_polymer = 12000`, no crosslinkers, no chewers, 5 min duration. Smoke parse/run passed locally; by 2 s the minus-x cap smoke grew from `0.1` to `0.893 um`, confirming filaments grow toward the 1 um cap while treadmilling turnover is active.
- [x] Analyzed the returned `job03` long treadmilling-only pilot in `ast/nuclear_aster_treadmilling_long_pilot/analysis/job03/`: 12/12 runs completed 5 min. The active aster cases all ended with 120 actin filaments at the 1.0 um length cap. The one-sided cap flips displacement sign as expected (`minus-x cap dx = -0.108 +/- 0.079 um`; `plus-x cap dx = +0.136 +/- 0.090 um`; SEM, n=3), while the symmetric dense aster remains near zero (`dx = +0.003 +/- 0.029 um`). This supports a weak geometry-dependent directional bias, but not yet a strong treadmilling-driven propulsion effect. The no-nucleator control wanders substantially as a freely diffusing confined solid and should not be used as the main propulsion baseline.
- [x] Built and smoke/quick-tested a mechanism-screen campaign in `ast/nuclear_propulsion_mechanism_mockups/`: 10 conditions x 3 production replicates plus 20 s quick-look configs. The 20 s one-replicate screen shows pure geometry remains weak (`plus-x cap no external dx = +0.053 um`), lab-frame polymerization-push cases are near zero, fixed non-motor binders can strongly pull a plus-x cap (`dx = +0.881 um`), fixed plus-end motors pull the cap in the opposite direction (`dx = -0.328 um`), and fixed minus-end motors give the strongest +x motion (`plus-x cap dx = +1.597 um`; symmetric aster dx = +1.374 um`). Immediate interpretation: this points toward external actin anchoring/motor traction and filament polarity as the dominant ingredients, not treadmilling-only propulsion.
- [ ] Test whether free actin forms the expected aster-like organization around the tethered microtubule structure.
- [ ] Sweep one actin polymerization speed and one reduced-crosslinking condition only after the baseline pilot runs.
- [ ] Keep the old actin-only smoke tests as background: `05_division_bound_asters`, `06_division_free_filaments`, and `07_division_perinuclear_free_turnover`.
- [ ] Replace the crude turnover-propulsion pilot with a more faithful nucleus-bound ACTIN10 cap: one-sided nucleation/turnover zone that moves with the nucleus, lower passive nucleus diffusion or realistic drag, matched polymerization-only versus turnover-only controls, and external/cortical actin traction controls informed by `ast/nuclear_propulsion_mechanism_mockups/`.
- [ ] Add a steric-on version only after smoke tests show promising qualitative behavior.
- [ ] For ACTIN8-to-ACTIN10 production sweeps, keep one fiber class and change polymerization speed instead of renucleating ACTIN10.
- [ ] Build analysis around nuclear spacing, pairwise distance change after division, actin density around nuclei, filament orientation relative to nuclei, and clump/bridge persistence.

Our next useful action:

- [x] Draft the smallest tethered-microtubule/free-actin config and run a local parse/smoke test. `ast/tethered_mt_free_actin_pilot/00_smoke_parse/config.cym` passed locally on 2026-06-16.
- [x] Run/watch the static steric-off and steric-on tethered-MT/free-actin pilots side by side enough to justify the fixed-vs-growing long follow-up.
- [x] Run quick local asymmetric nuclear-aster turnover propulsion test and save endpoint/time-course displacement summaries in `ast/nuclear_aster_turnover_propulsion_pilot/analysis/`.
- [ ] Watch representative movies from `analysis/job02_full_sweep/representative_movies.txt`: small fixed-vs-growing pairs `r0051` versus `r0059` and `r0061` versus `r0070`; large fixed-vs-growing pairs `r0016` versus `r0024` and `r0029` versus `r0035`.

## Current F-actin / Annulus / Motor-Cluster Project

Status: mini top-polymerization/bottom-chewing outputs analyzed; rough/mobile inner-wall branch remains paused while the turnover follow-up is chosen

Done recently:

- [x] Built initialized-minus-z slide deck.
- [x] Removed no-motor row from initialized-minus-z heatmaps/metrics because available no-motor controls were randomly initialized.
- [x] Verified fixed-global angular tolerance: `0.35` rad.
- [x] Implemented `rough_annulus` / imperfect rigid inner wall with seeded irregular radius variation and shaded wall display.
- [x] Added a generator for initialized-minus-z rough-wall campaigns.
- [x] Compiled and pushed cluster-safe source/helper updates; c80/m12, xlink ratio `1:8` rough-wall jobs are running on Bergamo-compatible nodes.
- [x] Analyzed first c80/m12 rough-wall outputs against smooth controls; roughness looked like a small perturbation rather than a strong behavioral change.
- [x] Generated a deeper/visibly rougher c80/m12, 1:8 anchored-motor campaign in `clu_rough_inner_wall_deep_c80_m12_init_minusz`.
- [x] Implemented sliding membrane-localized Picket anchors with `anchor_mode = slide`, preserving force transmission while allowing motor anchor diffusion on the wall.
- [x] Generated smooth/rough x rotatable/fixed-global mobile-anchor c80/m12, 1:8 configs in `clu_mobile_inner_wall_c80_m12_init_minusz`.
- [x] Corrected the canonical 1:8 run selection to use the actual `new 2048 crosslinker` configs: `r0002`, `r0005`, `r0008`, `r0011`, `r0014`, `r0017`, `r0020`, `r0023`, `r0026`, `r0029`.
- [x] 2026-06-09 daily check: identified the exact source/helper files that need staging before cluster rebuild: `CMakeLists.txt`, `src/sim/CMakeLists.txt`, `src/sim/space_prop.cc`, `src/sim/spaces/space_rough_annulus.cc`, `src/sim/spaces/space_rough_annulus.h`, `src/sim/single_prop.cc`, `src/sim/singles/picket.cc`, `analysis/generate_rough_inner_wall_campaign.py`, `analysis/generate_mobile_inner_wall_campaign.py`, and `python/run/compile_cytosim_cluster.sh`.
- [x] 2026-06-11 priority reset: top-biased polymerization and bottom-biased depolymerization are now the main F-actin mechanism to test.
- [x] 2026-06-11 mined existing bottom-chewer campaigns for full-height, top-2/3, and top-1/3 aligned conditions with crosslinkers and with/without motors.
- [x] 2026-06-11 checked bottom-chewer turnover: current source has `NEW_FIBER_END_CHEW = 0`, so chewer-mediated fiber-end depolymerization is compiled out. Existing bottom-chewer runs should be treated as growth plus bottom-localized binders, not true turnover.
- [x] 2026-06-11 enabled minus-end chewing with `NEW_FIBER_END_CHEW = 2` and added optional OpenMP support to the CMake build.
- [x] 2026-06-11 generated the mini minus-end-chewer campaign in `turnover_mini_minus_chew`: validation, short full-height aligned rotatable smoke config, six one-replicate planned cases, and a 45-config cluster production set.
- [x] 2026-06-11 local validation passed: chew-only total actin length dropped from `96.0` to `80.7` um in 30 s.
- [x] 2026-06-11 local full-height aligned rotatable smoke test completed: final report at frame 100 showed total actin length `903.9` um and cumulative chewed/off length `17.8` um. Longer local motor runs are slow and should be cluster jobs or shorter observation pilots.
- [x] 2026-06-11 local build `build_mini_turnover` produced `sim`, `report`, and `play`; `sim info` reports CytoSIM 3D built Jun 11 2026 20:02:32 with GCC 11.4.0.
- [x] 2026-06-11 archived partial local motor pilot under `turnover_mini_minus_chew/pilot_archive/full_height_aligned_rotatable_xlink_partial_20s/`; `fiber_length_partial.txt` shows actin length `320.0 -> 920.9` um and cumulative off length `0.0 -> 19.4` um by frame 109. The run ended with `killed 2` after heavy CPU time, so use it only as a growth/chewing sanity check, not as a completed motor trajectory.
- [x] 2026-06-11 added cluster helper `python/run/compile_cytosim_cluster.sh`: intended to run inside an allocated Bergamo compute shell, uses `CYTOSIM_NATIVE_ARCH=OFF`, `CYTOSIM_ENABLE_OPENMP=ON`, `CHEW_MODE=2`, fresh `build_cluster`, and an optional minus-chew sanity run.
- [x] 2026-06-11 added `turnover_mini_minus_chew/submit_cluster.sh`: submits all configs under `turnover_mini_minus_chew/cluster_runs` using `build_cluster/bin/sim`, defaults to `condo-sabel1`, account `ACF-UTK0049`, QoS `condo`, nodes `ber1528,ber1529`, 48 h, 8192 MB, 4 CPUs.
- [x] 2026-06-12 cluster compile attempt on `ber1528` failed before linking. Log `compile_diagnostics_20260612_003409.log` showed CMake used `/usr/bin/c++` GNU 8.5.0 even though `gcc/10.2.0` was loaded, and the compile failed at OpenMP pragmas placed over C++ range-for loops in `src/sim/meca.cc` and `src/sim/meca_precond.cc`.
- [x] 2026-06-15 fixed the cluster compile blocker locally by changing the OpenMP range-for loops in `src/sim/meca.cc` and `src/sim/meca_precond.cc` to index-based loops, and hardened `python/run/compile_cytosim_cluster.sh` to export/use the module `gcc`/`g++` paths explicitly via `CMAKE_C_COMPILER` and `CMAKE_CXX_COMPILER`.
- [x] 2026-06-15 local verification after the compile fix passed: `cmake --build build_mini_turnover --target sim report --parallel 4` rebuilt cleanly, `bash -n` passed for the compile/submit helpers, and the minus-chew validation total actin length dropped from `96.0` to `79.5` um by frame 30 with `16.5` um cumulative off length.
- [x] 2026-06-14/15 cluster compile succeeded after the CRLF/shebang and OpenMP/range-for fixes. The diagnostic log showed `targets: sim report`, `Built target sim`, `Built target report`, and `build_cluster/bin/sim help` ran successfully. The optional minus-chew sanity run was skipped only because the validation config was not present in that cluster copy.
- [x] 2026-06-15 user report: the mini-turnover production set is now running on the cluster. Treat this branch as outputs-pending until the jobs finish.
- [x] 2026-06-16 analyzed the returned mini-turnover job00 outputs in `analysis/results/turnover_mini_minus_chew/`. All 45 runs parsed successfully. Fixed-global motors gave the strongest post-motor axial transport, especially `top_two_thirds_aligned + fixed_global_xlink` with AUC(|v_z|) `9.386 +/- 0.516` um SEM and excess AUC `2.154 +/- 0.320` um above matched no-motor control. Rotatable motors produced weak bulk axial transport but best reduced final chewer-band mass, especially `top_two_thirds_aligned + rotatable_xlink` with final chewer-band fraction `0.048`. Exact `fiber:length` summaries confirm chewing is active but polymer length remains near `1012` um at the end, so this is not yet a balanced steady-turnover regime.
- [x] 2026-06-17 generated the scaled continuous top-supply follow-up in `turnover_continuous_supply_scaled/`: 20 cluster configs = 2 sizes x 2 conditions x 5 replicates. `small_long` uses `z = [-7.5, 7.5]`, 12 initial filaments, 3 top-supply filaments/min, 48 chewers, and 6 rotatable motor clusters when motors are present. `large_long` uses `z = [-15, 15]`, 96 initial filaments, 24 top-supply filaments/min, 384 chewers, and 23 rotatable motor clusters when motors are present. Filaments, polymer pool, crosslinkers, and chewers scale by annular volume; motor clusters scale by inner-wall surface area.
- [x] 2026-06-24 analyzed returned continuous-supply job09 in `analysis/results/turnover_continuous_supply_scaled_job09/`. All 20 runs parsed. Rotatable motors significantly increased axial transport in both sizes but also increased bottom/chewer-zone mass accumulation; mass retention remained >1.5 in all conditions, so the bottom chewers are still not strong enough to create balanced steady-state turnover.
- [x] 2026-06-24 generated aggressive-depolymerization follow-ups: `turnover_continuous_supply_scaled_chew_rate2x` keeps chewer count fixed but sets `chewing_speed = 1.6` and `max_chewing_speed = 3.0`; `turnover_continuous_supply_scaled_chewer_count2x` keeps per-chewer speed fixed but doubles diffuse bottom chewers; `turnover_continuous_supply_scaled_chewer_count2x_rate2x` doubles both chewer count and per-chewer chewing rate; `turnover_continuous_supply_scaled_sever_chew` adds true severing cutters in the bottom zone alongside baseline minus-end chewers. Each folder has 20 configs and its own `submit_cluster.sh`; the sever-chew config passed a shortened local smoke run.
- [x] 2026-06-25 analyzed returned aggressive turnover jobs: `job10 = chew_rate2x`, `job11 = chewer_count2x`, `job12 = chewer_count2x_rate2x`, and `job13 = sever_chew`. Jobs 10-12 parsed fully. The overkill count+rate case strongly reduced final chewer-band accumulation in rotatable runs (`large_long` 0.175 -> 0.036; `small_long` 0.461 -> 0.096) and increased off length, but also reduced AUC(|v_z|), especially in the small system. Job13 sever-chew exited with status 11 segmentation faults after only partial frames, so it should not be interpreted as a completed biological result.

Current mini-turnover campaign state:

- Source: `src/sim/fiber_prop.h` is currently configured with `NEW_FIBER_END_CHEW = 2`, enabling minus-end chewing.
- Campaign directory: `turnover_mini_minus_chew/`.
- Generation script: `turnover_mini_minus_chew/generate_turnover_mini_minus_chew.py`.
- Manifest: `turnover_mini_minus_chew/manifest.csv`.
- Validation config: `turnover_mini_minus_chew/validation/minus_end_chew/config.cym`.
- Short smoke config: `turnover_mini_minus_chew/pilot/full_height_aligned_rotatable_xlink_smoke/config.cym`, with only 10 s motor time.
- Production configs: `turnover_mini_minus_chew/cluster_runs/`, 45 configs total = 3 scenarios x 3 conditions x 5 replicates.
- Production scenarios: `full_height_aligned`, `top_two_thirds_aligned`, and exploratory `top_two_thirds_mixed_polarity`.
- Production conditions: `nomotor_xlink`, `rotatable_xlink`, and `fixed_global_xlink`.
- Timing in production configs: chewers plus growth for 120 s, crosslinkers for 60 s, motors for 600 s where present.
- Chewer setup: 256 diffuse bottom chewers in `z = [-10, -4]` um, `bind_only_end = minus_end, 0.5`, `hold_growing_end = 0, 1`, `chewing_speed = 0.8`, and fiber `max_chewing_speed = 1.5`.
- Motor setup: 15 wall clusters x 12 motors; fixed-global motors use `preferred_direction = 0, 0, 1` for aligned-up fibers in this mini campaign.
- Cluster output: `turnover_mini_minus_chew/job00/save/` contains 45 completed run directories, each with `objects.cmo`, `messages.cmo`, `properties.cmp`, `config.cym`, and `log.txt`.
- Analysis bundle: `analysis/results/turnover_mini_minus_chew/`, including replicate-level `run_metrics.csv`, condition-level `condition_metrics.csv`, exact `fiber_length_condition_summary.csv`, movement matrices, mean/SEM timecourses, medoid kymographs, and unwrapped line snapshots.
- Working interpretation: fixed-global motors can drive coherent downward transport into the bottom/chewer zone, but chewing capacity does not clear the transported material fast enough. Rotatable motors suppress bottom/chewer-zone accumulation but do not create strong directed transport in this setup.

Open modeling hypotheses and implementation order:

- [x] 1. Imperfect rigid inner wall: inner radius has fixed spatial variation along circumference and/or z. Most feasible; implemented first.
- [x] 2. Mobile or dispersible motor clusters: motors begin clustered but activity/rebinding allows redistribution. Intermediate difficulty; implemented as sliding membrane-localized Picket anchors.
- [ ] 3. Flexible/fluctuating inner wall: membrane boundary deforms over time. Highest difficulty; implement after the static and mobile-motor cases clarify what readouts matter.
- [ ] Matched initialized-minus-z no-motor controls only if collaborators ask for a formal baseline.

Our next useful action:

- [x] Use the existing bottom-chewer analysis to choose the first miniature polymerization/depolymerization pilot.
- [x] Before new turnover runs, enable or replace actual depolymerization. Options: set `NEW_FIBER_END_CHEW` for the relevant end(s), implement a local cut/shrink mechanism that works from the bottom zone, or use explicit minus-end shrink kinetics with correct polarity.
- [x] Prefer `sc02_top_two_thirds_aligned` with diffuse bottom chewers as the starting reference because it gave the clearest rotatable-motor gain over matched no-motor and fixed-global controls, while noting those old runs were not true chewing turnover.
- [x] Inspect local smoke/partial trajectory enough to confirm actual growth plus chewing and decide that longer motor readouts should run on the cluster.
- [x] Generate the reduced local/cluster pilot comparison with crosslinkers, no motors vs rotatable membrane-bound motors vs fixed-global motors, plus a mixed-polarity exploratory control.
- [x] Re-run `./python/run/compile_cytosim_cluster.sh` on the cluster after the OpenMP range-for and compiler-selection fixes; check that `build_cluster/bin/sim` and `report` build.
- [x] Submit/start `turnover_mini_minus_chew/cluster_runs` on the cluster.
- [x] When cluster jobs finish, copy/check outputs and run the first analysis pass.
- [x] After outputs return, analyze total actin length/off length, top-bottom actin density, filament polarity/orientation, bundle persistence, and motor/no-motor differences across the 45 production configs.
- [ ] Design the focused follow-up from the mini-turnover result: keep `top_two_thirds_aligned` as the primary scenario, carry `fixed_global_xlink` as the directional-transport positive case, carry `rotatable_xlink` as the clearance/low-chewer-band comparison, and test stronger or better-localized bottom depolymerization so transported material does not pile up.
- [ ] Decide whether the next follow-up should change chewer capacity first (`n_chewers`, `chewing_speed`, chewer-zone height, or `max_chewing_speed`) or add a lower-polymer-pool / shorter-growth control before scaling up.
- [ ] Submit or transfer `turnover_continuous_supply_scaled/` to the cluster. Because `.gitignore` ignores new files by default, either force-add it before pushing or copy the folder directly.
- [ ] After the scaled campaign returns, compare no-motor vs rotatable-motor runs for total polymer length, cumulative off length, top/bottom actin density, chewer-zone mass, axial velocity, and whether visible depletion is reduced.
- [ ] Keep deeper roughness and mobile-anchor campaigns on pause unless the polymerization/depolymerization branch fails to explain the collaborator phenotype.

## Experimental Actin Speed Measurement Support

Status: active next branch while mini-turnover cluster jobs run; first egg-cell velocity audit completed

Problem:

- Collaborators are currently estimating actin movement speed with kymographs, but the method is likely sensitive to line choice, projection, bundle crossing, and manual slope picking.

Proposed analysis direction:

- Start with movies or TIFF stacks from collaborators and define cell/annulus masks plus top-bottom axis orientation.
- Correct stage drift and intensity changes before measuring motion.
- Use dense optical flow or PIV-style local cross-correlation to estimate a velocity field over the actin channel.
- Summarize axial speed distributions rather than one hand-picked kymograph slope: median, interquartile range, fraction moving toward bottom/top, and region-specific speeds.
- Validate the workflow on simulation movies where the expected direction and approximate speeds are known.

Egg-cell F-actin velocity folder audit:

- Folder inspected: `/home/thekaka/project/cytosim/Kwaku and Berry`.
- Folder contents: `DN ROP9`, `CA ROP9`, `scar4-1`, and `Wild type` folders, plus `Egg cell mutant analysis.pptx`, `Protocol.jpeg`, `Velocity_Measurement_Tool-macro.txt`, and `Veleocity from single sample .csv`.
- File types: 22 Imaris `.ims` files, 5 Olympus `.oib` files, 1 PowerPoint, 1 protocol image, 1 ImageJ macro, and 1 example CSV.
- Audit note: `analysis/results/egg_cell_velocity_audit_2026-06-15/README.md`.
- Metadata inspector: `analysis/inspect_egg_cell_ims_metadata.py`.
- Metadata outputs: `analysis/results/egg_cell_velocity_audit_2026-06-15/ims_metadata/ims_metadata_summary.csv` and `ims_metadata.json`.
- First provisional network-scale flow prototype: `analysis/prototype_egg_cell_flow.py`, with WT/CA ROP9 outputs in `analysis/results/egg_cell_flow_prototype_2026-06-15/`. Treat this as QC only until masking, drift correction, axis definition, and frame interval are resolved.
- WT Fiji/kymograph reproduction pass: user-generated CSVs `c1a-d`, `c6a-d`, and `c7a-d` in `Kwaku and Berry/Wild type/`, summarized in `analysis/results/egg_cell_kymograph_reproduction_2026-06-15/`.
- WT reproduction result: manual final cumulative macro speeds were highly trace-dependent. Movie means were `1.595`, `2.389`, and `2.247` for `1_0001`, `1_0006`, and `1_0007`; the same movies had prototype flow medians of `0.829`, `0.414`, and `0.663` um/frame after rescaling with the Fiji-reported pixel sizes. Treat this as evidence of manual-trace sensitivity, not a final biological comparison.
- Updated network-scale flow workflow: `analysis/analyze_egg_cell_network_flow.py` adds foreground masking, optional global drift correction, local block cross-correlation vectors, per-movie summaries, and QC overlays. Pilot/all-IMS outputs are in `analysis/results/egg_cell_network_flow_pilot_2026-06-15/` and `analysis/results/egg_cell_network_flow_all_ims_2026-06-15_v2/`.
- First all-readable-IMS movie-level result: using direct `.ims` HDF5 reading and reporting movie median `um/frame`, WT `n=7` mean `0.407`, CA ROP9 `n=10` mean `0.395`, and DN ROP9 `n=5` mean `0.593`. Exact movie-level permutation checks gave WT vs CA p `0.826`, WT vs DN p `0.0227`, and CA vs DN p `0.0186`. Treat this as a method-development screen until Bio-Formats exported stacks confirm calibration and the true frame interval is resolved.
- scar4-1 `.oib` handling: added `analysis/ExportOibZProject.java` to use Fiji/Bio-Formats jars directly, max-project z slices 5-7, and write projected TIFF stacks plus metadata. Exported scar4-1 stacks are in `analysis/results/egg_cell_scar4_oib_exports_2026-06-16/`.
- First scar4-1 movie-level result: scar4-1 `n=5` projected TIFF movies had mean movie median speed `0.302 +/- 0.089` um/frame, median `0.291`, range `0.206-0.412`. Combined preliminary plot/table are in `analysis/results/egg_cell_network_flow_combined_with_scar4_2026-06-16/`. In the movie-level permutation screen, scar4-1 vs WT p `0.0707`, scar4-1 vs CA p `0.152`, and scar4-1 vs DN p `0.00794`.
- Confirmed 30 s/frame speed conversion: combined movie-level speed outputs are in `analysis/results/egg_cell_network_flow_combined_with_scar4_2026-06-16/combined_line_summary_30s.csv` and `movie_median_speed_um_s_30s_by_line_with_scar4.png`. Mean movie median speeds are WT `0.01357 +/- 0.00287` um/s, CA ROP9 `0.01318 +/- 0.00405` um/s, DN ROP9 `0.01975 +/- 0.00548` um/s, and scar4-1 `0.01007 +/- 0.00295` um/s. Linear 30 s conversion leaves permutation p-values unchanged.
- Standardized protocol: `analysis/standardized_actin_flow_protocol.md`; Windows collaborator protocol: `analysis/collaborator_actin_flow_protocol_windows.md`; manifest template: `analysis/egg_cell_actin_flow_manifest_template.csv`; manifest-aware analysis command is now supported by `analysis/analyze_egg_cell_network_flow.py --manifest <csv> --output <results_dir>`.
- Collaborator handoff package: `analysis/collaborator_actin_flow_package/` contains `README.md`, `analyze_egg_cell_network_flow.py`, `summarize_actin_flow_results.py`, `actin_flow_requirements.txt`, and `manifest_template.csv`.
- Direct `.ims` caveat: Python HDF5 reading can expose a different internal array size/calibration than the Bio-Formats/Fiji series selected manually. For final quantitative reporting, export standardized projected stacks from Fiji/Bio-Formats and rerun the network-flow script on those exported stacks with explicit pixel size/frame interval.
- Current local reader status: WSL Python did not have `h5py` initially, so `.ims` metadata was extracted with a temporary environment at `/tmp/cytosim_img_env`; `.oib` scar4-1 files need Fiji/Bio-Formats or OME-TIFF export.
- PowerPoint notes: Lifeact-Venus egg-cell movies, mostly 60x lens, z-project slices 5-7, slice width 2-3, and stated `Time interval - 30 sec`. The slides describe drawing a contour line to make a kymograph, then drawing segmented lines on visible F-actin cable ridges for velocity.
- Protocol image notes: Z project, repeated auto brightness/contrast, Gaussian blur radius 0.6, rolling-ball background subtraction radius 3, then ImageJ kymograph/velocity macro. Handwritten notes emphasize selecting moving cable segments/segmented lines.
- Macro behavior: `Velocity_Measurement_Tool-macro.txt` computes `actual speed = abs(dx) / abs(dy)` and `average speed = cumulative sum_dx / cumulative sum_dy`; exact `dy = 0` is forced to `dy = 1`. It discards direction, is sensitive to axis calibration, and can generate large artificial speeds from tiny `dy`.
- Example CSV red flag: one sample includes extreme segment speeds such as `20`, `32`, and `35` where `dy_now` is very small, consistent with manual-slope/outlier inflation.
- Axis ambiguity: the slides say kymograph x-axis is time and y-axis is distance, but the macro's `dx/dy` formula is only a velocity if x is distance and y is time. Confirm the actual Fiji/ImageJ kymograph axis convention before trusting old values.
- Metadata findings from `.ims`: one channel, mostly 10 time points, common data shape `8 x 512 x 512` or `8 x 512 x 256` stored as `Z, Y, X`, recoverable physical extents, and inferred x/y pixel sizes roughly `0.05-0.13 um/pixel`.
- Frame interval resolved by collaborator/grad-student confirmation on 2026-06-16: use `30 s/frame` for biological speed conversion. Imaris `.ims` time stamps still do not show this cleanly, so keep noting that file metadata is unreliable for timing.
- Current interpretation: the collaborator concern is valid. The existing kymograph approach can plausibly disagree with visual inspection because of manual cable selection, mixed/crossing cables, projection/preprocessing choices, direction loss, slope outliers, replicate handling, and unresolved time/axis calibration.
- Recommended replacement pipeline: build a movie metadata manifest; standardize projection and drift correction; estimate velocity fields with PIV/local cross-correlation or optical flow; retain a calibrated kymograph branch only for continuity/QC; aggregate by movie/sample as the replicate unit; export vector overlays and kymograph overlays for visual validation.

Two-track plan for next session:

Track A: reproduce the collaborator's Fiji/ImageJ result first.

- Open a small representative set in Fiji/ImageJ: at minimum one WT `.ims`, one CA ROP9 or DN ROP9 `.ims`, and one scar4-1 `.oib` if Bio-Formats opens it cleanly.
- Record Fiji's imported calibration for each file: x/y pixel size, z step, frame interval, channel count, z slices, time points, and whether Bio-Formats reports the same timing as the Imaris metadata.
- Run the collaborator protocol as faithfully as possible: z project slices 5-7, brightness/contrast autos, Gaussian blur radius `0.6`, rolling-ball background radius `3`, then the provided `Velocity_Measurement_Tool-macro.txt`.
- Save all intermediate artifacts: projected movie, kymograph image, cell-contour ROI, segmented velocity ROI, ImageJ results table, and notes on manual choices.
- Recompute the exact quantities they used: segment `actual speed`, cumulative `average speed`, per-cable mean, and per-sample mean. Keep the outliers at first so the reproduction is honest.
- Run sensitivity checks on the same kymographs: 30 sec frame interval vs file metadata interval, x/y axis interpretation, with/without tiny-`dy` outliers, repeated manual tracing of the same cable, and different visible-cable choices.
- Deliverable: a reproduction table showing reported-like values, our reproduced values, and which choice changes the result most.

Track B: make our independent measurement.

- Convert or read movies into a common format, preferably OME-TIFF, with a metadata manifest beside the image data.
- Use a fixed preprocessing recipe: projection choice, drift correction, intensity normalization if needed, and an egg-cell mask/top-bottom axis where possible.
- Measure at network scale with PIV/local cross-correlation or optical flow rather than only a few hand-picked filaments.
- Summarize per movie/sample: median speed, IQR, signed axial velocity if an axis is known, fraction of regions moving each direction, and low-confidence/low-signal fraction.
- Keep a controlled cable/kymograph branch only as a validation or continuity comparison: predefined sampling rules, saved ROIs, signed slopes, minimum segment length, and rejection of tiny time-axis steps.
- Export visual QC for every movie: projection preview, drift-correction check, velocity-vector overlay, and any kymograph overlays used.
- Aggregate statistics by movie/sample, not by every cable segment, to avoid pseudo-replication.
- Deliverable: a side-by-side comparison of Fiji/manual kymograph results versus network-scale velocity summaries, with a recommendation on which metric is defensible.

Decision after both tracks:

- If manual kymographs are still scientifically useful, write a new Fiji macro that enforces calibration, saves ROIs, keeps direction, rejects tiny-`dy` artifacts, and outputs a clean per-movie summary.
- If network-scale flow better matches visual inspection and is more reproducible, use Fiji only for export/QC and put the main quantitative pipeline in Python.
- If neither gives stable results, report that the movies are not suitable for velocity quantification under the current acquisition/export conditions and define what new acquisition metadata/settings are needed.

Initial XIG TrackMate workbook:

- File inspected: `C:\Users\User\Downloads\20250812_#521_h.xlsx`.
- Reproducible parser/summary script: `analysis/analyze_xig_trackmate_summary.py`.
- Results folder: `analysis/results/xig_cluster_speed_2026-06-15/`.
- Meeting brief: `analysis/results/xig_cluster_speed_2026-06-15/meeting_brief_2026-06-16.md`.
- Workbook contents: two TrackMate track-summary sheets, `#521_WT` and `#h_abcc8`, plus one embedded image. This is not the full per-frame trajectory export.
- Important parsing note: collaborator-added `Sample` and `File` columns shift the TrackMate data two columns to the right; `File` must be forward-filled down the sheet before aggregating tracks by movie.
- Track-level TrackMate mean speed supports the collaborator's qualitative claim: WT `n=13`, mean `0.1145`; abcc8 `n=11`, mean `0.1898`; abcc8 is about `1.66x` WT by track-level means, exact permutation p about `0.00044`.
- Movie/file-level means are the more conservative replicate unit: WT `n=5` files, mean `0.1176`; abcc8 `n=5` files, mean `0.1743`; abcc8 is about `1.48x` WT, exact permutation p about `0.071` for mean difference.
- Interpretation: the direction is consistent with faster XIG cluster motion in abcc8, but track-level statistics likely overstate evidence because multiple tracks come from the same movie. Use movie/file-level statistics or a mixed model if more data arrive.
- Limitation: this workbook cannot support velocity autocorrelation, reversal frequency, MSD, directional persistence, or correlation with membrane/tonoplast fluctuations. Need TrackMate `Spots`, `Tracks`, and `Edges` exports or equivalent per-frame trajectories with track ID, frame/time, x, y, z if available, quality/intensity, pixel size, frame interval, projection status, and drift-correction status.

Immediate TODOs:

- [x] Inspect the collaborator image/movie folder once available: list file types, dimensions, channel count, time/z dimensions, pixel size, frame interval, and whether the data are raw stacks or projections.
- [x] Resolve the frame-interval conflict with collaborators: use `30 s/frame` for speed conversion; note that `.ims` file time stamps are unreliable for this dataset.
- [ ] Confirm the actual kymograph axis convention used by Fiji/ImageJ for their macro session.
- [ ] In Fiji/ImageJ, reproduce the collaborator workflow on one WT file and one mutant file, saving projected movies, kymographs, ROIs, and results tables.
- [x] In Fiji/ImageJ, reproduce the collaborator workflow on WT files `1_0001`, `1_0006`, and `1_0007`; next repeat on a mutant file using the same naming/saving scheme.
- [ ] Build a reproduction/sensitivity table comparing 30 sec timing vs file timing, axis convention, outlier handling, and manual ROI choices.
- [ ] Export or read one representative movie per condition as OME-TIFF, including scar4-1 `.oib`, so the automated pipeline can use a common format.
- [x] Prototype max-projection preview plus PIV/optical-flow velocity fields on one WT `.ims` and one mutant `.ims`; save vector overlays for visual QC. First provisional block-correlation pass saved in `analysis/results/egg_cell_flow_prototype_2026-06-15/`; next pass needs foreground masking, drift correction, and resolved frame timing before quantitative reporting.
- [x] Implement and run masked, drift-aware network-flow batch script on all readable `.ims` files. Current output: `analysis/results/egg_cell_network_flow_all_ims_2026-06-15_v2/`.
- [x] Export scar4-1 `.oib` movies with Bio-Formats/Fiji Java and run the masked network-flow batch on the projected TIFFs. Current output: `analysis/results/egg_cell_network_flow_scar4_2026-06-16/`.
- [ ] Export representative Fiji/Bio-Formats projected TIFF/OME-TIFF stacks and rerun `analysis/analyze_egg_cell_network_flow.py` with explicit pixel size, then compare against direct `.ims` results.
- [ ] Ask collaborators for the TrackMate trajectory export, not only summary statistics: `Spots`, `Tracks`, and `Edges` tables, or a CSV with track ID, frame/time, x, y, and z if available.
- [ ] Ask collaborators for any missing metadata: pixel size, frame interval, channel meaning, objective/magnification, projection status, biological top-bottom axis orientation, and whether drift correction was applied.
- [ ] Decide which simulation object should represent XIG puncta before comparing to Cytosim: membrane/tether clusters, cargo/organelle-like puncta, or actin-associated bundle/cluster centroids.
- [ ] Prototype a Python notebook/script that reads TIFF/CZI/ND2 stacks as needed, performs drift correction, estimates velocity fields, and exports per-frame/per-region speed summaries.

## Vesicle / Manuscript Project

Status: waiting on PI feedback

Current state:

- Paper manuscript is with PI for feedback.
- Main paper scope from `vesicles/manuscript_draft/00_paper_plan_and_decisions.md`: spherical-shell study of confined actin networks with passive crosslinkers and vesicle-bound motor clusters.
- Need to be precise in methods language for metrics and parameter definitions.
- New geometry-comparison request: rerun the same paper simulation matrix in a same-volume rectangular slab with periodic boundary conditions in `x` and `y`, and in a same-volume solid sphere.
- Slab implementation choice: Cytosim `strip` space in 3D, with periodic `x/y` and hard confinement in `z`.
- Solid-sphere implementation choice: Cytosim native `sphere` space with radius chosen to match the paper spherical-shell volume.
- 2026-06-03: staged `vesicles/same_volume_geometry_comparison/` with 880 configs total: 440 solid-sphere configs and 440 periodic-slab configs.
- 2026-06-03: local direct smoke checks confirmed both geometries parse and start integrating with `/home/thekaka/project/cytosim/sim`.
- 2026-06-03: added `vesicles/same_volume_geometry_comparison/canonical_submit/` with only the canonical pair: `1:8` crosslinkers, `40` cargoes, `16` motors per vesicle, replicate `0`.
- 2026-06-03: cluster blocker was `Illegal instruction`. Since previous clean configs also failed, this pointed to binary, node, CPU-instruction, or library/module mismatch rather than the new geometry configs.
- 2026-06-06: blocker resolved operationally by running on Bergamo-compatible nodes and compiling with the cluster helper using private CMake/OpenBLAS setup and `CYTOSIM_NATIVE_ARCH=OFF`.

When feedback arrives:

- [ ] Update central claim and crosslinker wording.
- [ ] Tighten Methods definitions for all metrics.
- [ ] Revise figure captions before changing plots.
- [ ] Decide whether annulus/disc/slab geometries stay background, supplement, or are omitted.

Geometry-comparison TODOs:

- [x] Generate same-volume solid-sphere paper comparison runs.
- [x] Generate same-volume x/y-periodic slab paper comparison runs.
- [x] Run local parser/smoke checks for at least one solid-sphere and one slab config before relying on the campaign.
- [x] Create a two-config canonical submit folder for quick cluster submission.
- [x] Resolve cluster `Illegal instruction` failure using node restriction / portable compile helper.
- [ ] Submit canonical pair on a SLURM-capable cluster login node.
- [ ] Submit both staged full campaigns only after canonical pair runs cleanly.
- [ ] Extend analysis scripts to compare shell, solid sphere, and slab after outputs return.

## Daily Check-in Template

Date:

- Item worked on:
- What changed:
- Blocker:
- Next concrete step:
