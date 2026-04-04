# Cytosim Analysis Package

This folder contains the reusable analysis code plus the curated figures and summary tables for the fixed-global versus rotatable transport story.

## Code layout
- `metrics.py`, `engine.py`, `reader.py`: run-level metric extraction helpers
- `run_speckle_kymograph_batch.py`, `run_speckle_timecourse_batch.py`, `transport_regime_batch.py`: batch analysis entry points
- `plot_velocity_with_control.py`, `plot_transport_presentation.py`, `plot_kymographs_with_control.py`, `plot_vz_presentation.py`: figure generation scripts used for the current story

## Results layout
- `results/fixed_global_story/presentation_plots`: slide-ready summary plots
- `results/fixed_global_story/transport_regime`: run/condition tables plus transport-regime figures
- `results/fixed_global_story/kymographs`: cleaned kymograph panels and scale assets

Temporary preview/test figures and helper scripts were removed during cleanup so the tracked analysis package only contains the current useful assets.
