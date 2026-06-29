# Standardized F-Actin Flow Measurement Protocol

Last updated: 2026-06-16

This protocol measures movie-level F-actin network motion from time-lapse movies.
It is intended to replace hand-picked kymograph cable slopes for genotype/line
comparisons.

The main output is one value per movie: the median local actin speed in `um/s`.
For the current egg-cell dataset, use `30 s/frame`.

## What This Measures

- Local frame-to-frame displacement of actin-rich image regions.
- Movie-level summary: median speed across all accepted local vectors.
- Replicate unit: one movie/sample, not one cable segment and not one vector.

This does not yet measure signed top-to-bottom transport. It measures speed
magnitude. Signed axial velocity should be added only after the biological axis
is defined for every movie.

## Required Inputs

- Raw microscopy movie files, opened through Fiji/Bio-Formats.
- A projected time stack saved as TIFF or OME-TIFF for each movie.
- A manifest CSV with one row per movie.
- Confirmed frame interval. For the current egg-cell dataset: `30 s/frame`.

## Software

- Fiji with Bio-Formats.
- Python image-analysis environment on WSL:

```bash
cd /home/thekaka/project/cytosim
/tmp/cytosim_img_env/bin/python analysis/analyze_egg_cell_network_flow.py --help
```

## Folder Setup

Use this structure for a new analysis batch:

```text
analysis/input/projected_tiffs/<line>/
analysis/manifests/
analysis/results/<analysis_name>/
```

Example analysis name:

```text
central_cell_actin_flow_2026-06-16
```

## Part 1: Export Standardized Projected Movies In Fiji

Do this once per raw movie.

1. Open Fiji.

2. Open the raw movie:
   - `File > Import > Bio-Formats`
   - Select the raw `.ims`, `.oib`, or other microscope file.
   - If Bio-Formats shows multiple series/resolution levels, choose the
     highest-resolution real movie series. For Imaris pyramids, this is usually
     `Series 1`, not the lower-resolution duplicate.

3. Confirm the imported dimensions:
   - `Image > Properties`
   - Record:
     - channels
     - z slices
     - time frames
     - pixel width in `um`
     - pixel height in `um`
     - voxel depth
     - frame interval

4. If Fiji reports `Frame interval = 0 sec` but acquisition is known, use the
   confirmed acquisition interval in the manifest. For the current egg-cell
   dataset, use `30`.

5. Make the analysis projection:
   - `Image > Stacks > Z Project...`
   - Start slice: `5`
   - Stop slice: `7`
   - Projection type: `Max Intensity`
   - Use all time frames.

6. If the movie has fewer than 7 z slices, use the last three available z
   slices and record the actual `z_start` and `z_stop` in the manifest.

7. Save the projected time stack:
   - `File > Save As > Tiff...`
   - Recommended naming:

```text
<line>_<movie_id>_z<start>-<stop>_maxproj.tif
```

Example:

```text
WT_1_0001_z5-7_maxproj.tif
scar4_Image0001_z5-7_maxproj.tif
```

8. Do not apply brightness/contrast, blur, rolling-ball subtraction, or manual
   contrast enhancement to the analysis copy. Those display steps can be used
   for visualization, but not for the quantitative input stack.

## Part 2: Fill The Manifest

Start from:

```text
analysis/egg_cell_actin_flow_manifest_template.csv
```

Required columns:

```text
line,movie_id,projected_tiff,pixel_size_um,frame_interval_s,z_start,z_stop,raw_file,frame_interval_source,notes
```

Column meanings:

- `line`: genotype/line label, e.g. `Wild type`, `CA ROP9`, `DN ROP9`, `scar4-1`.
- `movie_id`: short unique ID for the movie.
- `projected_tiff`: path to the projected TIFF stack.
- `pixel_size_um`: Fiji/Bio-Formats pixel width in microns. Use this only when
  pixel width and height are effectively equal, as in the current data.
- `frame_interval_s`: confirmed seconds per frame. Current egg-cell data: `30`.
- `z_start`, `z_stop`: z slices used for max projection.
- `raw_file`: original microscope file path.
- `frame_interval_source`: e.g. `Bio-Formats`, `collaborator confirmed`.
- `notes`: optional.

Important: each movie gets one row. Do not add one row per cable or one row per
manual trace.

## Part 3: Run The Flow Analysis

From PowerShell:

```powershell
wsl bash -lc 'cd /home/thekaka/project/cytosim && /tmp/cytosim_img_env/bin/python analysis/analyze_egg_cell_network_flow.py --manifest analysis/manifests/my_actin_flow_manifest.csv --output analysis/results/my_actin_flow_run'
```

From a WSL shell:

```bash
cd /home/thekaka/project/cytosim
/tmp/cytosim_img_env/bin/python analysis/analyze_egg_cell_network_flow.py \
  --manifest analysis/manifests/my_actin_flow_manifest.csv \
  --output analysis/results/my_actin_flow_run
```

The script will write:

```text
analysis/results/my_actin_flow_run/movie_summary.csv
analysis/results/my_actin_flow_run/vectors/
analysis/results/my_actin_flow_run/drift/
analysis/results/my_actin_flow_run/qc/
analysis/results/my_actin_flow_run/README.md
```

## Part 4: QC Each Movie

Open the QC images before trusting the numbers.

Check these files:

```text
qc/*_preview.png
qc/*_mask.png
qc/*_flow.png
qc/*_drift.png
```

Pass criteria:

- The mask covers the cell/actin-rich region and excludes most empty background.
- Flow vectors sit mostly on visible actin structures.
- The vector field is not dominated by one obvious imaging artifact.
- Drift is small or smoothly varying.
- The projected movie is not saturated, blank, badly cropped, or mis-imported.

Flag a movie if:

- the wrong series/resolution was exported;
- the z projection missed the visible actin signal;
- the mask mostly covers background;
- vectors are mostly on noise;
- the movie has too few time frames;
- the biological axis or frame interval is uncertain.

## Part 5: Summarize By Movie, Then By Line

Use `movie_summary.csv` as the main analysis table.

Primary movie metric:

```text
median_speed_um_per_s_assumed_dt
```

For the current egg-cell dataset, this is the median local actin speed using
`30 s/frame`.

For line/genotype summaries:

- average the movie-level medians;
- report SD or SEM across movies;
- show each movie as a point;
- use movies/samples as the replicate unit.

Do not use individual local vectors as biological replicates.

## Part 6: Statistical Comparison

Recommended first-pass comparison:

- movie-level permutation test on mean movie median speed;
- optionally add a non-parametric rank test as a sensitivity check;
- avoid ANOVA on cable segments unless there is a mixed model that nests cables
  inside movies/samples.

Report clearly:

```text
Metric: movie median network-scale F-actin speed
Frame interval: 30 s
Replicate unit: movie/sample
Projection: max projection, z slices 5-7
```

## Part 7: Demonstration Script For Collaborators

Use this order when showing the protocol:

1. Show why the old kymograph method is trace-sensitive.
2. Open one raw movie in Fiji/Bio-Formats.
3. Choose the highest-resolution series.
4. Show `Image > Properties` and record pixel size/frame count.
5. Z-project slices 5-7 with max intensity.
6. Save the projected TIFF.
7. Add one manifest row.
8. Run the Python command.
9. Open the QC mask and vector overlay.
10. Open `movie_summary.csv` and point to the movie-level speed.
11. Show the final genotype/line plot with one point per movie.

## Current Egg-Cell Reference Values

Using this network-flow approach and `30 s/frame`, the current egg-cell movie
summaries are:

```text
Wild type: 0.01357 +/- 0.00287 um/s, n = 7 movies
CA ROP9:   0.01318 +/- 0.00405 um/s, n = 10 movies
DN ROP9:   0.01975 +/- 0.00548 um/s, n = 5 movies
scar4-1:   0.01007 +/- 0.00295 um/s, n = 5 movies
```

These are movie-level means of movie median speeds.

## Central-Cell Extension

For central-cell movies, use the same export and flow protocol first. Then add
a second analysis branch for signed axial motion:

- define the top-bottom axis for each movie;
- project each local vector onto that axis;
- report signed axial velocity in addition to speed magnitude;
- compare this directly with simulation `v_z` only after the axis convention is
  documented.

Until then, compare simulation and experiment cautiously:

- experiment: projected movie-level speed magnitude;
- simulation plot: mean axial velocity.

They can be compared for order-of-magnitude sanity, but not as identical
metrics.
