# Collaborator Protocol: Movie-Level F-Actin Flow Measurement On Windows

Last updated: 2026-06-16

This protocol measures F-actin movement from time-lapse microscopy movies in a
more reproducible way than manually drawing slopes on kymographs. It is written
for experimental users working on Windows with Fiji and Python.

The main result is one speed value per movie:

```text
movie median F-actin speed, in um/s
```

For the current egg-cell dataset, use:

```text
frame interval = 30 seconds
```

## Why We Are Using This

The old kymograph approach depends strongly on which contour, cable, and slope
the user selects. This workflow instead estimates local motion across the
actin-rich parts of the whole movie, then summarizes each movie with one median
speed. This makes genotype/line comparisons less sensitive to manual cable
choice.

## What You Will Install

### 1. Fiji

Install Fiji from:

```text
https://fiji.sc/
```

Fiji includes Bio-Formats, which is needed to open `.ims`, `.oib`, `.czi`,
`.nd2`, and other microscope formats.

### 2. Python

Install Python 3.10 or newer from:

```text
https://www.python.org/downloads/
```

During installation, check:

```text
Add python.exe to PATH
```

After installation, open Command Prompt and test:

```bat
python --version
```

## Files We Will Provide

Put these files in a folder called `scripts`:

```text
analyze_egg_cell_network_flow.py
summarize_actin_flow_results.py
actin_flow_requirements.txt
egg_cell_actin_flow_manifest_template.csv
BatchExportActinFlow.java
batch_export_projected_tiffs.ijm
setup_actin_flow_env.bat
run_full_actin_flow_pipeline.bat
```

Optional documentation:

```text
collaborator_actin_flow_protocol_windows.md
```

## Recommended Folder Layout

Create a working folder anywhere convenient, for example:

```text
C:\ActinFlowAnalysis
```

Inside it, make these folders:

```text
C:\ActinFlowAnalysis
  raw_movies
  projected_tiffs
  manifests
  results
  scripts
```

Use the folders like this:

- `raw_movies`: original microscope files from the microscope or Imaris.
- `projected_tiffs`: standardized Fiji exports used for analysis.
- `manifests`: CSV files that describe the movies.
- `results`: output from the Python analysis.
- `scripts`: Python scripts and requirements file.

## One-Time Python Setup

Open Command Prompt.

Move into the working folder:

```bat
cd C:\ActinFlowAnalysis
```

Create a Python environment:

```bat
python -m venv actin_flow_env
```

Activate it:

```bat
actin_flow_env\Scripts\activate
```

Install packages:

```bat
python -m pip install --upgrade pip
python -m pip install -r scripts\actin_flow_requirements.txt
```

Test that the analysis script opens:

```bat
python scripts\analyze_egg_cell_network_flow.py --help
```

You only need to do this setup once on a computer.

You can also run the setup helper:

```bat
scripts\setup_actin_flow_env.bat
```

## Automated Path: Export, Manifest, Analysis

Put raw movies into line/genotype folders:

```text
C:\ActinFlowAnalysis\raw_movies\Wild type
C:\ActinFlowAnalysis\raw_movies\CA ROP9
C:\ActinFlowAnalysis\raw_movies\DN ROP9
C:\ActinFlowAnalysis\raw_movies\scar4-1
```

Then run:

```bat
scripts\run_full_actin_flow_pipeline.bat
```

This runs Fiji/Bio-Formats from the command line, exports max-projected TIFF
time stacks, writes the manifest automatically, runs the Python flow analysis,
and creates the line summary plot.

Defaults:

```text
frame interval = 30 s
z projection = slices 5-7, max intensity
Bio-Formats series = 1
actin channel = 1
```

Override defaults before running if needed:

```bat
set FRAME_INTERVAL_S=30
set Z_START=5
set Z_STOP=7
set BIOFORMATS_SERIES=1
set ACTIN_CHANNEL=1
set FIJI_EXE=C:\Users\YourName\Downloads\Fiji\fiji-windows-x64.exe
scripts\run_full_actin_flow_pipeline.bat
```

Automated outputs:

```text
projected_tiffs\<line>\*_maxproj.tif
manifests\egg_cell_actin_flow_manifest.csv
manifests\fiji_export_metadata.csv
manifests\fiji_export_log.txt
results\egg_cell_actin_flow\movie_summary.csv
results\egg_cell_actin_flow\qc
results\egg_cell_actin_flow\summary_by_line
```

## Manual Fallback: Export Standardized Movies From Fiji

Do this for each raw movie.

### Step 1: Open The Movie

Open Fiji.

Use:

```text
File > Import > Bio-Formats
```

Select the raw movie file.

If Bio-Formats shows multiple series or resolution levels, choose the highest
resolution real movie series. For Imaris `.ims` files, this is usually:

```text
Series 1
```

Do not choose the smaller downsampled series.

### Step 2: Record Image Properties

In Fiji, go to:

```text
Image > Properties
```

Write down:

- number of channels
- number of z slices
- number of time frames
- pixel width in microns
- pixel height in microns
- voxel depth
- frame interval

For the current egg-cell dataset, use `30 seconds` even if Fiji or Imaris shows
`0 sec` or another incorrect value.

### Step 3: Make A Z Projection

Use:

```text
Image > Stacks > Z Project...
```

Settings:

```text
Start slice: 5
Stop slice: 7
Projection type: Max Intensity
```

Make sure all time frames are included.

If the movie has fewer than 7 z slices, use the last three z slices instead and
write down the actual start/stop slices.

### Step 4: Save The Projected Time Stack

Save the projected movie as a TIFF:

```text
File > Save As > Tiff...
```

Save it under:

```text
C:\ActinFlowAnalysis\projected_tiffs\<line>
```

Example filenames:

```text
WT_1_0001_z5-7_maxproj.tif
CA_ROP9_Image0004_z5-7_maxproj.tif
scar4_Image0001_z5-7_maxproj.tif
```

Important: do not apply brightness/contrast, blur, background subtraction, or
manual enhancement to the TIFF used for analysis. Those are fine for display
figures, but not for the quantitative input.

## Manual Fallback: Fill In The Manifest CSV

The automated Fiji export writes this file for you:

```text
manifests\egg_cell_actin_flow_manifest.csv
```

Only fill it manually if the automated Fiji export cannot open a file or if you
choose to export the projected TIFFs by hand.

Copy this file:

```text
scripts\egg_cell_actin_flow_manifest_template.csv
```

Paste it into:

```text
C:\ActinFlowAnalysis\manifests
```

Rename it, for example:

```text
egg_cell_actin_flow_manifest.csv
```

Each row is one movie.

Required columns:

```text
line,movie_id,projected_tiff,pixel_size_um,frame_interval_s,z_start,z_stop,raw_file,frame_interval_source,notes
```

Example row:

```csv
Wild type,WT_1_0001,projected_tiffs/Wild_type/WT_1_0001_z5-7_maxproj.tif,0.10358,30,5,7,raw_movies/Wild_type/1_0001_5ADVMLE.ims,collaborator confirmed,
```

Column meanings:

- `line`: genotype or line, such as `Wild type`, `CA ROP9`, `DN ROP9`, `scar4-1`.
- `movie_id`: short unique name for the movie.
- `projected_tiff`: path to the projected TIFF stack.
- `pixel_size_um`: pixel size from Fiji, in microns per pixel.
- `frame_interval_s`: seconds per frame. Current egg-cell data: `30`.
- `z_start`: first z slice used in the projection.
- `z_stop`: last z slice used in the projection.
- `raw_file`: original raw movie path.
- `frame_interval_source`: where the time interval came from.
- `notes`: optional.

Important: do not make one row per cable. The row is for the whole movie.

## Part 3: Run The Movie Flow Analysis

Open Command Prompt.

Move into the working folder:

```bat
cd C:\ActinFlowAnalysis
```

Activate the Python environment:

```bat
actin_flow_env\Scripts\activate
```

Run the analysis:

```bat
python scripts\analyze_egg_cell_network_flow.py --manifest manifests\egg_cell_actin_flow_manifest.csv --output results\egg_cell_actin_flow
```

The output folder will contain:

```text
results\egg_cell_actin_flow\movie_summary.csv
results\egg_cell_actin_flow\vectors
results\egg_cell_actin_flow\drift
results\egg_cell_actin_flow\qc
results\egg_cell_actin_flow\README.md
```

## Part 4: Check The QC Images

Before trusting the numbers, open the QC images in:

```text
results\egg_cell_actin_flow\qc
```

Check these files for each movie:

```text
*_preview.png
*_mask.png
*_flow.png
*_drift.png
```

What to look for:

- `preview`: projected movie is visible and not blank.
- `mask`: cyan outline covers the cell/actin-rich region, not mostly background.
- `flow`: arrows are mostly on actin structures.
- `drift`: whole-movie drift is not extreme or erratic.

If a movie fails QC, do not include it in genotype comparison until the import,
projection, or metadata are corrected.

## Part 5: Summarize By Line

Run:

```bat
python scripts\summarize_actin_flow_results.py results\egg_cell_actin_flow\movie_summary.csv --output results\egg_cell_actin_flow\summary_by_line
```

This writes:

```text
results\egg_cell_actin_flow\summary_by_line\line_summary.csv
results\egg_cell_actin_flow\summary_by_line\pairwise_permutation_tests.csv
results\egg_cell_actin_flow\summary_by_line\line_summary_plot.png
```

The plot uses one point per movie.

## Main Output Columns

In `movie_summary.csv`, the main value is:

```text
median_speed_um_per_s_assumed_dt
```

This is the median local F-actin speed for that movie, converted to `um/s`
using the frame interval in the manifest.

Other useful columns:

- `vectors`: number of local vector measurements accepted.
- `mask_fraction`: fraction of image included in the foreground mask.
- `median_speed_um_per_frame`: speed before time conversion.
- `median_peak_z`: rough confidence/quality of local correlations.

## What The Code Does

The script:

1. Reads each projected TIFF movie.
2. Normalizes image intensity frame by frame.
3. Builds a foreground mask to focus on actin-rich/cell regions.
4. Estimates any whole-movie drift.
5. Splits the image into small windows.
6. For each window, compares frame `t` to frame `t+1` using local
   cross-correlation.
7. Converts local pixel shifts into microns using `pixel_size_um`.
8. Converts per-frame displacement into speed using `frame_interval_s`.
9. Summarizes the movie by the median speed across accepted local vectors.
10. Saves QC overlays so the measurement can be visually checked.

## Reporting Rules

Use this wording:

```text
F-actin movement was quantified from max-projected time-lapse movies using
local cross-correlation on actin-rich regions. Each movie was summarized by the
median local actin speed, and movies were treated as the replicate unit.
```

Always report:

- projection method and z slices;
- frame interval;
- pixel size;
- number of movies per line;
- movie-level mean and SD;
- statistical test used.

Do not report:

- every local vector as an independent biological replicate;
- every cable/kymograph slope as an independent biological replicate;
- results from movies that failed QC.

## Current Egg-Cell Reference Values

Using this protocol with `30 s/frame`, current movie-level means are:

```text
Wild type: 0.01357 +/- 0.00287 um/s, n = 7 movies
CA ROP9:   0.01318 +/- 0.00405 um/s, n = 10 movies
DN ROP9:   0.01975 +/- 0.00548 um/s, n = 5 movies
scar4-1:   0.01007 +/- 0.00295 um/s, n = 5 movies
```

These values are means and SDs of movie median speeds.

## Troubleshooting

### Python is not recognized

Reinstall Python and check:

```text
Add python.exe to PATH
```

### Package installation fails

Try:

```bat
python -m pip install --upgrade pip
python -m pip install numpy scipy matplotlib pillow
```

### The script cannot find a TIFF file

Check the `projected_tiff` path in the manifest. Paths should be relative to:

```text
C:\ActinFlowAnalysis
```

or absolute Windows paths.

### The output speed is blank

Check that `frame_interval_s` is filled in the manifest.

### The QC mask is wrong

Possible causes:

- wrong Bio-Formats series selected;
- wrong z slices projected;
- movie is too dim or saturated;
- projected TIFF is not a time stack.

Fix the Fiji export first, then rerun the Python script.

### Fiji shows frame interval as 0 sec

Use the known acquisition interval from the microscope settings or collaborator
notes. For the current egg-cell dataset, use:

```text
30
```

in the `frame_interval_s` column.
