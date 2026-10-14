# RpEGEN

**Current version: v2.0.0, the MATLAB GUI workflow.**

For the previous script-only workflow, download the
**[Pre-GUI Script Archive](https://github.com/burchfisher/RpEGEN/releases/tag/pre-gui-archive)**.
It preserves the April 21, 2023 repository snapshot, including the original-style
scripts, example dataset, and dependencies. This snapshot includes post-publication
fixes; earlier revisions remain in Git history. Keep old and new versions in
separate folders. See the [v2.0.0 release notes](docs/RELEASE_NOTES.md) for the update.

RpEGEN quantifies and visualizes retinal pigment epithelium (RPE) regeneration
from images and regions of interest (ROIs). The original image-analysis workflow
is described in:

Leach, L. L., Fisher, G. B., and Gross, J. M. (2022).
"Nitroreductase/metronidazole-mediated ablation and a MATLAB platform (RpEGEN)
for studying regeneration of the zebrafish retinal pigment epithelium."
*Journal of Visualized Experiments*, 181.
[doi:10.3791/63658](https://doi.org/10.3791/63658).

The JoVE article includes a detailed protocol and video demonstration.
Please cite it when using RpEGEN and identify the software version used.

You may download and modify the code under the terms of the [license](LICENSE).

## Documentation

- **[User Guide](docs/USER_GUIDE.md):** installation, example walkthrough,
  input preparation, every GUI tab, plot interpretation, exports, and troubleshooting.
- **[Changes and Migration](docs/CHANGES_AND_MIGRATION.md):** moving from the
  published workflow, retained calculations, new methods, changed defaults,
  compatibility, verification, and known limitations.

These guides describe **v2.0.0**, released on **2026-10-14**.
Existing compatible group MAT files can be opened without rerunning extraction.
The GUI uses built-in MATLAB plotting. The older `RpEGEN_PermPlot.m` script
and its GRAMM dependency are available in the pre-GUI archive, not the current
checkout. For a smaller download without Git history, use GitHub's **Code >
Download ZIP** or clone with `git clone --depth 1 https://github.com/burchfisher/RpEGEN.git`.

## MATLAB GUI

In MATLAB, set the Current Folder to this repository and run:

```matlab
RpEGEN_GUI
```

The interface originally targeted MATLAB R2024b with the Image Processing,
Curve Fitting, and Statistics and Machine Learning toolboxes. Recent
GUI verification used R2026b. Parallel Computing Toolbox is optional for faster
permutations. See the User Guide for the tested scope and dependencies.

Follow the tabs in order: [Setup](#setup) > [Quality Control](#quality-control) >
[Results](#results) > [Compare Groups](#compare-groups).
For already-analyzed data, load a group MAT in Setup to review QC and Results,
or load group MATs directly in Compare Groups to make comparisons without extraction.
The **?** button explains the active tab or comparison view.

### Setup

Select the ROI, image, and output folders, enter a group name such as
`DMSO_4dpi`, and select the TIFF page containing 8-bit grayscale brightfield
data. **Example Data** selects the supplied DMSO dataset and a separate
`GUI_Output` folder. **Check Inputs** checks names, ROI geometry, image type,
and the selected page before **Run Analysis** starts the existing pipeline.

**Example Data does not run the analysis.** Click **Run Analysis** to create the
results and populate QC. Each run analyzes one treatment group and creates one
group MAT file. Repeat for other treatments using their own inputs and group names.

To skip extraction, use **Open Existing Results** in Setup. For example, load
`JoVE_dataset/Output/DMSO_4dpi.mat`, or open it directly from MATLAB:

```matlab
RpEGEN_GUI('ResultsFile', fullfile(pwd, 'JoVE_dataset', 'Output', 'DMSO_4dpi.mat'))
```

Images must follow the published dorsal-up, distal-left orientation and pointed
ROI endpoint convention. Wait for extraction to finish before closing the GUI.
A failed run may leave partial outputs; correct the reported issue before rerunning.

### Quality Control

After analysis or loading existing results in Setup, **Quality Control** displays
the image/ROI overlay, centerline angles, median intensities, and each sample's
intensity profile. Check the ROI and centerline against the image before comparing groups.

In Quality Control, use **Left Arrow / Right Arrow** or the previous/next buttons
to move through samples. Page Up / Page Down also work. The counter
shows the current sample and total. The dropdown remains available for jumping
to a named sample. Navigation stops at the first and last sample, and keyboard
navigation is active on the Quality Control tab.

On macOS, **Show QC PDF in Finder** selects the sample's original PDF; double-click
it in Finder to open it. This avoids the unreliable MATLAB-to-Preview direct-open
handoff. Other platforms retain **Open QC PDF**. No file permissions are changed.

**Export PDFs...** exports matching QC figures for every eye in the current group
from saved results, without rerunning extraction.

### Results

Inspect the current group's pooled intensity tables at 1-, 5-, or 10-degree
resolution. **Export Binned CSVs** exports all three tables, and **Open Output
Folder** opens the folder containing the group MAT file.

The `<group>.mat` file is the full result to load for QC and comparisons.
`<group>_run.mat` records settings, the input manifest, and processing output;
it is not a replacement group result. Each eye also has a QC PDF.

Saved group MAT files also contain images, ROI polygons, and centerlines that
[EyeLayerWidth](https://github.com/burchfisher/EyeLayerWidth) can use for width
measurements. Those files let you reuse the existing geometry without redrawing ROIs.

### Compare Groups

Open **Compare Groups** and use **Add Groups...** to select multiple saved group
MAT files, **Add Current Results** for the current QC dataset, or **Example
Groups** for the four supplied groups. Loading an existing file twice does not
duplicate it. The table lets you rename plot labels, select which groups appear,
and remove selected rows from the comparison without deleting files.

**Example Groups** loads four already-analyzed groups; it does not run extraction
or populate QC and Results. If you just ran Example Data, avoid also adding that
newly generated result alongside the four supplied groups: it duplicates the
DMSO example eyes under a different file path. Wait for the loading dialog to
finish before continuing.

Use the views below in order to inspect distributions, compare median curves,
and then run a permutation comparison. Loading groups does not run the test.

#### Heatmaps

The **Heatmaps** view shows pooled pixel distributions against dorsal-to-ventral angle.
Default bins span 5 degrees and 5 grayscale intensity units. All selected groups
share fixed edges (0-180 degrees and 0-255 intensity) and the same color limits.
These explicit shared edges can differ from the original GRAMM plot's automatic
bin placement and color limits. Endpoint values of 180 and 255 enter the final
bins. The heatmap view scrolls when many groups are selected.

**Pixel count** displays raw counts. **Percent per angle** normalizes each
angular column to 100 percent, allowing comparison of intensity distributions
when groups contain different numbers of pixels. It still pools pixels and does
not give each animal equal weight. Empty angular columns remain zero.
Use **Export Figure...** to save the heatmaps as a PDF.

#### Median Profiles

The **Median Profiles** view overlays the saved pooled-pixel medians for the selected
1-, 5-, or 10-degree bins without spline smoothing. Both interval methods share
one plot with one median curve per group. Independent **Pixel-based median CI (nominal
95%)** and **Whole-eye bootstrap CI (95%)** checkboxes show either, both, or neither
band. The original saved `CI95L`/`CI95U` use darker shading; whole-eye bootstrap
intervals use lighter shading underneath. Toggling either band preserves your
zoom. **Export Figure...** respects both checkboxes and exports the full angular
profile to vector PDF without modifying source MATs. Figure exports use a white
background and preserve vector plot elements and text for publication editing;
eye images and heatmap image layers remain raster within the PDF.

Both bands ask how precisely the group median is estimated. The pixel-based CI
treats each pixel as independent and can be very narrow. The whole-eye bootstrap
asks how much the same pooled median could move with a different sample of eyes.
Neither band shows the spread of pixel values. Toggle them without rerunning
extraction; see the [User Guide](docs/USER_GUIDE.md#median-profiles) for examples
and the resampling details.

#### Permutation Test

Load two or more group MAT files (or use **Example Groups**), then select
**Permutation Test** in Compare Groups. Pick different **Group A** and **Group B**
files. These selectors include all loaded groups, independently of the Show
checkboxes used by the descriptive plots. The group table and heatmap settings
are hidden on this tab; the test has its own 1-, 5-, or 10-degree bin setting.

**Auto** defaults to a five-minute soft budget and 10,000-50,000 permutations per
bin. Its approximate 95% simulation half-width target is 0.002 near p=0.05,
corresponding to 45,618 draws when the timing estimate permits.
**Manual** uses the chosen count in every bin. **Parallel CPU (up to 4)** is
initially checked when the Parallel Computing Toolbox is licensed. **Run Test**
displays completed-bin progress and elapsed time. **Cancel** discards the partial
run and restores controls; it does not modify input MATs or shut down an existing
pool. Cancellation during pool startup may wait until pool startup finishes.

The **P-value Plot** shows permutation p-values as steps spanning each angular
bin, the 0.05 reference, and dark gray, capped pointwise 95% simulation-error
bars at bin centers. **Show simulation intervals** is off by default; toggling it
preserves zoom and does not rerun the test. Figure exports follow this checkbox;
the results table and CSV/MAT exports always retain the interval limits.
Hollow orange markers remain visible when the bars are hidden and indicate intervals that
span 0.05, not an additional significance test. **Focus near p = 0.05** changes
only the plotted y-range to 0-0.12. **Shade p < 0.05**, on by default, adds pale
green background ribbons spanning bins below the unadjusted threshold. Adjacent
bins share a ribbon; the shading does not indicate effect size or direction.
Shading can be toggled without rerunning, preserves zoom, and is included in
figure exports when enabled. P-values are not adjusted for multiple testing;
zero still means no simulated exceedances were observed. The **Results Table** lists
actual counts, exceedances, and simulation limits for every bin.

Simulation intervals describe uncertainty in the calculated p-value from using
a finite number of shuffles. They help show whether more simulations might
change which side of 0.05 a result falls on; they do not measure biological variability.

**Export Results...** writes CSV or MAT; MAT also includes group/file identities,
settings, seed, execution details, and completion time. **Export Figure...** saves
the current p-value view to vector PDF. Changes to inputs or calculation settings
clear the previous result and disable export until a new run completes. Plot-only
changes do not recalculate. Analysis and group editing are locked during testing.

For a full walkthrough, including file meanings, statistical explanations, and
troubleshooting, see the [User Guide](docs/USER_GUIDE.md).

## Technical Reference

The sections below cover calculations and programmatic use, rather than steps
needed to operate the GUI. The published intensity calculations and pixel-level
permutation statistic are retained.

### Confidence Intervals

Bootstrap intervals use 5,000 whole-eye resamples with replacement within each
group and a reproducible local random stream (seed 20261012). Each resample uses
the same eye selection at every angle, pooling only that angle's pixels. The
2.5th and 97.5th percentiles of the resampled pooled medians give pointwise intervals,
not simultaneous confidence bands and not a median of per-eye medians. Original
pixel weights are retained. These intervals assume independent eyes; shared
wells or batches are not modeled. Results are cached by group and bin width for
the current GUI session. Sparse bins with fewer than two contributing eyes or
empty bootstrap draws have no bootstrap band, with a status message. Missing or
inconsistent raw eye IDs disable that group's bootstrap band, not its pixel-based CI plot.

The original calculation uses sorted pixel intensities at rounded ranks
`n/2 +/- 1.96*sqrt(n/4)`. It is a normal approximation, not an exact 95%
binomial interval. Exact symmetric order-statistic limits use ranks `r` and
`n-r+1`, with `r` chosen from binomial tail probabilities. The current extraction
also indexes zero for bins containing fewer than six pixels, including empty
bins; it needs small-bin handling before use on sparse datasets.

The audit reproduced every saved median and both limits in all 936 bins across
the four supplied groups (1-, 5-, and 10-degree tables), plus the 234 bins in
`GUI_Output/DMSO_4dpi.mat`. Compared with exact symmetric binomial ranks,
130 of the 936 supplied intervals would change, by at most two intensity units.
This is a numerical comparison only: it does not validate the independence
assumption. Pixels within the same eye are correlated, larger pixel counts give
eyes unequal weight, and these pointwise intervals are neither between-animal
confidence intervals nor simultaneous 95% bands over the entire profile.

For treatment inference, preserve each animal's profile and account for any
shared treatment well or experimental batch. The independent treatment unit
must be established before selecting animal-level or hierarchical resampling.
The legacy extraction and saved MAT files are unchanged by this display update.
References: [NIST median confidence limits](https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/mediancl.htm)
and [Lord et al., SuperPlots](https://pmc.ncbi.nlm.nih.gov/articles/PMC7265319/).

### Batch Extraction

For batch processing, the same adapter can be called directly:

```matlab
config = struct('ROIdir', '/path/to/ROIs', ...
    'IMGdir', '/path/to/TIFs', 'Outputdir', '/path/to/output', ...
    'GroupName', 'DMSO_4dpi', 'Channel', 3, 'Overwrite', false);
[data, outputFile] = rpegen_run(config);
```

### Permutation Runtime and Precision

The interactive parallel path uses `parfeval` futures for progress and cancellation;
script calls without callbacks retain `parfor`. Both execute the same per-bin
`permsim` calculation with the same seeds, and fixed-count results must match.

`permsim` still shuffles **pixels**, uses the absolute difference in medians,
and returns the original `exceedances / iterations` p-value. Fixed-count runs
preserve the original random draw sequence and output exactly for the same RNG
state. Shuffled pixels are now processed in bounded batches rather than allocating
all permutations at once. Empty, nonfinite, and non-column inputs are rejected.

The GUI's Auto mode uses a **five-minute soft budget**
for a whole comparison and 10,000-50,000 permutations per bin. Auto times a short
calibration, discards its results, restores the RNG, and selects a fixed count
before testing. It never uses significance or observed p-value changes to decide
when to stop. Budget is divided across bins. The 10,000 minimum takes precedence,
so a slow run can exceed the estimate; this is not a hard timeout. More demanding
precision or a larger maximum can be explicitly requested through the helpers.

```matlab
% One angular bin, with a time/precision compromise:
options = struct('TimeBudgetSeconds',300/180, 'HalfWidth',0.002, ...
    'MinIterations',10000, 'MaxIterations',50000);
[p, info] = permsim(pixelsA, pixelsB, 'auto', options);

% Explicit counts remain available:
[p, info] = permsim(pixelsA, pixelsB, 5000);
info.MonteCarloCI95
info.NearThreshold
```

Near p=0.05, approximate 95% Monte Carlo half-widths are 0.0043 at 10,000
permutations and 0.0019 at 50,000. Four times as many draws only halves this
uncertainty. Repetitions improve numerical precision, not biological power or
the validity of treating pixels as exchangeable. `rpegen_permutation_plan`
calculates an uncapped, approximate fixed count for a requested half-width;
the time-budgeted auto mode can deliberately choose fewer draws.

The GUI enables up to four CPU workers through
`rpegen_permutation_profile`. It reuses an existing pool or starts a four-thread
pool and falls back to serial execution if unavailable. It does not install or
purchase a toolbox. Time allocation accounts for worker count and startup.
This wrapper uses reproducible per-bin seeds (`20261012 + bin index`), so fixed
counts give identical serial and parallel results, although they need not match
an old unseeded script run. Existing pools are not closed or resized. Direct
three-argument `permsim` calls retain their original RNG behavior.

To run a comparison from MATLAB code without opening the GUI:

```matlab
[permutation_diagnostics, runInfo] = rpegen_permutation_profile( ...
    DMSO_4dpi.RAW_data, IWR1_4dpi.RAW_data, 1, 'auto', ...
    struct('TimeBudgetSeconds',300,'UseParallel',true,'NumWorkers',4));
```

Apple GPUs are not supported for MATLAB computation on macOS; CPU parallelism
is the supported acceleration used here. See
[MathWorks Mac requirements](https://www.mathworks.com/support/requirements/matlab-mac.html).

The optional second output contains exceedance counts, actual repetitions,
resolution, and exact binomial 95% Monte Carlo limits. The code example stores these
in `permutation_diagnostics`, flagging bins whose interval spans 0.05. These
are pointwise simulation-error intervals, not biological CIs or a multiple-testing
correction. A permutation p-value of zero means no exceedance was observed, not that
the true tail probability is zero. The recommended `(b+1)/(B+1)` correction is
**not applied**, pending approval to change the published p-value calculation.
See [MathWorks binofit](https://www.mathworks.com/help/stats/binofit.html) and
[Phipson and Smyth](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf).

## Verification

### Example-Data Benchmarks

With the current 10,000-50,000 Auto limits and 0.002 precision target on MATLAB
R2026b, the full DMSO 4 dpi versus IWR1 4 dpi GUI comparison (180 one-degree
bins) finished in 149 seconds with four CPU workers. Auto selected 45,618
permutations for every bin within the five-minute budget. All GUI, numerical
parity, cancellation, and threshold-ribbon checks passed.

With the previous 5,000-20,000 Auto limits on this M4 Max, the full DMSO 4 dpi
versus IWR1 4 dpi comparison (180 one-degree bins) finished in 69 seconds
including four-thread pool startup. Auto chose
20,000 permutations in every bin within the five-minute budget. None of its
pointwise Monte Carlo intervals crossed 0.05. A separate serial 5,000-draw audit
took 86 seconds and had one interval crossing 0.05. These timings are machine-
and workload-specific, not performance guarantees or validation of pixel independence.

A controlled fixed-count benchmark measured 2.76x faster execution with a warm
four-thread pool; serial and parallel p-values were identical with the same
per-bin seeds. Batching reduced shuffled-data storage about 20-fold in the
representative 5,000-draw bin but was slower in isolation; CPU parallelism,
not batching alone, supplies the demonstrated speedup. The scripts in
`comparison_preview` and `tests/verify_rpegen_update.m` regenerate comparison
reports and figures. Generated outputs are excluded from the repository.

### Regression Checks

Run the regression checks in MATLAB with:

```matlab
addpath('tests');
test_rpegen_extraction;
test_rpegen_resampling;
test_rpegen_confidence;
test_rpegen_permutation_gui;
test_rpegen_qc_navigation;
test_rpegen_qc_and_batch_exports;
test_rpegen_pdf_exports;
test_rpegen_group_loading;
test_rpegen_help;
```

For a full example-data Auto permutation run, also run
`verify_rpegen_permutation_integration`. It can take several minutes. Tests use
temporary folders or ignored output paths and do not overwrite the bundled data.
