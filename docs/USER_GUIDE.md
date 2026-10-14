# RpEGEN User Guide

This guide covers the MATLAB GUI for image analysis, quality control, group
intensity comparisons in **v2.0.0 (2026-10-14)**. Width analysis is provided
separately by [EyeLayerWidth](https://github.com/burchfisher/EyeLayerWidth).

The previous script-only workflow is available in the
[Pre-GUI Script Archive](https://github.com/burchfisher/RpEGEN/releases/tag/pre-gui-archive).
Download old and new versions into separate folders rather than mixing files.

New users can start with the example walkthrough below. Users moving from the
published scripts should also read [Changes and Migration](CHANGES_AND_MIGRATION.md).
The [Leach et al. JoVE protocol](https://pmc.ncbi.nlm.nih.gov/articles/PMC9036407/)
remains the reference for specimen preparation, imaging, orientation, and ROI
construction. This guide does not replace that laboratory protocol.

Click **?** at the top right of the GUI for help about the current tab, including
plot meanings and table-column definitions. In Compare Groups, help
follows the selected subtab. Hover over controls for brief definitions. After a
code update, close and relaunch `RpEGEN_GUI` to load the new interface.

## Contents

- [Installation and launch](#installation-and-launch)
- [Try the example data](#try-the-example-data)
- [Prepare your own inputs](#prepare-your-own-inputs)
- [Run an analysis](#run-an-analysis)
- [Review quality control](#review-quality-control)
- [Inspect and export results](#inspect-and-export-results)
- [Load groups and make heatmaps](#load-groups-and-make-heatmaps)
- [Read median profiles and confidence intervals](#read-median-profiles-and-confidence-intervals)
- [Run a permutation comparison](#run-a-permutation-comparison)
- [Width analysis has moved](#width-analysis-has-moved)
- [Save a reproducible analysis](#save-a-reproducible-analysis)
- [Troubleshooting](#troubleshooting)

## Installation and Launch

Keep the project folders together, including `ReadImageJROI-master` and
`JoVE_dataset`. Use MATLAB itself, not the macOS Terminal, for the commands below.
Set MATLAB's Current Folder to the project folder, then run:

```matlab
RpEGEN_GUI
```

For a project installed on your Desktop as `RpEGEN-main`, a Mac example is:

```matlab
cd('~/Desktop/RpEGEN-main')
RpEGEN_GUI
```

Use your actual folder location if different. On Windows, select the folder in
MATLAB's Current Folder browser before running the launch command.

| Component | Requirement |
| --- | --- |
| MATLAB | GUI originally targeted R2024b; recent verification used R2026b. Earlier releases have not been established as compatible. |
| Image Processing Toolbox | Required for the extraction and geometry workflow. |
| Curve Fitting Toolbox | Required by the extraction runner. |
| Statistics and Machine Learning Toolbox | Required by the extraction runner and statistical workflow. |
| ReadImageJROI | Included; the GUI adds its folder to the path. |
| Parallel Computing Toolbox | Optional for CPU acceleration of permutations; serial execution is available. |
| GRAMM | Not required by the GUI. Still required by the plotting sections of `RpEGEN_PermPlot.m`. |

An upgrade is not required solely to launch the original R2024b-targeted GUI.
The newest GUI features have been verified on R2026b; compatibility with every
release and operating system has not been tested. Use `ver` in MATLAB to list
installed products. A missing license is different from a missing installation.

After installing code updates, finish running calculations, close the RpEGEN
window, and launch it again. An already-open GUI can retain old callbacks.

## Try the Example Data

Start with the full workflow below to learn how to analyze images. If you only
want to explore the plots, the shortcut afterward uses results that have already
been analyzed and does not require **Run Analysis**.

### Run the example analysis from start to finish

1. Launch `RpEGEN_GUI`.
2. In **Setup**, click **Example Data**. This fills in the DMSO 4 dpi image and
   ROI folders, group name `DMSO_4dpi`, brightfield page 3, and output folder
   `GUI_Output`. It also checks the inputs, but **does not run the analysis**.
3. Confirm that the status reports **11 matched samples**. If you change any
   input settings, click **Check Inputs** again before proceeding.
4. Click **Run Analysis** in Setup and wait for it to finish. If results already
   exist, choose **Replace** only if you intend to overwrite that group's results;
   otherwise cancel and choose a different output folder before running again.
5. When analysis completes, the GUI opens **Quality Control**. Browse the eyes
   using the arrows and inspect the ROI, centerline, and intensity profile.
6. Open **Results** to inspect the current group's binned intensity tables.
   This run has produced **one group file**, `GUI_Output/DMSO_4dpi.mat`, not all
   four example groups.
7. In **Compare Groups**, click **Example Groups** to load the four supplied,
   already-analyzed group files from `JoVE_dataset/Output`. For this walkthrough,
   do not also click Add Current Results: that would add your newly generated
   DMSO result as a fifth entry containing the same example eyes.
8. View **Heatmaps**, then **Median Profiles**. Toggle the two CI bands separately.
9. Open **Permutation Test**, choose `DMSO_4dpi` and `IWR1_4dpi`, leave **Auto**
   and the five-minute budget selected, and click **Run Test** when ready.

**Example Data** prepares one group's inputs for extraction. **Example Groups**
loads four existing result files for comparisons. They perform different tasks;
Example Groups does not analyze the other three groups from their images.

The supplied groups contain 11 DMSO 4 dpi, 12 IWR1 4 dpi, 11 DMSO 9 dpf, and
10 IWR1 9 dpf eyes. **Example Groups** adds to the current list; it does not
clear other groups. Loading the same file path twice is prevented, but a copy
of the same dataset at another path can appear as another group.

### Shortcut using already analyzed results

To skip extraction, launch the GUI and select **Setup > Open Existing Results**,
then choose `JoVE_dataset/Output/DMSO_4dpi.mat`. **Do not click Run Analysis for
this route:** the file already contains the analysis. Continue at step 5 above
to review QC, inspect Results, and explore intensity comparisons.

Alternatively, open that saved example directly from MATLAB:

```matlab
RpEGEN_GUI('ResultsFile', ...
    fullfile(pwd,'JoVE_dataset','Output','DMSO_4dpi.mat'))
```

## Prepare Your Own Inputs

Process one treatment/time-point group at a time. Keep its TIFF images and
matching FIJI ROIs in separate folders. One example layout is:

```text
study/
  control/
    ROIs/eye01.roi
    TIFs/eye01.tif
  treated/
    ROIs/eye02.roi
    TIFs/eye02.tif
  outputs/
    control/
    treated/
```

Use one independent eye per animal for the independent-eye bootstrap described
here. The software counts image IDs, not verified independent biological units.
Multiple sections from one animal, paired eyes, or shared treatment wells/batches
require an analysis design that accounts for that dependence.

Before extraction:

- Orient each eye **dorsal up and distal/lens side left**; preserve the pointed
  dorsal and ventral ROI ends specified in the published protocol.
- Use matching `.tif` and `.roi` basenames, including case. The internal ROI
  name saved by FIJI must also match the image basename.
- Use individual ROI files, not an ROI ZIP archive. The validator accepts
  polygon, freehand, and traced area ROIs with at least three vertices.
- Select an **8-bit grayscale** brightfield TIFF page. A page is a position in
  the TIFF stack, not an RGB color channel. Page 3 is the example setting, not
  a universal setting.
- Use lowercase `.tif` and `.roi` extensions for the current file search. The
  folders should contain only the matched samples intended for that group.
- Choose a valid MATLAB identifier for **Group name**, such as `Control_4dpi`:
  start with a letter, use letters/digits/underscores, and omit spaces or hyphens.
- Use separate output folders when groups share sample basenames. QC PDF names
  contain the sample name, not the group name, so shared names can overwrite PDFs.

Apply consistent imaging and ROI criteria across groups. Check Inputs cannot
decide whether an ROI correctly represents the tissue or whether image exposure
is comparable across samples.

## Run an Analysis

1. Set **ROI folder**, **Image folder**, **Output folder**, and **Group name**.
2. Set **Brightfield page** to the correct page number.
3. Click **Check Inputs** and inspect the sample count and dimensions.
4. Click **Run Analysis**. Wait for completion before closing the interface.
5. Review the results in **Quality Control** before comparing groups.
6. Repeat for each additional group with its own name and inputs.

**Check Inputs** checks folder existence, matching filenames and ordering,
internal ROI names, selected TIFF page, image type, ROI type, nonzero polygon
area, finite coordinates, and bounds. It does not run extraction, approve the
anatomy, detect every problematic polygon, or verify output-folder writability.
The runner checks the output folder when analysis starts.

Extraction is synchronous and has no GUI Cancel button. If it fails, inspect
the message and MATLAB Command Window, correct the cause, and rerun. A failed
run can leave partial files; a MAT file alone does not prove the entire run and
all QC exports completed.

## Review Quality Control

The QC view shows the image/ROI/centerline, centerline angular coordinates,
centerline-associated median intensities, and the selected eye's intensity
profile. Angles run from 0 degrees dorsally to 180 degrees ventrally; intensity
runs from black at 0 to white at 255.

Use **Left Arrow / Right Arrow**, the previous/next buttons, or the sample dropdown.
Page Up / Page Down remain available. Arrow keys retain normal editing and selection
behavior in input controls; click an image or navigation button to resume eye navigation. Navigation stops
at the first and last eye. On macOS, **Show QC PDF in Finder** selects the
sample's saved PDF in Finder; double-click it there to open it. This avoids a
MATLAB-to-Preview handoff that can report a permission error even when the file
opens normally from Finder. On other platforms, **Open QC PDF** opens the PDF
directly. Initially the PDF is located beside the loaded group MAT.

**Export PDFs...** regenerates PDFs for every eye in this QC group from saved
results, without rerunning extraction or modifying the MAT. Select a folder;
replacing existing PDFs requires confirmation. The open/reveal button then uses
that export folder for this session. Newly generated PDFs match the GUI: three
equally sized brightfield overlays with color scales and the intensity profile
below. Colored dots are retained. Existing PDFs do not change until regenerated.

Check that the ROI follows the intended tissue, the centerline follows the RPE
without erroneous branches or jumps, the endpoints are sensible, and intensity
changes correspond to the image. The individual QC profile is not the pooled
group curve shown in Compare Groups.

There is no ROI editor, pass/fail annotation system, or per-eye exclusion button
in this GUI. Correct ROIs in FIJI and rerun the affected group. Record exclusions
outside the app using predefined criteria; hiding a group in a plot does not
exclude an eye from its saved data.

## Inspect and Export Results

**Results** displays the current group's saved 1-, 5-, or 10-degree intensity
tables. **Export Binned CSVs** writes all three tables, not just the selected one.
**Open Output Folder** opens the folder containing the loaded group MAT.

| File | Contents and use |
| --- | --- |
| `<group>.mat` | Full group result: images, geometry, raw pixels, centerline measurements, and binned intensity tables. Load this for QC and comparisons. |
| `<group>_run.mat` | Settings, input manifest, and processing log from the GUI/batch runner. This is not a group result file. |
| `<sample>_QC_plots.pdf` | Per-sample QC figure. |
| `<group>_BIN_1_deg_all.csv` | Pooled 1-degree intensity summaries. Equivalent files are exported for 5 and 10 degrees. |

In binned tables, `Numel` is a **pixel count**, not the animal count; `STD`
describes pixel variability, `SE` is the saved standard error of the mean, and
`CI95L`/`CI95U` are the saved pixel-based median limits. None of these columns
is the new whole-eye bootstrap interval.

**One analysis run produces one group MAT.** Run each treatment group separately,
then load those group MATs together in Compare Groups. No extraction rerun is
needed just to make plots, compute bootstrap bands, or run permutations.
Compatible full results can also be opened separately in EyeLayerWidth.

## Load Groups and Make Heatmaps

In **Compare Groups**, use **Add Groups...** to select multiple group MATs or
**Add Current Results** to add the group loaded in Setup. Do not select
`_run.mat` files or exported permutation/width MATs.

A **Loading Groups** dialog appears while files, plots, and bootstrap confidence
bands are prepared. Loading buttons are disabled until it closes. Click once;
files already in the comparison are skipped without rebuilding the plots.

The group table lets you change plot labels and toggle **Show**. **Remove
Selected** removes table entries from the comparison, not files from disk.
The names and visibility choices are session settings, not edits to the source
MAT. Loading files does not merge their animals into a new combined group.

**Heatmaps** is available immediately; no permutation test is required.

- **Angle bin**: 1, 5, or 10 degrees; default 5.
- **Intensity bin**: 1, 5, or 15 intensity units; default 5.
- **Pixel count**: total pixels in each angle/intensity cell.
- **Percent per angle**: each nonempty angular column totals 100 percent.

All shown groups share histogram edges and color limits. Changing the displayed
groups can therefore change the common color scale. Percent normalization
facilitates distribution comparisons but still pools pixels; it does not give
each eye equal weight. A zero/empty column is not evidence of zero intensity.
Use **Export Figure...** for a vector PDF with editable plot elements and text.
All publication figure exports use PDF with a white background. Source eye images
and heatmap image layers remain raster within the PDF rather than becoming
individual editable vector objects. Keep the source MAT files for
replotting data. Editing support depends on the PDF editor and available fonts.

## Read Median Profiles and Confidence Intervals

**Median Profiles** shows one pooled-pixel median curve per group, using the
saved bins without the original plotting script's spline interpolation. The
**Angle bin** control sets the displayed bin width. Intensity-bin and color
normalization settings affect heatmaps, not the median calculation.

Curve colors follow the order in which groups are loaded.

| Checkbox | What the band represents |
| --- | --- |
| Pixel-based median CI (nominal 95%) | Darker shading from the original saved median limits. Treats pixels as independent for the calculation; often very narrow. |
| Whole-eye bootstrap CI (95%) | Lighter shading from resampling independent eyes, then recalculating the pooled median at each angle. |

Both describe uncertainty in a median, not the range containing 95 percent of
pixel values. Both are **pointwise** intervals, not a simultaneous 95 percent
band for the entire eye. Overlap or non-overlap of the bands is not the
permutation test.

The bootstrap uses 5,000 resamples with replacement and seed `20261012`.
The original interval effectively asks how precisely the median is estimated
when each pixel is independent. Many pixels can make it very narrow, even with
few animals. The whole-eye bootstrap asks how much that same pooled median could
move with a different sample of independent eyes. It uses between-eye differences
to estimate uncertainty in the group median, rather than simply displaying the
spread of individual eyes. Pixels from the same eye may be related, which the
pixel-based interval does not account for. Bootstrap bands are often wider, but
not necessarily in every bin; neither method adds biological replication.

Each resample selects whole eyes, using the same selected eyes at every angle.
For a given angle bin, only pixels originally in that bin are pooled. Pixels
are never moved between angles and no composite eyes are assembled. Selecting
an eye twice repeats its contribution twice. This is **not a median of per-eye
medians**: eyes with more pixels in a bin still contribute more pixels.

Toggle either band without changing the median or rerunning extraction.
Bootstrap results are cached for the GUI session; their first calculation can
take longer. Missing/inconsistent eye IDs, fewer than two contributing eyes in
a bin, or empty bootstrap draws can prevent a band from being displayed. Read
the status line; a missing band is not a zero-width interval.

**Export Figure...** respects the CI checkboxes and exports the full angular
profile. The GUI does not currently export bootstrap numbers through the Results
CSV button. Programmatic export is available from the helper:

```matlab
loaded = load(fullfile('JoVE_dataset','Output','DMSO_4dpi.mat'));
group = loaded.DMSO_4dpi;
bootstrap = rpegen_bootstrap_ci(group.RAW_data,group.IMGname,5);
save('DMSO_4dpi_bootstrap_5deg.mat','bootstrap')
```

## Run a Permutation Comparison

1. Load at least two distinct groups in Compare Groups.
2. Open **Permutation Test** and select **Group A** and **Group B**.
3. Choose its **Angle bin** independently of the descriptive plots; default 1 degree.
4. Use **Auto** with **Budget (min)** of 5 for the standard runtime compromise.
5. Leave **Parallel CPU (up to 4)** enabled if a licensed toolbox is available.
6. Click **Run Test**, review both the plot and **Results Table**, and export.

The two selectors include all loaded groups, regardless of their Show checkboxes.
Group-table and heatmap controls are hidden while this subtab is selected.

### What is tested

Within each angle bin, the test pools the two groups' pixels and shuffles them
between groups while preserving their pixel counts. It compares the absolute
difference in medians, counting shuffled differences at least as large as the
observed one. The p-value is `Exceedances / Iterations`.

This remains the published **pixel-level** test, not a shuffle of whole animals
or a paired test. Within-eye correlation can make pixel-level evidence look more
precise than independent-animal evidence. The new bootstrap display does not
change or validate that permutation assumption. P-values are unadjusted across
angle bins and across separate comparisons; no multiple-testing correction is
applied.

### Auto and Manual

Auto targets an approximate 95 percent simulation half-width of 0.002 near
p = 0.05, subject to **10,000 to 50,000 draws per bin** and a soft runtime budget.
The target corresponds to 45,618 draws when time permits. A short timing
calibration selects each bin's count before its test, without using the observed
p-value to decide when to stop. Counts can differ between bins and machines.
The minimum takes precedence over time, so the budget is not a hard deadline.

Manual uses exactly the entered number in every bin. Its initial field value
is 5,000; select a count deliberately. The 50,000 cap belongs to the default
Auto configuration, not to Manual. More simulations improve Monte Carlo
precision, not biological sample size or power.

The parallel GUI uses `parfeval` for progress and cancellation; noninteractive
parallel script calls use `parfor`. The implementation uses CPU workers, not a
GPU. Failed/unavailable parallel startup can fall back to serial execution.
Fixed-count runs use reproducible per-bin seeds. Auto counts can vary with
timing even when the underlying seed is unchanged.

### Read the plot

| Element | Meaning |
| --- | --- |
| Blue step | Permutation p-value across the full angular bin, not a measurement at only its center. |
| Dashed 0.05 line | Unadjusted reference threshold. |
| Pale green ribbons | Bins with p < 0.05; **Shade p < 0.05** is initially on. Not effect size or direction. |
| Gray capped bars | Pointwise 95 percent Monte Carlo simulation intervals; **Show simulation intervals** is initially off. |
| Open orange markers | Simulation interval includes 0.05. Threshold classification is numerically uncertain. |
| Focus near p = 0.05 | Restricts the displayed y-range to 0-0.12; does not recalculate anything. |

If a bin is flagged uncertain, inspect its interval and the median profiles.
A larger, predefined Manual run can clarify numerical precision, but do not
repeat runs until a preferred significance label appears. Report a result near
the threshold as borderline when appropriate. A p-value of zero means **no
simulated exceedances**, not an exactly zero probability.

**Export Results...** writes CSV or MAT with per-bin counts and interval limits.
MAT also retains structured run metadata. The p-value column is still called
`LegacyP` internally for compatibility; the plot calls it **Permutation p-value**.
**Export Figure...** saves a vector PDF with the display options applied. Export plots
are redrawn rather than saved as exact screenshots of arbitrary manual zooms.

**Cancel** discards partial results. Cancellation during pool startup may wait
for startup to finish. Changing test settings or loaded group identities clears
the previous result; export before changing them. Exported permutation MATs are
for records/programmatic use, not reloadable comparison sessions in the GUI.

## Width Analysis Has Moved

Width analysis is now in **EyeLayerWidth**, a separate MATLAB GUI and repository.
Launch `EyeLayerWidth_GUI` from that project folder. Use Add Groups to open your
existing full group MATs, or Example Groups for the independent test-data copies.
No RpEGEN extraction rerun is needed. EyeLayerWidth retains per-eye profiles and
geometry inspection, Eye Summary, IQR, group medians with both CI methods,
contextual help, and all-eye exports. Consult its own User Guide for instructions.

## Save a Reproducible Analysis

Keep the input images/ROIs, full group MATs, `_run.mat` files, QC PDFs, exported
tables/figures, and a copy or release identifier of the code used. GUI sessions
are not saved automatically: loaded group lists, renamed plot labels, cached
bootstrap calculations, and permutation results do not survive closing
unless separately recorded/exported.

Record at least:

- MATLAB release, toolbox versions, analysis date, code snapshot, group
  definitions, independent-animal counts, and any exclusions with reasons.
- Image orientation, TIFF page, acquisition/calibration information, and ROI criteria.
- Median/heatmap bin widths, heatmap normalization, and which CI bands were shown.
- Bootstrap count/seed and the independent-eye assumption when using that band.
- Permutation pair, bin width, mode, actual per-bin counts, seed, threshold,
  pixel-level unit, and absence of multiple-testing correction.

Export permutation results to MAT as well as CSV when you need their
full metadata. Do not overwrite input group MATs with these exports.

## Troubleshooting

| Symptom | What to check or do |
| --- | --- |
| `RpEGEN_GUI` is not found | Set MATLAB's Current Folder to the project root. Use `which RpEGEN_GUI -all` to identify duplicate copies on the path. |
| A callback cannot be found after updating files | Finish calculations, close the old GUI window, and relaunch from the intended project folder. Restart MATLAB if stale callbacks persist. |
| Cannot acquire a MATLAB license | Resolve MATLAB sign-in/license availability before continuing; restarting MATLAB helped in development but is not a guaranteed license fix. |
| Check Inputs reports mismatched names | Match image filename, ROI filename, internal FIJI ROI name, and case; check for extra unmatched files. |
| Selected page is absent or not grayscale | Inspect the TIFF stack in FIJI and choose an 8-bit grayscale brightfield page. Do not simply assume page 3. |
| Index error during intensity CI calculation | The retained extraction code is unsafe for sparse bins (fewer than six pixels). Stop and investigate; Check Inputs does not detect this. Changing a display bin cannot fix a failed extraction. |
| Only one group after analysis | Expected: each extraction produces one group MAT. Run/load the other group files separately. |
| Five groups after Example Groups | Another group was already loaded, possibly a newly generated copy of an example. Remove the extra table entry; no disk deletion is needed. |
| MAT file is rejected | Choose a full `<group>.mat`, not `<group>_run.mat`, a width export, or a permutation export. |
| Very narrow pixel-based CI | It may lie beneath the median line. Zoom or toggle the bootstrap band off to inspect it; narrowness is not proof of independent-animal precision. |
| Bootstrap CI is missing | Read the status line for sparse bins, unmatched eye IDs, or inconsistency between raw and saved medians. |
| Run Test is disabled | Load at least two distinct group files, select different A/B files, and wait for other calculations to finish. |
| Permutation reports an empty bin | Both groups must contribute pixels to every tested bin. Review inputs and scientifically appropriate binning; empty bins are not silently skipped. |
| Simulation bars are invisible | Enable Show simulation intervals and Focus near p = 0.05. With many draws the bars may genuinely be small. |
| QC PDF is missing | Load the full group MAT in Setup and use Quality Control > Export PDFs to regenerate all-eye PDFs without rerunning extraction. |
| Preview reports a permission error when launched from MATLAB | Close and relaunch the updated GUI. On macOS the button now reads Show QC PDF in Finder; double-click the selected file there. No permission changes or analysis rerun are needed if the file opens from Finder. |

For support, provide the full MATLAB error, release/toolboxes, active GUI tab,
settings, code snapshot, and a shareable example that reproduces the issue.

## Citation

Leach LL, Fisher GB, Gross JM (2022). Nitroreductase/Metronidazole-Mediated
Ablation and a MATLAB Platform (RpEGEN) for Studying Regeneration of the Zebrafish
Retinal Pigment Epithelium. *Journal of Visualized Experiments*, 181.
[doi:10.3791/63658](https://doi.org/10.3791/63658).

Cite the original method and identify this GUI code snapshot separately. The
new bootstrap and runtime controls should not be described as
features validated by the 2022 publication.
