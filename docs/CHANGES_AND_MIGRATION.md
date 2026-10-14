# RpEGEN Changes and Migration

This document compares the published 2022 script workflow with the current GUI
release **v2.0.0 (2026-10-14)**. It is a migration and methods record, not a
claim that every new feature was part of, or validated by, the original paper.
For operating instructions, see the [User Guide](USER_GUIDE.md).

The main transition is from editing scripts to using a GUI. The pooled-pixel
intensity calculation and pixel-level permutation statistic are retained.
Whole-eye bootstrap bands are an addition. Width analysis has moved to the
separate EyeLayerWidth repository and is no longer part of the RpEGEN GUI.
Plots and permutation runtime defaults have changed, so visual or stochastic
output need not exactly match a previously exported figure.

## Published Baseline

The baseline is Leach, Fisher, and Gross (2022),
[doi:10.3791/63658](https://doi.org/10.3791/63658), especially protocol sections
6-8 and Figures 4, 9, and 10. The workflow used FIJI-prepared images/ROIs,
`RpEGEN.m` for extraction, and `RpEGEN_PermPlot.m` plus GRAMM for comparisons and
figures. The published example used 20,000 permutations per one-degree bin.
[Full article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9036407/).

The files with those original script names are still present but have been adapted.
Download the [Pre-GUI Script Archive](https://github.com/burchfisher/RpEGEN/releases/tag/pre-gui-archive)
for the unchanged April 21, 2023 snapshot at commit
`1e8d10bb879b4db61ed6b9f7f6a65d3e65cc9bc8`. This is the last pre-GUI version,
including fixes made after the paper; it is not claimed to be the exact paper-era
code. Earlier revisions remain available in Git history. Keep the archived code
in a separate folder for reproduction of analyses that used that snapshot.

## Move an Existing Project

1. Back up the original code, input images/ROIs, group MATs, and published figures.
2. Place the new code in a separate folder; do not mix individual files from
   old and new versions. Keep the bundled ROI reader with it.
3. Finish calculations and close any GUI instance using an older code copy.
4. In MATLAB, change to the new project folder and run `RpEGEN_GUI`.
5. Use **Setup > Open Existing Results** for an existing full group MAT.
   Inspect representative eyes and the Results tables.
6. Use **Compare Groups > Add Groups...** to load each group's MAT. Do not
   select `_run.mat` files, exported statistics, or multiple copies of one group.
7. Check heatmap binning/normalization and CI display settings before exporting.
8. Run new permutation comparisons deliberately: use Manual 20,000 for the
   published example's count, or record the current Auto settings and actual counts.
9. For width measurements, open existing geometry in the separate EyeLayerWidth
   application; follow its User Guide and save separate exports.
10. Archive the code snapshot and outputs used for the new analysis.

Existing compatible results do **not** need intensity extraction repeated just
to use the new interface. Re-extract when images, ROIs, group membership,
orientation, or the selected TIFF page change, or when required saved fields
are unavailable. A MAT file is not upgraded or overwritten merely by opening it.

## Workflow Comparison

Contextual help is available through the top-right **?** button on every main
tab and follows the selected Compare Groups subtab.

QC PDF generation shares the GUI drawing code: all three panels overlay the
brightfield image with identical extents, plus the intensity profile below.
All three GUI panels reserve matching space for colorbars. **Export PDFs** in QC
regenerates the whole group's PDFs from saved results without extraction.

On 2026-10-14, the Widths tab and its geometry helpers/tests were separated into
EyeLayerWidth. Its independent GUI retains linked profiles, contextual help,
per-eye and all-eye exports, and pooled group curves with point-based and
whole-eye bootstrap confidence bands. RpEGEN's intensity calculations and
permutation behavior were not changed by this separation.

| Task | Published script workflow | Current GUI workflow |
| --- | --- | --- |
| Setup | Edit paths, output name, and TIFF-page variable in `RpEGEN.m`. | Setup fields, Example Data, Check Inputs, and Run Analysis. |
| Input checks | Script name matching and existing checks. | Additional preflight validation of TIFF pages/type and ROI names/type/area/bounds. Anatomical QC remains manual. |
| Extraction | Execute script and inspect outputs. | Adapter calls the retained extraction calculations, reports progress, and isolates its workspace. |
| Output records | Group MAT and per-eye QC PDFs. | Same primary outputs, plus `<group>_run.mat` with settings, manifest, and log. |
| QC | Review individual exported plots. | Image/profile viewer with arrows, dropdown, sample counter, and keyboard navigation; PDFs remain available. |
| Group loading | Edit/load group variables in MATLAB. | Add Groups, Add Current Results, or Example Groups. One group per extraction is unchanged. |
| Tables | Inspect MATLAB structures. | Results tab and CSV exports of saved 1-, 5-, and 10-degree tables. |
| Heatmaps | GRAMM plotting code and chosen bin counts. | Built-in MATLAB plots, explicit shared edges/color limits, optional percent-per-angle normalization. |
| Median profiles | Supplied script spline-interpolates the example binned medians and CI limits. | Saved bin values connected directly, without spline interpolation; independent CI toggles. |
| Median uncertainty | Saved pixel-based median limits. | Original limits retained, with optional whole-eye bootstrap intervals overlaid. |
| Permutations | User-specified fixed count; published example 20,000. | Auto runtime/precision planning or exact Manual count, CPU parallelism, cancellation, and diagnostics. |
| P-value display | Supplied plotting script uses bars and example-specific anatomical shading. | Bin-spanning steps, optional simulation intervals, uncertain-bin markers, and optional unadjusted p < 0.05 ribbons. |
| Width analysis | Outside the published RpEGEN scope. | Separate EyeLayerWidth application; no Widths tab in RpEGEN. |
| Plot dependency | GRAMM for `RpEGEN_PermPlot.m`. | GUI uses built-in plotting; GRAMM remains bundled for the older script workflow. |

Publication figure exports now require PDF and explicitly request vector content
with a white background, including intensity group plots, permutation plots, and
new QC PDFs. Plot lines, bands, and text remain editable in compatible PDF
editors; eye images and heatmap image layers remain raster. The export helpers
no longer accept PNG filenames. Existing files are untouched; re-export figures
from loaded results rather than rerunning extraction. Existing QC PDFs are not
retroactively converted. Developer UI screenshots remain PNG.

## Calculations Retained

The GUI extraction adapter retains the existing ROI-mask/centerline workflow,
nearest-centerline pixel assignment, RpEGEN angular coordinates, pixel intensity
summaries, and pooled 1-, 5-, and 10-degree tables. The adapter adds validation,
configuration, progress reporting, and controlled QC export; it does not replace
the intensity statistic with an equal-eye-weighted statistic.

The group median remains the median of all contributing **pixels** in each
angular bin. Eyes with more pixels in a bin contribute more observations.
Neither the heatmap percentage option nor the new bootstrap changes that into a
median of per-eye medians.

The original saved `CI95L`/`CI95U` values remain available and are not overwritten
by bootstrap limits. Their calculation uses sorted pixel values at rounded
ranks `n/2 +/- 1.96*sqrt(n/4)`, a normal approximation rather than an exact
binomial 95 percent interval. The GUI label is now **Pixel-based median CI
(nominal 95%)**, replacing development-era "Legacy pixel CI" wording.

The permutation unit remains the **pixel within its angular bin**. The statistic
is the absolute between-group median difference, the exceedance comparison uses
`>=`, and p remains `b/B`. No whole-eye label shuffle, composite-eye procedure,
`(b+1)/(B+1)` correction, or multiple-testing adjustment has been introduced.

These retained assumptions matter: correlated pixels are not independent
animals. Agreement with the published calculation is numerical compatibility,
not evidence that the test provides animal-level treatment inference under all
experimental designs.

## New Statistical Display

The whole-eye bootstrap estimates pointwise uncertainty in the **pooled-pixel
median**, using 5,000 resamples with replacement and seed `20261012`. Each draw
selects whole eyes within a group and uses that same selection at every angular
position. Pixel values remain in their original angle bins. The 2.5th and 97.5th
percentiles form the displayed interval.

This adds an independent-eye resampling perspective without changing the curve
it surrounds. It assumes independent eye IDs and does not model shared wells,
batches, paired eyes, or repeated sections. Bins without adequate contributions
are not assigned an artificial zero-width band. These are pointwise intervals,
not simultaneous bands or a replacement for the permutation test.

Both CI types share one plot and can be shown separately. Pixel-based intervals
use darker shading; bootstrap intervals use lighter shading. The Results CSVs
still contain the saved pixel-based limits, not bootstrap values. The helper
`rpegen_bootstrap_ci` provides numeric bootstrap output for separate saving.

## Permutation Runtime Changes

| Setting | Current behavior |
| --- | --- |
| Default mode | Auto in the GUI and adapted plotting script. |
| Whole-comparison budget | Five minutes, soft estimate rather than timeout. |
| Auto minimum and maximum | 10,000 and 50,000 per bin. |
| Precision target | Approximate 95 percent simulation half-width 0.002 near p = 0.05; desired count 45,618. |
| Manual | Exactly the entered count in every bin; GUI field initially 5,000. |
| Count selection | Per-bin timing calibration before testing; not significance-dependent stopping. |
| Random seeds | Profile wrapper uses base seed `20261012` plus bin index. |
| Parallel path | Up to four CPU workers; GUI uses `parfeval`, noninteractive calls use `parfor`. Serial fallback is available. |
| Memory | Shuffled data processed in bounded batches instead of allocating all iterations together. |
| Diagnostics | Actual counts, exceedances, pointwise exact binomial simulation intervals, and threshold uncertainty flags. |

The 10,000 minimum can exceed the runtime estimate. Auto counts can differ with
bin sizes, workload, and timing; they are not a power calculation. More draws
reduce simulation error but cannot compensate for too few independent animals.

For fixed-count direct `permsim` calls, tests verify parity with the original
shuffle calculation and random-number state when inputs, count, and starting RNG
state are the same. The seeded profile wrapper need not reproduce an old
unseeded run exactly, even with the same count. Auto can choose different counts
between runs because it measures runtime. For exact repeatability, preserve
counts/seeds and use the corresponding fixed-count settings.

Simulation intervals measure Monte Carlo uncertainty, not biological variability.
Zero means no observed exceedances, and `NearThreshold` means the interval
includes 0.05. The exported field `LegacyP` retains its name for compatibility;
the user-facing label is **Permutation p-value**.

## Plot Differences to Expect

Heatmaps now use fixed edges from 0 to 180 degrees and 0 to 255 intensity units
and shared color limits across shown groups. This can look different from
GRAMM's automatic bin placement/scaling. Percent-per-angle is an additional
normalization, not a different extraction.

Median curves no longer use the supplied plotting script's spline interpolation.
Identical saved bin medians can therefore produce visibly different connecting
lines. The additional bootstrap shading is not a change to saved median values.

The supplied older plotting script hard-codes a central example region around
41-142 degrees. The GUI does not automatically identify that region. Green
ribbons in the new p-value plot instead mark bins with unadjusted p < 0.05;
they are not an anatomical injury-zone annotation or an effect-size display.

### Endpoint conventions

The retained intensity tables, whole-eye bootstrap, and pixel permutation test
use lower-inclusive, upper-exclusive angle intervals. Thus exactly 180-degree
pixels are excluded from their final intensity bin. GUI heatmaps include exactly
180 degrees and intensity 255 in the final histogram bins. EyeLayerWidth profiles
include exactly 180 degrees in their last bin. This difference is deliberate
compatibility behavior to record when comparing counts across views, not a
silent repair of the published intensity tables.

## Width Analysis Separation

EyeLayerWidth is a standalone successor to the exploratory width workflow. Existing
full RpEGEN group MATs can be opened there without conversion. It also bundles
geometry-only copies of all 44 example eyes. Its local-normal algorithm, total
25-unit window default, optional endpoint exclusion, IQR, and pooled group
confidence-band calculations are retained. See EyeLayerWidth's migration record
and User Guide for method details and limitations. This separation does not
validate the measurements for photoreceptor regeneration or change the original
RpEGEN image/intensity scope.

## Compatibility and Limitations

- Full original-format group MATs can be loaded without conversion. QC needs
  images/ROI/centerline fields; comparisons need raw pixels and all three binned
  tables. Bootstrap additionally requires consistent unique eye identities.
- Group MATs still contain one group per extraction run. A file with several
  candidate structures is not a supported multi-group import: loaders select
  the first compatible structure.
- `_run.mat`, width-export MATs, and permutation-export MATs are different formats,
  not replacement group results or saved GUI sessions.
- GRAMM is optional for GUI users. The retained `RpEGEN_PermPlot.m` still needs
  GRAMM, edited paths/group variables, and its script-specific plot configuration.
- The retained extraction's median-CI indexing can fail in bins with fewer than
  six pixels. Input preflight does not detect sparse angular bins; all three
  intensity resolutions are computed during extraction. This remains unresolved.
- The permutation workflow requires nonempty bins in both groups. Bootstrap
  cannot provide a band from fewer than two contributing eyes or empty resamples.
- Check Inputs is not an anatomical approval or complete polygon-topology check.
- QC PDFs are named by sample, so use separate output folders if sample names
  recur in different groups. Failed extraction may leave partial outputs.
- There is no GUI ROI editing, saved QC approval list, saved comparison session,
  or GPU permutation path.

## Verification Record

The development checks include reconstruction of saved medians and pixel-based
CI limits in all 936 bins of the four supplied groups, fixed-count permutation
and RNG parity, literal whole-eye bootstrap parity, serial/parallel agreement,
GUI navigation/toggles/cancellation/exports. Width-specific geometry and GUI
tests now live in EyeLayerWidth. These checks establish specific implementation
properties, not universal biological validity.

Recent GUI verification used MATLAB R2026b. A measured example
permutation comparison completed in about 149 seconds with four CPU workers
and 45,618 draws in each of 180 bins on the development machine. Runtime is not
a guarantee for other computers or datasets.

From the project root, developers can run:

```matlab
addpath('tests')
test_rpegen_resampling
test_rpegen_confidence
test_rpegen_permutation_gui
```

Some GUI verification tests open windows and write preview/export files. Run them
with no active analysis and preserve any previews you need beforehand. See
[tests](../tests) for the additional integration checks.

## Reporting the Updated Workflow

Identify the original method and the exact updated code used. A reporting
checklist is provided in the [User Guide](USER_GUIDE.md#save-a-reproducible-analysis).
State which additions were actually used rather than naming all available
features. In particular, distinguish:

- the pooled-pixel median from a median of animal summaries;
- pixel-based nominal median intervals from whole-eye bootstrap intervals;
- Monte Carlo intervals from biological confidence intervals;
- pixel-level permutation evidence from independent-animal inference.

Continue citing Leach et al. (2022) for RpEGEN. Do not attribute the new bootstrap
implementation or EyeLayerWidth methods to the JoVE publication. Preserve the current [license](../LICENSE) and bundled
dependency attribution when distributing the updated project.
