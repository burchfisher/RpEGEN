# RpEGEN v2.0.0

Released October 14, 2026. This release adds a MATLAB GUI to the image and
intensity-analysis workflow, with documentation and regression tests.

## Archived Scripts

The [Pre-GUI Script Archive](https://github.com/burchfisher/RpEGEN/releases/tag/pre-gui-archive)
preserves commit `1e8d10bb879b4db61ed6b9f7f6a65d3e65cc9bc8` from April 21, 2023.
It contains the unchanged pre-GUI scripts, dataset, and dependencies. It includes
post-publication fixes and is not labeled as the exact version used in the 2022
paper. Older commits remain accessible in repository history.

## Highlights

- Setup controls, input checking, example-data loading, extraction progress,
  and saved run settings.
- Image-based QC with adjacent navigation buttons, left/right keyboard arrows,
  and group-wide PDF exports matching the GUI overlays.
- Saved intensity tables and CSV exports.
- Group heatmaps and median profiles with independently selectable pixel-based
  and whole-eye bootstrap confidence bands.
- GUI permutation comparisons with Manual counts or Auto runtime/precision
  planning, CPU parallelism, cancellation, and simulation-precision diagnostics.
- Editable PDF figure exports, contextual help, a User Guide, and a detailed
  changes-and-migration record.

The existing public example dataset and bundled dependencies are unchanged.
Generated local analysis outputs are not part of this release.

## Start

Open the project folder in MATLAB and run `RpEGEN_GUI`. Follow the
[User Guide](USER_GUIDE.md) for the example workflow. Existing compatible group
MATs can be loaded without rerunning extraction.

Verified on MATLAB R2026b. Image Processing, Curve Fitting, and Statistics and
Machine Learning toolboxes are required; Parallel Computing Toolbox is optional.
The GUI uses built-in plots. GRAMM remains bundled for the script plotting workflow.

## Analysis Scope and Limitations

Pooled-pixel medians and the pixel-level permutation statistic are retained.
The whole-eye bootstrap bands are a separate addition; they do not change the
permutation resampling unit. Auto plans a count before the final test, targeting
45,618 draws when time permits, with a 10,000 minimum and 50,000 ceiling.

Pixels within an eye are related, and pixel-level uncertainty and p-values do
not account for that relationship. P-values are unadjusted for multiple testing.
The retained p-value is exceedances divided by draws, without a plus-one adjustment.
The original intensity confidence-limit calculation still requires sufficiently
populated bins; fewer than six pixels can cause extraction to fail. See
[Changes and Migration](CHANGES_AND_MIGRATION.md) for these and other assumptions.

For reproducibility, record the software version, source data, settings, and
exported results. Keep each archived software version in its own folder.

## Release Verification

All nine regression tests and the real-data GUI integration check passed on
MATLAB R2026b from a clean checkout:

- Fresh extraction of 11 DMSO 4 dpi eyes reproduced the saved raw intensity data
  and all three binned tables exactly, with 11 QC PDFs and a saved run record.
- The confidence tests reproduced all 936 saved example medians and limits.
  Seeded permutation parity, whole-eye bootstrap calculations, serial/parallel
  agreement, cancellation, loading, navigation, help, and PDF exports passed.
- The full DMSO 4 dpi versus IWR1 4 dpi Auto comparison completed 180 one-degree
  bins in 263.4 seconds using four CPU workers and 10,872-45,618 draws per bin.
  Counts varied with the timing estimates. This is one verification run, not a
  runtime guarantee or biological validation.

These checks were run locally, not by GitHub CI. See the README for commands.
