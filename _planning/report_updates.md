# REPORT.Rmd — Required Updates

**Source review:** `pipeline_refactor_log.md` (Phases 1–6), `app.R`, Scripts 6a, 6b, 7, 8, 9.
**Status:** All six pipeline refactor phases are complete; `REPORT.Rmd` has not been touched and still reflects the pre-refactor state.

---

## Summary of required changes

The report was written against three methods (`davenport`, `adaptive`, `quartile`). The pipeline now has three different methods (`depth_thresh`, `frac_thresh_3class`, `frac_thresh_2class`). The `quartile` method has been removed entirely. The `adaptive` method (GMM + SVM bivariate boundaries) has been replaced by the `frac_thresh` family (univariate rain_frac threshold + population-median split). Across the report this affects: data loading, shared style constants, and five sections of content (Introduction, Sections 4–5, and Results).

---

## 1. Setup chunk — data loading and shared constants

### 1a. Data loading (lines 33–38)

The three `readRDS` calls at the top of the setup chunk load files that no longer exist as named:

| Current | Required |
|---|---|
| `readRDS("../_data/pot_peaks_classified_davenport.rds")` → `peaks_dav` | → `readRDS("../_data/pot_peaks_classified_depth_thresh.rds")` → rename object `peaks_dt` (or keep `peaks_dav` if preferred, but the filename must change) |
| `readRDS("../_data/pot_peaks_classified_adaptive.rds")` → `peaks_adap` | → `readRDS("../_data/pot_peaks_classified_frac_thresh_3class.rds")` → rename object `peaks_frac` |
| `readRDS("../_data/pot_peaks_classified_quartile.rds")` → `peaks_qrt` | Remove entirely — quartile method no longer exists |

The FFA data files (`ffa_q`, `ffa_p`) can stay but the method keys used to filter them (`"adaptive"`, `"davenport"`, `"quartile"`) must be updated to `"frac_thresh_3class"`, `"frac_thresh_2class"`, and `"depth_thresh"` wherever they appear in code chunks.

### 1b. `TYPE_COLORS` constant (lines 44–50)

Remove the `"ROS"` entry — `"ROS"` was standardized to `"rain_on_snow"` in Phase 3. The constant should read:

```r
TYPE_COLORS <- c(
  "Rainfall"      = "#b2182b",
  "Snowmelt"      = "#2166ac",
  "Rain-on-Snow"  = "#9970ab",
  "Unclassified"  = "#cccccc"
)
```

---

## 2. Section 2 — Classification Approaches Overview (lines 136–166)

The overview paragraph currently describes three methods as:
1. Davenport et al. (2020) fixed threshold
2. Adaptive rain-fraction method (GMM + SVM bivariate)
3. Quartile method

This must be updated to describe the actual current methods:
1. **Depth Threshold** (`depth_thresh` / Script 6a) — the Davenport fixed-depth rules. Description is broadly correct but the method name must be updated throughout.
2. **Fractional Threshold 3-class** (`frac_thresh_3class` / Script 6b) — univariate `rain_frac` threshold: `rain_frac >= 0.99` → Rainfall; remaining events split at the population median (`MIXED_SPLIT = 0.2979`) into Snowmelt (lower half) and Rain-on-Snow (upper half). **No GMM, no SVM, no SWE feature.** This is a substantial rewrite — the current description is entirely wrong for the new method.
3. **Fractional Threshold 2-class** (`frac_thresh_2class`) — a Script 7 downstream remap collapsing Snowmelt and Rain-on-Snow into Non-Rainfall. Not a standalone script output.

The final sentence about the quartile method should be removed.

---

## 3. Section 4 — Adaptive Classification (lines 353–826)

This entire section describes the GMM+SVM bivariate adaptive classifier that no longer exists in the pipeline. It must be substantially rewritten to document the new `frac_thresh_3class` method. Key changes:

### 3a. Section 4.1 Method Description (lines 356–379)

- Remove all references to SWE as a classification feature — the new Script 6b uses `rain_frac` only.
- Remove the `swe_mean_mm` / `log10(swe_mean_mm + 1)` bivariate feature discussion.
- Replace with the actual classification logic: (1) Unclassified if window coverage < 90% of `search_radius` or `rain_frac` is NA; (2) Rainfall if `rain_frac >= RAIN_FRAC_THRESH` (0.99); (3) Snowmelt if `rain_frac < MIXED_SPLIT` (population median of mixed events = 0.2979); (4) Rain-on-Snow if `rain_frac >= MIXED_SPLIT`.
- Add the operational-label caveat that appears in the Script 6b header: "Snowmelt" and "rain_on_snow" labels are operationally defined as below/above the population median of mixed events, not direct physical process attribution.
- Document the `RAIN_FRAC_THRESH = 0.99` tolerance rationale (absorbs floating-point artifacts and trace-melt days just above the noise floor).
- Document `NOISE_FLOOR_MM = 0.1` mm/day.

### 3b. Section 4.2 Decision Tree flowchart (lines 383–440)

The `adaptive-flowchart` code chunk renders a GMM+SVM decision tree (B1/B2 boundaries in 2D space). Replace with a new flowchart matching the current Script 6b logic:

```
POT peak → Complete window? (n_days >= 0.90 × search_radius)
  → No  → Unclassified
  → Yes → Compute rain_frac over search_radius window
         → rain_frac >= 0.99?
           → Yes → Rainfall
           → No  → rain_frac < MIXED_SPLIT (0.2979)?
                   → Yes → Snowmelt
                   → No  → Rain-on-Snow
```

### 3c. `fig-gmm-fit` code chunk (lines 442–622)

Remove entirely. This chunk re-runs the GMM and SVM fitting (`mclust`, `e1071`) to regenerate Figure X (GMM component plot with B1/B2 boundaries). The GMM method is no longer in the pipeline; re-running it would be inconsistent with the current outputs. The associated figure caption (lines 625–625) should also be removed.

Remove the `suppressPackageStartupMessages(library(mclust))` and `suppressPackageStartupMessages(library(e1071))` calls.

### 3d. `fig-adaptive-space` code chunk (lines 627–701)

This chunk plots classification in `rain_frac × log10_swe` space using `peaks_adap` and the B1/B2 SVM boundaries. Replace with a new visualization appropriate for the frac_thresh method: a 1-D density or histogram of `rain_frac` across all data-sufficient peaks, with vertical lines at `RAIN_FRAC_THRESH = 0.99` and `MIXED_SPLIT = 0.2979`, colored by the resulting flood type. Data source becomes `peaks_frac` (the new `pot_peaks_classified_frac_thresh_3class.rds` object).

The associated Figure X caption (line 703) must be rewritten to describe the new plot.

### 3e. Sections 4.3–4.6 (Feature Selection, GMM Boundary Derivation, Decision Tree Ordering, Noise Floor)

- **Section 4.3** (lines 705–735): Discusses `ros_flux_frac` removal and `SNOW_MULT` window. These decisions belong to the now-archived GMM/adaptive version. Remove or replace with a brief rationale for why `rain_frac` alone (without SWE) was sufficient for the fractional threshold method.
- **Section 4.4** (lines 737–776): GMM boundary derivation. Remove entirely — there is no GMM in the current method.
- **Section 4.5** (lines 778–797): Decision tree ordering (Snowmelt before ROS). This discussion applies to the old bivariate boundary approach. Remove or rewrite as a note about the sequential evaluation order in Script 6b (Rainfall first, then MIXED_SPLIT split).
- **Section 4.6** (lines 799–812): Noise floor decision. Partially still relevant (NOISE_FLOOR_MM = 0.1 is still used in Script 6b for rain_frac computation), but the specific diagnostic references to the flux-floor sensitivity sweep and the WLDAS mismatch analysis should be reviewed for whether they still apply.
- **Section 4.7** (lines 814–826): Unclassified events. The coverage criterion is unchanged (< 90% of `search_radius` window), so this section is mostly still correct but references to `swe_mean_mm` / `data_sufficient` should be updated to match the Script 6b column names (`n_days_data`, `rain_frac` NA check).

---

## 4. Section 5 — Quartile Classification (lines 829–909)

**Remove entirely.** The quartile method was deleted in Phase 1. This includes:
- Section 5.1 (Motivation), lines 831–849
- Section 5.2 (Middle-Bin Sparsity) including the `quartile-middle-bin` code chunk, lines 851–899
- Section 5.3 (Retained as Comparison Method), lines 901–909

The `peaks_qrt` object referenced in this section (and loaded in the setup chunk) no longer exists.

---

## 5. Results section — Section 8

### 5a. Section 8.1 Classification Composition figure (lines 998–1069)

The `composition-figure` chunk builds a comparison bar chart across three methods: `"Davenport"`, `"Adaptive"`, `"Quartile"`. All three data objects, method labels, and the `label_qrt` block must be updated:

- Rename `label_dav` → sourced from `peaks_dt` (depth_thresh), method label `"Depth Threshold"`.
- Rename `label_adap` → sourced from `peaks_frac` (frac_thresh_3class), method label `"Frac Thresh (3-class)"`. The `flood_type_raw` factor now uses `"rain_on_snow"` (not `"ROS"`).
- Remove `label_qrt` block entirely.
- Update `comp_data` binding to remove `label_qrt`.
- Update `comp_colors` to remove `"Q1 (most melt)"`, `"Q2"`, `"Q3"`, `"Q4 (most rain)"` entries.
- Update `method` factor levels from `c("Davenport", "Adaptive", "Quartile")` to `c("Depth Threshold", "Frac Thresh (3-class)")`.
- Optionally add a third bar for `frac_thresh_2class` if desired.

### 5b. Section 8.2 Cross-Method 100-Year Flood Comparison (lines 1072–1138)

The `cross-method-scatter` chunk filters `ffa_q` for `method == "adaptive"` and `method == "davenport"`. Both method keys no longer exist in the FFA outputs. Update to:

- `method == "frac_thresh_3class"` in place of `"adaptive"`
- `method == "depth_thresh"` in place of `"davenport"`
- Update axis labels: `"Depth Threshold — 100-yr mixture quantile"` and `"Frac Thresh (3-class) — 100-yr mixture quantile"`.
- Update `dominant_type` derivation: `flood_type_raw` for `peaks_frac` uses `"rain_on_snow"` not `"ROS"`.
- Update `dominant` case_when: `"ROS"` → `"rain_on_snow"`.
- Update chart title and subtitle.

### 5c. Section 8.3 Dominant Flood Type by Site (lines 1140–1150)

Placeholder text references "adaptive classification". Update to `"frac_thresh_3class"` or `"depth_thresh"` per the intended map method. Both are now available in the app.

### 5d. Section 8.4 GPD Shape Parameter Distributions (lines 1151–1168)

The code stub filters `ffa_p` by `curve_type` and plots `xi` grouped by `flood_type`. The stub will work as-is once the data object is updated (`ffa_p` still comes from the same file). However, the `facet_wrap(~method)` will now show `depth_thresh`, `frac_thresh_3class`, `frac_thresh_2class` rather than the old method names. The stub itself can stay but the surrounding placeholder text should reference the new method names.

---

## 6. Section 6 — Minimum Sample Size Thresholds (lines 912–940)

Lines 933–939 describe a correction to single-sample pool ordering that was applied across "all three methods" including Script 6c (quartile). Since Script 6c has been deleted, revise the reference to "all three methods" to "both fixed-threshold methods" or similar. The correction itself is still historically accurate and worth retaining.

---

## 7. Section title and section numbering

If Section 5 (Quartile) is removed, renumber:
- Current Section 6 → Section 5
- Current Section 7 → Section 6
- Current Section 8 → Section 7
- etc.

Update all cross-references within the document (e.g., "Section 4.4" in Section 4.6 becomes irrelevant once Section 4.4 is removed).

---

## 8. `library()` calls in setup chunk (lines 22–30)

With the GMM figure chunk removed, the following libraries are no longer needed:
- `library(arrow)` — not used elsewhere in the report (verify before removing)
- `mclust` and `e1071` — called inline in the removed `fig-gmm-fit` chunk, not loaded in the setup chunk but referenced via `suppressPackageStartupMessages` inside the chunk

The setup chunk `library(DiagrammeR)`, `library(DiagrammeRsvg)`, `library(rsvg)` are still needed for the Davenport flowchart in Section 3.2. Retain them unless that figure is also being replaced.

---

## Priority order for edits

1. **Setup chunk** — fix data loading (broken file paths will prevent the document from knitting at all)
2. **Section 4** — full rewrite of method description and figure chunks
3. **Section 5** — delete quartile section
4. **Results Section 8** — update method keys in all code chunks
5. **Section 2** — update overview prose
6. **Section 6 / numbering** — minor cleanup
