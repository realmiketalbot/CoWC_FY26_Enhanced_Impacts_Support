# Logic Review — Findings Report

Scope: `_scripts/6a_classify_usgs_unit_data_peaks_davenport.R`,
`_scripts/6b_classify_usgs_unit_data_peaks_adaptive.R`,
`_scripts/6c_classify_usgs_unit_data_peaks_quartiles.R`,
`_scripts/7_compile_ffa_results.R`, `functions.R`, `app.R`.

**Status:** initial findings produced without edits; items 1, 2 (with the
6a/6c-aligned-to-6b decision), and 3 from the action list below have since
been applied. See "Resolutions" at the end of each affected section.

## 1. `flood_type` vs `flood_type_raw` consistency

**No issues found.** The pre/post-demotion column convention is applied correctly in every location traced:

- **Script 6a** (davenport): `flood_type_raw` set from thresholds (lines 122–128); `flood_type` overwritten by demotion (lines 182–190).
- **Script 6b** (adaptive): decision-tree result assigned to `flood_type` (lines 190–199), copied to `flood_type_raw` at line 202 *before* the factor conversion, then `flood_type` overwritten by demotion (lines 285–295).
- **Script 6c** (quartile): `flood_type_raw` set from quartiles (lines 76–82); `flood_type` overwritten by demotion (lines 139–148).
- **Script 7** `remap_2class` correctly reads from `flood_type_raw` (line 109) so peaks demoted in 3-class are still eligible for the pooled non-rainfall type — the documented intent.
- **Script 7** `class_rows` (lines 270, 275) correctly uses `flood_type_raw` so the per-method composition reflects the decision tree, not the MIN_TYPE_SAMPLE-filtered view.
- **Script 7** `ros_rows` (line 303) uses `flood_type` (post-demotion). This is a deliberate choice — only "fittable" ROS peaks count toward the temporal trend — and is internally consistent.
- **app.R** `site_peaks_labeled` 3-class branch (line 417) uses `flood_type` (post-demotion); 2-class branch (lines 422–428) uses `flood_type_raw`. Both correct.
- **app.R** `site_emp_peaks` 2-class non-rainfall pool (line 454) uses `flood_type_raw`, mirroring `remap_2class`.

## 2. Fit-to-plot consistency

**One real inconsistency for the single-sample annotation; per-type and mixture fits are consistent.**

### 2a. Per-type curves (davenport, adaptive 3-class, quartile)
- **Script 7** `get_sample(typ)` for non-pooled types (line 327): `filter(flood_type == typ, in_type_sample)`.
- **app.R** `site_emp_peaks` non-2class branch (line 448): `filter(in_type_sample)`, then narrowed in the FFA render at line 545 by `flood_type %in% type_curve_types`.

These produce the same set: by construction in scripts 6a/6b/6c, `in_type_sample == TRUE` implies `flood_type` was *not* demoted (line 174 in 6a, 268 in 6b, 131 in 6c all gate on `n() >= MIN_TYPE_SAMPLE` before the flag is assigned), and `type_curve_types` is exactly the set of types whose fits succeeded.

The previously surfaced **ROS-below-MIN_TYPE_SAMPLE failure cannot recur**: when ROS has < 20 peaks, `in_type_sample == FALSE` for every ROS peak at that site, no ROS fit is attempted, no ROS curve appears in `curves_df`, and the render's `flood_type %in% type_curve_types` filter excludes those points. Confirmed for adaptive, davenport, and quartile.

### 2b. Per-type curves (adaptive 2-class)
- **Script 7** `get_sample("non_rainfall")` (lines 316–325): post-`remap_2class`, filters `flood_type == "non_rainfall"`, then `slice_max` per water year, then top-n at MIN_TYPE_SAMPLE.
- **app.R** `site_emp_peaks` 2-class branch (lines 451–465): filters `flood_type_raw %in% non_rainfall_types`, then `slice_max` per water year, then top-n at 20.

Equivalent — `flood_type == "non_rainfall"` after `remap_2class` is exactly `flood_type_raw %in% NON_RAINFALL_TYPES`.

For 2-class rainfall, both sides filter `flood_type_raw == "rainfall" & in_type_sample`.

### 2c. Single-sample empirical annotation — **RESOLVED**
- **Script 7** uses `peaks_single <- site_peaks %>% filter(in_single_sample)` (line 332), where `in_single_sample` carries the script-6 logic: `flood_type != "unclassified"` filter **and** `slice_max(peak_flow, n=1)` per water year **before** the top-`m` rank.
- **Original app.R** `ss_ids` derivation re-derived top-n from `site_peaks()` without either of those steps — no `flood_type != "unclassified"` filter, no per-water-year `slice_max`. The open-circle "In Single Sample" annotation could therefore flag peaks that were *not* in the fitted sample.

**Resolution (applied):** `ss_ids` (app.R lines 558–560) now reads the `in_single_sample` flag from the RDS directly:

```r
ss_ids <- site_peaks() %>%
  dplyr::filter(in_single_sample) %>%
  dplyr::pull(peak_flow)
```

Annotation now matches script 7's `peaks_single` by construction.

## 3. `in_type_sample` / `in_single_sample` flag logic

### 3a. `in_type_sample`
**No issues.** Assigned consistently across 6a/6b/6c — top-n per (site, type, water-year) gated on `n() >= MIN_TYPE_SAMPLE`. Used identically by script 7 (line 327) and app.R (line 448).

### 3b. `in_single_sample` — **RESOLVED (6a/6c aligned to 6b)**
- **Original 6a** (lines 192–203): single-sample selection ran *after* demotion. The `flood_type != "unclassified"` filter saw the demoted column, so demoted peaks were excluded from the single-sample pool.
- **Original 6c** (lines 138–161): same ordering — demote first, then single-sample.
- **6b** (lines 245–256, then 285–295): single-sample selection runs *before* demotion. Treated as the canonical Yu et al. (2022) ordering.

**Consequence (before fix):** at sites where a 3-class type was demoted under adaptive (e.g., ROS = 10 peaks at site X), those ROS peaks were eligible for adaptive single-sample selection but the equivalent peaks under davenport/quartile were not. Script 7's `peaks_single` inherited whichever convention the source RDS used — so the single-sample fit at that site differed across methods for reasons unrelated to the classification itself.

**Resolution (applied):** in both 6a and 6c, the `single_sample_ids` block was moved ahead of the per-type/per-quartile demotion block. The single-sample pool is now drawn from all peaks irrespective of flood type, before any per-type MIN_TYPE_SAMPLE filtering — matching 6b and Yu et al. (2022).
- 6a: `in_single_sample` assigned at line 177; `flood_type` overwritten to `"unclassified"` at line 203.
- 6c: `in_single_sample` assigned at line 134; `flood_type` overwritten to `"unclassified"` at line 160.

The inline comments at 6a:163 and 6c:120 (which previously justified the demote-first ordering) have been replaced with comments stating that single-sample selection precedes demotion, following Yu et al. (2022).

**Required pipeline rerun:** scripts 6a, 6c, and 7 must be re-executed to regenerate `pot_peaks_classified_davenport.rds`, `pot_peaks_classified_quartile.rds`, and the downstream `ffa_*` outputs with the new `in_single_sample` semantics.

### 3c. App's `ss_ids` — **RESOLVED**
Already covered in 2c. After the fix, the app reads `in_single_sample` directly from the RDS, inheriting whichever convention scripts 6a/6b/6c wrote — now uniformly the Yu et al. (2022) ordering.

## 4. Silent failure points

- **functions.R does not `library(lmomco)`.** ~~The script-7 fix (line 58) is in place, but the same risk exists for any other entry point that sources `functions.R` without first loading lmomco — `pargpa`, `lmoms`, `parpe3`, `vec2par`, `cdfpe3`, `quape3` would fail.~~ **RESOLVED:** `library(lmomco)` is now line 1 of `functions.R`. Any caller that sources the file (including future entry points beyond app.R and script 7) gets lmomco loaded automatically.
- **Script 7 lines 201, 211, 221:** `tryCatch(dist_curve / dist_mixture_curve, error = rep(NA, length(ri_grid)))`. A failed evaluation produces an all-NA row in `curves_df`. The row is still written to `ffa_curves_*.rds`, so a fully-failed curve is detectable downstream as an all-NA group, but a *partially* failed grid (some NAs) silently shows up as gaps in the plotted curve. Not visible in the script-7 console diagnostics.
- **Script 7 lines 344, 353:** `tryCatch(fit_dist, error = NULL)`. A NULL fit causes the type to be skipped silently. The diagnostic block (lines 415–420) only counts non-NULL `class_rows` — not non-NULL fits per type — so a site whose rainfall-fit fails consistently across families would show normal class_rows but missing curves. Not surfaced.
- **Script 7 line 297:** `ros_rows <- if (!is.null(ros_type) && "water_year" %in% names(site_peaks))`. If `water_year` were ever missing from a future RDS, ROS temporal output would silently be NULL for that method. Low risk given current scripts always add it, but worth noting.
- **Script 7 lines 246–247:** empty `site_peaks` returns all-NULL early. Detectable from the per-method count diagnostic.
- **functions.R line 27 (`fit_gpd`), line 119 (`fit_lp3`)**: `if (anyNA(pars$para) || ...) return(NULL)` — a degenerate L-moments fit returns NULL silently. Visible only as a missing curve.
- **functions.R lines 84–88 (`mixture_curve`), lines 175–179 (`lp3_mixture_curve`)**: `tryCatch(uniroot, error = NA_real_)`. Silent NA in the mixture curve grid. Same gap-in-plot symptom as above.

No additional missing-library risks beyond the `lmomco`-in-`functions.R` note.

## 5. Namespace collisions

- **`dplyr::count` and `dplyr::case_when`** are explicitly qualified at the per-method summary blocks (6a:138/219, 6b:213/225/302, 6c:95/177, script 7:269) but most other dplyr verbs are unqualified — `select`, `filter`, `mutate`, `group_by`, `slice_max`, `summarise`, `arrange`, `bind_rows`, `left_join`, `inner_join`, `if_else`, `replace_na`, `pull`, `distinct`, `rename`, `desc`, `row_number`, `n`, `n_distinct`.
- **`MASS::select`** masks `dplyr::select` whenever MASS is loaded after dplyr. None of the six files loads MASS directly, but the project's GMM analysis cited in 6b's comments uses mclust, and mclust transitively loads MASS in some configurations. Any future addition that brings MASS into the session would silently break every unqualified `select(...)` call — including `select(-any_of(c("in_single_sample", "in_type_sample", "water_year")))` at 6c:117, the regulated-column joins at 6a:120, 6b:176, script-7:86, and the `select(site_no, peak_id, ...)` chains at 6a:176/202/206, 6b:255/270/275, 6c:133/160/165, script-7:286/309.
- **`stats::filter`** is masked by `library(tidyverse)` immediately after load, so `filter(...)` is safe in current execution. If anything re-attaches stats above tidyverse, this protection lapses silently.
- **`purrr::map`** is qualified at all five usage sites in script 7 (lines 200, 374, 418, 433, 438, 443). Should be the standard everywhere.
- **app.R** uses `dplyr::filter` / `dplyr::mutate` / `dplyr::select` only in the download handler. Everywhere else (e.g., reactives at 410, 415, 444, 453, 468, 477, 484), they are unqualified. App load order does not bring in MASS, so no current breakage, but the inconsistency is a maintenance hazard.

**Recommendation:** if mclust/MASS will ever be loaded in the same session as these scripts (e.g., adding GMM diagnostics into the pipeline), every `dplyr` verb should be qualified. As a minimum, the `select(-any_of(...))` calls and any `select` immediately after `group_by` are worth qualifying preemptively.

---

## Summary of action items

1. **app.R single-sample annotation (area 2c / 3c).** ✅ Applied — `ss_ids` now reads `in_single_sample` from the RDS (app.R lines 558–560).
2. **Single-sample ordering (area 3b).** ✅ Applied — 6a and 6c moved single-sample selection ahead of demotion to match 6b and Yu et al. (2022). Pipeline rerun of 6a, 6c, and 7 is required to regenerate the affected RDS outputs.
3. **functions.R lmomco import (area 4).** ✅ Applied — `library(lmomco)` is now line 1 of `functions.R`.
4. **Namespace qualification (area 5).** ⏳ Not yet applied. Optional hardening — qualify `dplyr::select`/`filter` in scripts and app.R, especially the `select(-any_of(...))` patterns. Recommended before any addition that brings MASS or mclust into the session.

Areas 1 and 2 (per-type and 2-class fit-to-plot, including ROS-below-threshold) were clean and required no action.
