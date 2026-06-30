# FFA Pipeline Refactor Log

Tracks results and notes for each phase of the six-phase FFA classification framework refactor.
Goal: simplify from four methods (`davenport`, `adaptive`, `adaptive_2class`, `quartile`) down to three new methods.

---

## Phase 1 — Archive old artifacts and remove the quartile method

**Date completed:** 2026-05-18

### Files moved to `_data/_archive/`
- `pot_peaks_classified_davenport.rds`
- `pot_peaks_classified_adaptive.rds`
- `pot_peaks_classified_quartile.rds`
- `6b_classify_usgs_unit_data_peaks_adaptive_GMM_VERSION.R` (copy from `_scripts/`; GMM-derived classification logic preserved for future PhD chapter work)

### Files deleted
- `_scripts/6c_classify_usgs_unit_data_peaks_quartiles.R`

### Files edited

| File | What was removed |
|---|---|
| `_scripts/7_compile_ffa_results.R` | `quartile` header comment; `pot_peaks_classified_quartile.rds` INPUTS line; `peaks_qrt <- readRDS(...)` load; Quartile log message; entire `quartile = list(...)` block in `methods`; `unique(peaks_qrt$site_no)` from `all_sites`; `quartile = NULL` from `ros_type` switch |
| `app.R` | `quartile` comment block in file header; `6a-6c` → `6a-6b`; Q1–Q4 entries from `type_labels`; Q1–Q4 color entries from `type_colors`; `quartile` entry from `method_labels`; `"Quartiles" = "quartile"` from `classification_mode` selectInput; quartile choice + selected from `class_methods`, `ffa_methods`, and `ffa_active_method` inputs; `"quartile"` from both `method_label = factor(...)` level vectors; Q1–Q4 label pairs from `stack_order`; `"Adaptive (Quartile)"` entry from `legend_schemes` |

### Final grep verification

| Pattern | Matches outside `_data/_archive/` | Status |
|---|---|---|
| `quartile` | `_diagnostics/` only (elevation/area quartile stats — unrelated to FFA method) | Clean |
| `peaks_qrt` | None | Clean |
| `pot_peaks_classified_quartile` | None | Clean |
| `\bQ[1-4]\b` | `_diagnostics/` only (drainage area / elevation bin labels) | Clean |

### Unexpected findings
None. Scripts 8, 9, `ui.R`, and `functions.R` had no FFA quartile method references and were left untouched.

### Exit criteria
All met. Live `_data/` has no `pot_peaks_classified_*.rds` files. Archive contains all four expected files. Live `6a` and `6b` scripts are unchanged.

### Notes

---

## Phase 2 — Rename `davenport` → `depth_thresh`

**Date completed:** 2026-05-21

### File renamed

`_scripts/6a_classify_usgs_unit_data_peaks_davenport.R` → `_scripts/6a_classify_usgs_unit_data_peaks_depth_thresh.R`  
(via `git mv`; history preserved)

### Files edited

| File | What changed |
|---|---|
| `_scripts/6a_classify_usgs_unit_data_peaks_depth_thresh.R` | Header description: "Davenport Fixed Thresholds" → "Depth Threshold (Davenport et al. 2020 fixed-depth rules)"; OUTPUT docstring filename → `pot_peaks_classified_depth_thresh.rds`; `output_file` constant → `_data/pot_peaks_classified_depth_thresh.rds`. Citation references in messages (lines 94, 116) left unchanged. |
| `_scripts/7_compile_ffa_results.R` | Header method list comment; INPUTS docstring filename; `readRDS` path → `pot_peaks_classified_depth_thresh.rds`; log message "Davenport" → "Depth Threshold"; `methods` list key `davenport` → `depth_thresh`; `ros_type` switch entry `davenport` → `depth_thresh`. Local variable `peaks_dav` left unchanged. |
| `app.R` | Header docstring; `method_labels` entry (`davenport = "Fixed"` → `depth_thresh = "Depth Threshold"`); `classification_mode` selectInput choices; `class_methods` checkboxGroupInput choices + selected; `ffa_methods` checkboxGroupInput choices + selected; `ffa_active_method` selectInput choices; both `method_labels[c(...)]` factor level vectors; `legend_schemes` key `"Fixed"` → `"Depth Threshold"`. |

### Final grep verification

| Pattern | Matches outside `_archive/` and `_planning/` (excl. "et al") | Status |
|---|---|---|
| `davenport` (method key) | None in live pipeline scripts | Clean |
| `pot_peaks_classified_davenport` | None | Clean |
| `depth_thresh` | 17 occurrences across Script 6a, Script 7, and `app.R` — all correct | Clean |
| `ls _scripts/6a_*.R` | `6a_classify_usgs_unit_data_peaks_depth_thresh.R` only | Clean |
| `_data/pot_peaks_classified_*.rds` | None (Phase 4 will rerun scripts) | Clean |

### Unexpected findings

- **`/davenport-scratch/` filesystem paths** in Scripts 0, 1, 2, 3, 4, and two utility scripts. These are server/scratch-disk paths, not method keys — left untouched.
- **`_diagnostics/`** scripts use `davenport` as local variable names and citation references. Out-of-scope diagnostic scripts — left untouched.
- **`"Fixed"` in `_diagnostics/threshold_evaluation/gmm_rain_frac_swe/b2_threshold_comparison.R`** — refers to a "fixed `rain_frac` threshold" in a GMM comparison diagnostic; unrelated to the FFA method label. No change made.

### Exit criteria
All met. `6a_classify_usgs_unit_data_peaks_depth_thresh.R` exists; davenport version does not. Script 7 uses `depth_thresh` key with correct input path. `app.R` uses `depth_thresh` key and `"Depth Threshold"` label throughout. No live file references `davenport` as a method key. Live `_data/` has no `pot_peaks_classified_*.rds` files.

### Notes
Scripts 8, 9, `ui.R`, and `functions.R` had no `davenport` method-key references and were left untouched (confirmed by fresh grep).

---

## Phase 3 — Write new Script 6b and standardize `rain_on_snow` factor level

**Date completed:** 2026-05-21

### Task 1: New Script 6b created

`_scripts/6b_classify_usgs_unit_data_peaks_frac_thresh.R` — fractional-threshold 3-class classifier using `rain_frac` over the basin-adaptive `search_radius` window. No SWE. `MIXED_SPLIT` is computed at run time as the population median of mixed-event `rain_frac` and saved as an attribute on the output. Author line: Mike Talbot, PhD Student, Colorado State University (no author line was present in either the old live 6b or the archived GMM version).

Old live `_scripts/6b_classify_usgs_unit_data_peaks_adaptive.R` deleted (it was supposed to be removed in Phase 1 but was missed; removed now via `git rm`).

### Task 2: `"ROS"` → `"rain_on_snow"` standardization

| File | What changed |
|---|---|
| `_scripts/7_compile_ffa_results.R` | `NON_RAINFALL_TYPES`: `"ROS"` → `"rain_on_snow"`; `adaptive` method `types` vector: `"ROS"` → `"rain_on_snow"`; `ros_type` switch `adaptive` entry: `"ROS"` → `"rain_on_snow"`. |
| `app.R` | `type_labels`: removed duplicate `ROS = "Rain-on-Snow"` entry (kept `rain_on_snow = "Rain-on-Snow"`); both `non_rainfall_types` vectors in `site_peaks_labeled` and `site_emp_peaks` reactives: `"ROS"` → `"rain_on_snow"`. |
| `_scripts/9_analyze_flood_type_temporal_trends.R` | `FLOOD_TYPES` vector; `TYPE_COLORS` and `TYPE_LABELS` keys; `annual_ts` column `n_ROS` → `n_rain_on_snow`; `prop_ROS` → `prop_rain_on_snow`; `site_record_stats` column `n_ROS_yrs` → `n_rain_on_snow_yrs`; both `pivot_longer` `cols` vectors; three `flood_type == "ROS"` filters; `TYPE_COLORS["ROS"]` lookup. |

Script 6a was already clean from Phase 2 — no changes needed.

### Final grep verification

| Check | Result |
|---|---|
| `ls _scripts/6b_*.R` | `6b_classify_usgs_unit_data_peaks_frac_thresh.R` only — correct |
| `"ROS"` as string literal in `_scripts/` and `app.R` | Zero matches — clean |
| `"rain_on_snow"` occurrences | 26 occurrences across Scripts 6a, 6b, 7, 9, and `app.R` — all correct |
| `pot_peaks_classified_adaptive` in live files | Remains in Scripts 7, 8, 9 — expected; replaced in Phase 4 |
| `_data/pot_peaks_classified_*.rds` | None — clean |

### Unexpected findings

- **Old `6b_classify_usgs_unit_data_peaks_adaptive.R` still existed** in `_scripts/` — Phase 1 archived a GMM copy but did not delete the original. Removed in this phase via `git rm`.
- **No author line** in the old 6b (live or archive). Used the project-standard author line from Script 7.

### Exit criteria
All met. New `6b_classify_usgs_unit_data_peaks_frac_thresh.R` exists and matches spec. No live script uses `"ROS"` as a factor level string value. `app.R` `type_labels` has exactly one entry mapping to `"Rain-on-Snow"` (keyed to `rain_on_snow`). `_data/` has no `pot_peaks_classified_*.rds` files.

---

## Phase 4 — Integrate new methods into Script 7 and rerun

**Date completed:** 2026-05-21

### Task 1: Script 7 edits (five sections)

| Section | Old → New |
|---|---|
| Header method list | `adaptive (3-class), adaptive_2class` → `frac_thresh_3class, frac_thresh_2class` |
| INPUTS comment | `pot_peaks_classified_adaptive.rds` → `pot_peaks_classified_frac_thresh_3class.rds` |
| LOAD DATA | `peaks_adap <- readRDS("..._adaptive.rds")` + "Adaptive" log msg → `peaks_frac <- readRDS("..._frac_thresh_3class.rds")` + "Frac Thresh 3cls" msg |
| `methods` list | Replaced `adaptive` and `adaptive_2class` entries with `frac_thresh_3class` and `frac_thresh_2class`; `peaks_adap` → `peaks_frac`; removed "merged pool" comment on pooled_types |
| `all_sites` | `peaks_adap$site_no` → `peaks_frac$site_no` |
| `ros_type` switch | `adaptive = "rain_on_snow", adaptive_2class = NULL` → `frac_thresh_3class = "rain_on_snow", frac_thresh_2class = NULL` |
| Two internal comments | `adaptive_2class` → `frac_thresh_2class`; `ROS ->` → `rain_on_snow ->` |

### Task 2: Scripts run

Scripts 6a, 6b, and 7 all completed without errors. Key console output:

**Script 6a (tail):**
```
Sites with >= 20 rain_on_snow peaks  : 16
Saved 21211 classified peaks to _data/pot_peaks_classified_depth_thresh.rds
  MIN_TYPE_SAMPLE  : 20
```

**Script 6b (tail):**
```
MIXED_SPLIT       : 0.2979   rain_frac range: [0.0000, 0.9900]
Sites with >= 20 rainfall peaks      : 156
Sites with >= 20 snowmelt peaks      : 111
Sites with >= 20 rain-on-snow peaks  : 96
Saved 21211 classified peaks to _data/pot_peaks_classified_frac_thresh_3class.rds
  RAIN_FRAC_THRESH : 0.99  MIXED_SPLIT : 0.2979  MIN_TYPE_SAMPLE : 20
```

**Script 7 (tail):**
```
206 sites x 3 methods x 2 families = 1236 total fits
GPD curves    : 896808 rows | 204 sites | 6 curve types
GPD quantiles : 13904  rows | 8 return periods  |  GPD params : 1449 rows
LP3 curves    : 895776 rows | 204 sites | 6 curve types
LP3 quantiles : 13888  rows | 8 return periods  |  LP3 params : 1448 rows
ffa_classification_summary : 2043 rows | 206 sites | methods: depth_thresh, frac_thresh_2class, frac_thresh_3class
ffa_ros_temporal           : 13055 rows | water years 1986-2024
Script 7 complete.
```

### Output file sizes

| File | Size | Timestamp |
|---|---|---|
| `pot_peaks_classified_depth_thresh.rds` | 370 KB | 2026-05-21 10:35 |
| `pot_peaks_classified_frac_thresh_3class.rds` | 701 KB | 2026-05-21 10:35 |
| `ffa_curves_gpd.rds` | 5.1 MB | 2026-05-21 10:37 |
| `ffa_curves_lp3.rds` | 5.1 MB | 2026-05-21 10:37 |
| `ffa_quantiles_gpd.rds` | 89 KB | 2026-05-21 10:37 |
| `ffa_quantiles_lp3.rds` | 89 KB | 2026-05-21 10:37 |
| `ffa_params_gpd.rds` | 22 KB | 2026-05-21 10:37 |
| `ffa_params_lp3.rds` | 30 KB | 2026-05-21 10:37 |
| `ffa_classification_summary.rds` | 13 KB | 2026-05-21 10:37 |
| `ffa_ros_temporal.rds` | 14 KB | 2026-05-21 10:37 |

### Spot-checks

All six checks passed:

1. Methods in `ffa_curves_gpd`: `"depth_thresh"` `"frac_thresh_2class"` `"frac_thresh_3class"` — correct
2. No legacy method keys — OK
3. Curve types: `depth_thresh` and `frac_thresh_3class` each have mixture, rain_on_snow, rainfall, single_sample, snowmelt; `frac_thresh_2class` has mixture, non_rainfall, rainfall, single_sample — correct
4. No legacy strings in output columns — OK
5. `MIXED_SPLIT = 0.2979`, `RAIN_FRAC_THRESH = 0.99` — both present and numeric
6. Sites in `ffa_curves_gpd`: 204 (206 sites total; 2 lacked sufficient data for any fit)

### Verification grep
`grep -n "adaptive|quartile|davenport|peaks_adap|peaks_qrt" _scripts/7_compile_ffa_results.R` → zero matches

### Unexpected findings
None. Scripts 8, 9, and `app.R` were not touched.

### Exit criteria
All met.

---

## Phase 5 — Repoint Scripts 8 and 9 to the new classified file

**Date completed:** 2026-05-21

### Files edited

**`_scripts/8_analyze_climate_forcing_trends.R`** — 3 lines changed:

| Line | Old | New |
|---|---|---|
| 37 | `_data/pot_peaks_classified_adaptive.rds  (Script 6b)` | `_data/pot_peaks_classified_frac_thresh_3class.rds (Script 6b)` |
| 78 | `PEAKS_FILE <- "_data/pot_peaks_classified_adaptive.rds"` | `PEAKS_FILE <- "_data/pot_peaks_classified_frac_thresh_3class.rds"` |
| 151 | comment: `in pot_peaks_classified_adaptive.rds` | `in pot_peaks_classified_frac_thresh_3class.rds` |

**`_scripts/9_analyze_flood_type_temporal_trends.R`** — 6 sections changed:

| Lines | Change |
|---|---|
| 10–12 | Header description: "adaptive 3-class classification scheme" → "fractional threshold 3-class classification (frac_thresh_3class) from Script 6b" with rationale note |
| 62 | INPUTS comment: `pot_peaks_classified_adaptive.rds` → `pot_peaks_classified_frac_thresh_3class.rds` |
| 68–69 | OUTPUTS comment: both adaptive output filenames → `frac_thresh_3class` equivalents |
| 125 | `readRDS("_data/pot_peaks_classified_adaptive.rds")` → `readRDS("_data/pot_peaks_classified_frac_thresh_3class.rds")` |
| 561–567 | Both `saveRDS` calls + `message()` lines: `temporal_trends_adaptive.rds` and `annual_type_counts_adaptive.rds` → `frac_thresh_3class` equivalents |
| 811–812 | Final summary message block: both output filenames updated |

### Task 3: Archive verification

`_data/_archive/pot_peaks_classified_adaptive.rds` — present (archived in Phase 1). Not in live `_data/` — correct.

### Final grep verification

| Check | Result |
|---|---|
| `pot_peaks_classified_adaptive` in live scripts | Zero matches — clean (line 151 comment in Script 8 also updated) |
| `frac_thresh_3class` in Script 8 | Lines 37 and 78 — correct |
| `frac_thresh_3class` in Script 9 | Lines 11, 62, 68, 69, 125, 563–564, 566–567, 811–812 — all correct |
| `"adaptive"` as string literal in `_scripts/` | Zero matches — clean |
| `"adaptive"` in `app.R` | Remains (UI method key, display choices, reactives) — expected; Phase 6 scope |
| `adaptive.rds` in Script 9 | Zero matches — clean |

### Unexpected findings

- Script 8 had a third `pot_peaks_classified_adaptive.rds` reference at line 151 (an ASSUMPTION comment), not mentioned in the spec. Updated along with the two specified changes.
- `app.R` still uses `"adaptive"` as a method key throughout the UI — expected and correct; that is Phase 6's scope.

### Exit criteria
All met. Scripts 6a, 6b, 7, and `app.R` were not touched.

---

## Phase 6 — Update `app.R`: final refactor phase

**Date completed:** 2026-05-21

### Changes made to `app.R`

| Task | Lines (approx.) | Change |
|---|---|---|
| 1 — Header | 9–18 | Classification methods block: replaced `adaptive`/`adaptive_2class` descriptions with `frac_thresh_3class`/`frac_thresh_2class` with full rationale text |
| 2 — `method_labels` | 112–117 | Replaced `adaptive = "Adaptive (3-class)"` and `adaptive_2class = "Adaptive (2-class)"` with `frac_thresh_3class = "Fractional Threshold (3-class)"` and `frac_thresh_2class = "Fractional Threshold (2-class)"` |
| 3 — Map tab `selectInput` | 194–205 | Replaced 2-choice + `conditionalPanel` (adaptive + granularity radio) with 3-choice `selectInput`; `selected = "depth_thresh"` |
| 4 — `class_methods` | 251–259 | Replaced adaptive choices/selected with frac_thresh equivalents |
| 5 — `ffa_methods` + `ffa_active_method` | 284–311 | Replaced adaptive choices/selected in both inputs; `ffa_active_method` default → `"depth_thresh"` |
| 6 — `site_peaks_labeled` | 407–423 | Simplified from 16-line conditional remap to 3-line passthrough — `frac_thresh_2class` has `flood_type` pre-remapped by Script 7 |
| 7 — `active_method` reactive | 425–431 | Simplified from granularity conditional to `input$classification_mode` directly |
| 7 — Map observer `method_key` | 690–702 | Replaced 4-line granularity conditional with `method_key <- input$classification_mode`; updated comment |
| 8 — `site_emp_peaks` | 433–439 | Updated comment and `is_2class` condition: granularity test → `input$classification_mode == "frac_thresh_2class"` |
| 9 & 10 — Factor levels | 848–849, 872–873 | Both `method_labels[c("depth_thresh", "adaptive", "adaptive_2class")]` → `method_labels[c("depth_thresh", "frac_thresh_3class", "frac_thresh_2class")]` (replace_all) |
| 11 — `legend_schemes` | 926–930 | `"Adaptive (3-class)"` and `"Adaptive (2-class)"` keys → `"Fractional Threshold (3-class)"` and `"Fractional Threshold (2-class)"` |
| 12 — Two fallback defaults | 1070, 1134 | Both `"adaptive"` fallbacks → `"depth_thresh"` (second required a separate edit — different whitespace) |

### Verification greps

| Check | Result |
|---|---|
| `"adaptive"` as method key | Zero matches — clean |
| `"adaptive_2class"` | Zero matches — clean |
| `adaptive_granularity` | Zero matches — clean |
| `frac_thresh_3class\|frac_thresh_2class` | 20 occurrences — header, `method_labels`, all UI inputs, reactives, factor levels — all correct |

### Source check
`Rscript -e "source('app.R')"` — no errors. One unrelated bslib deprecation warning: `page_navbar()` `underline` argument consolidation (pre-existing, not introduced by this refactor).

### Unexpected findings
- The map symbology observer (lines 690–702) had its own copy of the `active_method()` conditional logic that the spec didn't explicitly list — updated along with Task 7 as it mirrors the same pattern.
- The second `"adaptive"` fallback default had double whitespace (`"adaptive"` instead of `"adaptive"`) so `replace_all` on the single-space form left it unmatched. Caught by the verification grep and fixed with a targeted edit.

### Exit criteria
All met. `app.R` sources without errors. No `"adaptive"`, `"adaptive_2class"`, or `"adaptive_granularity"` strings remain. All three UI tabs offer the correct three methods with correct display labels. `method_labels`, factor levels, and `legend_schemes` are internally consistent. Pipeline scripts are untouched.
