# Shiny App Reorganization Plan

> **Current tab order (as of 2026-06-22):** About → FFA → Classification → Methods → Trends
> - "FFA" was formerly "Map"
> - "Methods" was formerly "FFA Comparison"

## Step 1: Resolve the antecedent window question before anything else

This is a methodological prerequisite because the answer affects both the pipeline outputs
and what the app displays. The question has two parts: (a) should `frac_thresh` use a
fixed 8-day window instead of the adaptive Livneh et al. `search_radius`? And (b) should
`depth_thresh` adopt the Livneh et al. window instead of its current fixed 8-day window?

**Decision: keep methods as-is.** `depth_thresh` retains the fixed 8-day window faithful
to Davenport et al. (2020). `frac_thresh` retains the adaptive `search_radius` window,
internally consistent with the Livneh et al. peak independence framework. The download
schema discrepancy between the two classified output files is a naming/documentation fix
to be addressed in Step 10.

### Key findings

**Search radius distribution (across 206 sites):**
- The Livneh et al. equation `5 + ln(1.6092 × A)` produces a median `search_radius` of
  12 days for the Colorado GAGES-II network, ranging from 7 to 17 days. Only ~2% of
  sites have a `search_radius` at or below 8 days. Forcing a fixed 8-day window for
  `frac_thresh` would mean using a systematically too-short antecedent period for nearly
  every site in the network.

**rain_frac sensitivity to window length:**
- Substituting a fixed 8-day window (W1) for the adaptive window (W2) changes the
  classification of 6.7% of peaks overall. The change rate is highly non-uniform: rainfall
  is essentially unaffected (<1%), snowmelt is moderately sensitive (3–7%), and
  rain_on_snow is most sensitive (15–19%). The aggregate `rain_frac` distribution is
  nearly identical between W1 and W2 (the density curves are nearly indistinguishable),
  but this masks real individual-peak disagreement that has downstream FFA consequences.

**Operational impact — mixture-FFA site viability:**
- This is the smoking gun. Switching to a fixed 8-day window for `frac_thresh` would
  cause 27 of 96 currently viable rain_on_snow sites (28%) to fall below the
  `MIN_TYPE_SAMPLE = 20` threshold and lose mixture-FFA viability. Rainfall loses 1 site
  (1%) and snowmelt loses 5 sites (5%). Rain_on_snow is already the type with the
  smallest viable site count, and losing 28% of those sites would directly reduce the
  spatial coverage of the mixture FFA in the dashboard.

**Mechanism — why rain_on_snow is most affected:**
- The reclassification to unclassified under W1 is driven almost entirely by NA
  `rain_frac` (zero total active flux above the 0.1 mm/day noise floor in just 8 days),
  not by data insufficiency. The seasonal pattern is concentrated in January and December
  — cold-season events where neither rain nor melt is active within the shortened window.
  The cumulative flux profiles confirm the physical mechanism: the melt signal for
  snowmelt and rain_on_snow events is still actively accumulating past day 8 and has not
  plateaued by the W1 cutoff. Rainfall events, by contrast, are episodic and front-loaded,
  which is why an 8-day window is adequate for that type.

**Why the Livneh window cannot simply be applied to depth_thresh instead:**
- The Davenport et al. (2020) depth thresholds (rain > 10 mm, melt > 10 mm over the
  antecedent window) are absolute depths, not ratios. They were derived and validated
  specifically at an 8-day window. Extending the window to 12 days accumulates more
  total rain and melt by construction, so applying the same thresholds to a longer window
  mechanically inflates rainfall and snowmelt classification rates — not because the
  physical signal is stronger, but because more days are being summed. A proportional
  threshold scaling approach (e.g., `thresh_scaled = thresh_8d × window / 8`) was
  considered as an exploratory option, but has no external validation and would no longer
  be attributable to Davenport et al. The diagnostic script for this question was scoped
  but not run, as the methodological argument against it is sufficient without empirical
  confirmation.

### Diagnostic scripts

- `_diagnostics/antecedent_window/diag_antecedent_window.R` — Part 1: characterizes the
  `search_radius` distribution, computes `rain_frac` under W1 (fixed 8d), W2 (adaptive),
  and W3 (fixed 12d, rounded median), reports classification change fraction (6.7%
  overall), and documents the column schema discrepancy between the two classified output
  files. Outputs: `fig_search_radius_distribution.png`, `fig_search_radius_vs_drain_area.png`,
  `fig_rainfrac_W1_vs_W2.png`, `fig_rainfrac_W1_vs_W3.png`,
  `fig_rainfrac_density_W1_vs_W2.png`, `fig_rainfrac_sensitivity_vs_window.png`,
  `summary.txt`.

- `_diagnostics/antecedent_window/diag_antecedent_window_extended.R` — Part 2: assesses
  the operational impact of W1 substitution for `frac_thresh`. Reports NA `rain_frac`
  rates by flood type and month, reclassification-to-unclassified counts by season and
  reason (data insufficiency vs. zero active flux), per-site-type mixture-FFA viability
  loss under W1, and cumulative flux profiles by flood type showing whether the
  hydrometeorological signal is still accumulating past day 8. Appends results to
  `summary.txt`. Outputs: `fig_na_rainfrac_rate.png`,
  `fig_reclassified_unclassified_by_season.png`, `fig_site_type_viability.png`,
  `fig_cumulative_flux_profiles.png`.

- `_diagnostics/antecedent_window/diag_depth_thresh_window_sensitivity.R` — Part 3
  (scoped, not run): would have quantified how `depth_thresh` classification changes when
  the Davenport thresholds are applied at the adaptive `search_radius` window rather than
  the fixed 8-day window, under both unscaled and proportionally scaled threshold regimes.
  Not run because the methodological argument against the Livneh window for `depth_thresh`
  is sufficient on analytical grounds: the Davenport thresholds are absolute depths tied
  to an 8-day window by design, and applying them at a longer window without recalibration
  is not a neutral change.

---

## Step 2: Structural reorganization — split Classification from FFA ✓

The current structure conflates classification results with FFA results. These need to be separated into two distinct tabs, each with its own sub-structure:

- A **Classification tab** focused on how peaks are classified at each site — including the fraction of peaks in each class, spatial patterns in dominant flood type, and the ability to compare 2-class vs. 3-class methods.
- An **FFA tab** focused on flood frequency results — quantile estimates, spatial patterns in return-period flood magnitudes, and comparisons between mixture and single-sample approaches.

The FFA Comparison → Spatial Patterns content should be merged into the Maps tab (add a return period dropdown there), freeing the FFA Comparison tab to focus on regional-level summaries and method comparisons.

### Changes made (2026-06-22)

- Removed the "Classification Composition" sub-tab from the FFA Comparison tab.
- Added a new top-level **Classification** tab (positioned between FFA Comparison and Trends) with its own sidebar containing independent `class_methods` and `class_reg_filter` inputs.
- The Classification Composition plot and legend (`summary_composition_plot`, `summary_composition_legend`) now live in the Classification tab; `summary_class_data()` was updated to read from the new inputs.
- Removed the "Classification Composition" branch from the FFA sub-tab input gating observer.
- FFA Comparison retains its three remaining sub-tabs (Mixture vs. Single-Sample, Method Intercomparison, Spatial Patterns) and its existing sidebar inputs unchanged.
- Tab order: About → Map → FFA Comparison → **Classification** → Trends.

---

## Step 3: Reorganize the Maps tab ✓

Remove the "100-year dominant flood" layer — this is a derived FFA result and belongs in the FFA tab, not the Maps tab. Consider whether Reference Class and Regulation Status layers are better placed as static maps on an About page rather than interactive layers here. Add a return period dropdown (absorbing the Spatial Patterns content from FFA Comparison). Add hover pop-ups to the map plots. Stack the FFA and Time Series plots vertically on this tab.

### Changes made (2026-06-22)

- Removed "Reference class" and "Regulation" from the `map_symbology` radio buttons (earmarked for the About page — see Step 11).
- Renamed "100-yr dominant flood" → "Dominant flood type" in the radio buttons; removed the hardcoded return period of 100.
- Added a conditional `map_return_period` dropdown to the Maps tab sidebar (visible only when "Dominant flood type" is selected); the map symbology observer now uses this value instead of the hardcoded 100.
- Removed the `site_regulation` join and `ref_label`/`reg_label` columns from the map symbology observer; updated popup HTML accordingly.
- Added `label = ~hover_label` to `addCircleMarkers` on the station map for hover tooltips (site ID + name on mouse-over).
- Removed the "Spatial Patterns" sub-tab from FFA Comparison; removed the `ffa_dominant_data` reactive and `ffa_dominant_map` output that exclusively served it.
- Removed the "Spatial Patterns" branch from the FFA sub-tab input gating observer.
- Merged the "Analysis" and "Data" pill panels in the site analysis container: FFA plot and timeline are now stacked vertically in the "Analysis" panel, with the download button below. "Parameters" remains a separate pill.
- Increased FFA plot height to 750px and timeline plot height to 600px (~50% increase from original).
- Made the analysis card independently scrollable: both cards in the Map tab layout are fixed at 660px via explicit `style = "height: 660px;"` (replacing `heights_equal = "row"`); the analysis card body has `overflow-y: auto` and `min-height: 0` (the latter is required to override the default flex `min-height: auto` that would otherwise allow the card body to grow past the card boundary). Switched `navset_card_pill` to `navset_pill` in `analysis_container` to eliminate the nested card element that was also blocking scroll.

---

## Step 4: Build out the Classification tab ✓

Organize into two sub-tabs. General rule: if a plot does not change based on user selections or site clicks, it belongs in the technical report rather than the app.

**Sub-tab 1 — "Classification Map":** A Leaflet map where sites are colored by either (a) dominant flood type under a selected method and return period (uses `ffa_quantiles`), or (b) the fraction of peaks in a selected type under a selected method (uses `ffa_classification_summary`). Driven by a radio button ("Dominant flood type" | "Peak classification") with conditional dropdowns. Clicking a site reveals a site-level detail panel showing classification composition (fraction of peaks per type) across all methods for that site — this replaces the former regional `summary_composition_plot`.

**Sub-tab 2 — "Method Comparison":** A Leaflet map showing classification-level agreement between any two user-selected methods. Agreement is computed as the fraction of peaks assigned the same type by both methods. Precomputed at startup for all three method pairs. This is purely a classification comparison — not FFA magnitudes.

**Implementation notes:**
- Remove the regional `summary_composition_plot` and `summary_composition_legend` outputs entirely; replace with site-level `class_site_composition_plot`.
- Remove the old `class_methods` / `class_reg_filter` inputs from the sidebar (replaced by sub-tab–specific controls).
- Precompute `classification_agreement` at startup by joining peak-level classified files for all three method pairs; store per-site agreement fraction and N peaks.
- For frac_thresh_2class comparisons, collapse snowmelt and rain_on_snow to non_rainfall before comparing.
- Update `class_type` dropdown choices dynamically based on selected method (2-class has only rainfall / non_rainfall).
- Prevent Method A == Method B in Sub-tab 2 by dynamically excluding the selected Method A from Method B choices.

### Changes made (2026-06-22)

- Replaced the old single-card Classification tab with a `navset_card_tab` containing two sub-tabs.
- **Sub-tab 1 ("Classification Map"):** Leaflet map with sidebar radio ("Dominant flood type" | "Peak classification"); conditional dropdowns for RP + distribution (dominant mode) or flood type (proportion mode); method selector; regulation filter. Map markers update via `leafletProxy`. Site detail panel shows a stacked bar chart of classification composition across all methods for the selected site (`class_site_composition_plot`), replacing the former regional `summary_composition_plot`. Site selection is bidirectional (dropdown ↔ map click).
- **Sub-tab 2 ("Method Comparison"):** Leaflet map colored by per-site agreement fraction between two user-selected methods. Method B choices dynamically exclude Method A. `class_type` choices update dynamically when method changes (2-class has only rainfall / non_rainfall).
- Added startup precomputation of `classification_agreement` (using `local()`) — loads depth_thresh and frac_thresh_3class peak files, derives frac_thresh_2class by remapping, joins all three pairs on `site_no + peak_date`, computes per-site agreement fraction and N peaks.
- Removed `summary_class_data`, `summary_composition_plot`, `summary_composition_legend`, and the old `class_methods` / `class_reg_filter` sidebar inputs entirely.

---

## Step 5: Build out the FFA Comparison tab ✓

Replace or extend the Method Intercomparison sub-tab with a spatial Leaflet map showing the quantile ratio between any two user-selected methods at a selected return period, with point size scaled to ratio magnitude. The user should be able to select any two methods as the comparison pair (not fixed to Depth as a reference axis). This follows the Yu et al. (2022) framing — sites where mixture and single-sample estimates diverge most should be visually prominent.

### Changes made (2026-06-22)

- Renamed "Active threshold method" sidebar label → "Method A (intercomparison)"; added `ffa_compare_method_b` ("Method B") selector with dynamic exclusion of Method A (same pattern as Classification tab).
- Replaced the "Method Intercomparison" sub-tab content: removed `ffa_pairwise_plot` scatter and `ffa_pairwise_tooltip`; added `ffa_intercomp_map` Leaflet map.
- `ffa_intercomp_data` reactive: joins mixture Q for the two selected methods at the selected RP, computes ratio = Q_A / Q_B, log_ratio, and point radius scaled to `4 + 8 * min(|log_ratio|, log(5)) / log(5)`.
- Map uses a reversed RdBu color scale on log_ratio (red = Method A larger, blue = Method B larger); legend labels show ratio values (exp of log_ratio).
- Updated FFA sub-tab gating observer: "Mixture vs. Single-Sample" disables Method A/B and enables `ffa_methods`; "Method Intercomparison" disables `ffa_methods` and enables Method A/B.
- Removed `ffa_pairwise_data` reactive and `ffa_pairwise_plot` / `ffa_pairwise_tooltip` outputs entirely.

---

## Step 6: Revise the Trends tab ✓

Shift the analysis from POT peaks to annual maximum peaks only. Implement two complementary views:

- Rain fraction as a continuous trend over time (Mann-Kendall on the annual series).
- Flood type frequency as a shift comparison between the first and last half (or thirds) of the record — essentially, what caused the largest event each year and is the distribution of flood types changing?

Pull the Kampf (2015) paper Frances sent and check whether it informs the framing or methodology for the frequency-shift analysis.

### Changes made (2026-06-22)

- Scripts 8 and 9 updated to produce AM series outputs before this step was implemented:
  - Script 8: added `am_rain_frac` as a sixth variable to `VAR_KEYS` (orange, `#f28e2b`); derives the AM peak per site-year via `slice_max(peak_flow, n=1)` on classified peaks; joins `am_rain_frac` and `am_flood_type` onto `annual_climate_forcing_summary.rds`; runs Mann-Kendall / Sen slope for `am_rain_frac` through the existing trend pipeline; masking rule matches `total_precip_mm` (NA when no POT peaks that year).
  - Script 9: added a `frac_thresh_3class`-only AM section (inside the method loop) that produces `annual_am_type_frac_thresh_3class.rds` (one row per site-year), `am_type_shift_frac_thresh_3class.rds` (first- vs. second-half AM type frequency per site × type), and `temporal_trends_am_frac_thresh_3class.rds` (named list with `am_rain_frac_trends` MK results and `am_type_shift`).
- Replaced the Trends tab UI: old `navset_tab(Climate Forcing, Flood Frequency)` replaced by `layout_sidebar(shared sidebar, navset_card_tab(Rain Fraction Trend, Flood Type Shift))`. Shared sidebar has a single `trends_reg_filter` checkbox driving both sub-tabs.
- Added five new startup data loads: `am_annual`, `am_trends_obj`, `am_type_shift`, `am_rf_trends`, `clim_forcing`.
- **Sub-tab 1 — "Rain Fraction Trend":** Leaflet map colors sites by AM rain fraction MK trend direction and significance (sig. positive red, sig. negative blue, not significant grey, insufficient record light grey); 4-category `addLegend`. Site detail panel shows a scatter of `am_rain_frac` vs. `water_year`, points colored by `am_flood_type` via `type_colors`, with Theil-Sen segment overlaid and tau/slope/p/n in subtitle.
- **Sub-tab 2 — "Flood Type Shift":** Leaflet map colors each site by the dominant AM flood type in the second half of the record (`slice_max(n_second_half)`), using `type_colors`. Site detail panel shows a dodged bar chart (`geom_col`, `position_dodge`) of first-half vs. second-half year counts per flood type, with split-year annotation.
- Both maps share `trends_selected_id` reactiveVal (already defined for the former Climate Forcing sub-tab). Map click observers update it; detail panels react to it.
- All existing server-side code for the former Climate Forcing and Flood Frequency sub-tabs is retained in place (dormant but not removed).

### Revised (2026-06-22)

**Script 9 — binary AM type indicator trend tests (Step 3b):**
- Added Step 3b inside the `frac_thresh_3class` block (after half-period split, before the combined save): runs Mann-Kendall + Sen slope on a 0/1 annual indicator series for each eligible site × flood type (`rainfall`, `snowmelt`, `rain_on_snow`). NA `am_flood_type` is treated as 0 (not a type occurrence; no per-type minimum count applied — the series is well-defined even if the type never generated the AM peak).
- Output tibble columns: `site_no, flood_type, regulated, n_years, n_type_years, tau, p_value, sen_slope, significant, sig_positive, sig_negative`. Direction flags pre-computed in Script 9 to avoid downstream inference in the app.
- Updated `saveRDS` to include `am_type_indicator_trends` as a third list element alongside `am_rain_frac_trends` and `am_type_shift`.

**App data loads:**
- Renamed `am_trends_obj` → `am_trends`; added `am_type_indicator_trends <- am_trends[["am_type_indicator_trends"]]`; removed `clim_forcing` (not used by the Trends tab). Previous bullet "Added five new startup data loads" is superseded.

**Flood Type Shift map — revised to indicator-based directional coloring:**
- `trends_shift_map_data` reactive now driven by `am_type_indicator_trends` (not `am_type_shift`). Per-site dominant shift: filter to `significant == TRUE`, pick highest `|tau|` via `slice_max`. Eight categories: 3 flood types × 2 directions (increasing/decreasing) + "No significant shift" (#888888) + "Insufficient record" (#dddddd). Left-join from all stations so absent sites always appear on the map.
- `fill_color` and `shift_label` computed via `case_when` before passing to Leaflet — no named vector subsetting in `~formula` expressions (known failure mode in this project).
- Popup shows τ/p/n for significant sites; plain label only for no-data and insufficient-record sites.
- Map legend updated to static 8-category color/label pairs.
- Map click observer updated to also call `updateSelectizeInput` to sync the sidebar site dropdown on marker click.
- `req(!is.null(am_type_indicator_trends))` guard added to reactive and to `trends_shift_indicator_plot` to prevent NULL crash while Script 9 has not yet been re-run.

**Flood Type Shift detail panel — replaced single plot with two stacked plots:**
- `renderUI` now returns `tagList(plotOutput("trends_shift_indicator_plot", height="280px"), plotOutput("trends_shift_halfperiod_plot", height="280px"))`.
- `trends_shift_indicator_plot` (upper): bar chart of `frac_am = n_type_years / n_years`, colored by `FLOOD_TYPE_COLORS`, with τ/p annotation above each bar (* appended for significant results). Y-axis formatted as percent.
- `trends_shift_halfperiod_plot` (lower): half-period grouped bar chart (first half vs. second half year counts per flood type), sourced from `am_type_shift`, `base_size = 11`. Replaces the former `trends_shift_site_plot` output.

**Bug fix:** The initial Flood Type Shift map rendered all markers black due to named vector subsetting inside `dplyr::mutate()` returning named vectors that corrupted Leaflet `~formula` evaluation in dplyr 1.2.0. The full reactive replacement above resolved this structurally.

**Other:** Methods tab icon changed from `fa-water` to `fa-magnifying-glass-chart`.

---

## Step 7: Decide what belongs in the app vs. the technical report

Audit the existing plots with the rule: anything that isn't spatially dynamic (i.e., doesn't change based on user selections) is probably better suited to the technical report. Candidates for removal from the app include the regional-level bar charts (classification summary) and the scatterplots. Removing these will simplify the UI considerably and let the app focus on what it does better than a static report — interactive, site-level exploration.

---

## Step 8: Terminology and labeling cleanup ✓

Update "Threshold method" → "Classification method" throughout. Update method labels to "Depth threshold" and "Rain fraction (2-class)" / "Rain fraction (3-class)". Confirm `"rain_on_snow"` factor level is displayed consistently in all UI-facing labels.

### Changes made (2026-06-22)

- Replaced the two `"Threshold method"` / `"Threshold methods"` UI labels with `"Classification method"` / `"Classification methods"` (lines 366 and 579).
- Replaced all eight UI-facing `"Depth"` method labels with `"Depth threshold"`: `method_labels` lookup vector (line 187), the `classification_mode` selectInput (line 368), `class_method` selectInput (line 440), `class_compare_method_a` selectInput (line 508), `ffa_active_method` selectInput (line 552), `ffa_methods` checkboxGroupInput (line 581), and the two `all_methods` vectors in the FFA and Classification comparison observers (lines 1188 and 1253). Internal value `"depth_thresh"` is unchanged throughout.
- Audited all UI-facing rain-on-snow display labels (`type_labels`, `FLOOD_TYPE_LABELS`, `type_colors`, `FLOOD_TYPE_COLORS`, selectInput choices, `stack_order`, `case_when` shift labels, `scale_fill_manual` value names). All are already consistently `"Rain-on-Snow"`. No changes required.

---

## Step 9: UI polish — tooltips and hover behavior ✓

Add tooltips to help users interpret classification methods, FFA outputs, and trend statistics. Add hover pop-ups on map plots (flagged in Step 3 but implement as a unified pass across all tabs at this stage, after the structural reorganization is stable).

### Changes made (2026-06-22)

- Added `library(bsicons)` (line 51); `bslib` was already loaded. All tooltips use `bslib::tooltip(bsicons::bs_icon("info-circle"), ...)` as the trigger pattern.
- Defined `TT_REGULATION` string constant (before the UI definition) for shared regulation status tooltip text; referenced at all four `checkboxGroupInput` label locations rather than hardcoding the text four times.
- **FFA tab — `classification_mode`** (line 371): label wrapped in `tags$span` with tooltip explaining depth-threshold vs. rain-fraction classification approaches.
- **FFA tab — `map_symbology`** (line 389): converted `choices` to `choiceNames`/`choiceValues`; tooltip on "Data sufficiency" explains per-type and total-peak minimums; tooltip on "Dominant flood type" explains return-period quantile coloring.
- **Classification Map — `class_map_mode`** (line 455): converted `choices` to `choiceNames`/`choiceValues`; tooltip on "Dominant flood type" explains most-frequent-type coloring; tooltip on "Peak classification" explains peak-fraction coloring.
- **Classification Map — `class_map_reg_filter`** (line 507): label uses `TT_REGULATION`.
- **Classification → Method Comparison — `class_compare_reg_filter`** (line 567): label uses `TT_REGULATION`.
- **Trends — `trends_reg_filter`** (line 674): label uses `TT_REGULATION`.

Nine `bslib::tooltip` instances total; `grep -n "bslib::tooltip" app.R` returns exactly nine lines.

- **Methods (FFA Comparison) — `ffa_reg_filter`** (line 634): label uses `TT_REGULATION`. Added in correction pass after initial implementation.

---

## Step 10: Data download

Standardize headers across all output files (resolve the antecedent window naming inconsistency decided in Step 1). Add an option to download all data at once rather than requiring three separate downloads.

---

## Step 11: Build out the About tab

Content earmarked here as it is removed from other tabs during reorganization:

- **Reference class map** — static choropleth or dot map showing GAGES-II reference vs. non-reference classification for all sites. Removed from Maps tab `map_symbology` in Step 3.
- **Regulation status map** — static dot map showing regulated vs. unregulated sites. Removed from Maps tab `map_symbology` in Step 3.
- Colorado Water Center logo (`_includes/Colorado-Water-Center.png`).
- Project description, data sources, and methodology overview.

---

## Step 12: Catalog figures for the Technical Report and About page

This step does not involve app changes. It produces the static figures needed to explain the methodology and results to Colorado Water Center staff and watershed managers who may use the app. The audience is professional but not academic — figures should answer "what does this tool show, which sites have it, and what does Colorado look like overall?" rather than justifying methodological choices. All figures should be generated as standalone R scripts in `_diagnostics/technical_report/`, sourcing `functions.R` and loading existing `.rds` outputs — no recomputation of intermediate results. Use `ggplot2` with the CSU brand palette (primary `#1E4D2B`, secondary `#C8C372`). Target 300 dpi PNG or PDF. Script names should match their figure output (e.g., `fig_network_classification_composition.R`).

---

### 12.1 About page static maps

Both maps belong in the About page (Step 11), not the body of the Technical Report. They were removed from the Maps tab `map_symbology` radio buttons in Step 3 and earmarked here. Step 12 is the implementation point.

- **Reference class map** — static dot map showing GAGES-II reference vs. non-reference classification for all 206 sites. Orients the reader to which sites meet the reference standard for unimpaired hydrology.
- **Regulation status map** — static dot map showing regulated (143) vs. unregulated (63) sites. Provides spatial context for why the unregulated subset is the primary analysis network.

---

### 12.2 Network classification composition

The regional `summary_composition_plot` was removed from the app in Step 4 and replaced by the site-level detail panel. A network-level version belongs in the Technical Report as a summary of what the classification actually produces across Colorado.

- **Network classification composition** — grouped or faceted bar chart showing the fraction of peaks assigned to each flood type (rainfall / snowmelt / rain_on_snow) for the `frac_thresh_3class` method, faceted by regulation status (unregulated / regulated). Answers "what types of floods drive peak flows in Colorado, and does that differ between regulated and unregulated basins?" Source: `pot_peaks_classified_frac_thresh_3class.rds`.

---

### 12.3 Spatial FFA coverage map

This figure has no existing output and requires a new standalone script. It is the clearest single-figure answer to "which sites have full mixture FFA results and which don't?"

- **Spatial FFA coverage map** — dot map showing, for each site, the level of mixture-FFA support: single-sample only / partial mixture (one or two flood types viable) / full mixture (all three types viable), where viability is determined by the `MIN_TYPE_SAMPLE = 20` threshold. Source: `ffa_quantiles_{gpd|lp3}.rds`, filtered by which type-specific quantiles are populated per site.

---

### 12.4 Network-level trend summaries

The app Trends tab is site-driven; a network-level summary belongs in the Technical Report to give watershed managers a Colorado-wide picture of how flood-generating conditions are shifting.

- **Network MK trend summary** — for each of the six Script 8 variables (`total_precip_mm`, `max_swe_mm`, `melt_season_length`, `peak_swe_mm`, `total_melt_mm`, `am_rain_frac`), show the count and fraction of sites with significant positive, significant negative, and non-significant trends, for unregulated sites. Suggested form: grouped bar chart or dot plot with significance indicated by fill. Source: MK outputs from Script 8.
- **Flood type shift summary** — grouped bar chart of first-half vs. second-half annual maximum type frequency, aggregated across unregulated sites. The network-level complement to the app's Flood Type Shift site detail panel — answers "across Colorado, is the mix of flood-generating processes changing?" Source: `am_type_shift_frac_thresh_3class.rds` from Script 9.