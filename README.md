# CoWC FY26 Enhanced Impacts Support

Analysis pipeline for Colorado flood type classification and trend analysis. The pipeline classifies peaks-over-threshold (POT) streamflow events at Colorado GAGES-II sites by dominant flood-generating process (rainfall, snowmelt, rain-on-snow), runs mixture-model flood frequency analysis (FFA), and tests for temporal trends in flood type composition and climate forcing variables.

## Repository structure

```
<<<<<<< HEAD
_scripts/       Pipeline scripts (0–9) and supporting utilities
_data/          Intermediate and final data products (mostly .rds and .parquet)
_diagnostics/   Diagnostic figures organized by analysis stage
functions.R     Shared GPD/LP3 fitting and curve evaluation functions
renv/           Package environment (renv lockfile)
```

## Pipeline

Scripts are numbered in execution order. Each script reads from `_data/` and writes its outputs back to `_data/`. Run from the repository root so relative paths resolve correctly.
=======
_scripts/         Pipeline scripts (0–9) and supporting utilities
_data/            Intermediate and final data products (.rds, one .csv). WLDAS
                   raw/processed files, including Script 4's .parquet output,
                   live under _data/WLDAS/, a symlink to external storage —
                   not part of this git repository.
_documentation/   report_figures.R and its rendered output (_documentation/report_figures/),
                   which build the figures used in the project's technical report
functions.R       Shared GPD/LP3 fitting, curve evaluation, and trend-test
                   (Mann-Kendall/Theil-Sen, Pettitt) helper functions
renv/             Package environment (renv lockfile)
```

Scripts 8 and 9 also write diagnostic figures to a local `_diagnostics/` directory (paths given below). That directory is generated on demand and is not tracked in this repository — the project's exploratory/diagnostic work now lives in the companion [`CoWC_FY26_Enhanced_Impacts_Diagnostics`](https://github.com/realmiketalbot/CoWC_FY26_Enhanced_Impacts_Diagnostics) repository.

## Pipeline

Each script reads from `_data/` and writes its outputs back to `_data/`. Run from the repository root so relative paths resolve correctly.

Script numbering reflects logical grouping, not a strict run order: Script 0 (`site_metadata.rds`) reads `_data/pot_peaks_all.rds` to get its final site list, so it must be run *after* Script 3, not before it. Script 1 likewise depends on `_data/selected_unit_stations.csv` from Script 2. A working run order is: **2 → 3 → 0 → 1 → 4 → 5 → 6a/6b → 7 → 8 → 9**.
>>>>>>> ed86e0dbeba2b3cceef4c6c9e94a61400a872d51

---

### Script 0 — `0_extract_site_metadata.R`
<<<<<<< HEAD
Assembles `_data/site_metadata.rds` from four GAGES-II attribute tables (basin ID, topography, dams, classification). Derives a regulation flag using the Davenport et al. (2020) upstream storage threshold of 1,000 acre-feet (1.233 mcm).
=======
Assembles `_data/site_metadata.rds` from four GAGES-II attribute tables (basin ID, topography, dams, classification). Derives a regulation flag using the Davenport et al. (2020) upstream storage threshold of 1,000 acre-feet (1.2335 mcm).
>>>>>>> ed86e0dbeba2b3cceef4c6c9e94a61400a872d51

**Inputs:** GAGES-II spreadsheets (`conterm_basinid.txt`, `conterm_topo.txt`, `conterm_hydromod_dams.txt`, `conterm_bas_classif.txt`), `_data/pot_peaks_all.rds`
**Output:** `_data/site_metadata.rds`

---

### Script 1 — `1_create_spatial_data_layers.R`
Builds an `sf` point object (`_data/station_spatial.rds`) from selected stations and GAGES-II metadata (location, topography, reference/non-reference classification). Used by the Shiny dashboard map tab and by downstream trend analysis scripts.

**Inputs:** `_data/selected_unit_stations.csv`, GAGES-II spreadsheets
**Output:** `_data/station_spatial.rds`

---

### Script 2 — `2_evaluate_usgs_unit_data_completeness.R`
Evaluates the hourly (unit) data completeness of all Colorado GAGES-II sites. Processes year-by-year files to manage memory, summarizes completeness by site-year, and writes diagnostic figures. Applies a station selection criterion (100% complete data for 2022 and 20+ consecutive qualifying years) to produce the final station list.

Supports optional chunked execution (`Rscript 2_... <chunk_index> <num_chunks>`) for HPC use.

**Inputs:** USGS unit discharge files (`/path/<site_no>/<year>.rds`), GAGES-II metadata
<<<<<<< HEAD
**Outputs:** `_data/selected_unit_stations.csv`, completeness heatmap and threshold summary figures
=======
**Outputs:** `_data/selected_unit_stations.csv`; completeness heatmap and threshold summary figures written to `_scripts/`
>>>>>>> ed86e0dbeba2b3cceef4c6c9e94a61400a872d51

---

### Script 3 — `3_extract_usgs_unit_data_peaks.R`
Extracts independent peaks-over-threshold (POT) for all selected sites. Applies the Livneh et al. (2024) independence criteria: minimum time separation `θ = 5 + ln(1.6092 × A)` days, and an inter-peak recession condition (flow between peaks must drop below 75% of the lesser peak). The per-site POT threshold is set at the 50th percentile of daily maximum flow (following Yu et al. 2022) and stored as `pot_threshold_cfs` for use as the GPD location parameter `μ` in downstream fitting. Discharge is also converted to mm/day for mass balance diagnostics against WLDAS precipitation.

**Inputs:** `_data/selected_unit_stations.csv`, USGS unit discharge files, GAGES-II metadata
**Output:** `_data/pot_peaks_all.rds` (one row per site-day within each peak window)

---

### Script 4 — `4_extract_wldas_daily_to_gages-ii_basin_mean.R`
Extracts basin-averaged daily WLDAS (NOAHMP001-DA1) variables for all Colorado GAGES-II watersheds. Reads each NetCDF file once and extracts all watersheds simultaneously using `exactextractr`. Full spatial statistics (max, standard deviation, coefficient of variation, wetted fraction) are computed for rainfall, snowfall, and snowmelt; basin mean only is computed for all other variables. Results are written as Parquet files (one per date) for efficient downstream querying.

Supports chunked execution and resumption via a checkpoint file.

**Inputs:** WLDAS `.nc` files, GAGES-II shapefiles, GAGES-II metadata
**Outputs:** `_data/WLDAS/processed/daily_by_date/YYYYMMDD_wldas.parquet` (one per date)

---

### Script 5 — `5_build_accumulation_curves.R`
Constructs 60-day antecedent accumulation curves for each POT peak by joining WLDAS daily Parquet files to the peak windows. For each (peak, lookback day) row the script stores daily rainfall and snowmelt rates and their cumulative sums back from the peak date. Also identifies end-member peaks — those with dominant (top 50th percentile) rainfall or snowmelt and negligible opposing signal (bottom 25th percentile) — for use in classification validation.

**Inputs:** `_data/pot_peaks_all.rds`, WLDAS Parquet files, `_data/site_metadata.rds`
**Outputs:** `_data/accum_curves_60d.rds`, `_data/endmembers_60d.rds`

---

### Script 6a — `6a_classify_usgs_unit_data_peaks_depth_thresh.R`
<<<<<<< HEAD
Classifies each POT peak by dominant flood type using the fixed antecedent-window depth thresholds of Davenport et al. (2020): **Rainfall** (rain > 10 mm, melt < 5 mm), **Snowmelt** (melt > 10 mm, rain < 5 mm), **Rain-on-Snow** (both > 5 mm), over the 8-day window preceding each peak. Applies Yu et al. (2022) top-*n* sample selection (single-sample pool: ≥25 events; per-type pool: ≥20 events, at most one per water year). Types below the minimum sample size are demoted to "unclassified".
=======
Classifies each POT peak by dominant flood type using the fixed antecedent-window depth thresholds of Davenport et al. (2020): **Rainfall** (rain > 10 mm, melt < 5 mm), **Snowmelt** (melt > 10 mm, rain < 5 mm), **Rain-on-Snow** (both > 5 mm), over the 8-day window preceding each peak. Applies Yu et al. (2022) top-*n* sample selection: the largest *m* = max(25, n_water_years) independent events for the single-sample pool, and the largest *n* = max(20, n_water_years) independent events per type, both ranked directly by peak flow with no per-water-year cap. Types below the minimum sample size are demoted to "unclassified".
>>>>>>> ed86e0dbeba2b3cceef4c6c9e94a61400a872d51

**Inputs:** `_data/accum_curves_60d.rds`, `_data/pot_peaks_all.rds`, `_data/site_metadata.rds`
**Output:** `_data/pot_peaks_classified_depth_thresh.rds`

---

### Script 6b — `6b_classify_usgs_unit_data_peaks_frac_thresh.R`
Classifies each POT peak using a fractional approach based on `rain_frac` — the weighted fraction of total active flux (rain + snowmelt above a 0.1 mm/day noise floor) that is attributable to rainfall, computed over each site's basin-adaptive antecedent window (`search_radius` days). Classification rules:

1. **Rainfall:** `rain_frac ≥ 0.99`
2. **Snowmelt:** `rain_frac < MIXED_SPLIT` (population median of mixed events)
3. **Rain-on-Snow:** `rain_frac ≥ MIXED_SPLIT`
4. **Unclassified:** insufficient data

`MIXED_SPLIT` is computed once from all classified non-rainfall peaks pooled across sites, and stored as an attribute on the output. The same Yu et al. (2022) top-*n* selection as Script 6a is applied. This is the primary classification used in Scripts 7–9.

**Inputs:** `_data/accum_curves_60d.rds`, `_data/pot_peaks_all.rds`, `_data/site_metadata.rds`
**Output:** `_data/pot_peaks_classified_frac_thresh_3class.rds`

---

### Script 7 — `7_compile_ffa_results.R`
Pre-computes mixture-model flood frequency analysis (Yu et al. 2022) for all sites and three classification methods:

| Method | Types |
|---|---|
| `depth_thresh` | Rainfall / Snowmelt / Rain-on-Snow |
| `frac_thresh_3class` | Rainfall / Snowmelt / Rain-on-Snow |
<<<<<<< HEAD
| `frac_thresh_2class` | Rainfall / Non-Rainfall (snowmelt + ROS pooled) |
=======
| `frac_thresh_2class` | Rainfall / Snow-Influenced (snowmelt + rain-on-snow pooled) |
>>>>>>> ed86e0dbeba2b3cceef4c6c9e94a61400a872d51

For each (site, method) combination, per-type and mixture return period curves are fitted using both GPD and LP3 distributions, evaluated on a dense RI grid for map tab rendering and a sparse grid (2–500 years) for the summary scatter tab. Parallelised with `furrr` over site × method jobs.

**Inputs:** `_data/pot_peaks_classified_depth_thresh.rds`, `_data/pot_peaks_classified_frac_thresh_3class.rds`, `_data/site_metadata.rds`, `functions.R`
**Outputs:**
- `_data/ffa_curves_{gpd,lp3}.rds` — dense RI curves for map tab
- `_data/ffa_quantiles_{gpd,lp3}.rds` — sparse quantiles for summary tab
- `_data/ffa_params_{gpd,lp3}.rds` — fitted distribution parameters
- `_data/ffa_classification_summary.rds` — peak composition by site and method
- `_data/ffa_ros_temporal.rds` — annual ROS peak counts for temporal tab

---

### Script 8 — `8_analyze_climate_forcing_trends.R`
Derives five POT-based annual time series per site and tests each for monotonic trends using Mann-Kendall / Theil-Sen: total precipitation, peak SWE, annual maximum peak discharge, POT event rate, median rain fraction, and AM rain fraction (rain fraction of the single highest classified peak per water year). Generates five diagnostic figures including a record coverage heatmap, regional aggregate time series, trend maps, slope distribution violins, and an aggregate bar chart tested against a 5% binomial null.

**Inputs:** WLDAS Parquet files, `_data/pot_peaks_classified_frac_thresh_3class.rds`, `_data/site_metadata.rds`, `_data/station_spatial.rds`
<<<<<<< HEAD
**Outputs:** `_data/annual_climate_forcing_summary.rds`, `_data/climate_forcing_trends.rds`, figures in `_diagnostics/climate_forcing_trends/`
=======
**Outputs:** `_data/annual_climate_forcing_summary.rds`, `_data/climate_forcing_trends.rds`, figures in local `_diagnostics/climate_forcing_trends/` (generated on run, not tracked in this repo)
>>>>>>> ed86e0dbeba2b3cceef4c6c9e94a61400a872d51

---

### Script 9 — `9_analyze_flood_type_temporal_trends.R`
Tests for monotonic trends in annual flood type frequency (count and proportion) at all Colorado GAGES-II sites, across all three classification methods. Uses the pre-demotion `flood_type_raw` label so that every classified event is included regardless of per-type sample size eligibility. Also performs:

- **AM indicator tests:** binary Mann-Kendall on whether each flood type generated the annual maximum peak in a given water year
- **Half-period comparison:** first-half vs. second-half AM type frequency
- **AM rain fraction trends** (frac_thresh_3class only)

<<<<<<< HEAD
Produces per-method diagnostic figures (site × year heatmap, regional stacked bar, record-length histograms, trend maps, aggregate summary heatmap, top-6 proportion time series) and saves all trend results and annual count series to `_data/`.

**Inputs:** `_data/pot_peaks_classified_depth_thresh.rds`, `_data/pot_peaks_classified_frac_thresh_3class.rds`, `_data/site_metadata.rds`, `_data/station_spatial.rds`
**Outputs (per method):**
- `_data/temporal_trends_{method}.rds`
- `_data/annual_type_counts_{method}.rds`
- `_data/temporal_trends_am_{method}.rds`
- `_data/annual_am_type_{method}.rds`
- `_data/am_type_shift_{method}.rds`
- Figures in `_diagnostics/temporal_trends/{method}/`
=======
Also computes annual maximum (AM) peak-flow stationarity once, independent of flood-type classification (`_data/annual_am_peakflow.rds`); the resulting trend is duplicated into every `temporal_trends_am_{method}.rds` file below so the app can load it consistently.

Produces per-method diagnostic figures (site × year heatmap, regional stacked bar, record-length histograms, trend maps, aggregate summary heatmap, top-6 proportion time series) and saves all trend results and annual count series to `_data/`.

**Inputs:** `_data/pot_peaks_classified_depth_thresh.rds`, `_data/pot_peaks_classified_frac_thresh_3class.rds`, `_data/pot_peaks_all.rds`, `_data/site_metadata.rds`, `_data/station_spatial.rds`
**Outputs (per method):**
- `_data/temporal_trends_{method}.rds`
- `_data/annual_type_counts_{method}.rds`
- `_data/temporal_trends_am_{method}.rds` — list(am_type_shift, am_type_indicator_trends, am_peak_flow_trends); `frac_thresh_3class` also includes am_rain_frac_trends
- `_data/annual_am_type_{method}.rds`
- `_data/am_type_shift_{method}.rds`
- Figures in local `_diagnostics/temporal_trends/{method}/` (generated on run, not tracked in this repo)

**Output (computed once, classification-independent):** `_data/annual_am_peakflow.rds`
>>>>>>> ed86e0dbeba2b3cceef4c6c9e94a61400a872d51

---

### Supporting scripts

| Script | Purpose |
|---|---|
| `assess_candidate_colorado_usgs_gages.R` | Preliminary screening of candidate Colorado USGS gages prior to the main pipeline |
| `identify_mixture_viable_sites.R` | Identifies sites meeting the Yu et al. (2022) minimum threshold (≥20 peaks per type) for the full three-component mixture model |

---

## Shared utilities

<<<<<<< HEAD
**`functions.R`** — GPD and LP3 distribution fitting and curve evaluation functions shared by Scripts 7–9. Includes:
- `fit_gpd()` / `fit_lp3()` — parameter estimation via L-moments, with `μ` fixed at the site-wide POT threshold (Yu et al. 2022)
- `gpd_sf()`, `gpd_quantile()` — GPD survival function and quantile inversion
- `mixture_sf()`, `mixture_curve()` — mixture return period curve using the Yu et al. (2022) CDF product formulation
- `dist_curve()`, `dist_mixture_curve()` — unified curve evaluation interface used by Script 7
=======
**`functions.R`** — shared by Scripts 7–9 and by `_documentation/report_figures.R`. Includes:
- `fit_gpd()` / `fit_lp3()` — parameter estimation via L-moments, with `μ` fixed at the site-wide POT threshold (Yu et al. 2022)
- `gpd_sf()`, `gpd_quantile()` — GPD survival function and quantile inversion
- `mixture_sf()`, `mixture_curve()` — mixture return period curve using the Yu et al. (2022) CDF product formulation
- `fit_dist()`, `dist_curve()`, `dist_mixture_curve()` — family-agnostic (GPD/LP3) fitting and curve evaluation interface used by Script 7
- `water_year()` — converts a date to USGS water year
- `mk_sens_test()` — Mann-Kendall trend test with Theil-Sen slope, used by Scripts 8–9
- `pettitt_shift_test()` — Pettitt change-point test, used by Script 9
- `empirical_ri()` — empirical (Weibull plotting-position) recurrence interval, used by `report_figures.R`

---

## Report figures

**`_documentation/report_figures.R`** — Builds the figures used in the project's technical report/manuscript (study area map, classification composition, return period and ratio plots, rain-fraction trend maps, etc.) from the pipeline's `_data/` outputs. Rendered PNGs are saved to `_documentation/report_figures/`.
>>>>>>> ed86e0dbeba2b3cceef4c6c9e94a61400a872d51

---

## Data dependencies

| Dataset | Location | Used by |
|---|---|---|
| GAGES-II attribute tables | `/davenport-scratch/DATA/GAGES-II/.../spreadsheets-in-csv-format/` | Scripts 0, 1, 2, 3, 4 |
| GAGES-II watershed shapefiles | `/davenport-scratch/DATA/GAGES-II/.../boundaries-shapefiles-by-aggeco/` | Script 4 |
| USGS unit discharge files | `/davenport-scratch/DATA/USGS/raw/Discharge/Unit/<site_no>/<year>.rds` | Scripts 2, 3 |
| WLDAS NetCDF files | `_data/WLDAS/raw/WLDAS_NOAHMP001_DA1_YYYYMMDD.D10.nc` | Script 4 |

---

## Key references

- **Yu et al. (2022)** — Mixture model FFA framework; GPD fitting with fixed μ; per-type and mixture return period curves
- **Davenport et al. (2020)** — Fixed-depth classification thresholds (Script 6a); upstream storage regulation screen
- **Livneh et al. (2024)** — POT peak independence criteria: minimum separation and inter-peak recession conditions

---

<<<<<<< HEAD
## Dashboard

The analysis outputs produced by this pipeline (Scripts 7–9) are consumed by the companion Shiny dashboard:

**[`CoWC_FY26_Enhanced_Impacts_App`](github.com/realmiketalbot/CoWC_FY26_Enhanced_Impacts_App)** — Interactive R Shiny dashboard for exploring flood type composition, mixture FFA curves, and AM rain fraction trends across the gage network. Deployed at https://019ef057-7c20-639c-982f-e6c1d717ecd0.share.connect.posit.cloud.

The app reads all `.rds` files from `_data/` directly; no additional processing is required between running the pipeline and serving the dashboard.
=======
## Companion repositories

| Repo | Role |
|---|---|
| [`CoWC_FY26_Enhanced_Impacts_App`](https://github.com/realmiketalbot/CoWC_FY26_Enhanced_Impacts_App) | Interactive R Shiny dashboard for exploring flood type composition, mixture FFA curves, and AM rain fraction trends across the gage network. Consumes this pipeline's `.rds` outputs (Scripts 7–9) directly from `_data/`, with no additional processing required. Deployed at https://019ef057-7c20-639c-982f-e6c1d717ecd0.share.connect.posit.cloud. |
| [`CoWC_FY26_Enhanced_Impacts_Diagnostics`](https://github.com/realmiketalbot/CoWC_FY26_Enhanced_Impacts_Diagnostics) | Exploratory analysis, sensitivity studies, and decision documentation behind the methodology implemented in this pipeline. |
>>>>>>> ed86e0dbeba2b3cceef4c6c9e94a61400a872d51

---

## Environment

Package environment is managed with `renv`. To restore:

```r
renv::restore()
```

R project file: `CoWC_FY26_Enhanced_Impacts.Rproj`
