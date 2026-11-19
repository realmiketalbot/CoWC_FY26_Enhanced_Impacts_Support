# AGU 2026 poster: content outline

> Status: content outline, 2026-11-19. Session H024 (Advancing flood characterization, modeling and
> communication). Poster only. Figures come from `_documentation/poster_figures.R` and
> `_documentation/poster_showcase.R`; numbers are from the M1 Poisson–GPD outputs (Script 7).
> Related: `agu2026_interim_plan.md` (decisions D-A … D-F, the QR pin).

## Header

- **Title:** A Regional Mixed-Population Flood Frequency Analysis: Insights from Colorado
  (the abstract title without "Across 206 Watersheds").
- **Authors:** Michael Talbot and Frances V. Davenport, Civil and Environmental Engineering,
  Colorado State University.
- **Acknowledgement:** FY26 CoWC Enhanced Impacts Initiative.
- **QR code:** to the AGU version of the app. PINNED until the app is deployed (plan D-D).
- Label: "Preliminary results".

## 1. Motivation

- Flood frequency analysis (FFA) assumes peaks come from one identically distributed population.
- In Colorado, floods come from rainfall, snowmelt and rain-on-snow, which may follow different
  distributions.
- **Question:** does a mixed-population FFA (one distribution per flood type, then combined) change
  design-flood estimates, at which return periods, and where?
- Visual: the "identically distributed" schematic from the Hydrology Days talk, or a three-curves-into-
  a-mixture sketch.

## 2. Data and methods

- **Streamgages:** USGS instantaneous discharge, water years 1986–2023 (unit records start in WY1986).
  A water year qualifies if May–Sep is ≥ 90% observed and ≥ 95% approved; sites need ≥ 20 qualifying
  years. **215 gages**, regulated and unregulated.
- **Peaks:** independent peaks declustered with the USWRC/Lang criterion (separation grows with basin
  area).
- **Classification:** WLDAS (Noah-MP) rainfall and snowmelt fluxes summed over an antecedent window
  that scales with basin size.
  - *Adaptive*: antecedent rain fraction. Rainfall ≥ 0.99; the remaining "mixed" events split at
    the network median into rain-on-snow and snowmelt (3-class), or kept together as non-rainfall
    (2-class).
  - *Fixed*: the antecedent depth thresholds of earlier work, for comparison.
- **Frequency model:** Poisson–GPD per flood type (L-moments; threshold ≈ one exceedance per
  qualifying year, ≥ 20 exceedances). The mixture multiplies the types' annual-maximum distributions.
  It is compared with a single Poisson–GPD fitted to all peaks. Return periods are annual.
- **Visuals:** the antecedent-window plot and decision tree from the talk (check they don't say
  206), plus `poster_figures/02_composition_by_method.png`.
  - Callout: *Adaptive classification labels ~98% of peaks, vs ~55% for Fixed.*

## 3. Example site: Tomichi Creek at Sargents (09115500)

- **Figure:** `poster_figures/showcase/09115500_A_fixed_vs_adaptive3.png`.
  - Optional: `09115500_B_adaptive3_vs_adaptive2.png`.
- **Message:**
  - Fixed misses many peaks, so its mixture falls below observed common floods.
  - Adaptive tracks the observations.
  - The mixture rises above the single sample for rarer floods: 100-yr mixture/single 1.30
    (Adaptive 3-class) vs 1.16 (Fixed).

## 4. Network results (the core)

1. **Headline:** `poster_figures/04_ratio_by_return_period.png`.
   - *Mixture and single-sample estimates agree for common floods but diverge for rare ones.*
   - At 100 yr the mixture is higher at 85% of sites (median +10%, Adaptive 3-class).
   - At 2 yr it is higher at about 27% of sites (median −2%).
2. **Elevation:** NEW panel, to build. Q100 ratio vs mean basin elevation, or a box plot by band.
   - Median ratio 1.03 below 2,500 m, 1.11 at 2,500–3,000 m, 1.26 above 3,000 m.
   - Spearman ρ = 0.47, p < 10⁻¹¹.
   - *Largest differences in high-elevation basins where flood drivers compete.*
3. **Flood type in the tail:** `poster_figures/09_dominant_type_by_rp.png`.
   - Rain-on-snow becomes the largest flood type at more sites as return periods lengthen:
     ~15% of sites at 2 yr, ~35% at 100 yr.
4. **Map:** `poster_figures/06_ratio_map_q100.png` (Tomichi circled).
5. **Robustness:** one line or a small inset; detail in `14_`/`15_` as backup.
   - The divergence holds across fitting thresholds (median Q100 ratio 1.05–1.28) and rain-fraction
     cutoffs (1.10–1.16).
   - It depends on classification: Fixed 1.16, Adaptive 3-class 1.10, Adaptive 2-class 1.01. So
     separating rain-on-snow from snowmelt drives most of it.
- **Context (text):** 91% of sites (196/215) have at least two flood types each making up ≥ 10% of
  their peaks.

## 5. Takeaways

1. Floods at nearly all Colorado gages come from more than one process. Adaptive rain-fraction
   classification labels nearly every peak.
2. Mixed-population and single-sample estimates agree for common floods and diverge for rare ones.
   The mixture is higher at most sites by the 100-yr flood, most strongly above 3,000 m.
3. Separating rain-on-snow from snowmelt matters most, which raises the question of whether that
   split is physically distinct.

## 6. Next steps

- Threshold selection by EQD (Murphy et al. 2025), replacing the interim rule.
- Double-bootstrap confidence intervals: are the differences significant?
- Rain-on-snow vs snowmelt diagnostic.
- Goodness-of-fit comparison across classification methods (LP3 was dropped from this version).
- **Trends:** preliminary tests show no detectable change in the flood type of annual maxima,
  1986–2023 (5.7% of sites significant, consistent with chance). Smaller peaks suggest a shift
  toward snowmelt that we are testing for record-coverage and spatial-correlation artefacts.
- Interactive dashboard for Colorado water managers (QR code).

## Backup (keep at hand, not on the poster)

- Dispersion check (`13_`): 83–89% of fits are consistent with Poisson.
- Shape-parameter distribution (`12_`).
- Sensitivity panels (`14_`, `15_`).
- Other showcase candidates: Ohio Creek (09113980), East River (09112500), San Juan at Pagosa
  (09342500).

## To do before layout

- [ ] Poster size and orientation; then re-export figures at print size with larger text.
- [ ] Build the elevation panel (section 4.2) in the poster style.
- [ ] Statement-style figure titles (e.g., "Estimates diverge only for rare floods").
- [ ] Check the talk's methods visuals for site-count references.
- [ ] Phase A app deploy, then the new QR code (PIN).
