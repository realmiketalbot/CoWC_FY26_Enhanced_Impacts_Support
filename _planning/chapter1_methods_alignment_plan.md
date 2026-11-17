# Chapter 1 Methods Alignment Plan

> **Created:** 2026-10-03
> **Source:** `_documentation/Preliminary_Proposal.pdf`, Chapter 1 (§1.1–1.4) and Appendix A
> **Status:** Planning only — no code has been changed

## Purpose

The current pipeline produces the proposal's **initial results** (§1.3): the Yu et al. (2022)
simplification, which fits a plain GPD to the top-*n* POT peaks as a proxy for an annual maximum
series. The proposal commits to a different **final** method. This document lists every gap
between the code and §1.2, then groups the fixes into workstreams.

Chapter 1 research questions the revised pipeline has to answer:

1. Do the data support three distinct populations (rainfall, snowmelt, rain-on-snow)?
2. Which watersheds show mixtures of these types, and at what relative frequencies?
3. How sensitive are mixed-population quantiles to the classification thresholds?
4. Where, by how much, and at which return periods do mixture and single-sample quantiles differ?

---

## 1. Gap summary

| Area | Proposal (§1.2) | Current code | Status |
|---|---|---|---|
| Streamflow source | 206 gages via `dataRetrieval`, completeness screened, daily max | Pre-downloaded unit files on `/davenport-scratch`; Script 2 screens; Script 3 builds daily max | ✅ Aligned (wording only) |
| Forcing | WLDAS 1979–2023, area-weighted catchment means via `exactextractr` | Script 4 | ✅ Aligned |
| Candidate peak detection | `findpeaks`, jitter, descending pad at gaps, above per-site median | Script 3 Stages 1–2 | ✅ Aligned |
| Declustering | Eq. 1.1: θ > 5 + ln(0.3861·A), X_min < 0.75·min(Q1,Q2) | Script 3 `find_independent_peaks()` | ⚠️ `≥` vs `>` boundary (see §2.1) |
| Candidate vs flood threshold | Median = *candidate* threshold only; separate per-type *model-fitting* threshold u_i | `pot_threshold_cfs` (median) doubles as GPD μ | ❌ Must decouple |
| Depth-threshold classification | 8-day window; 10/5 mm rules | Script 6a | ✅ Aligned |
| Rain-fraction classification | α_r = 0.99; α_s = median α of non-rainfall peaks | Script 6b | ⚠️ Window length and noise floor not in proposal |
| Sample selection | No top-*n*; flood peaks are those above the EQD-selected u_i | Yu top-*n* (`in_type_sample`, `in_single_sample`) and demotion in 6a/6b | ❌ Remove from classification |
| Minimum type count | Swept as a sensitivity parameter | Fixed `MIN_TYPE_SAMPLE = 20`, `MIN_SINGLE_SAMPLE = 25` | ❌ Make configurable and sweep |
| Distribution | Poisson–GPD (Eq. 1.3–1.4), L-moments via `lmomco` | Plain GPD on PoT-RI axis; LP3 alongside it | ❌ Replace |
| Mixture | Product of Poisson–GPD annual-max CDFs (Eq. 1.5) | Product of conditional GPD CDFs | ❌ Replace |
| Threshold selection | EQD (Murphy et al. 2025) with L-moments in place of MLE | None (μ = median) | ❌ New |
| Threshold uncertainty | Double bootstrap (Murphy et al. 2025) → CIs on return levels | None | ❌ New |
| L-moment EQD validation | Limited simulation study at type-specific sample sizes | None | ❌ New |
| Inter-type independence | Co-occurrence and magnitude association tests, network distribution | None | ❌ New |
| Poisson assumption | Dispersion index with χ² bounds (App. A.2), at selected and candidate thresholds | None | ❌ New |
| 3-class vs 2-class support | Diagnostics on whether a snowmelt / rain-on-snow split is supported | 2-class remap exists in Script 7; no diagnostic | ⚠️ Diagnostic needed |
| Classification sensitivity | Vary thresholds (depth 5/10 mm; α_r, α_s) | Single fixed setting per method | ❌ New |

---

## 2. Workstreams

### WS1 — Peak extraction (Script 3)

**1.1 Decouple the candidate threshold from the fitting threshold.**
Rename `pot_threshold_cfs` to `candidate_threshold_cfs` (or similar). It filters local maxima
and nothing else. The GPD location now comes from WS3. Every consumer of `pot_threshold_cfs`
must change: 6a, 6b, 7, `functions.R` (`fit_gpd(mu_threshold = ...)`), and the App repo.

**1.2 Check the declustering inequality.**
Eq. 1.1 says θ **>** 5 + ln(0.3861·A). `is_dependent()` treats `< min_sep_days` as dependent,
so a separation exactly equal to θ counts as independent. Because θ is non-integer and dates are
integers, this probably never binds. Confirm that, then fix either the code or the equation so they match.

**1.3 Record effective years of record per site.**
The Poisson rate λ_i needs a proper denominator. Compute the number of complete (or
sufficiently complete) water years per site from the daily series and save it in the output.
`fit_gpd()` currently uses `diff(range(peak dates))`, which is wrong for λ (see WS3.2).

**1.4 Housekeeping.**
- Remove or gate the "compare against previous methodology" block at the end of Script 3. It
  reads an archived `.rds` and is a one-off diagnostic.
- Make the header citation match the proposal: USWRC (1975/1976), Lang et al. (1999).

### WS2 — Classification (Scripts 5, 6a, 6b)

**2.1 Remove top-*n* selection and demotion from 6a/6b.**
Classification should only label peaks. Drop `in_single_sample`, `in_type_sample`, the
`MIN_*` constants, and the demotion of under-sized types. Keep `flood_type_raw` as the one
label; `flood_type` becomes redundant. Magnitude filtering moves to the u_i threshold in WS3.
*Downstream impact:* `remap_2class()` and `get_sample()` in Script 7, `report_figures.R`
(lines ~900–908), Script 9 comments, and the App repo.

**2.2 Settle the rain-fraction window and noise floor (decision, see §4).**
6b uses the basin-adaptive `search_radius` (7–17 days) and `NOISE_FLOOR_MM = 0.1`. Section 1.2.2 of the proposal
doesn't give the α window or a noise floor. Either keep both and document them in the
manuscript, or change the code. The 2026-06-22 development plan already argues for the adaptive window.

**2.3 Parameterize thresholds for sensitivity runs (RQ3).**
Expose the depth thresholds (5/10 mm), α_r, and the α_s rule as arguments or a config list so
the classification can run over a grid. Store the threshold settings as attributes or columns on the output.

**2.4 Diagnostic: is a snowmelt / rain-on-snow split supported? (RQ1)**
The proposal says no dividing line is apparent below α_r. Candidate checks: multimodality tests
on α for α < α_r (pooled and by region or elevation), mixture-model BIC (1 vs 2 components), and
whether 3-class and 2-class quantiles differ materially. The result decides whether
`frac_thresh_3class` stays a primary method.

**2.5 Diagnostic: mixture prevalence (RQ2).**
Per-site type fractions, the number of types meeting the minimum count, and a map. Most of
this already exists in `ffa_classification_summary.rds`. Rebuild it on unfiltered labels.

**2.6 Optional: validate against the Hamshaw et al. (2026) benchmark** wherever it overlaps
Colorado sites.

### WS3 — Flood quantile estimation (`functions.R`, Script 7)

**3.1 Poisson–GPD fitting functions.**
New `fit_pgpd(x, u, n_years)` returning `u`, `sigma`, `xi`, `lambda`, `n_exc`, with:
- `H(y) = 1 − (1 + ξy/σ̃)^(−1/ξ)` fit to excesses y = x − u by L-moments
- `F(x) = exp(−λ[1 − H(x − u)])` for x ≥ u
- `T(x) = 1 / (1 − F(x))`, a true **annual** return period. The current code uses PoT RI,
  so every downstream axis changes meaning.

⚠️ **Verify the `lmomco::pargpa` call.** The existing `fit_gpd()` calls
`pargpa(lmoms(y))` without fixing the lower bound. If `pargpa` takes an optional known lower bound
(`xi =`) and otherwise estimates it, then the current fits estimate a location that the code
then discards and replaces with μ. The new function should pass the threshold as a known lower
bound (`xi = 0` on excesses). I couldn't confirm this locally because renv points at the HPC library path.

**3.2 Exceedance rate λ_i.**
λ_i = (number of type-i flood peaks > u_i) / (effective years of record), using the WS1.3
denominator. Use the same denominator for every type at a site.

**3.3 Mixture CDF.**
`F_mix(x) = ∏ F_i(x) = exp(−Σ λ_i[1 − H_i(x − u_i)])`, inverted numerically.
*Design issue:* F_i is defined only for x ≥ u_i, and each type gets its own u_i. Choose and
document the treatment below u_i: F_i(x) = exp(−λ_i) as a floor, or report mixture quantiles
only for x ≥ max(u_i). This sets the lowest reportable return period per site.

**3.4 EQD threshold selection with L-moments.**
For each (site, type) and for the single sample:
1. Build a candidate threshold grid from the flood-peak distribution, e.g. empirical quantiles,
   subject to a minimum number of exceedances.
2. For each candidate, bootstrap the excesses, fit the GPD by L-moments, and compute the expected
   quantile discrepancy against the empirical exceedances.
3. Select the candidate that minimizes the EQD.
Implementation details such as grid spacing, bootstrap count, and the discrepancy metric must
follow Murphy et al. (2025). Read the paper or its code before implementing. Nothing here comes from it yet.

**3.5 Double bootstrap for threshold uncertainty.**
Follow the Murphy et al. (2025) procedure: an outer bootstrap re-runs threshold selection and an
inner one fits parameters, giving return-level CIs that include threshold uncertainty. New outputs: lower and upper bounds per
(site, method, curve_type, return_period). This is the most expensive step. Plan for HPC
(`furrr` is already used) and save the selected thresholds and the bootstrap draws.

**3.6 Replace the Yu-style sample construction in Script 7.**
Remove `MIN_TYPE_SAMPLE`/`MIN_SINGLE_SAMPLE` gating and `get_sample()` top-*n* logic. Replace
them with a configurable `MIN_TYPE_COUNT` applied to exceedances of u_i (WS4.1 sweeps it).

**3.7 Single-sample definition (decision, see §4).**
State which peaks the single-sample fit uses. Options: all candidate peaks, including ones
unclassified by the depth rules, or only classified peaks. The choice changes the depth-threshold
comparison noticeably, since many peaks are unclassified under that method.

**3.8 LP3 (decision, see §4).** Chapter 1 doesn't use LP3. Retire it, or keep it as clearly
labelled supplementary output.

**3.9 Outputs.** Rework `ffa_curves_*`, `ffa_quantiles_*`, and `ffa_params_*` around Poisson–GPD
(add `u`, `lambda`, `n_exc`, CI columns, selected-threshold metadata). Update `empirical_ri()`
to annual plotting positions consistent with λ.

### WS4 — Sensitivity analysis and diagnostics (new scripts)

**4.1 Minimum-count sweep.** Re-run WS3 across a range of `MIN_TYPE_COUNT`. Report the
number and identity of analyzable sites and the mixture/single-sample quantile ratio by return period.

**4.2 Classification sensitivity (RQ3).** Cross the WS2.3 threshold grid with both methods,
then report how quantiles and analyzable sites respond.

**4.3 Inter-type independence.** Per site:
- *Co-occurrence:* whether types occur in the same water years more than chance predicts
  (e.g. Fisher's exact test on annual presence/absence)
- *Magnitude association:* Kendall's τ between type-specific annual maxima in years where both occur
Report the network distribution of the statistics, because power at individual sites is low.

**4.4 Dispersion index (App. A.2).** DI = s² / λ̂ of annual exceedance counts. Accept if
χ²_{α/2, M−1}/(M−1) < DI < χ²_{1−α/2, M−1}/(M−1). Compute it per (site, type) at the
selected u_i and across the EQD candidate grid. Label the direction: overdispersed points to
declustering failure, underdispersed to regular snowmelt.

**4.5 Mixture vs single-sample comparison (RQ4).** Ratio by return period and site, with
fraction ≥ 1, maps, and stratification by regulation, region, and drainage area. Use the
WS3.5 CIs to judge which differences exceed uncertainty.

### WS5 — Simulation study for L-moment EQD

A standalone script or folder, independent of the observed data:
- Follow the Murphy et al. (2025) design: synthetic series with a known threshold and known GPD tail
- Use sample sizes typical of type-specific series. Take them from the observed distribution of
  per-type candidate counts.
- Compare L-moment EQD with MLE EQD on threshold recovery and quantile bias/RMSE
- Gate: WS3.4 results are provisional until this is done.

### WS6 — Documentation and downstream consumers

- Update the `README.md` script descriptions, which still cite Livneh et al. 2024 with the 1.6092 constant, the Yu top-*n* selection, and PoT RI.
- Update `report_figures.R` for the new outputs and the annual RI axis.
- Update the App repo (`CoWC_FY26_Enhanced_Impacts_App`) for schema changes. Track that work in its own plan.
- Archive the initial-results `_data/ffa_*` files before regenerating them, so the §1.3 figures stay reproducible.

---

## 3. Proposed sequence

1. **WS1**: decouple thresholds and add effective years. Small, and unblocks everything else.
2. **WS2.1–2.3**: strip selection from classification and parameterize thresholds.
3. **WS3.1–3.3**: Poisson–GPD and mixture functions, tested on a few sites with a fixed u.
4. **WS5**: simulation study, which can run in parallel with step 3.
5. **WS3.4–3.9**: EQD, double bootstrap, and Script 7 rewrite.
6. **WS4**: diagnostics and sensitivity sweeps.
7. **WS2.4–2.6**: classification diagnostics. These can start any time after step 2.
8. **WS6**: documentation and consumers.

---

## 4. Open decisions

| # | Decision | Options | Notes |
|---|---|---|---|
| D1 | α window for rain-fraction method | Adaptive `search_radius` (current) / fixed 8-day | Proposal is silent; 2026-06-22 plan favors adaptive |
| D2 | Noise floor in 6b | Keep 0.1 mm and document / remove | Not mentioned in proposal |
| D3 | Peaks in the single-sample fit | All candidates / classified only | Matters most for depth-threshold method |
| D4 | Mixture below type thresholds | Floor F_i(x<u_i)=exp(−λ_i) / truncate at max(u_i) | Sets lowest reportable return period |
| D5 | LP3 | Retire / keep as supplementary | Not in Chapter 1 |
| D6 | 3-class vs 2-class as primary | Pending WS2.4 | Directly answers RQ1 |
| D7 | `MIN_TYPE_COUNT` sweep range | e.g. 10–30 | Minimum must also suit EQD candidate grid |
| D8 | Bootstrap sizes (EQD, double bootstrap) | Per Murphy et al. (2025) | Drives HPC cost |

## 5. Assumptions and items to verify

- Murphy et al. (2025) procedural details (EQD metric, candidate grid, bootstrap design) are
  not yet captured here. Read the source before WS3.4/3.5.
- `lmomco::pargpa` lower-bound behaviour (WS3.1 warning) has not been verified in this environment.
- I assumed the proposal's "per-catchment median flow" means the median of daily maximum flow, as
  Script 3 computes it.
- I didn't review Scripts 8–9 (trend analysis) in detail. Chapter 1 doesn't describe trend
  methods, so only the schema changes in WS2.1 should affect them.
