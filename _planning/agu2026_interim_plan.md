# AGU 2026 Interim Plan: Milestone M1 ("Poisson–GPD at interim thresholds")

> **Created:** 2026-11-17 · **Revised:** 2026-11-17 (upstream work moved ahead of Script 7; full rerun)
> **Parent plan:** `_planning/chapter1_methods_alignment_plan.md` (workstream IDs such as WS3.1 refer to it)
> **Source:** `_documentation/Preliminary_Proposal.pdf`, Ch. 1 §1.2–1.3, App. A
> **Deadline:** poster print in ~2 weeks (target ≈ 2026-12-01)
> **Status:** Phase 0 done (tags). Phase U not started. D-E open.

## 1. The stopping point

The full alignment adds four layers to the initial results: (a) the right *distribution*,
(b) the right *sample*, (c) the right *threshold*, and (d) the *uncertainty and diagnostics* around them.
Layers (a) and (b) are cheap and mostly deterministic. Layers (c) and (d) carry the real research
risk: the EQD procedure, the L-moment substitution and its simulation study, and the double bootstrap.

**M1 = upstream data contract finalized, layers (a) and (b) done, (c) held by a documented interim rule, (d) deferred.**

At M1:

- The upstream scripts (3, 6a, 6b) produce their **final-shaped outputs**. The candidate threshold has
  its own name, effective record length is computed from the daily series, and classification only labels
  peaks: no top-*n* sampling and no demotion.
- Every quantile is a **Poisson–GPD annual** quantile (Eq. 1.3–1.5), not a PoT recurrence interval.
- Each type is fit to its exceedances of a type-specific threshold u_i, set by a simple, stated
  **interim rule** (§6, D-A). EQD replaces it later without touching anything else.
- LP3 is retired from the primary pipeline.
- The App shows these results, labels them as preliminary, and its schema already holds `u`, `λ`, `n_exc`
  and empty CI columns. Later steps fill in columns instead of changing structure.

### Why stop here and not one step further

EQD (WS3.4) is the next step. The proposal says outright that the L-moment variant is untested,
and makes WS5 (the simulation study) its gate. Shipping EQD thresholds on a poster before
that gate means presenting a method you've told your committee you will validate first. An interim rule
that's clearly labelled is easier to defend than an unvalidated version of the final method.

## 2. A defect to fix regardless (found while planning)

`fit_gpd()` in `functions.R` calls `lmomco::pargpa(lmoms(y))` without `xi =`, so it **estimates a
location** from the excesses. The code then discards that location and substitutes μ. The scale and
shape it reports therefore come from a 3-parameter fit, while the quantiles assume a 2-parameter
fit. I checked this locally (lmomco, `set.seed(1)`, 40 exponential excesses):

| Call | ξ_loc | α (scale) | κ (= −ξ) |
|---|---|---|---|
| `pargpa(lmoms(y))` | 1.30 | 16.7 | −0.078 |
| `pargpa(lmoms(y), xi = 0)` | 0 (fixed) | 20.6 | +0.060 |

The sign of the shape flips. This resolves the ⚠️ in parent WS3.1. **Fix: pass `xi = 0` on excesses.**
The initial results in the proposal (§1.3, Fig. 1.4) have this defect. That's worth knowing, but
it doesn't need to be raised on the poster, because M1 supersedes those numbers anyway.

There is a related issue. The current single μ is the site **median** daily max, but the fitted sample is the
top-*n* peaks, which mostly sit well above the median. The excesses therefore don't start near 0, and
a GPD with a fixed lower bound of 0 fits them poorly. Per-type u_i (item F3) fixes this.

## 3. Pre-AGU scope (do now)

The work runs in three phases. Phase U is upstream and runs on the HPC. Phase F (fitting) and Phase A
(App) run locally.

### Phase 0: Freeze

| # | Item | Effort |
|---|---|---|
| Z1 | ✅ **Done 2026-11-17.** Annotated tag `initial-results-prelim-proposal` on Support `04b4f7c` and App `8c40def`. Both working trees were clean, and `_data/` is tracked in both repos, so the tag captures every data file exactly. **No copy into `_data/_archive/` was made:** it would have added ~90 MB of duplicates per repo on a nearly full disk. Read any baseline file with `git show initial-results-prelim-proposal:_data/<file>.rds > <scratch path>`. The tags are local until pushed. | 15 min |

### Phase U: Upstream pass and full rerun (HPC)

Make every change in U1–U3, then rerun Scripts 3 → 5 → 6a/6b → 8 → 9 once. Check the run against
the invariants in U4 before any Script 7 work.

| # | Item | Parent WS | Notes |
|---|---|---|---|
| U1 | **Script 3** | WS1 | |
| | • Rename `pot_threshold_cfs` → `candidate_threshold_cfs` | WS1.1 | Ripple is small: Scripts 3 (7 uses), 6a (5), 6b (1), and 7 (1, which is rewritten anyway). The App never references it. Six old `_diagnostics/` scripts use it; leave them as historical and add a note to the folder's README. |
| | • **Effective record length** per site from the daily-max series: `n_days_obs`, `eff_years = n_days_obs / 365.25`, and `n_wy_complete` (water years ≥ 90% covered). Save `_data/site_record_length.rds`. | WS1.3 | The data is already in memory in Script 3 (`all_flows`). See D-E. |
| | • Declustering check: assert or log that no gap equals θ exactly | WS1.2 | Expected to be a no-op. Records the confirmation. |
| | • Remove the "compare against previous methodology" block; fix the header citation (USWRC 1975/76, Lang et al. 1999) | WS1.4 | |
| U2 | **Scripts 6a/6b** | WS2.1, 2.3 | |
| | • Remove the top-*n* selection, `in_single_sample`, `in_type_sample`, `MIN_*`, and demotion. `flood_type_raw` becomes the only label. Either drop `flood_type`, or keep it as an exact copy for one release so the App doesn't break mid-transition. | WS2.1 | Recommend keeping the copy. Remove it in M2. |
| | • Add an explicit **eligibility flag** `fit_eligible`: 6a uses `complete_window`; 6b uses `data_sufficient & flux_sufficient & !is.na(rain_frac)`. | supports D3 | This is the logic the single-sample selection uses today, made explicit so it isn't lost when that selection block is removed. |
| | • Move thresholds into a config list at the top (6a: 5/10 mm; 6b: α_r = 0.99, α_s rule, noise floor), with **defaults equal to today's values**. Write the settings to output attributes. | WS2.3 | No sweep yet. This is only so S4 and M4 don't need another edit. |
| U3 | Copy `remap_2class()` into `functions.R` so Script 7, the App and `report_figures.R` share one definition | housekeeping | Prevents the three copies drifting apart. |
| U4 | **Invariant check** script (`_scripts/checks/check_m1_invariants.R`, not `_diagnostics/`, which `.gitignore` excludes), comparing against files extracted from the Z1 tag: | — | Any failure: stop and explain it before moving on. |
| | • Script 3: identical `(site_no, peak_id, peak_date, peak_flow)` and threshold values (only the column name changes) | | Confirms the rename and housekeeping didn't change the peak set |
| | • Script 5: identical accumulation curves and endmembers | | |
| | • 6a/6b: identical `flood_type_raw`, `rain_frac`, and flux columns | | Confirms removing top-*n* and parameterizing changed only what was meant to change |
| | • Scripts 8/9: identical outputs (both read `flood_type_raw`) | | |
| U5 | Sync regenerated `_data/` to the local Support repo (and the subset the App needs to the App repo, after Phase A) | — | |

**Why one combined run works here:** every Phase U change is either a refactor (rename, config,
housekeeping), purely additive (record length, `fit_eligible`), or a removal of downstream-only columns
(top-*n* flags). None should change which peaks exist or how they're labelled. The U4 checks
confirm that directly, so any difference can be traced without needing separate runs.

### Phase F: Fitting (local; Script 7 + `functions.R`)

| # | Item | Parent WS | Effort |
|---|---|---|---|
| F1 | Fix `pargpa` lower bound (`xi = 0`). | WS3.1 | trivial |
| F2 | New `fit_pgpd(x, u, eff_years)` → `u, sigma, xi, lambda, n_exc`; `pgpd_cdf`, `pgpd_quantile(T)` with T = 1/(1−F) **annual**; `pgpd_mixture_curve` via F_mix = exp(−Σ λ_i[1−H_i]). Below u_i, set F_i = exp(−λ_i) (D4 floor) and report quantiles only where x ≥ max(u_i). | WS3.1–3.3 | ½–1 day |
| F3 | **Interim threshold rule** (D-A) per (site, method, type) and for the single sample. Store `u`, `u_rule`, `n_exc`, `lambda`. | stand-in for WS3.4 | ½ day |
| F4 | Build samples from `flood_type_raw` (after the 2-class remap) and `fit_eligible`. Replace `MIN_TYPE_SAMPLE` with `MIN_TYPE_COUNT` applied to `n_exc`. λ uses `eff_years` from `site_record_length.rds`. | WS3.2, 3.6 | ½ day |
| F5 | **Retire LP3** (`FAMILIES <- "pgpd"`). Old LP3 files stay in the archive. | D5 | trivial |
| F6 | Outputs `ffa_{curves,quantiles,params}_pgpd.rds` with the new columns, plus empty `lo`/`hi` CI columns. | WS3.9 | ½ day |
| F7 | Spot-check 3–5 sites by eye (curves vs annualized empirical points), then rerun all 206 sites. | — | ½ day |

### Cheap diagnostics worth one poster panel each (pick ≤ 2)

| # | Item | Parent WS | Why it's cheap now |
|---|---|---|---|
| S1 | **Dispersion index** at the interim u_i, per (site, type), with χ² bounds (App. A.2), over the complete water years (`n_wy_complete`). Map or histogram showing over- and under-dispersion. | WS4.4 | Only needs annual counts above u_i |
| S2 | **Old vs new** Q100 mixture/single ratio (Yu-GPD vs Poisson–GPD). This is a direct "sensitivity to model specification" result and fits the chapter title. | RQ4 framing | Both outputs exist after Z1 and F7 |
| S3 | **Interim-threshold sensitivity**: rerun F3 at 2–3 rule settings and show how the ratio by return period moves. It shows *why* EQD is needed. | WS4.1-lite | It's a loop over Script 7 |
| S4 | **Rain-fraction α sensitivity**: rerun 6b at 2–3 (α_r, α_s) settings through the U2 config. | WS4.2-lite | The config is already in place after U2 |

### Phase A: App (App repo; light touch, schema-forward)

| # | Item | Notes |
|---|---|---|
| A1 | Load `ffa_*_pgpd.rds`. **Remove the Distribution toggle** (GPD/LP3). | Touches the ~10 `input$dist_family` references; replace with a constant `"pgpd"` |
| A2 | Replace the `in_type_sample`/`in_single_sample` logic in `site_emp_peaks()` and the `in_ss` flag with "peak > u_i", using `u` from the params table. The 2-class top-*n* re-derivation goes away entirely. Use the shared `remap_2class()`. | ~app.R:1011–1050, 1148–1155 |
| A3 | **Annual empirical plotting positions** for exceedances: T = (M+1)/rank with M = `eff_years`, which matches λ. Update `empirical_ri()` in `functions.R`. | Makes the points line up with the curves |
| A4 | Parameters table: show `u`, `λ`, `n_exc`, `σ`, `ξ`. | ~app.R:1381 |
| A5 | **Text**: change "Yu et al. (2022) framework" (~app.R:5, 329, 2748) to Poisson–GPD (Beguería 2005). Add a one-line "Preliminary: interim thresholds; EQD threshold selection and uncertainty bounds in progress" banner. Remove the LP3 description. The return-period labels already say "% annual chance", and M1 makes that statement true. | |
| A6 | Copy the regenerated classified-peak, trend and FFA files into the App's `_data/`. | |
| A7 | Optional: CWC logo on About (open TODO). | 10 min |

Also update the `report_figures.R` lines ~900–908 the same way as A2. That's the only part of it that
depends on the old top-*n* flags.

The App deliberately **doesn't** get CI ribbons, a threshold-diagnostic view, a sensitivity explorer,
or new tabs. Those come with M2 (§5).

## 4. Hold off (post-AGU)

| Item | Parent WS | Reason to wait |
|---|---|---|
| EQD threshold selection (L-moment) | WS3.4 | Needs Murphy et al. (2025) details. Gated on WS5. |
| Simulation study | WS5 | Gate for EQD. Weeks of work, not days. |
| Double bootstrap CIs | WS3.5 | Depends on EQD. Heavy HPC. |
| Full classification-sensitivity grid incl. the depth method | WS4.2 | U2's config makes this a loop later. S4 covers the rain-fraction slice. |
| `MIN_TYPE_COUNT` sweep as a formal analysis | WS4.1 | Belongs with EQD, since minimum count and threshold choice interact. |
| Inter-type independence tests | WS4.3 | Cheap, but adds a claim that needs careful interpretation. Good early-M2 item. |
| 3-class vs 2-class diagnostic | WS2.4 | The poster can show both without arbitrating. |
| Rain-fraction window and noise floor | WS2.2 | A decision, not code. The June plan already chose to keep the adaptive window. Write it up in the methods text. |
| Hamshaw benchmark | WS2.6 | Optional. |
| Trends: exclude incomplete first and last years (open TODO) | — | Unrelated to Ch. 1 methods. Phase U shouldn't change trend outputs, and U4 checks that. |
| README rewrite, full `report_figures.R` overhaul | WS6 | Update only what the poster needs. |
| Drop the `flood_type` copy | WS2.1 | M2, once the App no longer reads it. |

## 5. What comes after M1

- **M2: threshold.** WS5, then WS3.4. This swaps the F3 rule for EQD and fills `u_rule = "eqd"`.
  WS4.1 and WS4.4 run across the candidate grid. App adds a threshold-diagnostic plot.
- **M3: uncertainty.** WS3.5 fills the `lo`/`hi` columns. App adds CI ribbons, and RQ4 claims
  move to "difference exceeds uncertainty".
- **M4: sensitivity and diagnostics.** WS4.2 (now a config loop), WS2.4, WS4.3, WS6.

## 6. Decisions

**Settled (2026-11-17):** D-A, D3, D4, D5, D-B, D-C and D-D as recommended below.
Upstream work happens before Script 7, with a full rerun.

| # | Decision | Choice | Alternatives |
|---|---|---|---|
| **D-A** | Interim threshold rule | **Rate-based:** u_i = the k-th largest eligible type-*i* peak with k = max(`MIN_TYPE_COUNT`, round(r·eff_years)), r = 1 exceedance/yr. This is the closest analogue to the current top-*n*, so old-vs-new differences reflect the *model* change rather than the sample change. | Fixed quantile of type-specific candidates (e.g., 0.75) |
| D3 | Single-sample peaks | **All peaks eligible for classification** (`fit_eligible`), with the same rule applied. | Classified-only |
| D4 | Mixture below u_i | **Floor F_i = exp(−λ_i)**, and report only x ≥ max(u_i) | Truncate only |
| D5 | LP3 | **Retire from the pipeline and the App.** Archive the old outputs. | Keep as supplementary |
| D-B | `MIN_TYPE_COUNT` | **20** (unchanged) | 10–15 would admit more sites |
| D-C | Poster's primary method | **Rain fraction (3-class)**, with 2-class alongside. Depth threshold as context only. | — |
| D-D | Redeploy to the proposal's dashboard URL? | **Deploy M1 under a new app name**, keeping the old deployment as a frozen snapshot | Overwrite |
| **D-E** *(settled 2026-11-17)* | Effective record length definition | **λ denominator = days with daily-max data / 365.25.** It needs no completeness cutoff and matches how λ is used. **Dispersion index M = water years ≥ 90% covered**, since DI needs whole-year counts. | Complete water years for both (simpler, but discards partial-year peaks from λ's denominator and not its numerator) |

## 7. Schedule

| Days | Work |
|---|---|
| 1 | Z1. U1–U3 code changes. Run Script 3 on the HPC. |
| 2 | Run Scripts 5 → 6a/6b → 8 → 9. U4 invariant checks. U5 sync. |
| 3–4 | F1–F6, plus spot-checks on 3–5 sites |
| 5 | F7 full run. Regenerate the Fig. 1.4-style panels. |
| 6–7 | S1 + S2 (or S3/S4). Poster figures. |
| 8–9 | A1–A6, local test, deploy (per D-D) |
| 10 | Poster layout. **This is the only buffer left.** If time runs short, cut the S items, then A5 polish, in that order. |

Phase U has used most of the original buffer. **If the U4 checks fail on day 2 and the cause isn't
clear within half a day,** fall back: keep the frozen upstream data and point Script 7 at it with
compatibility code. Debug Phase U after AGU.

**Done means:** both repos are tagged `agu2026-m1`; U4 passes; the App runs on Poisson–GPD outputs
with the preliminary banner; the poster figures regenerate from committed scripts; and this file's
status is updated with what slipped.

## 8. Assumptions and open items

- The rate-based rule (D-A) is a placeholder. The poster should say "fixed exceedance-rate
  threshold", not imply it is optimal.
- I haven't run any script. Script run times come from the user ("don't take that long").
- U4 assumes Script 5 is deterministic for the same inputs. If it uses parallel reads or any
  randomness, compare with a tolerance rather than exact equality.
- Script 3's "jitter" (proposal §1.2.2) is deterministic: it subtracts a monotonic offset,
  `eps * seq_len(n)` with `eps = 1e-6 × max`, so there's no RNG and no seed needed. Rerun peaks should
  match exactly. The proposal's word "jitter" suggests randomness. Consider "a monotonic tie-breaking
  perturbation" in the manuscript.
- The App runtime cost should drop (one family instead of two). I haven't measured it.
- The Support `.gitignore` contains committed merge-conflict markers (`<<<<<<< HEAD` … `>>>>>>>`).
  They're harmless, because git reads them as patterns that match nothing, and `_diagnostics/` is
  still ignored. Clean them up whenever convenient.
