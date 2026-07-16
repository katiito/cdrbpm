# Intervention Analysis for HIV Cluster Detection Strategies

## Quick Start

```r
# Load the pipeline
source("src/run_analysis.R")

# Option 1: Run everything fresh (interventions + paired comparison + mechanism analysis)
run_full_analysis()

# Option 2: Reuse cached intervention results, regenerate all four plots
run_full_analysis(use_cached = TRUE)

# Option 3: Reuse cache but skip the slow mechanism analysis (faster)
run_full_analysis(use_cached = TRUE, run_mechanism = FALSE)

# Option 4: Reuse cache, efficiency plot only (fastest)
run_full_analysis(use_cached = TRUE, run_mechanism = FALSE, run_paired = FALSE)

# Option 5: Regenerate plots from one *specific* saved run
#   (timestamp = a YYYYMMDD_HHMMSS from intervention-results/; omit the arg for most recent)
generate_plots_from_cache(load_cached_results("YYYYMMDD_HHMMSS"))

# Option 6: Merge every run from a given day (e.g. several seeds) into one combined figure
plot_grant_comparison(load_todays_cached_results())  # defaults to today
```

`run_full_analysis()` runs all eight strategies, prints a unified summary table, saves
per-strategy CSVs, and generates four figures (efficiency distributions, paired
comparison vs random, probability of beating random, and the mechanism analysis).
Cluster-size histograms are written during the intervention run itself.

Helper functions:

- `load_cached_results(timestamp = NULL)` – load a single run (most recent by default).
- `load_todays_cached_results(date = today)` – load every run from a given day, e.g. to
  merge results across multiple seeds before plotting.

## The eight intervention strategies

The pipeline evaluates and compares eight strategies (see [src/intervention0.R](src/intervention0.R)):

| # | Strategy | Trigger | Key parameters |
|---|----------|---------|----------------|
| 1 | **Distance–size (k=5)** | 5 cases within a genetic-distance cluster | `cluster_size_5`, `distance_threshold` |
| 2 | **Distance–size (k=2)** | 2 cases within a genetic-distance cluster | `cluster_size_2`, `distance_threshold` |
| 3 | **Growth-rate (k=5)** | 5 cases sequenced within a sliding time window | `cluster_size_5`, `lookback_window_months`, `growth_distance_threshold` |
| 4 | **Growth-rate (k=2)** | 2 cases sequenced within a sliding time window | `cluster_size_2`, `lookback_window_months`, `growth_distance_threshold` |
| 5 | **Random allocation** | Baseline: randomly selected diagnosed individuals | `random_coverage` |
| 6 | **RITA** | Recently-infected cases (Recent Infection Testing Algorithm) | `rita_window_months` |
| 7 | **Network degree** | Well-connected individuals (many contacts) | `network_degree_threshold` |
| 8 | **RITA + Secondary** | RITA-positive index cases plus 2-degree contact tracing | `rita_window_months` |

**Growth-rate vs distance–size:** distance–size triggers on the *k*-th sequenced case in a
cluster regardless of time span; growth-rate triggers only when *k* cases are sequenced
within `lookback_window_months`, so it detects *rapidly growing* clusters rather than simply
large ones. The growth strategies use a separate (typically larger) distance threshold,
`growth_distance_threshold`.

**RITA + Secondary:** identifies RITA-positive index cases, then traces their transmission
network backward and forward (primary and secondary links) within the notification window,
and notifies the contacts of every traced individual.

**Two subnetwork assumptions** are computed for every cluster-based strategy:

- **`small` (dense)** – cluster members are assumed to be highly connected to each other, so
  overlapping contacts are counted once (fewer total contacts).
- **`large` (sparse)** – cluster members are assumed to have minimal internal connections, so
  contacts are summed (more total contacts).

Individual-based strategies (Random, RITA, Network) report a single row (subnetwork = `-`).

## Output Files

Plots are written to `intervention-plots/` and metrics to `intervention-results/`. Most files
carry a `_YYYYMMDD_HHMMSS` timestamp so runs don't overwrite each other.

| Output | Location | Produced by |
|--------|----------|-------------|
| **IDA efficiency violin plot** | `intervention-plots/efficiency_distributions_<TS>.png` | `plot_efficiency_distributions()` |
| **Paired comparison vs random** | `intervention-plots/paired_comparison_percent_<TS>.png` | `plot_paired_comparisons()` |
| **Probability of beating random** | `intervention-plots/prob_beats_random_<TS>.png` | `plot_prob_beats_random()` |
| **Mechanism analysis plot** | `intervention-plots/mechanism_analysis_<TS>.png` | `run_mechanism_analysis()` |
| Cluster-size histograms | `intervention-plots/cluster_sizes_<TS>.png` | `run_intervention_analysis()` |
| Grant-figure comparison | `intervention-plots/grant_comparison.{png,svg,pdf}` | `plot_grant_comparison()` (run manually) |
| Summary counts (units + contacts) | `intervention-results/counts_<TS>.csv` | `run_intervention_analysis()` |
| Per-strategy detail rows | `intervention-results/details_<strategy>_<TS>.csv` | `run_intervention_analysis()` |
| Run parameters | `intervention-results/parameters_<TS>.csv` | `run_intervention_analysis()` |

Per-strategy detail files are written for each of the eight strategies: `distsize5`,
`distsize2`, `growth5`, `growth2`, `random`, `rita`, `network`, `ritasecondary`.

## Key parameters & defaults

Defaults for `run_intervention_analysis()` (see [src/intervention0.R](src/intervention0.R#L106)):

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `cluster_size_5` | 5 | Larger cluster/growth trigger size |
| `cluster_size_2` | 2 | Smaller cluster/growth trigger size |
| `distance_threshold` | 0.005 | Genetic distance for distance–size clusters |
| `growth_distance_threshold` | 0.01 | Genetic distance for growth-rate clusters |
| `lookback_window_months` | 6 | Sliding window for the growth-rate trigger |
| `network_degree_threshold` | 4 | Contacts required to be a network-targeting unit |
| `partner_notification_window_months` | 6 | Contact-tracing lookback (3 or 6 months → 90d/180d) |
| `random_coverage` | 0.30 | Proportion of the eligible population sampled at random |
| `rita_window_months` | 6 | Mean RITA detection window |
| `analysis_delay_days` | 14 | Cluster-analysis delay (cluster-based strategies only) |
| `implementation_delay_days` | 14 | Intervention-implementation delay (all strategies) |
| `seed` | `NULL` | Random seed; defaults to system time and is printed |

## Notes

- For mechanism analysis, use `n_sims = NULL` (default) to include all simulations.
- With smaller `n_sims`, the size-5 growth-cluster data can be sparse; the size-2 growth
  trigger fires far more often and is much less sensitive to `n_sims`.
- RITA status is assigned **once, up front** (in `run_intervention_analysis()`, stored on
  `Gall$rita`) and reused by every strategy, so all strategies see the same recent-infection
  labelling and the comparison is fair.
- `intervention0.R` uses lookup tables for the per-individual transmission/contact
  computations, which substantially speeds up full runs.
- Cached results load from the most recent files in `intervention-results/`.

## Comparing strategies against random

Two complementary analyses quantify how each strategy performs relative to the random
baseline.

### Paired comparison (`plot_paired_comparisons`)

For each strategy, matches outbreaks to the random baseline by simulation ID and reports the
**relative impact** (percent difference in per-contact efficiency vs random) on the large
(sparse) subnetwork. This is included in `run_full_analysis()` by default (toggle with
`run_paired`) and saved as `paired_comparison_percent_<TS>.png`.

### Probability of beating random (`plot_prob_beats_random`)

Estimates the probability that a strategy beats random allocation, at two levels:

- **Per outbreak** – the fraction of shared simulations in which the strategy's per-contact
  efficiency exceeds random.
- **Population level (pooled bootstrap)** – restricts to simulation IDs present in *both* the
  strategy and random, pools IDA and contacts by summing across outbreaks, and bootstraps over
  the shared simulation IDs (default 1,000 replicates) to estimate
  `P(pooled strategy efficiency > pooled random efficiency)`.

Saved as `prob_beats_random_<TS>.png`. Plot subtitles annotate the number of runs and the
number of simulations with onward transmission.

---

# Intervention model – functional changes (New vs Old)

This repository is a **fork** of the original CDR/BPM code by Erik Volz. Throughout this
section:

- **"Old"** = the original upstream code as first pulled — the initial April 2024 commit, in
  which `src/intervention0.R` was ~257 lines and used top-level scripts with parameters like
  `thdist`, `thsize`, and a stochastic `ritdist()` delay function.
- **"New"** = the current rewritten implementation in this repository (`intervention0.R` is now
  ~1,800 lines, wrapped in `run_intervention_analysis()`).

The subsections below document how the rewrite differs from that original, function by
function; the "New" bullets describe current behaviour.

---
## 0. Outcome terminology (PUTA → IDA)

The two reported outcomes were **relabelled**; the underlying calculations did not change.

| Old label | New label | What it measures | Definition (identical old & new) |
|-----------|-----------|------------------|----------------------------------|
| **PUTA** — potential undiagnosed time averted | **IDA** — infectious days averted | Days people remain undiagnosed (and thus potentially infectious) after the intervention, which earlier diagnosis averts | `sum(timediagnosed - IT)` over people infected before `IT` but diagnosed after it |
| **PIA** — potential infections averted | **PIA** — potential infections averted *(unchanged)* | Onward infections that occur after the intervention | `sum(timeinfected > IT)` over the connected population |

`PUTA` was renamed to `IDA` because the old acronym is an offensive word in Spanish; the change
is cosmetic — `IDA` and the old `PUTA` are the same quantity. The rename is applied throughout
the current code: the summary table and `summary_<TS>.csv` use `Total_IDA` / `Median_IDA` /
`Low_IDA` / `High_IDA`, and the per-strategy `details_<strategy>_<TS>.csv` files use an `ida`
column. Cached results saved before mid-2026 store these values in a `puta` column instead; the
plotting code reads either name for backward compatibility.

---
## 1. Overall workflow

* **New**
    * Single entry point: `run_intervention_analysis(...)` loads data, runs all eight strategies, prints a unified table, and saves per-strategy CSVs.
    * Robust error handling: each strategy wrapped in `tryCatch`; failures return zeroed summaries instead of stopping execution.
* **Old**
    * Top-level script blocks executed sequentially; final table structure differed.
    * Uncaught errors could halt the run.

---
## 2. Randomness & seeds

* **New**
    * Optional `seed` argument (defaults to current system time; printed).
    * Seed set explicitly via `set.seed(seed)` for reproducibility.
    * RITA status is drawn once up front and reused by every strategy for a fair comparison.
* **Old**
    * No explicit seed parameter; stochastic elements driven directly by `rexp()` / `ritdist()`.

---
## 3. Inputs & parameterization

* **File I/O**
    * New: Data loaded inside `run_intervention_analysis()` via `d_file`, `g_file`.
    * Old: Data loaded at script top level.
* **Renamed / normalized parameters**
    * `thdist` → `distance_threshold`
    * `thsize` → `cluster_size` (now split into `cluster_size_5` and `cluster_size_2`)
    * `thdegree` → `network_degree_threshold`
    * `ritdist()` (function) → `intervention_rate` (scalar passed to `rexp()`)
* **Additional new options**
    * `subnetwork = "small" | "large"` (dense vs sparse internal cluster contacts)
    * `random_coverage` (proportion sampled; default 0.30)
    * `rita_window_months`
    * `growth_distance_threshold`, `lookback_window_months` (growth-rate trigger)

---
## 4. Cluster definition & intervention time (distance–size strategy)

**Cluster definition (shared by New and Old).** Starting from index case `"0"`, follow
transmission links with genetic distance ≤ `distance_threshold` to build the connected cluster
`G1`, then drop the last generation (not yet observable). Members are ordered by sequencing
time.

**Intervention timing — New.** Let `t_k = G1$timesequenced[cluster_size]` be the time the
*k*-th cluster member is sequenced — the moment the cluster reaches the trigger size. If the
cluster has at least `cluster_size` members, three events follow, each separated by a fixed
(deterministic) delay:

| Event | Time | With defaults |
|-------|------|---------------|
| **Detection** — *k*-th case sequenced, cluster reaches trigger size | `t_k` | — |
| **Analysis complete** — cluster investigated | `t_k + analysis_delay_days` | `t_k + 14 d` |
| **Intervention (`IT`)** — intervention deployed | `t_k + analysis_delay_days + implementation_delay_days` | `t_k + 28 d` |

So `IT = timesequenced[cluster_size] + analysis_delay_days + implementation_delay_days` — i.e.
**28 days after detection** with the defaults.

**Which members get traced — New (knowledge cutoff).** Sequence data only become visible after
the analysis lag, so at the intervention moment `IT` the team knows only about cases sequenced
up to `IT - analysis_delay_days`. The cluster that is actually contact-traced is therefore
restricted to members with `timesequenced < IT - analysis_delay_days` (which equals
`t_k + implementation_delay_days`). Members sequenced after that cutoff exist but aren't yet
visible, so they are excluded — they don't count toward `nc` (the number of traced contacts),
and this filtered `G1` is also what defines the PIA/IDA population (Section 5).

**Old (original forked code).**
* Trigger: loop through the cluster until it reaches `thsize`, then
  `IT = timesequenced[thsize] + ritdist()`, where `ritdist()` draws a **stochastic** delay from
  `rexp(1, rate = 1/90)` (mean 90 days) — not fixed 14 + 14-day delays.
* Membership: excluded only cases with `timesequenced ≥ IT`; there was no separate
  analysis-lag knowledge cutoff.

---
## 5. Scope of PIA / IDA (behavioural change)

* **New**
    * `piapids = recipients of cluster members' transmissions ∪ cluster members ∪ donors who transmitted to cluster members` (one outward step both directions).
    * Broader set increases potential PIA/IDA.
* **Old**
    * `piapids = recipients of cluster donors ∪ cluster members` (donors who transmitted to cluster members omitted).

---
## 6. Contact modelling

**The headline change:** the new code introduces a *per-contact denominator for every strategy*
by explicitly counting time-windowed contacts and summing them. The old code had **no summed
contact count anywhere** — the only weighted-degree calculation lived in the network strategy
and was used solely to *select* targets, never to measure per-contact efficiency.

**New.** Contacts are actual (unweighted) counts within the partner-notification window
(default 6 months = 180 days), taken from the window before diagnosis rather than degrees at
infection:

* `total_contacts = Fcontacts_180d + Gcontacts_180d + Hcontacts_180d` (or the `_90d` columns for a 3-month window).
* Cluster-based strategies compute two subnetwork estimates:
    * Dense (`small`): `contacts = nrow(G1) + sum(pmax(total_contacts - (nrow(G1)-1), 0))` — members assumed fully connected internally.
    * Sparse (`large`): `contacts = sum(total_contacts) - (nrow(G1) - 2)` — minimal internal connections.
* Individual-based strategies: `contacts = total_contacts + 1` (self + partners).
* Efficiency is reported per contact: `sum(IDA) / sum(contacts)` and `sum(PIA) / sum(contacts)`.

**Old.** No per-contact denominator anywhere; handling differed by strategy:

| Strategy | Contact / degree handling | Used for |
|----------|---------------------------|----------|
| Distance–size | `nc = nrow(G1)` (cluster member count) | denominator of per-member metrics (`puta/nc`, `pia/nc`) |
| Random | none | per-person `puta`, `pia` (no contacts) |
| RITA | none | per-person `puta`, `pia` (no contacts) |
| Network | weighted degree `Fdegree + (7/2)·Gdegree + (7·30)·Hdegree` | **selection threshold only** (`degree ≥ thdegree`, default 30); metrics still per person |

The weighted degree therefore existed in just one strategy (network) and only chose *whom* to
target — it was never a count of contacts that got summed.

---
## 7. Strategy-specific changes

### A. Distance–size (`distsize_intervention`)
* **New**
    * Run at two trigger sizes (k=5 and k=2).
    * Reports totals + per-contact efficiency + quantiles.
    * Drops rows with `total_contacts == 0`; returns all-zeros if nothing remains.
* **Old**
    * Means/variances per member; Wald CIs.
* **Both**
    * `propintervened = sum(nc) / sum(generation > 0 & generation < last)`.

### B. Random selection (`random_intervention`)
* **New**
    * Samples `random_coverage` proportion of eligible cases (default 0.30 = 30% of eligible population, gens 1 .. last-1).
    * Sample size computed as: `random_sample_size = round(random_coverage × eligible_population)`
    * Contacts = `total_contacts + 1` where `total_contacts = Fcontacts_XXd + Gcontacts_XXd + Hcontacts_XXd` (time-windowed).
    * Summaries: total IDA, IDA/contact, PIA quantiles (10th / 90th percentiles).
* **Old**
    * Uses all eligible cases.
    * No contact notion; PIA/IDA mean & variance per person.
* **Prop intervened**: Not meaningful / `NA` in both.

### C. RITA-based (`rita_intervention`)
* **New**
    * Configurable detection window: `rita_window_months` (default 6); RITA status pre-assigned once via `rexp(1/(months*30))`.
    * Contacts = `total_contacts + 1` where `total_contacts = Fcontacts_XXd + Gcontacts_XXd + Hcontacts_XXd` (time-windowed).
    * Reports totals + per-contact efficiency + quantiles (10th / 90th percentiles); returns sum of contacts.
* **Old**
    * Fixed 6-month window.
    * Per-person mean / variance only.
* **Prop intervened**
    * New: `NA`
    * Old: `nrow(G1) / nrow(Gall)`

### D. Network-degree targeting (`network_intervention`)
* **New**
    * Targets individuals with `total_contacts ≥ network_degree_threshold` (default 4).
    * Contacts = `Fcontacts_XXd + Gcontacts_XXd + Hcontacts_XXd` (time-windowed, within partner notification window).
    * Outputs: per-contact IDA/PIA efficiency + totals; quantiles (10th / 90th percentiles).
* **Old**
    * Weighted degree ( `F + (7/2)*G + (7*30)*H` ) approximating partners over 7 months.
    * Reports per-person metrics and proportion intervened.
* **Parameter name**
    * New: `network_degree_threshold`
    * Old: `thdegree`

### E. Growth-rate clusters (`growthrate_intervention`) — new
* Triggers when `cluster_size` cases are sequenced within a sliding `lookback_window_months`
  window (run at k=5 and k=2), rather than on the k-th case at any time span.
* Builds the genetic cluster with `growth_distance_threshold` (default 0.01), then applies the
  same fixed `analysis_delay_days` + `implementation_delay_days` delays as distance–size.
* Reports the same cluster-based outputs (small/large subnetwork totals, per-contact
  efficiency, quantiles).

### F. RITA + secondary contact tracing (`rita_secondary_intervention`) — new
* Starts from RITA-positive index cases (gens 1..last-1), then uses the transmission tree to
  identify **primary** traced individuals (linked to/from the index) and **secondary** traced
  individuals (linked to/from primary), all within the notification window.
* Transmission links are used to *identify* individuals; contacts are then counted from
  `Fcontacts + Gcontacts + Hcontacts`. Small subnetwork sums contact counts (no overlap);
  large subnetwork takes the max (complete overlap).
* Reported as a cluster-based strategy (small/large rows).

---
## 8. Aggregation & summary table

* **New**
    * Detail columns: `Strategy`, `Subnetwork`, `Contacts`, total IDA, `IDA/contact`, median/low/high IDA, total PIA, `PIA/contact`, median/low/high PIA.
    * Cluster-based strategies emit two rows each (`small`/`large`); individual-based strategies emit one row (`subnetwork = "-"`).
    * Rows cover all eight strategies, including Growth (k=5/k=2) and RITA+Secondary.
    * A companion `counts` table reports `Units`, `Contacts_Small`, `Contacts_Large` per strategy.
    * Rounded; printed via `knitr::kable()` if available, else plain `data.frame`.
* **Old**
    * Columns: `Proportion intervened`, `Mean`, `Variance`, `lb`, `ub` (CI around mean IDA per member).

---
## 9. Edge cases & defaults

* **New**
    * Empty selections return structured zeroed summaries (`total_contacts = 0`).
    * Uses empirical quantiles: 10th / 90th percentiles uniformly for all strategies and metrics (IDA and PIA efficiency).
* **Old**
    * Empty sets could yield undefined means/variances.
    * Confidence intervals via normal approximation.

---
## 10. Backward-incompatible changes (attention)

* Shift from per-member means/variances to per-contact efficiency + totals.
* Broader PIA/IDA population in distance–size (adds donors of cluster recipients) → larger effects possible.
* Network targeting degree definition changed (unweighted vs weighted) → threshold semantics differ.
* Random strategy now samples a fixed proportion (`random_coverage`, default 30%) rather than all.
* RITA proportion intervened removed (now `NA`).
* Strategy set expanded from four to eight (added Growth k=5/k=2 and RITA+Secondary).
* Summary table shape & column names changed, and per-strategy CSVs are timestamped; any downstream parsing must be updated.

---

## Re-running the Julia simulations

The Julia code generates the synthetic datasets consumed by the R analysis. You can regenerate them from the project root on macOS/zsh.

- One-time (or after pulling changes): install the exact dependencies pinned by `Manifest.toml`.

```zsh
julia --project -e 'using Pkg; Pkg.instantiate()'
```

- Quick rerun with the defaults in `src/run_generate_experiment.jl` (10 sims, 5 generations, prefix "experiment1"):

```zsh
julia --project src/run_generate_experiment.jl
```

This writes CSVs into `src/` (now including the generation count in the name):
- `src/experiment1-N10-gens5-G.csv`
- `src/experiment1-N10-gens5-D.csv`

- Custom run via CLI (choose how many simulations and generations):

```zsh
# Usage: julia --project src/generate_experiment.jl <n_sims> [out_prefix] [maxgenerations]
julia --project src/generate_experiment.jl 100 experiment1 5
```

Outputs will be placed in `src/` with names like:
- `src/experiment1-N100-gens5-G.csv`
- `src/experiment1-N100-gens5-D.csv`

- From the Julia REPL (optional):

```zsh
julia --project
```

Then in the REPL:

```julia
using Pkg; Pkg.instantiate()
include("src/generate_experiment.jl")
generate_experiment(100; out_prefix="experiment1", maxgenerations=5)
```

Notes:
- Filenames now include both the number of simulations (N) and the number of generations (gens), e.g., `experiment1-N100-gens5-*.csv`. Update any downstream scripts to point at the new names.
- The R pipeline defaults to `experiment1-N10000-gens7-*.csv`. A larger `experiment1-N50000-gens7-*.csv` dataset is also available; pass it via `d_file` / `g_file` (or `D_path` / `G_path`) for higher-precision runs.
- You can switch the initial contact type by passing `initialcontact = :G` (default), `:F`, or `:H` to `generate_experiment(...)`.
- Threading is optional; if you want to enable threads for future parallel runs, you can prefix commands with `JULIA_NUM_THREADS=auto` (the current scripts themselves do not require it).

### HIV-TRACE-like pairwise distances (`compute_pairwise_distances.jl`)

The main simulation only stores genetic distances between direct donor–recipient pairs.
HIV-TRACE, by contrast, computes distances between **all** sequenced individuals — including
pairs not directly linked by transmission (e.g. "siblings" infected by the same source).
[src/compute_pairwise_distances.jl](src/compute_pairwise_distances.jl) reconstructs those
full pairwise distances from the transmission tree within each simulation, producing
HIV-TRACE-style output for downstream clustering comparisons.

```zsh
julia --project src/compute_pairwise_distances.jl
```

---

## Plotting & Analysis Functions

### Source Files

- **`src/intervention0.R`** – Core intervention analysis functions that compute IDA, PIA, and efficiency metrics for each strategy
- **`src/plot_interventions.R`** – Visualization functions for intervention results
- **`src/run_analysis.R`** – Master pipeline (`run_full_analysis`, `generate_plots_from_cache`, cache loaders)

### Running the Mechanism Analysis

The mechanism analysis figure explains *why* different intervention strategies have different effectiveness. It is run automatically by `run_full_analysis()`, or standalone:

```r
source("src/plot_interventions.R")
run_mechanism_analysis(
  D_path = "src/experiment1-N10000-gens7-D.csv",
  G_path = "src/experiment1-N10000-gens7-G.csv",
  save_path = "intervention-plots/mechanism_analysis.png",
  n_sims = NULL  # Use all simulations (NULL) or specify a number for faster testing
)
```

**Note:** Running with all 10,000 simulations takes several minutes. Progress updates are printed every 10%. (When called from the pipeline, the plot is written to `intervention-plots/`.)

### Mechanism Analysis Figure (5 Panels)

The figure `intervention-plots/mechanism_analysis_<TS>.png` contains:

**Panel A: Targeting active transmission chains**
- Shows % of intervention targets whose donor is still undiagnosed at intervention time
- Higher % = strategy catches people in active transmission chains
- RITA performs best because recently-infected people are more likely to have undiagnosed donors

**Panel B: RITA intervenes early in transmission**
- Among people who transmitted, shows % of transmissions that would be prevented by intervention
- Higher % = intervention happens earlier in transmission career
- RITA catches people before they've done most of their transmitting

**Panel C: Growth clusters identify higher transmitters**
- Mean number of transmissions per intervention target
- Growth cluster members have higher transmission counts than random population
- This is due to survivorship bias (Panel D), not inherent risk

**Panel D: Survivorship bias in growth clusters**
- Compares transmission rates of growth cluster offspring vs general population by generation
- Key insight: Enrichment is ~2x in generation 1, ~1.7x in generation 2, but disappears by generations 3-4
- This demonstrates that growth clusters identify "lucky" epidemic branches, not inherently high-risk subpopulations

**Panel E: Growth cluster delay components**
- Breaks down the time from diagnosis to intervention for growth cluster strategy:
  - **Dx to Sequencing** (sequencing delay): Time from diagnosis to sequence availability
  - **Sequencing to Trigger** (analysis delay): Time from sequencing until growth trigger condition met
  - **Trigger to Intervention** (implementation delay): Fixed delay after trigger (28 days = 14 analysis + 14 implementation)

### Parameters

The mechanism analysis uses parameters matching `intervention0.R` defaults:
- `implementation_delay_days` (default: 14 days)
- `analysis_delay_days` (default: 14 days)
- `rita_window_months` (default: 6 months, exponential distribution with mean 180 days)
- `distance_threshold_growth` (default: 0.01, for growth clusters)
- `distance_threshold_distsize` (default: 0.005, for size-based clusters)
- `lookback_window_months` (default: 6 months, for growth trigger)
- `partner_notification_window_months` (default: 6 months = 180 days, determines contact window)
- `network_degree_threshold` (default: 4 contacts, for network-based targeting)

When called from the pipeline, `partner_notification_window_months` and
`network_degree_threshold` are read from the saved `parameters_<TS>.csv` and passed through so
the figure matches the intervention run; the remaining values use the mechanism function's
internal defaults.

### Efficiency Distributions

To plot efficiency distributions across strategies:

```r
source("src/intervention0.R")
source("src/plot_interventions.R")

# Run interventions first
results <- run_intervention_analysis()

# Plot efficiency distributions (violin)
plot_efficiency_distributions(results)
```

`plot_efficiency_distributions(results, title_prefix, timestamp, n_sims, n_sims_total)` returns
a ggplot object; `n_sims` / `n_sims_total` are used only to annotate the subtitle.

## Exploratory / diagnostic scripts

These standalone scripts in `src/` are not part of the main pipeline; they were used to
investigate specific questions and can be sourced individually:

- [src/analyze_paired_comparison.R](src/analyze_paired_comparison.R) – earlier standalone paired comparison (now folded into `run_full_analysis()`).
- [src/analyze_fair_rita_growth_comparison.R](src/analyze_fair_rita_growth_comparison.R) – decomposes RITA's apparent delay advantage into strategy effects vs population (infection→diagnosis) effects.
- [src/analyze_growth_rita_efficiency.R](src/analyze_growth_rita_efficiency.R) – whether the number and timing of RITA-positive cases inside a growth cluster affect IDA/contact.
- [src/analyze_rita_position_in_growth.R](src/analyze_rita_position_in_growth.R) – effect of a RITA-positive case's position within the growth trigger window.
- [src/analyze_trigger_window_only.R](src/analyze_trigger_window_only.R) – fair RITA-vs-growth comparison restricted to the trigger-window members only.
- [src/rerun_interventions_with_simid.R](src/rerun_interventions_with_simid.R) – re-runs the full analysis preserving simulation IDs for paired comparison.
