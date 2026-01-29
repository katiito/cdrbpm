# Intervention Analysis for HIV Cluster Detection Strategies

## Quick Start

```r
# Load the pipeline
source("src/run_analysis.R")

# Option 1: Run everything fresh (~2 min interventions + ~5 min mechanism analysis)
run_full_analysis()

# Option 2: Use cached results (fast - only regenerates plots)
run_full_analysis(use_cached = TRUE)

# Option 3: Skip mechanism analysis (fastest)
run_full_analysis(use_cached = TRUE, run_mechanism = FALSE)

# Option 4: Just regenerate plots from most recent saved results
generate_plots_from_cache()
```

## Output Files

| Output | Location |
|--------|----------|
| **IDA efficiency violin plot** | `intervention-plots/efficiency_distributions_violin.pdf` |
| **Mechanism analysis plot** | `intervention-plots/mechanism_analysis.png` |
| Cluster size histograms | `intervention-plots/cluster_sizes_*.png` |
| Intervention metrics (CSV) | `intervention-results/*.csv` |

## Notes

- For mechanism analysis, use `n_sims = NULL` (default) to include all simulations
- With smaller `n_sims`, growth cluster data may be sparse (~0.14% of simulations trigger)
- Cached results load from the most recent files in `intervention-results/`

---

# Intervention model – functional changes (New vs Old)

---
## 1. Overall workflow

* **New**
    * Single entry point: `run_intervention_analysis(...)` loads data, runs all strategies, prints a unified table.
    * Robust error handling: each strategy wrapped in `tryCatch`; failures return zeroed summaries instead of stopping execution.
* **Old**
    * Top-level script blocks executed sequentially; final table structure differed.
    * Uncaught errors could halt the run.

---
## 2. Randomness & seeds

* **New**
    * Optional `seed` argument (defaults to current system time; printed).
    * Seed set explicitly via `set.seed(seed)` for reproducibility.
* **Old**
    * No explicit seed parameter; stochastic elements driven directly by `rexp()` / `ritdist()`.

---
## 3. Inputs & parameterization

* **File I/O**
    * New: Data loaded inside `run_intervention_analysis()` via `d_file`, `g_file`.
    * Old: Data loaded at script top level.
* **Renamed / normalized parameters**
    * `thdist` → `distance_threshold`
    * `thsize` → `cluster_size`
    * `thdegree` → `network_degree_threshold`
    * `ritdist()` (function) → `intervention_rate` (scalar passed to `rexp()`)
* **Additional new options**
    * `subnetwork = "small" | "large"` (dense vs sparse internal cluster contacts)
    * `random_sample_size`
    * `rita_window_months`

---
## 4. Cluster definition & intervention time (distance–size strategy)

* **Shared concept**
    * Grow from index `"0"` along edges ≤ threshold; drop last generation.
* **Intervention trigger**
    * New: If `nrow(G1) ≥ cluster_size`, set `IT = G1$timesequenced[cluster_size] + analysis_delay_days + implementation_delay_days`.
    * Two-stage delay system:
        * `analysis_delay_days` (default 14): Time to analyze cluster after k-th sequence arrives
        * `implementation_delay_days` (default 14): Time to deploy intervention after analysis complete
        * Total delay: 28 days from detection to intervention
    * Old: Loop until size reaches threshold then `timesequenced[i] + ritdist()`.
    * → New uses deterministic fixed delays; old used stochastic delays.
* **Post-intervention filtering**
    * New: Excludes cases with `timesequenced ≥ (IT - analysis_delay_days)` (knowledge cutoff - accounts for data processing lag at intervention time).
    * Old: Excluded cases with `timesequenced ≥ IT`.

---
## 5. Scope of PIA / IDA (behavioural change)

* **New**
    * `piapids = recipients of cluster members' transmissions ∪ cluster members ∪ donors who transmitted to cluster members` (one outward step both directions).
    * Broader set increases potential PIA/IDA.
* **Old**
    * `piapids = recipients of cluster donors ∪ cluster members` (donors who transmitted to cluster members omitted).

---
## 6. Contact modelling

* **New**
    * Computes time-windowed contact counts from partner notification window (default 6 months = 180 days):
        * `total_contacts = Fcontacts_180d + Gcontacts_180d + Hcontacts_180d` (or 90d for 3-month window)
        * These are actual contact counts in the window before diagnosis, not degrees at infection
    * Cluster-based strategies compute two subnetwork contact estimates:
        * Dense (`small`) subnetwork: `contacts = nrow(G1) + sum(pmax(total_contacts - (nrow(G1)-1), 0))`
            * Assumes cluster members are all connected to each other internally
        * Sparse (`large`) subnetwork: `contacts = sum(total_contacts) - (nrow(G1) - 2)`
            * Assumes minimal internal connections among cluster members
    * Individual-based strategies: `contacts = total_contacts + 1` (self + partners)
    * Outputs emphasize per-contact efficiencies.
* **Old**
    * No contact accounting; metrics mostly per member (`/ nc`).

---
## 7. Strategy-specific changes

### A. Distance–size (`distsize_intervention`)
* **New**
    * Reports totals + per-contact efficiency + quantiles.
    * Drops rows with `total_contacts == 0`; returns all-zeros if nothing remains.
* **Old**
    * Means/variances per member; Wald CIs.
* **Both**
    * `propintervened = sum(nc) / sum(generation > 0 & generation < last)`.

### B. Random selection (`random_intervention`)
* **New**
    * Samples `random_coverage` proportion of eligible cases (default 0.10 = 10% of eligible population, gens 1 .. last-1).
    * Sample size computed as: `random_sample_size = round(random_coverage × eligible_population)`
    * Contacts = `total_contacts + 1` where `total_contacts = Fcontacts_XXd + Gcontacts_XXd + Hcontacts_XXd` (time-windowed).
    * Summaries: total IDA, IDA/contact, PIA quantiles (10th / 90th percentiles).
* **Old**
    * Uses all eligible cases.
    * No contact notion; PIA/IDA mean & variance per person.
* **Prop intervened**: Not meaningful / `NA` in both.

### C. RITA-based (`rita_intervention`)
* **New**
    * Configurable detection window: `rita_window_months` (default 6); test simulated via `rexp(1/(months*30))`.
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

---
## 8. Aggregation & summary table

* **New**
    * Columns: `Contacted Total`, `Total IDA`, `IDA/contacted`, `Low`, `High`, `PIA`, `Low`, `High`.
    * Rows: cluster size / threshold labels plus `Random`, `RITA`, `Network, partners>threshold`.
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
* Random strategy now samples a fixed number (`random_sample_size`) rather than all.
* RITA proportion intervened removed (now `NA`).
* Summary table shape & column names changed; any downstream parsing of old columns must be updated.

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
- You can switch the initial contact type by passing `initialcontact = :G` (default), `:F`, or `:H` to `generate_experiment(...)`.
- Threading is optional; if you want to enable threads for future parallel runs, you can prefix commands with `JULIA_NUM_THREADS=auto` (the current scripts themselves do not require it).

---

## Plotting & Analysis Functions

### Source Files

- **`src/intervention0.R`** – Core intervention analysis functions that compute IDA, PIA, and efficiency metrics for each strategy
- **`src/plot_interventions.R`** – Visualization functions for intervention results

### Running the Mechanism Analysis

The mechanism analysis figure explains *why* different intervention strategies have different effectiveness. To generate it:

```r
source("src/plot_interventions.R")
run_mechanism_analysis(
  D_path = "src/experiment1-N10000-gens7-D.csv",
  G_path = "src/experiment1-N10000-gens7-G.csv",
  save_path = "intervention-results/mechanism_analysis.png",
  n_sims = NULL  # Use all simulations (NULL) or specify a number for faster testing
)
```

**Note:** Running with all 10,000 simulations takes several minutes. Progress updates are printed every 10%.

### Mechanism Analysis Figure (5 Panels)

The figure `intervention-results/mechanism_analysis.png` contains:

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

These parameters can be passed to `run_mechanism_analysis()` to match specific intervention analysis settings.

### Efficiency Distributions

To plot efficiency distributions across strategies:

```r
source("src/intervention0.R")
source("src/plot_interventions.R")

# Run interventions first
results <- run_interventions(Dall, Gall)

# Plot efficiency distributions
plot_efficiency_distributions(results, plot_type = "violin")  # or "density", "boxplot"
```
