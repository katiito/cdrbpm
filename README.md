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
    * New: If `nrow(G1) ≥ cluster_size`, set `IT = G1$timesequenced[cluster_size] + rexp(1, intervention_rate)`.
    * Old: Loop until size reaches threshold then `timesequenced[i] + ritdist()`.
    * → Functionally equivalent; new is simpler & explicit.
* **Post-intervention filtering**
    * Both exclude cases with `timesequenced ≥ IT`.

---
## 5. Scope of PIA / PUTA (behavioural change)

* **New**
    * `piapids = recipients of cluster donors ∪ cluster members ∪ donors of cluster recipients` (one outward step both directions).
    * Broader set increases potential PIA/PUTA.
* **Old**
    * `piapids = recipients of cluster donors ∪ cluster members` (donors-of-recipients omitted).

---
## 6. Contact modelling

* **New**
    * Computes contact counts where needed.
    * Per member degree = `Fdegree + Gdegree + Hdegree`.
    * Dense (`small`) subnetwork: `total_contacts = nrow(G1) + sum(pmax(degree - (nrow(G1)-1), 0))`.
    * Sparse (`large`) subnetwork: `total_contacts = total_degree - (nrow(G1) - 2)`.
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
    * Samples up to `random_sample_size` (default 30) from eligible cases (gens 1 .. last-1).
    * Degree = `F + G + H`; contacts per individual = `degree + 1`.
    * Summaries: total PUTA, PUTA/contact, PIA quantiles (1% / 99%).
* **Old**
    * Uses all eligible cases.
    * No contact notion; PIA/PUTA mean & variance per person.
* **Prop intervened**: Not meaningful / `NA` in both.

### C. RITA-based (`rita_intervention`)
* **New**
    * Configurable detection window: `rita_window_months` (default 6); test simulated via `rexp(1/(months*30))`.
    * Adds degree-based contacts; reports totals + per-contact efficiency + quantiles; returns sum of contacts.
* **Old**
    * Fixed 6-month window.
    * Per-person mean / variance only.
* **Prop intervened**
    * New: `NA`
    * Old: `nrow(G1) / nrow(Gall)`

### D. Network-degree targeting (`network_intervention`)
* **New**
    * Unweighted degree = `F + G + H` (changes qualifying set).
    * Outputs: per-contact PUTA efficiency + totals; (comment notes mixed scaling in PIA quantiles).
* **Old**
    * Weighted degree ( `F + (7/2)*G + (7*30)*H` ) approximating partners over 7 months.
    * Reports per-person metrics and proportion intervened.
* **Parameter name**
    * New: `network_degree_threshold`
    * Old: `thdegree`

---
## 8. Aggregation & summary table

* **New**
    * Columns: `Contacted Total`, `Total PUTA`, `PUTA/contacted`, `Low`, `High`, `PIA`, `Low`, `High`.
    * Rows: cluster size / threshold labels plus `Random`, `RITA`, `Network, partners>threshold`.
    * Rounded; printed via `knitr::kable()` if available, else plain `data.frame`.
* **Old**
    * Columns: `Proportion intervened`, `Mean`, `Variance`, `lb`, `ub` (CI around mean PUTA per member).

---
## 9. Edge cases & defaults

* **New**
    * Empty selections return structured zeroed summaries (`total_contacts = 0`).
    * Uses empirical quantiles (10% / 90% for PUTA/contact; 10% / 90% or 1% / 99% for PIA depending on strategy).
* **Old**
    * Empty sets could yield undefined means/variances.
    * Confidence intervals via normal approximation.

---
## 10. Backward-incompatible changes (attention)

* Shift from per-member means/variances to per-contact efficiency + totals.
* Broader PIA/PUTA population in distance–size (adds donors of cluster recipients) → larger effects possible.
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

- **`src/intervention0.R`** – Core intervention analysis functions that compute PUTA, PIA, and efficiency metrics for each strategy
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

The mechanism analysis uses parameters matching `intervention0.R`:
- `implementation_delay_days = 14`
- `analysis_delay_days = 14`
- `rita_window_months = 6` (exponential distribution with mean 180 days)
- `distance_threshold_growth = 0.01` (for growth clusters)
- `distance_threshold_distsize = 0.005` (for size-based clusters)
- `lookback_window_months = 6` (for growth trigger)

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
