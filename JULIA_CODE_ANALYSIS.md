# Julia Simulation Code Analysis

## Executive Summary

The Julia code (`j0.jl`, `generate_experiment.jl`) simulates HIV transmission dynamics using a branching process model with detailed contact networks, disease progression, and diagnostic/sequencing delays. It generates two CSV files (G.csv and D.csv) consumed by downstream R analysis.

**Status**: Code matches README documentation. No critical issues found, but see recommendations below.

---

## What is Being Computed and Saved

### Core Simulation: `simbp()` in j0.jl

The simulation models HIV transmission through a **branching process** with:

1. **Infection dynamics** (Stochastic jump process)
   - Disease progression: EHI → Chronic → AIDS → Death
   - Dynamic contact networks with 3 partner types
   - Time-varying transmission probability (high during acute infection)
   - Diagnosis, sequencing, and care transitions

2. **Contact network types**:
   - **F links** (frequent partners): Stable partnerships, gain/loss rate = 1/7/30 per day
   - **G links** (group partners): Semi-stable, gain/loss rate = 1/2/30 per day
   - **H links** (one-off contacts): Gamma-distributed contact rate
   - Degrees sampled from empirical distributions; excess degree used for infected individuals

3. **Transmission mechanics**:
   - Transmission probability: τ(t) = T × ρ × [τ₁ + (τ₀-τ₁) × exp(-ν×t)]
   - Higher during acute infection (τ₀ = 0.001171), decays exponentially (half-life ~6 months)
   - Reduced by diagnosis (ρ = 0.5) and suppression (ρ = 0)

4. **Genetic distance**:
   - Simulated using relaxed molecular clock: NegativeBinomial(μ×t/ω, 1/(1+ω))
   - μ = 0.001/365 (mean clock rate per day)
   - ω = 0.5 (variance inflation)
   - Distance = genetic divergence between donor's sequence time and recipient's sequence time

5. **Contact tracing windows**:
   - Tracks contacts in 90-day and 180-day windows before diagnosis
   - Separates F, G, H contact counts for network interventions

### Output Files

#### **G.csv** (Individual-level data)
Each row = one infection

| Column | Description | Example |
|--------|-------------|---------|
| `pid` | Person ID (hierarchical: "0", "0.1", "0.1.1") | "0.1" |
| `timesequenced` | Days until sequence available | 4045.33 |
| `timediagnosed` | Days until diagnosis | 4044.56 |
| `timeinfected` | Days when infected (0 for index) | 595.79 |
| `generation` | Transmission generation (0=index) | 1 |
| `Fdegree` | Number of F partners at infection | 0 |
| `Gdegree` | Number of G partners at infection | 3 |
| `Hdegree` | One-off contact rate (continuous) | 0.00236 |
| `Fcontacts_90d` | F contacts in 90 days before diagnosis | 0 |
| `Gcontacts_90d` | G contacts in 90 days before diagnosis | 9 |
| `Hcontacts_90d` | H contacts in 90 days before diagnosis | 0 |
| `Fcontacts_180d` | F contacts in 180 days before diagnosis | 0 |
| `Gcontacts_180d` | G contacts in 180 days before diagnosis | 12 |
| `Hcontacts_180d` | H contacts in 180 days before diagnosis | 1 |
| `simid` | Simulation ID | "1" |

**Note**: `timesequenced = Inf` if never sequenced (e.g., died/lost before diagnosis or sequencing)

#### **D.csv** (Transmission pairs)
Each row = one transmission event

| Column | Description | Example |
|--------|-------------|---------|
| `donor` | Donor person ID | "0" |
| `recipient` | Recipient person ID | "0.1" |
| `distance` | Genetic distance | Inf (if donor/recipient never sequenced) |
| `timetransmission` | Days when transmission occurred | 595.79 |
| `contacttype` | Link type that transmitted (:F, :G, :H) | :G |
| `simid` | Simulation ID | "1" |

**Important**: `distance = Inf` if either donor or recipient never gets sequenced. This is critical for cluster detection (only sequenced pairs have finite distances).

---

## README Consistency Check

### ✅ What Matches

1. **File naming convention**: ✅ Correct
   - README states: `experiment1-N<n_sims>-gens<maxgenerations>-G.csv` and `-D.csv`
   - Code produces: `experiment1-N10-gens5-G.csv` (line 62 in generate_experiment.jl)

2. **CLI usage**: ✅ Correct
   - README: `julia --project src/generate_experiment.jl 100 experiment1 5`
   - Code: Parses ARGS[1]=n_sims, ARGS[2]=prefix, ARGS[3]=maxgenerations (lines 77-82)

3. **Quick rerun command**: ✅ Correct
   - README: `julia --project src/run_generate_experiment.jl`
   - File exists and calls `generate_experiment(10; out_prefix="experiment1", maxgenerations=5)`

4. **Function signature**: ✅ Correct
   - README mentions `initialcontact` parameter
   - Code supports `:G`, `:F`, `:H` options (line 26, used in line 402)

5. **Output location**: ✅ Correct
   - README: Files written to `src/`
   - Code: `out_dir = @__DIR__` writes to src/ (line 63)

6. **Default parameters**: ✅ Correct
   - README: Default maxgenerations=5
   - Code: Line 80 defaults to 5

### ⚠️ Discrepancies / Clarifications Needed

**NONE FOUND** - README accurately describes the Julia code behavior!

---

## Potential Issues and Recommendations

### 🟡 Minor Issues

#### 1. **Contact window calculation may over-count initial partners**
**Location**: `j0.jl` lines 240-241, 244-245

**Issue**: Initial partners (u₀[1], u₀[2]) are added to contact counts unconditionally:
```julia
fcontacts_90d = u₀[1] + count(t -> t >= t_90 && t <= int.t, fgain_times)
```

**Problem**: Initial partners present at infection (time=tinf) might have been lost by the time of diagnosis (especially if diagnosis is >6 months later, given partner turnover rates). The code assumes all initial partners are still "active" for contact tracing purposes.

**Impact**: May slightly overestimate contact counts for contact tracing, particularly for F/G contacts in individuals with long delays to diagnosis.

**Recommendation**: Track partner loss times or use a stochastic approach to determine which initial partners are still present at diagnosis time.

---

#### 2. **Hdegree represents rate, not count**
**Location**: `j0.jl` line 84, output in G.csv

**Current behavior**: `Hdegree` column in G.csv is the continuous rate (e.g., 0.00236), not a count.

**Potential confusion**: Users might expect integer counts like Fdegree/Gdegree. The column name "degree" is somewhat misleading since it's actually a rate parameter.

**Impact**: Low - downstream R code correctly interprets this as a rate. However, could confuse users examining raw CSV files.

**Recommendation**:
- Rename to `Hrate` in output or add documentation
- OR compute actual expected H partners over some time window (e.g., 90 days) for consistency

---

#### 3. **Missing validation: maxgenerations vs SIMDURATION**
**Location**: `generate_experiment.jl` line 29, `j0.jl` line 29

**Issue**: No check that SIMDURATION (20 years) is sufficient for requested number of generations.

**Scenario**: If transmission chains are slow (long intervals between infections), late generations might not complete before SIMDURATION ends, leading to truncated chains.

**Impact**: With current parameters, unlikely to be a problem for ≤7 generations. Could affect very deep trees (>10 generations).

**Recommendation**: Add assertion or warning if maxgenerations > reasonable threshold (e.g., 10).

---

#### 4. **Sequential simid overwrite in generate_experiment.jl**
**Location**: `generate_experiment.jl` lines 49-50

**Current behavior**: Overwrites UUID simids from `simbp()` with sequential strings "1", "2", ..., "N"

**Code**:
```julia
G_df[!, :simid] = fill(string(i), nrow(G_df))
D_df[!, :simid] = fill(string(i), nrow(D_df))
```

**Reasoning**: Human-readable IDs are better for analysis than UUIDs.

**Potential issue**: If batches are merged from different runs, simids could collide.

**Impact**: Low - as long as users don't manually merge CSVs from different runs.

**Recommendation**: Consider adding run timestamp or prefix to simid if merging is ever needed.

---

#### 5. **No seed control in generate_experiment**
**Location**: `generate_experiment.jl` does not set Random.seed!

**Issue**: Simulations are not reproducible unless Julia is launched with same RNG state.

**Impact**: Cannot exactly reproduce specific simulation batches, which is important for debugging or validation.

**Recommendation**: Add optional `seed` parameter:
```julia
function generate_experiment(n_sims::Int;
    seed::Union{Nothing,Int} = nothing,
    ...
)
    if !isnothing(seed)
        Random.seed!(seed)
        println("Set RNG seed: $seed")
    end
    ...
end
```

---

### 🟢 Good Practices Observed

1. **Progress reporting**: Every 1000 sims with ETA (lines 40-44)
2. **Error handling**: Try-catch for Pkg.activate (lines 11-16)
3. **Flexible CLI**: Supports positional args with sensible defaults
4. **Documented parameters**: Clear comments in P struct
5. **Output validation**: Prints row counts and file paths on completion

---

## Key Parameters Used in Simulations

Based on `j0.jl` lines 51-78:

| Parameter | Value | Meaning | Source |
|-----------|-------|---------|--------|
| **Transmission** |
| τ₀ | 0.001171 | Initial transmission prob per act (acute) | O4Y4epigen |
| τ₁ | 0.0002647 | Final transmission prob (chronic) | O4Y4epigen |
| T | 4.2997 | Transmission probability scale factor | - |
| ν | 1.386 yr⁻¹ | Decay rate (t½ = 6 months) | - |
| ρ_diagnosed | 0.50 | Transmission reduction when diagnosed | - |
| ρ_suppressed | 0.0 | Transmission reduction when suppressed | - |
| **Contact rates** |
| fcont | 1.5/7 d⁻¹ | F contact rate per partner | - |
| gcont | 1.0/7 d⁻¹ | G contact rate per partner | - |
| frate | 1/(7×30) d⁻¹ | F partner turnover rate (~weekly/monthly) | - |
| grate | 1/(2×30) d⁻¹ | G partner turnover rate (~bimonthly) | - |
| **Disease progression** |
| γ_ehi | 1/365 yr⁻¹ | EHI duration (~1 year) | Volz PMED 2013 |
| γ_chron | 2.1 yr | Chronic stage scale (Gamma) | Volz PMED 2013 |
| shape_chron | 3.0 | Chronic stage shape parameter | - |
| γ_aids | 1/(2×365) yr⁻¹ | AIDS → death (~2 years) | - |
| **Diagnosis & care** |
| δ | 0.00063 d⁻¹ | Diagnosis rate (~4 years) | O4Y4epigen |
| α | 0.0262 d⁻¹ | Sequencing rate (~38 days) | O4Y4epigen |
| κ | 0.004736 d⁻¹ | Care/suppression rate (~211 days) | O4Y4epigen |
| **Genetic evolution** |
| μ | 0.001/365 d⁻¹ | Mean substitution rate | - |
| ω | 0.5 | Variance inflation for relaxed clock | - |

---

## Usage Verification

### Current Analysis Uses:
Looking at R code in `src/run_analysis.R`:
- Default files: `src/experiment1-N10000-gens7-D.csv` and `-G.csv`
- ✅ These files exist (from `ls` output earlier)
- ✅ Analysis pipeline correctly references these paths

### README Default vs Actual:
- README quick-start example: Runs with 10 sims, 5 generations
- Actual analysis data: Uses 10,000 sims, 7 generations
- **This is fine** - README shows minimal example for testing; actual analysis uses larger batch

---

## Recommendations for README

### Suggested Additions:

1. **Document column meanings** more explicitly:
   ```markdown
   ### Output File Formats

   **G.csv columns**:
   - `Fdegree`, `Gdegree`: Partner counts at infection time
   - `Hdegree`: One-off contact *rate* (not count)
   - `timesequenced = Inf`: Individual never sequenced
   - `*contacts_90d`: Contacts in 90 days *before diagnosis*
   ```

2. **Add reproducibility note**:
   ```markdown
   ### Reproducibility

   Note: Simulations use Julia's default RNG and are not currently seeded.
   For reproducible results, set `Random.seed!(123)` before calling
   `generate_experiment()` in REPL mode.
   ```

3. **Clarify existing vs example data**:
   ```markdown
   ### Current Analysis Data

   The analysis currently uses:
   - `experiment1-N10000-gens7-*.csv` (10,000 simulations, 7 generations)

   The quick-start examples generate smaller test datasets (10 sims, 5 gens)
   for verification purposes.
   ```

---

## Summary

**Overall Assessment**: ✅ **Code is well-structured and matches documentation**

**Critical Issues**: None found

**Minor Issues**: 5 minor issues identified (see above), none blocking

**Alignment with README**: Excellent - all documented behaviors match implementation

**Recommendation**: Code is production-ready. Consider addressing minor issues (#1, #5) if:
- Exact reproducibility is required (add seed parameter)
- Contact tracing accuracy is critical (fix initial partner counting)

The simulation framework is sophisticated and appropriate for HIV cluster-based intervention analysis.
