# Plotting Code Fixes - Implementation Summary

## Date: 2026-01-15

All recommended changes from the plotting code analysis have been implemented to ensure parameter consistency between intervention analysis and mechanism visualization.

---

## Summary of Changes

### Files Modified
1. `src/plot_interventions.R` - Fixed parameter inconsistencies and added validation
2. `README.md` - Updated to reference parameter names instead of hardcoded values

### Issues Fixed
- ✅ Fix 1: Network degree threshold now uses parameter instead of hardcoded value
- ✅ Fix 2: Contact window now inherits from partner_notification_window_months
- ✅ Fix 3: Updated misleading comments
- ✅ Fix 4: Added input validation for plot_efficiency_distributions
- ✅ Fix 5: Updated README to use parameter names

---

## Detailed Changes

### Fix 1: Network Degree Threshold Parameter

**Location**: `src/plot_interventions.R`

**Issue**: Code used hardcoded threshold of 8, but intervention0.R default is 4

**Before** (line 444):
```r
# Network (high degree) strategy - top decile of contacts (>=12)
high_degree <- G1_eligible_sorted[G1_eligible_sorted$total_contacts >= 8, ]
```

**After**:
```r
# Network (high degree) strategy - individuals with >= network_degree_threshold contacts
high_degree <- G1_eligible_sorted[G1_eligible_sorted$total_contacts >= network_degree_threshold, ]
```

**Added parameter** (line 258-267):
```r
plot_mechanism_analysis <- function(D, G,
                                    distance_threshold_distsize = 0.005,
                                    distance_threshold_growth = 0.01,
                                    partner_notification_window_months = 6,
                                    network_degree_threshold = 4,
                                    n_sims = 100) {
```

**Impact**: Mechanism analysis now correctly uses the same threshold (4) as intervention analysis.

---

### Fix 2: Contact Window Parameter

**Location**: `src/plot_interventions.R`

**Issue**: Code used hardcoded 90-day window, but intervention0.R default is 180 days

**Before** (lines 441-443):
```r
G1_eligible_sorted$total_contacts <- G1_eligible_sorted$Fcontacts_90d +
                                      G1_eligible_sorted$Gcontacts_90d +
                                      G1_eligible_sorted$Hcontacts_90d
```

**After** (lines 287-288, 443-445):
```r
# Determine contact window based on partner notification window
contact_suffix <- if (partner_notification_window_months == 3) "90d" else "180d"

...

G1_eligible_sorted$total_contacts <- G1_eligible_sorted[[paste0("Fcontacts_", contact_suffix)]] +
                                      G1_eligible_sorted[[paste0("Gcontacts_", contact_suffix)]] +
                                      G1_eligible_sorted[[paste0("Hcontacts_", contact_suffix)]]
```

**Impact**: Contact counts now dynamically match the partner notification window used in intervention analysis (90d for 3-month, 180d for 6-month).

---

### Fix 3: Updated Comments

**Location**: `src/plot_interventions.R` line 439

**Before**:
```r
# Network (high degree) strategy - top decile of contacts (>=12)
```

**After**:
```r
# Network (high degree) strategy - individuals with >= network_degree_threshold contacts
```

**Impact**: Comment now accurately describes the code logic and references the parameter name.

---

### Fix 4: Input Validation

**Location**: `src/plot_interventions.R` lines 31-34

**Added**:
```r
# Validate input
if (is.null(results) || !is.list(results) || is.null(results$details)) {
  stop("Invalid results object. Must contain 'details' component from run_intervention_analysis().")
}
```

**Impact**: Provides clear error message if invalid results object is passed, preventing cryptic errors downstream.

---

### Fix 5: Updated run_mechanism_analysis

**Location**: `src/plot_interventions.R` lines 994-1019

**Before**:
```r
run_mechanism_analysis <- function(D_path = "src/experiment1-N10000-gens7-D.csv",
                                   G_path = "src/experiment1-N10000-gens7-G.csv",
                                   save_path = "intervention-plots/mechanism_analysis.png",
                                   n_sims = 100,
                                   width = 12,
                                   height = 12) {
  ...
  p <- plot_mechanism_analysis(D, G, n_sims = n_sims)
```

**After**:
```r
run_mechanism_analysis <- function(D_path = "src/experiment1-N10000-gens7-D.csv",
                                   G_path = "src/experiment1-N10000-gens7-G.csv",
                                   save_path = "intervention-plots/mechanism_analysis.png",
                                   partner_notification_window_months = 6,
                                   network_degree_threshold = 4,
                                   n_sims = 100,
                                   width = 12,
                                   height = 12) {
  ...
  p <- plot_mechanism_analysis(D, G,
                                partner_notification_window_months = partner_notification_window_months,
                                network_degree_threshold = network_degree_threshold,
                                n_sims = n_sims)
```

**Impact**: Convenience function now accepts and passes through the new parameters.

---

### Fix 6: Enhanced Documentation

**Location**: `src/plot_interventions.R` lines 253-260

**Before**:
```r
#' @param n_sims Number of simulations to analyze (default 100)
```

**After**:
```r
#' @param partner_notification_window_months Partner notification window (default 6, matching intervention0.R)
#' @param network_degree_threshold Network degree threshold (default 4, matching intervention0.R)
#' @param n_sims Number of simulations to analyze (default 100). For robust Panel D (survivorship),
#'   use n_sims >= 1000 or NULL (all simulations) since growth clusters trigger in ~0.14% of simulations.
```

**Impact**: Better documentation of parameters and expectations for Panel D (survivorship bias analysis).

---

### Fix 7: README Updates

**Location**: `README.md` lines 303-314

**Before**:
```markdown
### Parameters

The mechanism analysis uses parameters matching `intervention0.R`:
- `implementation_delay_days = 14`
- `analysis_delay_days = 14`
- `rita_window_months = 6` (exponential distribution with mean 180 days)
- `distance_threshold_growth = 0.01` (for growth clusters)
- `distance_threshold_distsize = 0.005` (for size-based clusters)
- `lookback_window_months = 6` (for growth trigger)
```

**After**:
```markdown
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
```

**Impact**: Documentation now emphasizes that values are defaults (not hardcoded) and can be customized.

---

## Verification

### Manual Testing Steps

1. **Test with default parameters** (should match intervention0.R defaults):
   ```r
   source("src/plot_interventions.R")
   run_mechanism_analysis()
   ```

2. **Test with custom parameters**:
   ```r
   run_mechanism_analysis(
     partner_notification_window_months = 3,  # Use 90-day window
     network_degree_threshold = 8,             # Use higher threshold
     n_sims = 100
   )
   ```

3. **Verify input validation**:
   ```r
   plot_efficiency_distributions(NULL)  # Should error with clear message
   plot_efficiency_distributions(list())  # Should error with clear message
   ```

### Expected Behavior

- ✅ Mechanism analysis uses same parameters as intervention analysis by default
- ✅ Custom parameters correctly propagate through function calls
- ✅ Contact window automatically switches between 90d and 180d based on partner_notification_window_months
- ✅ Network strategy uses consistent degree threshold
- ✅ Invalid inputs produce clear error messages

---

## Parameter Consistency Check

| Parameter | intervention0.R | plot_mechanism_analysis | Match? |
|-----------|-----------------|-------------------------|--------|
| `distance_threshold_distsize` | 0.005 | 0.005 (default) | ✅ |
| `distance_threshold_growth` | 0.01 | 0.01 (default) | ✅ |
| `implementation_delay_days` | 14 | 14 (hardcoded) | ✅ |
| `analysis_delay_days` | 14 | 14 (hardcoded) | ✅ |
| `lookback_window_months` | 6 | 6 (hardcoded) | ✅ |
| `rita_window_months` | 6 | 6 (hardcoded) | ✅ |
| `partner_notification_window_months` | 6 (180d) | 6 (parameter, default) | ✅ |
| `network_degree_threshold` | 4 | 4 (parameter, default) | ✅ |

**Status**: ✅ All parameters now consistent between intervention analysis and mechanism visualization.

---

## Remaining Minor Issues (Low Priority)

These were identified but not critical to fix:

1. **Progress reporting for small n_sims** (line 302):
   - Could set `progress_interval <- max(10, floor(length(simids) / 10))` to avoid console spam
   - Current behavior: reports every simulation for n_sims < 10

2. **Pseudo-log transformation scales** (lines 109, 117):
   - `puta_scale = 0.5` and `pia_scale = 10` are hardcoded
   - Could be function parameters, but current values work well for existing data

These are quality-of-life improvements that don't affect correctness.

---

## Testing Recommendation

Run full mechanism analysis to verify all changes:

```r
source("src/run_analysis.R")

# Test 1: Default parameters (should match intervention analysis)
run_full_analysis(use_cached = TRUE, run_mechanism = TRUE, n_sims = 100)

# Test 2: Custom parameters
source("src/plot_interventions.R")
run_mechanism_analysis(
  partner_notification_window_months = 3,
  network_degree_threshold = 6,
  n_sims = 100
)
```

Expected output:
- ✅ Mechanism analysis completes without errors
- ✅ Panel A (network strategy) reflects correct threshold
- ✅ Contact counts match the specified window (90d or 180d)
- ✅ All 5 panels render correctly

---

## Verification Checklist

- [x] Network degree threshold uses parameter (not hardcoded 8)
- [x] Contact window dynamically matches partner_notification_window_months
- [x] Comments accurately describe code behavior
- [x] Input validation provides clear error messages
- [x] run_mechanism_analysis passes new parameters
- [x] README references parameter names (not hardcoded values)
- [x] Documentation notes Panel D requires large n_sims
- [x] All parameters match intervention0.R defaults

**Status**: ✅ All critical fixes successfully implemented

---

## Impact Summary

**Before fixes**:
- Mechanism analysis used inconsistent parameters (8 vs 4 for network threshold, 90d vs 180d for contacts)
- Could show different network strategy behavior than actual intervention
- Hardcoded values in README could cause confusion if defaults changed

**After fixes**:
- Mechanism analysis perfectly mirrors intervention analysis parameters
- Parameters are explicit and customizable
- Documentation clearly indicates values are defaults, not hardcoded
- Users can ensure visualizations match their specific analysis settings

The plotting code now maintains perfect parameter consistency with intervention analysis, eliminating potential confusion and ensuring mechanism visualizations accurately represent the actual intervention strategies being analyzed.
