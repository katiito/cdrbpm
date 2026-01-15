# Parameter Inheritance Fix - Implementation Summary

## Date: 2026-01-15

Fixed critical issue where mechanism analysis did not inherit parameters from intervention analysis, causing potential mismatches in visualization.

---

## Problem Identified

The user correctly identified that `run_mechanism_analysis()` in the pipeline was not inheriting parameters from the intervention analysis that was run before it. This meant:

1. **Intervention analysis** could be run with custom parameters (e.g., `partner_notification_window_months = 3`, `network_degree_threshold = 8`)
2. **Mechanism analysis** would always use defaults (6 months, threshold 4)
3. **Result**: Mechanism plots would not accurately represent the intervention being analyzed

### Example of the Problem

```r
# Run intervention with 3-month window and threshold 8
results <- run_intervention_analysis(
  partner_notification_window_months = 3,
  network_degree_threshold = 8
)

# Generate mechanism plots - WOULD USE DEFAULTS (6 months, threshold 4)
run_mechanism_analysis(n_sims = 100)  # ❌ Parameter mismatch!
```

---

## Solution Implemented

### 1. Added Parameters to Return Value

**File**: `src/intervention0.R` lines 527-552

**Change**: Added `parameters` component to the list returned by `run_intervention_analysis()`

**Before**:
```r
invisible(list(
  summary = odf,
  counts = counts_df,
  details = list(
    distsize5 = ods5,
    distsize2 = ods2,
    growth = ogrowth,
    random = orand,
    rita = orita,
    network = onet
  )
))
```

**After**:
```r
invisible(list(
  summary = odf,
  counts = counts_df,
  details = list(
    distsize5 = ods5,
    distsize2 = ods2,
    growth = ogrowth,
    random = orand,
    rita = orita,
    network = onet
  ),
  parameters = list(
    partner_notification_window_months = partner_notification_window_months,
    network_degree_threshold = network_degree_threshold,
    distance_threshold = distance_threshold,
    growth_distance_threshold = growth_distance_threshold,
    cluster_size_5 = cluster_size_5,
    cluster_size_2 = cluster_size_2,
    random_sample_size = random_sample_size,
    rita_window_months = rita_window_months,
    lookback_window_months = lookback_window_months,
    analysis_delay_days = analysis_delay_days,
    implementation_delay_days = implementation_delay_days,
    seed = seed
  )
))
```

---

### 2. Updated Parameter Extraction in Plotting Pipeline

**File**: `src/run_analysis.R` lines 144-176

**Change**: `generate_plots_from_cache()` now extracts parameters from results and passes them to mechanism analysis

**Added Logic**:
```r
# Extract parameters from results if available, otherwise use defaults
partner_notification_window_months <- 6  # default
network_degree_threshold <- 4  # default

if (!is.null(results$parameters)) {
  if (!is.null(results$parameters$partner_notification_window_months)) {
    partner_notification_window_months <- results$parameters$partner_notification_window_months
  }
  if (!is.null(results$parameters$network_degree_threshold)) {
    network_degree_threshold <- results$parameters$network_degree_threshold
  }
  cat(sprintf("  Using parameters from intervention analysis: partner_notification_window_months=%d, network_degree_threshold=%d\n",
              partner_notification_window_months, network_degree_threshold))
} else {
  cat("  Warning: No parameter information found in results, using defaults\n")
}

p_mechanism <- run_mechanism_analysis(
  partner_notification_window_months = partner_notification_window_months,
  network_degree_threshold = network_degree_threshold,
  n_sims = n_sims
)
```

---

### 3. Updated Cached Results Loading

**File**: `src/run_analysis.R` lines 105-128

**Change**: `load_cached_results()` now loads parameters from the saved `parameters_*.csv` file

**Added**:
```r
# Load parameters if available
params_file <- find_most_recent("^parameters_")
parameters <- NULL
if (!is.null(params_file) && file.exists(params_file)) {
  params_df <- read.csv(params_file)
  # Convert to named list
  parameters <- setNames(as.list(params_df$value), params_df$parameter)
  # Convert numeric strings back to numeric where appropriate
  numeric_params <- c("cluster_size_5", "cluster_size_2", "distance_threshold",
                      "network_degree_threshold", "random_sample_size",
                      "rita_window_months", "lookback_window_months",
                      "growth_distance_threshold", "analysis_delay_days",
                      "implementation_delay_days", "partner_notification_window_months")
  for (param in numeric_params) {
    if (!is.null(parameters[[param]])) {
      parameters[[param]] <- as.numeric(parameters[[param]])
    }
  }
  cat(sprintf("  Loaded: %s\n", basename(params_file)))
} else {
  cat("  Warning: parameters file not found\n")
}
```

---

## How It Works Now

### Fresh Analysis Run

```r
# Run intervention with custom parameters
results <- run_intervention_analysis(
  partner_notification_window_months = 3,
  network_degree_threshold = 8
)

# Parameters are stored in results$parameters
# results$parameters$partner_notification_window_months = 3
# results$parameters$network_degree_threshold = 8

# Generate plots - automatically inherits parameters
run_full_analysis(use_cached = FALSE)
# ✅ Mechanism analysis uses 3-month window and threshold 8
```

### Cached Analysis Run

```r
# Load cached results
run_full_analysis(use_cached = TRUE)

# Parameters loaded from parameters_TIMESTAMP.csv
# Automatically extracted and passed to mechanism analysis
# ✅ Mechanism analysis matches the original intervention parameters
```

---

## Benefits

### 1. **Parameter Consistency**
- Mechanism analysis always uses same parameters as intervention analysis
- No more manual parameter passing required
- Eliminates risk of parameter mismatch

### 2. **Automatic Inheritance**
- Parameters automatically flow from intervention → mechanism analysis
- Works for both fresh runs and cached results
- User doesn't need to remember to pass parameters

### 3. **Clear Feedback**
- Console output confirms which parameters are being used
- Warns if parameters can't be loaded (fallback to defaults)
- Transparent parameter flow

### 4. **Backward Compatibility**
- If `results$parameters` is NULL, defaults are used
- Graceful degradation for old cached results
- No breaking changes to existing code

---

## Testing

### Test Case 1: Fresh Run with Custom Parameters

```r
source("src/run_analysis.R")

# Run with custom parameters
run_full_analysis(
  use_cached = FALSE,
  partner_notification_window_months = 3,
  network_degree_threshold = 8,
  n_sims = 100
)

# Expected: Console output shows:
# "Using parameters from intervention analysis: partner_notification_window_months=3, network_degree_threshold=8"
# Mechanism plot Panel A uses threshold 8
# Network strategy uses 90-day contact window
```

### Test Case 2: Cached Results

```r
source("src/run_analysis.R")

# Run fresh to create cache
run_full_analysis(
  use_cached = FALSE,
  partner_notification_window_months = 3,
  network_degree_threshold = 6
)

# Later: Load cached results
run_full_analysis(use_cached = TRUE, n_sims = 100)

# Expected: Parameters loaded from parameters_*.csv
# Mechanism analysis uses same parameters as original run
```

### Test Case 3: Old Cached Results (No Parameters)

```r
# Simulate old results without parameters component
results_old <- list(
  counts = data.frame(),
  details = list()
  # No parameters component
)

generate_plots_from_cache(results = results_old, run_mechanism = TRUE)

# Expected: Warning message + defaults used
# "Warning: No parameter information found in results, using defaults"
```

---

## Files Modified

1. **src/intervention0.R**
   - Lines 527-552: Added `parameters` component to return value

2. **src/run_analysis.R**
   - Lines 105-128: Updated `load_cached_results()` to load parameters from CSV
   - Lines 144-176: Updated `generate_plots_from_cache()` to extract and pass parameters

---

## Verification Checklist

- [x] Parameters included in `run_intervention_analysis()` return value
- [x] `load_cached_results()` loads parameters from CSV file
- [x] `generate_plots_from_cache()` extracts parameters from results
- [x] Parameters passed to `run_mechanism_analysis()`
- [x] Console output confirms parameter usage
- [x] Graceful fallback to defaults if parameters unavailable
- [x] Works for both fresh runs and cached results
- [x] Backward compatible with old cached results

---

## Impact

**Before**:
- Mechanism analysis always used defaults (6 months, threshold 4)
- Could show different network strategy behavior than actual intervention
- Required manual parameter passing to match intervention settings
- Risk of parameter mismatch and confusing results

**After**:
- Mechanism analysis automatically inherits intervention parameters
- Perfect consistency between intervention and visualization
- No manual parameter management needed
- Clear console feedback on parameter usage
- Cached results preserve original parameters

This fix ensures that mechanism analysis plots always accurately represent the intervention analysis being visualized, eliminating a major source of potential confusion and ensuring scientific accuracy.
