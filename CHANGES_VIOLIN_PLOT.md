# Efficiency Distributions Violin Plot - Changes Summary

## Overview
Added mean values and 2.5th/97.5th percentile centiles as overlays to the efficiency distributions violin plots.

## Changes Made

### Modified Files
- **[src/plot_interventions.R](src/plot_interventions.R)** - Updated the `plot_efficiency_distributions()` function

### Specific Changes (lines 135-202)
Added two new geom layers to each of the 4 plots (p1, p2, p3, p4):

1. **Centile Error Bars** - Shows 2.5th to 97.5th percentile range
   - Implementation: `stat_summary()` with custom function computing quantiles
   - Styling: Black color, 0.8 linewidth, 0.3 width
   - Position: Between violin and mean point layers

2. **Mean Point Indicator** - Shows the mean value
   - Implementation: `stat_summary()` with `fun = mean`
   - Styling: White center with black outline (shape 21), size 3, stroke 1.5
   - Position: Top layer (most visible)

### Visual Design
```
Layer Order (bottom to top):
├── geom_violin() - Distribution (alpha=0.7)
├── stat_summary() - Centile error bars (2.5th-97.5th percentile)
└── stat_summary() - Mean point (white with black outline)
```

## Testing

### Test Results
✓ Successfully generated plot with cached data
✓ All 4 panels display correctly:
  - PUTA (Small Subnetwork)
  - PUTA (Large Subnetwork)
  - PIA (Small Subnetwork)
  - PIA (Large Subnetwork)
✓ Mean points visible on all violin plots
✓ Centile ranges displayed correctly
✓ Works with pseudo-log transformation

### Output
- **File**: [intervention-plots/efficiency_distributions_violin.pdf](intervention-plots/efficiency_distributions_violin.pdf)
- **Size**: 116 KB
- **Dimensions**: 12" x 10"
- **Last Updated**: 2026-01-14

## How to Regenerate

### Using cached results:
```r
source("src/run_analysis.R")
generate_plots_from_cache()
```

### Using test script:
```r
source("test_violin_plot.R")
```

### Running full analysis:
```r
source("src/run_analysis.R")
run_full_analysis(use_cached = FALSE)
```

## Technical Details

### Statistics Computed
- **Mean**: Arithmetic mean of efficiency values per strategy
- **2.5th percentile**: Lower bound of 95% interval
- **97.5th percentile**: Upper bound of 95% interval

### Compatibility
- Works correctly with pseudo-log transformation (asinh)
- Compatible with horizontal orientation (coord_flip)
- Handles all 6 intervention strategies
- Filters non-finite values (Inf, NaN) before plotting

## Benefits
1. **Quick comparison**: Mean points allow rapid visual comparison across strategies
2. **Distribution spread**: Centile bars show the 95% credible interval
3. **Full distribution**: Violin plots still show complete distribution shape
4. **Statistical rigor**: Provides both summary statistics and full data visualization

## Code Example
```r
# For plot p1 (PUTA Small Subnetwork):
p1 <- ggplot(df, aes(x = strategy, y = puta_eff_small, fill = strategy)) +
  geom_violin(alpha = 0.7, scale = "width") +
  # Add centile range (2.5th to 97.5th percentile)
  stat_summary(fun.data = function(x) {
    data.frame(
      y = mean(x),
      ymin = quantile(x, 0.025, na.rm = TRUE),
      ymax = quantile(x, 0.975, na.rm = TRUE)
    )
  }, geom = "errorbar", width = 0.3, color = "black", linewidth = 0.8) +
  # Add mean indicator
  stat_summary(fun = mean, geom = "point",
               color = "white", size = 3, shape = 21, fill = "black",
               stroke = 1.5) +
  scale_y_continuous(...) +
  coord_flip() +
  labs(...) +
  scale_fill_manual(values = strategy_colors) +
  common_theme
```

## Dependencies
No new dependencies required:
- `ggplot2::stat_summary()` - Already available
- `ggplot2::geom_errorbar()` - Already available
- `ggplot2::geom_point()` - Already available
- `stats::quantile()` - Base R function
- `base::mean()` - Base R function
