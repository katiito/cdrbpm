# =============================================================================
# analyze_paired_comparison.R
# =============================================================================
# Standalone script for paired comparison of strategies vs random baseline
#
# This is a convenience wrapper that calls plot_paired_comparisons() from
# plot_interventions.R. For full analysis including all plots, use:
#   source("src/run_analysis.R")
#   run_full_analysis(use_cached = TRUE)
#
# Usage:
#   source("src/analyze_paired_comparison.R")
# =============================================================================

library(here)

# Source required scripts
source(here::here("src", "run_analysis.R"))

cat("Running paired comparison analysis (standalone)...\n\n")

# Load cached results and run paired comparison
results <- load_cached_results()
plots <- plot_paired_comparisons(results)

cat("\nPaired comparison analysis complete!\n")
cat("\nOutputs:\n")
cat("  - intervention-results/paired_comparison_results.csv\n")
cat("  - intervention-plots/paired_comparison_distributions.png\n")
cat("  - intervention-plots/paired_comparison_percent.png\n")
cat("\nNote: For all plots together, use run_full_analysis(use_cached = TRUE)\n")
