# Test script for updated efficiency distributions violin plot
# This script tests that the mean and centile overlays work correctly

cat("Testing updated efficiency_distributions_violin plot...\n\n")

# Source the required functions
source("src/intervention0.R")
source("src/plot_interventions.R")
source("src/run_analysis.R")

# Load cached results
cat("Loading cached results from most recent run...\n")
tryCatch({
  results <- load_cached_results()

  cat("\nGenerating efficiency distributions plot with mean and centile overlays...\n")
  p <- plot_efficiency_distributions(results)

  # Save the plot
  output_file <- "intervention-plots/efficiency_distributions_violin.pdf"
  if (!dir.exists("intervention-plots")) {
    dir.create("intervention-plots", recursive = TRUE)
  }

  cat(sprintf("Saving plot to %s...\n", output_file))
  ggsave(output_file, p, width = 12, height = 10)

  cat("\n=================================================\n")
  cat("SUCCESS! Plot saved with the following features:\n")
  cat("=================================================\n")
  cat("  ✓ Violin plots showing full efficiency distributions\n")
  cat("  ✓ Black error bars showing 2.5th to 97.5th percentiles\n")
  cat("  ✓ White points with black outline showing mean values\n")
  cat("\nPlease visually inspect the plot to verify:\n")
  cat("  1. Mean points are clearly visible on all 4 panels\n")
  cat("  2. Centile error bars span from 2.5th to 97.5th percentile\n")
  cat("  3. Overlays work correctly with the pseudo-log scale\n")
  cat("  4. All elements are visible against the violin colors\n")
  cat("\nPlot location: ", output_file, "\n")

}, error = function(e) {
  cat("\nERROR: Could not load cached results.\n")
  cat("Error message:", conditionMessage(e), "\n\n")
  cat("To generate new results, run:\n")
  cat("  source('src/run_analysis.R')\n")
  cat("  run_full_analysis(use_cached = FALSE)\n")
})
