# Re-run intervention analysis with simid preservation
# This script re-runs the full analysis after modifying intervention functions
# to preserve simulation IDs for paired comparison

library(here)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                                                                  ║\n")
cat("║   Re-running intervention analysis with simid preservation      ║\n")
cat("║                                                                  ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Source the analysis script
source(here::here("src", "run_analysis.R"))

# Run full analysis (not using cache)
cat("This will take several minutes...\n\n")
results <- run_full_analysis(use_cached = FALSE,
                             run_mechanism = FALSE)  # Skip mechanism plot for speed

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║   Analysis complete! Results saved with simid preserved.         ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")
cat("\n")
