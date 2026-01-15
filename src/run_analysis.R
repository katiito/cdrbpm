# =============================================================================
# Master Analysis Pipeline
# =============================================================================
#
# This script runs the complete intervention analysis pipeline, generating:
#   1. PUTA efficiency violin plot (efficiency_distributions_violin.pdf)
#   2. Mechanism analysis plot (mechanism_analysis.png)
#
# Usage:
#   source("src/run_analysis.R")
#   
#   # Option 1: Run everything fresh
#   run_full_analysis()
#   
#   # Option 2: Use most recent saved results (faster)
#   run_full_analysis(use_cached = TRUE)
#   
#   # Option 3: Just regenerate plots from cached results
#   generate_plots_from_cache()
#
# =============================================================================

library(here)

# Source required scripts
source(here::here("src", "intervention0.R"))
source(here::here("src", "plot_interventions.R"))

#' Find the most recent results file by pattern
#' 
#' @param pattern Glob pattern for files (e.g., "counts_*.csv")
#' @param dir Directory to search in
#' @return Path to most recent file, or NULL if none found
find_most_recent <- function(pattern, dir = here::here("intervention-results")) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  
  # Extract timestamps (YYYYMMDD_HHMMSS) and find most recent
  # Handle both formats: name_TIMESTAMP.csv and name_extra_TIMESTAMP.csv
  timestamps <- sapply(files, function(f) {
    # Extract the last occurrence of YYYYMMDD_HHMMSS pattern before .csv
    m <- regmatches(f, regexpr("\\d{8}_\\d{6}(?=\\.csv$)", f, perl = TRUE))
    if (length(m) == 0) NA_character_ else m
  })
  
  valid_idx <- !is.na(timestamps)
  if (!any(valid_idx)) return(NULL)
  
  most_recent_idx <- which(timestamps == max(timestamps[valid_idx], na.rm = TRUE))[1]
  files[most_recent_idx]
}

#' Load cached intervention results
#' 
#' @return List with results in the same format as run_intervention_analysis()
#' @details Loads the most recent CSV files from intervention-results/
load_cached_results <- function() {
  cat("Looking for cached results...\n")
  
  # Find most recent counts file to get the timestamp
  counts_file <- find_most_recent("^counts_")
  
  if (is.null(counts_file)) {
    stop("No cached results found. Run with use_cached = FALSE first.")
  }
  
  # Extract timestamp for display
  timestamp <- regmatches(counts_file, regexpr("\\d{8}_\\d{6}", counts_file))
  cat(sprintf("  Found results from: %s\n", timestamp))
  
  # Load counts
  counts_df <- read.csv(counts_file)
  cat(sprintf("  Loaded: %s\n", basename(counts_file)))
  
  # Load details from individual strategy files
  details <- list()
  
  # Strategy files and their keys
  strategy_files <- list(
    distsize5 = find_most_recent("^details_distsize5_"),
    distsize2 = find_most_recent("^details_distsize2_"),
    growth = find_most_recent("^details_growth_"),
    random = find_most_recent("^details_random_"),
    rita = find_most_recent("^details_rita_"),
    network = find_most_recent("^details_network_")
  )
  
  for (key in names(strategy_files)) {
    file <- strategy_files[[key]]
    if (!is.null(file) && file.exists(file)) {
      df <- read.csv(file)
      # Normalize column names for individual-based strategies
      # They have 'contacts' instead of 'contacts_small'/'contacts_large'
      if ("contacts" %in% names(df) && !"contacts_small" %in% names(df)) {
        df$contacts_small <- df$contacts
        df$contacts_large <- df$contacts
      }
      details[[key]] <- list(o = df)
      cat(sprintf("  Loaded: %s (%d rows)\n", basename(file), nrow(df)))
    } else {
      cat(sprintf("  Warning: %s details not found\n", key))
    }
  }
  
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

  cat("  Cache loaded successfully!\n\n")

  list(
    counts = counts_df,
    details = details,
    parameters = parameters,
    timestamp = timestamp
  )
}

#' Generate plots from cached results
#' 
#' @param results Optional: pre-loaded results. If NULL, loads from cache.
#' @param run_mechanism Whether to run mechanism analysis (default TRUE)
#' @param n_sims Number of simulations for mechanism analysis (NULL = all)
#' @return Invisible list of plot objects
generate_plots_from_cache <- function(results = NULL, 
                                       run_mechanism = TRUE,
                                       n_sims = NULL) {
  
  plot_dir <- here::here("intervention-plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # Load cached results if not provided
  if (is.null(results)) {
    results <- load_cached_results()
  }
  
  plots <- list()
  
  # -------------------------------------------------------------------------
  # Plot 1: Efficiency distributions (violin plot)
  # -------------------------------------------------------------------------
  cat("Generating efficiency distribution plot...\n")
  p_violin <- plot_efficiency_distributions(results)
  violin_path <- file.path(plot_dir, "efficiency_distributions_violin.pdf")
  ggsave(violin_path, p_violin, width = 12, height = 10)
  cat(sprintf("  Saved: %s\n", violin_path))
  plots$violin <- p_violin
  
  # -------------------------------------------------------------------------
  # Plot 2: Mechanism analysis
  # -------------------------------------------------------------------------
  if (run_mechanism) {
    cat("\nGenerating mechanism analysis plot...\n")
    if (is.null(n_sims)) {
      cat("  Using all simulations (this may take a few minutes)...\n")
    } else {
      cat(sprintf("  Using %d simulations...\n", n_sims))
    }

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
    plots$mechanism <- p_mechanism
  }
  
  cat("\nAll plots saved to: intervention-plots/\n")
  invisible(plots)
}

#' Run the full analysis pipeline
#' 
#' @param use_cached If TRUE, use most recent saved results instead of re-running
#' @param run_mechanism Whether to run mechanism analysis (default TRUE)
#' @param n_sims Number of simulations for mechanism analysis (NULL = all)
#' @param d_file Path to D (transmission) data file
#' @param g_file Path to G (individual) data file
#' @param ... Additional arguments passed to run_intervention_analysis()
#' @return Invisible list with results and plots
run_full_analysis <- function(use_cached = FALSE,
                               run_mechanism = TRUE,
                               n_sims = NULL,
                               d_file = "src/experiment1-N10000-gens7-D.csv",
                               g_file = "src/experiment1-N10000-gens7-G.csv",
                               ...) {
  
  cat("=================================================\n")
  cat("  Intervention Analysis Pipeline\n")
  cat("=================================================\n\n")
  
  # -------------------------------------------------------------------------
  # Step 1: Get intervention results (cached or fresh)
  # -------------------------------------------------------------------------
  if (use_cached) {
    cat("STEP 1: Loading cached intervention results\n")
    cat("-------------------------------------------------\n")
    results <- load_cached_results()
  } else {
    cat("STEP 1: Running intervention analysis (this takes ~2 minutes)\n")
    cat("-------------------------------------------------\n")
    results <- run_intervention_analysis(
      d_file = d_file,
      g_file = g_file,
      ...
    )
  }
  
  # -------------------------------------------------------------------------
  # Step 2: Generate plots
  # -------------------------------------------------------------------------
  cat("\nSTEP 2: Generating plots\n")
  cat("-------------------------------------------------\n")
  plots <- generate_plots_from_cache(
    results = results,
    run_mechanism = run_mechanism,
    n_sims = n_sims
  )
  
  cat("\n=================================================\n")
  cat("  Pipeline complete!\n")
  cat("=================================================\n")
  cat("\nOutputs:\n")
  cat("  - intervention-results/*.csv (intervention metrics)\n")
  cat("  - intervention-plots/efficiency_distributions_violin.pdf\n")
  if (run_mechanism) {
    cat("  - intervention-plots/mechanism_analysis.png\n")
  }
  
  invisible(list(results = results, plots = plots))
}

# =============================================================================
# Quick usage examples (uncomment to run)
# =============================================================================

# # Run everything fresh:
# run_full_analysis()

# # Use cached results (much faster):
# run_full_analysis(use_cached = TRUE)

# # Just regenerate plots from most recent results:
# generate_plots_from_cache()

# # Skip mechanism analysis (even faster):
# run_full_analysis(use_cached = TRUE, run_mechanism = FALSE)
