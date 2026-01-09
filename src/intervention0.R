# =============================================================================
# Intervention Analysis for HIV Cluster Detection Strategies
# =============================================================================
#
# This script evaluates the effectiveness of different intervention strategies
# for HIV cluster detection and response. It processes simulated epidemic data
# to compare how various cluster prioritization methods perform in terms of:
#   - Person-years of Untreated infection Averted (PUTA)
#   - Potential Infections Averted (PIA)
#   - Contact tracing efficiency
#
# The six intervention strategies implemented are:
#   1. Distance-size (k=5): Trigger when 5 cases within genetic distance threshold
#   2. Distance-size (k=2): Trigger when 2 cases within genetic distance threshold
#   3. Growth-rate: Trigger when k infections occur within a time window
#   4. Random allocation: Baseline comparator with randomly selected individuals
#   5. RITA (Recent Infection Testing Algorithm): Target recently infected cases
#   6. Network degree: Target high-degree (well-connected) individuals
#
# Input data:
#   - D file (transmissions): donor, recipient, distance, timetransmission, simid
#   - G file (individuals): pid, generation, timeinfected, timediagnosed, 
#                           timesequenced, Fdegree, Gdegree, Hdegree, simid
#
# Output metrics:
#   - PUTA: Person-years of untreated infection averted by intervention
#   - PUTA/contact: Efficiency metric (PUTA per contact traced)
#   - PIA: Number of potential infections averted
#   - Contacts: Number of individuals requiring follow-up (contact tracing burden)
#
# Two subnetwork assumptions are computed for cluster-based strategies:
#   - "small" (dense): Assumes cluster members are highly connected internally
#   - "large" (sparse): Assumes cluster members have minimal internal connections
#
# Usage:
#   source("intervention0.R")
#   run_intervention_analysis(seed = 123, output_dir = "results")
#
# =============================================================================

library(glue)
library(ggplot2)

# -----------------------------------------------------------------------------
# Path resolution using the 'here' package (required)
# Ensures file paths work correctly regardless of working directory
# -----------------------------------------------------------------------------
if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required to run this script. Install it with install.packages('here').")
}

#' Resolve input file paths relative to project root
#' 
#' @param file Character string: filename or relative path
#' @return Absolute path to the file
#' @details
#'   - Absolute paths (Unix or Windows) are returned unchanged
#'   - Bare filenames check src/ first, then project root
#'   - Relative paths with directories resolve from project root
resolve_input_path <- function(file) {
  # Absolute paths (Unix or Windows) are returned unchanged
  if (grepl("^/|^[A-Za-z]:", file)) {
    return(file)
  }
  # Prefer src/<file> for unqualified filenames commonly kept under src/
  if (!grepl("[\\/]", file)) {
    candidate <- here::here("src", file)
    if (file.exists(candidate)) return(candidate)
  }
  # Fallback: resolve (possibly nested relative) path from project root
  here::here(file)
}


# =============================================================================
# MAIN ANALYSIS FUNCTION
# =============================================================================

#' Run complete intervention analysis across all strategies
#'
#' @param d_file Path to transmission data CSV (D file)
#' @param g_file Path to individual data CSV (G file)
#' @param seed Random seed for reproducibility (NULL = use system time)
#' @param cluster_size_5 Cluster size threshold for strategy 1 (default: 5)
#' @param cluster_size_2 Cluster size threshold for strategy 2 (default: 2)
#' @param distance_threshold Genetic distance threshold for clustering (default: 0.005)
#' @param network_degree_threshold Minimum degree to trigger network intervention (default: 4)
#' @param partner_notification_window_months Lookback window for partner notification (default: 6 months; alternative: 3 months)
#' @param random_coverage Proportion of eligible population to randomly sample (default: 0.30 = 30%)
#' @param rita_window_months Average RITA detection window in months (default: 6)
#' @param lookback_window_months Growth-rate trigger window in months (default: 6)
#' @param growth_distance_threshold Distance threshold for growth-based clustering (default: 0.1)
#' @param analysis_delay_days Delay for cluster analysis after sequencing (default: 14 days = 2 weeks)
#' @param implementation_delay_days Delay for implementing intervention after trigger (default: 14 days = 2 weeks)
#' @param show_table Whether to display results table (default: TRUE)
#' @param output_dir Directory to save results (NULL = don't save)
#'
#' @return List containing summary, counts, and detailed results for each strategy
run_intervention_analysis <- function(
  d_file = 'experiment1-N10000-gens7-D.csv',
  g_file = 'experiment1-N10000-gens7-G.csv',
  seed = NULL,  # Default to random seed
  cluster_size_5 = 5,
  cluster_size_2 = 2,
  distance_threshold = 0.005,
  network_degree_threshold = 4,
  partner_notification_window_months = 6,  # 6 months (180d) or 3 months (90d) lookback for contact tracing
  random_coverage = 0.10,  # 10% of eligible population
  rita_window_months = 6, # average RITA detection window
  lookback_window_months = 6, # growth-rate trigger window
  growth_distance_threshold = 0.01, # separate D for growth-based trigger
  analysis_delay_days = 14, # 2 weeks for cluster analysis (cluster-based only)
  implementation_delay_days = 14, # 2 weeks to implement intervention (all strategies)
  show_table = TRUE,
  output_dir = "intervention-results"  # if non-NULL, save results to this directory
) {
    # -------------------------------------------------------------------------
    # Set random seed for reproducibility
    # -------------------------------------------------------------------------
    if (is.null(seed)) {
      seed <- as.numeric(Sys.time())
      cat("Using random seed:", seed, "\n")
    } else {
      cat("Using fixed seed:", seed, "\n")
    }
    set.seed(seed)

    # -------------------------------------------------------------------------
    # Load and prepare simulation data
    # -------------------------------------------------------------------------
    cat("Loading data...\n")
    Dall <- read.csv(resolve_input_path(d_file), stringsAs = FALSE)
    Gall <- read.csv(resolve_input_path(g_file), stringsAs = FALSE)

    # Split by simulation ID
    # CRITICAL: Use character keys to ensure proper name-based list subsetting
    # (integer keys would cause numeric indexing, leading to mismatched D/G pairs)
    simids <- as.character(unique(Dall$simid))
    Ds <- split(Dall, Dall$simid)[simids]
    Gs <- split(Gall, Gall$simid)[simids]
    cat("  Loaded", nrow(Dall), "D rows,", nrow(Gall), "G rows across", length(simids), "simulations\n")

    # -------------------------------------------------------------------------
    # Calculate random sample size based on coverage of eligible population
    # -------------------------------------------------------------------------
    # Eligible population: generations 1 to (max-1), excluding seed and last gen
    lastgen <- max(Gall$generation)
    eligible_pop <- sum(Gall$generation > 0 & Gall$generation < lastgen)
    random_sample_size <- round(random_coverage * eligible_pop)
    cat("  Random sample size:", random_sample_size, "(", round(100*random_coverage), "% of", eligible_pop, "eligible)\n")

    # -------------------------------------------------------------------------
    # Run all six intervention strategies
    # -------------------------------------------------------------------------
    cat("Running interventions (6 strategies)...\n")
    t_start <- Sys.time()

    # Strategy 1: Distance-size with cluster_size = 5
    cat("  [1/6] Distance-size (size=", cluster_size_5, ", D=", distance_threshold, ")...", sep = "")
    ods5 <- tryCatch(
      distsize_intervention(
        Ds = Ds, Gs = Gs, Gall = Gall,
        distance_threshold = distance_threshold,
        cluster_size = cluster_size_5,
        analysis_delay_days = analysis_delay_days,
        implementation_delay_days = implementation_delay_days
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta_small = c(0, 0, 0, 0, 0), puta_large = c(0, 0, 0, 0, 0),
             pia_small = c(0, 0, 0, 0, 0), pia_large = c(0, 0, 0, 0, 0),
             total_contacts_small = 0, total_contacts_large = 0)
      }
    )
    cat(" done (", ods5$n_units, " units)\n", sep = "")

    # Strategy 2: Distance-size with cluster_size = 2
    cat("  [2/6] Distance-size (size=", cluster_size_2, ", D=", distance_threshold, ")...", sep = "")
    ods2 <- tryCatch(
      distsize_intervention(
        Ds = Ds, Gs = Gs, Gall = Gall,
        distance_threshold = distance_threshold,
        cluster_size = cluster_size_2,
        analysis_delay_days = analysis_delay_days,
        implementation_delay_days = implementation_delay_days
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta_small = c(0, 0, 0, 0, 0), puta_large = c(0, 0, 0, 0, 0),
             pia_small = c(0, 0, 0, 0, 0), pia_large = c(0, 0, 0, 0, 0),
             total_contacts_small = 0, total_contacts_large = 0)
      }
    )
    cat(" done (", ods2$n_units, " units)\n", sep = "")

    # Strategy 3: Growth-rate based intervention
    cat("  [3/6] Growth-rate (size=", cluster_size_5, ", W=", lookback_window_months, "mo, D=", growth_distance_threshold, ")...", sep = "")
    ogrowth <- tryCatch(
      growthrate_intervention(
        Ds = Ds, Gs = Gs, Gall = Gall,
        growth_distance_threshold = growth_distance_threshold,
        cluster_size = cluster_size_5,
        lookback_window_months = lookback_window_months,
        analysis_delay_days = analysis_delay_days,
        implementation_delay_days = implementation_delay_days
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta_small = c(0, 0, 0, 0, 0), puta_large = c(0, 0, 0, 0, 0),
             pia_small = c(0, 0, 0, 0, 0), pia_large = c(0, 0, 0, 0, 0),
             total_contacts_small = 0, total_contacts_large = 0)
      }
    )
    cat(" done (", ogrowth$n_units, " units)\n", sep = "")

    # Strategy 4: Random allocation (baseline comparator)
    cat("  [4/6] Random allocation (", round(100*random_coverage), "% coverage, n=", random_sample_size, ")...", sep = "")
    orand <- tryCatch(
      random_intervention(
        Dall = Dall, Gall = Gall,
        random_sample_size = random_sample_size,
        implementation_delay_days = implementation_delay_days
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta = c(0, 0, 0, 0, 0), pia = c(0, 0, 0, 0, 0), total_contacts = 0)
      }
    )
    cat(" done (", orand$n_units, " units)\n", sep = "")

    # Strategy 5: RITA (Recent Infection Testing Algorithm)
    cat("  [5/6] RITA (window=", rita_window_months, "mo)...", sep = "")
    orita <- tryCatch(
      rita_intervention(
        Dall = Dall, Gall = Gall,
        rita_window_months = rita_window_months,
        implementation_delay_days = implementation_delay_days
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta = c(0, 0, 0, 0, 0), pia = c(0, 0, 0, 0, 0), total_contacts = 0)
      }
    )
    cat(" done (", orita$n_units, " units)\n", sep = "")

    # Strategy 6: Network degree-based intervention
    cat("  [6/6] Network degree (threshold=", network_degree_threshold, ", window=", partner_notification_window_months, "mo)...", sep = "")
    onet <- tryCatch(
      network_intervention(
        Dall = Dall, Gall = Gall,
        network_degree_threshold = network_degree_threshold,
        implementation_delay_days = implementation_delay_days,
        partner_notification_window_months = partner_notification_window_months
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta = c(0, 0, 0, 0, 0), pia = c(0, 0, 0, 0, 0), total_contacts = 0)
      }
    )
    cat(" done (", onet$n_units, " units)\n", sep = "")

    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat("All interventions completed in", round(elapsed, 1), "seconds\n\n")

    # -------------------------------------------------------------------------
    # Compile summary table
    # -------------------------------------------------------------------------
    # Cluster-based strategies get 2 rows each (small/large subnetwork)
    # Individual-based strategies get 1 row (subnetwork = "-")
    
    # Helper to build a row for each strategy+subnetwork combination
    # puta_vec and pia_vec both have 5 elements: Total, Mean/contact, Median, Low, High
    build_row <- function(strategy_name, subnetwork, contacts, puta_vec, pia_vec) {
      c(Strategy = strategy_name, Subnetwork = subnetwork,
        Contacts = contacts, 
        Total_PUTA = puta_vec[1], 
        `PUTA/contact` = puta_vec[2],
        Median_PUTA = puta_vec[3], 
        Low_PUTA = puta_vec[4], 
        High_PUTA = puta_vec[5],
        Total_PIA = pia_vec[1], 
        `PIA/contact` = pia_vec[2],
        Median_PIA = pia_vec[3],
        Low_PIA = pia_vec[4], 
        High_PIA = pia_vec[5])
    }
    
    odf <- rbind(
      # Size=5: small and large rows
      build_row(paste0('Size=', cluster_size_5, ',D=', distance_threshold), "small",
                ods5$total_contacts_small, ods5$puta_small, ods5$pia_small),
      build_row(paste0('Size=', cluster_size_5, ',D=', distance_threshold), "large",
                ods5$total_contacts_large, ods5$puta_large, ods5$pia_large),
      # Size=2: small and large rows
      build_row(paste0('Size=', cluster_size_2, ',D=', distance_threshold), "small",
                ods2$total_contacts_small, ods2$puta_small, ods2$pia_small),
      build_row(paste0('Size=', cluster_size_2, ',D=', distance_threshold), "large",
                ods2$total_contacts_large, ods2$puta_large, ods2$pia_large),
      # Growth: small and large rows
      build_row(paste0('Growth,size=', cluster_size_5, ',W=', lookback_window_months, 'mo,D=', growth_distance_threshold), "small",
                ogrowth$total_contacts_small, ogrowth$puta_small, ogrowth$pia_small),
      build_row(paste0('Growth,size=', cluster_size_5, ',W=', lookback_window_months, 'mo,D=', growth_distance_threshold), "large",
                ogrowth$total_contacts_large, ogrowth$puta_large, ogrowth$pia_large),
      # Individual-based strategies: single row each (subnetwork = NA)
      build_row(paste0('Random,', round(100*random_coverage), '%'), "-", orand$total_contacts, orand$puta, orand$pia),
      build_row('RITA', "-", orita$total_contacts, orita$puta, orita$pia),
      build_row(paste0('Network,partners>', network_degree_threshold), "-", onet$total_contacts, onet$puta, onet$pia)
    ) |> as.data.frame()
    
    # Convert numeric columns from character
    num_cols <- c("Contacts", "Total_PUTA", "PUTA/contact", "Median_PUTA", "Low_PUTA", "High_PUTA",
                  "Total_PIA", "PIA/contact", "Median_PIA", "Low_PIA", "High_PIA")
    odf[num_cols] <- lapply(odf[num_cols], as.numeric)
    odf1 <- odf
    # Round PUTA metrics to 0dp, PIA per-contact metrics to 2dp
    cols_0dp <- c("Contacts", "Total_PUTA", "PUTA/contact", "Median_PUTA", "Low_PUTA", "High_PUTA", "Total_PIA")
    cols_2dp <- c("PIA/contact", "Median_PIA", "Low_PIA", "High_PIA")
    odf1[cols_0dp] <- lapply(odf1[cols_0dp], function(x) round(x, 0))
    odf1[cols_2dp] <- lapply(odf1[cols_2dp], function(x) round(x, 2))

    # Secondary counts table: units and contacts per strategy
    counts_df <- data.frame(
      Strategy = c(
        paste0('Size=', cluster_size_5, ',D=', distance_threshold),
        paste0('Size=', cluster_size_2, ',D=', distance_threshold),
        paste0('Growth,size=', cluster_size_5, ',W=', lookback_window_months, 'mo,D=', growth_distance_threshold),
        paste0('Random,', round(100*random_coverage), '%'),
        'RITA',
        paste0('Network,partners>', network_degree_threshold)
      ),
      Units = c(ods5$n_units, ods2$n_units, ogrowth$n_units, orand$n_units, orita$n_units, onet$n_units),
      Contacts_Small = round(c(ods5$total_contacts_small, ods2$total_contacts_small, ogrowth$total_contacts_small,
                               orand$total_contacts, orita$total_contacts, onet$total_contacts), 0),
      Contacts_Large = round(c(ods5$total_contacts_large, ods2$total_contacts_large, ogrowth$total_contacts_large,
                               orand$total_contacts, orita$total_contacts, onet$total_contacts), 0)
    )

    # -------------------------------------------------------------------------
    # Save outputs to files (if output_dir specified)
    # -------------------------------------------------------------------------
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        cat("Created output directory:", output_dir, "\n")
      }
      
      # Generate timestamp for unique filenames
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      
      # Save summary table
      summary_path <- file.path(output_dir, paste0("summary_", timestamp, ".csv"))
      write.csv(odf, summary_path, row.names = FALSE)
      
      # Save counts table
      counts_path <- file.path(output_dir, paste0("counts_", timestamp, ".csv"))
      write.csv(counts_df, counts_path, row.names = FALSE)
      
      # Save per-simulation details for each strategy (useful for further analysis)
      # Distance-size 5
      if (!is.null(ods5$o) && nrow(ods5$o) > 0) {
        ods5$o$strategy <- paste0("Size=", cluster_size_5, ",D=", distance_threshold)
        write.csv(ods5$o, file.path(output_dir, paste0("details_distsize5_", timestamp, ".csv")), row.names = FALSE)
      }
      # Distance-size 2
      if (!is.null(ods2$o) && nrow(ods2$o) > 0) {
        ods2$o$strategy <- paste0("Size=", cluster_size_2, ",D=", distance_threshold)
        write.csv(ods2$o, file.path(output_dir, paste0("details_distsize2_", timestamp, ".csv")), row.names = FALSE)
      }
      # Growth-rate
      if (!is.null(ogrowth$o) && nrow(ogrowth$o) > 0) {
        ogrowth$o$strategy <- paste0("Growth,size=", cluster_size_5, ",W=", lookback_window_months, "mo,D=", growth_distance_threshold)
        write.csv(ogrowth$o, file.path(output_dir, paste0("details_growth_", timestamp, ".csv")), row.names = FALSE)
      }
      # Random
      if (!is.null(orand$o) && nrow(orand$o) > 0) {
        orand$o$strategy <- "Random"
        write.csv(orand$o, file.path(output_dir, paste0("details_random_", timestamp, ".csv")), row.names = FALSE)
      }
      # RITA
      if (!is.null(orita$o) && nrow(orita$o) > 0) {
        orita$o$strategy <- "RITA"
        write.csv(orita$o, file.path(output_dir, paste0("details_rita_", timestamp, ".csv")), row.names = FALSE)
      }
      # Network
      if (!is.null(onet$o) && nrow(onet$o) > 0) {
        onet$o$strategy <- paste0("Network,partners>", network_degree_threshold)
        write.csv(onet$o, file.path(output_dir, paste0("details_network_", timestamp, ".csv")), row.names = FALSE)
      }
      
      # Also save a combined long-format details file for easy comparison
      all_details <- list()
      if (!is.null(ods5$o) && nrow(ods5$o) > 0) all_details$distsize5 <- ods5$o
      if (!is.null(ods2$o) && nrow(ods2$o) > 0) all_details$distsize2 <- ods2$o
      if (!is.null(ogrowth$o) && nrow(ogrowth$o) > 0) all_details$growth <- ogrowth$o
      # Random/RITA/Network have different column structures; combine cluster-based ones only
      if (length(all_details) > 0) {
        combined_clusters <- do.call(rbind, all_details)
        write.csv(combined_clusters, file.path(output_dir, paste0("details_all_clusters_", timestamp, ".csv")), row.names = FALSE)
      }
      
      # Save run parameters for reproducibility
      params <- data.frame(
        parameter = c("d_file", "g_file", "seed", "cluster_size_5", "cluster_size_2",
                      "distance_threshold", "network_degree_threshold",
                      "random_sample_size", "rita_window_months", "lookback_window_months",
                      "growth_distance_threshold", "analysis_delay_days", "implementation_delay_days", "n_simulations"),
        value = c(d_file, g_file, as.character(seed), cluster_size_5, cluster_size_2,
                  distance_threshold, network_degree_threshold,
                  random_sample_size, rita_window_months, lookback_window_months,
                  growth_distance_threshold, analysis_delay_days, implementation_delay_days, length(simids))
      )
      write.csv(params, file.path(output_dir, paste0("parameters_", timestamp, ".csv")), row.names = FALSE)
      
      cat("\nOutputs saved to:", output_dir, "\n")
      cat("  - summary_", timestamp, ".csv (aggregated metrics)\n", sep = "")
      cat("  - counts_", timestamp, ".csv (units and contacts per strategy)\n", sep = "")
      cat("  - details_*_", timestamp, ".csv (per-simulation/unit details)\n", sep = "")
      cat("  - parameters_", timestamp, ".csv (run parameters for reproducibility)\n", sep = "")
    }

    # -------------------------------------------------------------------------
    # Display results
    # -------------------------------------------------------------------------
    if (show_table) {
      if (require(knitr, quietly = TRUE)) {
        print(knitr::kable(odf1))
        cat("\nSample sizes and totals (for context):\n")
        print(knitr::kable(counts_df))
      } else {
        print(odf1)
        cat("\nSample sizes and totals (for context):\n")
        print(counts_df)
      }
    }
    
    # -------------------------------------------------------------------------
    # Plot cluster size distributions
    # -------------------------------------------------------------------------
    if (!is.null(output_dir) || show_table) {
      # Extract cluster sizes for cluster-based interventions
      distsize5_nc <- if (!is.null(ods5$o) && nrow(ods5$o) > 0) ods5$o$nc else numeric(0)
      distsize2_nc <- if (!is.null(ods2$o) && nrow(ods2$o) > 0) ods2$o$nc else numeric(0)
      growth_nc <- if (!is.null(ogrowth$o) && nrow(ogrowth$o) > 0) ogrowth$o$nc else numeric(0)
      
      if (length(distsize5_nc) > 0 || length(distsize2_nc) > 0 || length(growth_nc) > 0) {
        # Combine into data frame for plotting
        plot_df <- rbind(
          if (length(distsize5_nc) > 0) 
            data.frame(strategy = paste0('Distance-size (k=', cluster_size_5, ')'), nc = distsize5_nc) 
          else NULL,
          if (length(distsize2_nc) > 0) 
            data.frame(strategy = paste0('Distance-size (k=', cluster_size_2, ')'), nc = distsize2_nc) 
          else NULL,
          if (length(growth_nc) > 0) 
            data.frame(strategy = paste0('Growth-rate (k=', cluster_size_5, ', W=', lookback_window_months, 'mo)'), nc = growth_nc) 
          else NULL
        )
        
        if (nrow(plot_df) > 0) {
          p <- ggplot(plot_df, aes(x = nc, fill = strategy)) +
            geom_histogram(binwidth = 1, color = 'black', alpha = 0.7) +
            facet_wrap(~strategy, scales = 'free_y', ncol = 1) +
            labs(
              title = 'Distribution of Cluster Sizes at Intervention Time',
              x = 'Number of cluster members (nc) at intervention',
              y = 'Count'
            ) +
            theme_bw() +
            theme(legend.position = 'none') +
            scale_x_continuous(breaks = seq(0, max(plot_df$nc, na.rm = TRUE), by = 5))
          
          # Save plot if output_dir specified
          if (!is.null(output_dir)) {
            plot_path <- file.path(output_dir, paste0("cluster_sizes_", timestamp, ".png"))
            ggsave(plot_path, p, width = 10, height = 8, dpi = 150)
            cat("  - cluster_sizes_", timestamp, ".png (cluster size distributions)\n", sep = "")
          }
          
          # Display plot if show_table is TRUE
          if (show_table) {
            print(p)
          }
        }
      }
    }
    
    # Return results invisibly for programmatic use
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
  }

# =============================================================================
# INTERVENTION STRATEGY FUNCTIONS
# =============================================================================
# Each function implements a specific cluster detection/intervention strategy.
# Cluster-based strategies (proc_cluster, distsize, growthrate) compute both
# "small" and "large" subnetwork contact estimates.
# Individual-based strategies (random, rita, network) compute single contacts.
# =============================================================================


# -----------------------------------------------------------------------------
# Core cluster processing function
# -----------------------------------------------------------------------------

#' Process a single simulation to find and evaluate a cluster
#'
#' This function identifies a genetic cluster within a single simulation,
#' determines when intervention would be triggered, and calculates outcomes.
#'
#' Cluster formation logic:
#'   1. Filter transmissions to those within distance_threshold
#'   2. Build cluster by traversing from patient 0 through connected cases
#'   3. Exclude the final generation (not yet observable)
#'   4. Trigger intervention when cluster_size cases are sequenced
#'   5. Add exponential delay for intervention implementation
#'
#' Contact network estimation (two assumptions):
#'   - Small/Dense: contacts = n + sum(max(degree - (n-1), 0))
#'     Assumes cluster members are all connected to each other internally
#'   - Large/Sparse: contacts = sum(degrees) - (n - 2)
#'     Assumes minimal internal connections among cluster members
#'
#' @param D Data frame of transmissions for one simulation
#' @param G Data frame of individuals for one simulation
#' @param distance_threshold Genetic distance threshold for clustering
#' @param cluster_size Number of sequenced cases to trigger intervention
#' @param analysis_delay_days Fixed delay (days) for cluster analysis after detection
#' @param implementation_delay_days Fixed delay (days) for intervention implementation
#'
#' @return Named vector: pia, puta, interventiontime, nc, sum_degrees,
#'         sum_excess, contacts_small, contacts_large
proc_cluster <- function(
    D, G, 
    distance_threshold,
    cluster_size, 
    analysis_delay_days,
    implementation_delay_days
) 
{
  # Sort data by relevant time variables for proper temporal ordering
  G <- G[ order(G$timesequenced), ]
  D <- D[ order(D$timetransmission),]
  
  lastgeneration <- max( G$generation )
  
  # -------------------------------------------------------------------------
  # Step 1: Build genetic cluster from patient 0
  # -------------------------------------------------------------------------
  # Filter to transmissions within distance threshold
  D1 <- D[ D$distance <= distance_threshold , ]
  
  # Traverse transmission network starting from patient 0
  # Build set of all reachable patients (connected component)
  keeppids <- "0" 
  addpids  <- D1$recipient[ D1$donor %in% keeppids ]
  
  while( length( addpids ) > 0 ){
    keeppids <- union( addpids, keeppids )
    addpids  <- setdiff( D1$recipient[ D1$donor %in% keeppids ], keeppids )
  }
  
  # Filter to cluster members, excluding last generation (not yet observable)
  G1 <- G[ G$pid %in% keeppids , ]
  G1 <- G1[ G1$generation != lastgeneration, ]
  D1 <- D1[ D1$donor %in% G1$pid & D1$recipient %in% G1$pid, ]
  
  # -------------------------------------------------------------------------
  # Step 2: Determine intervention time (IT)
  # -------------------------------------------------------------------------
  # Trigger: when cluster_size cases are detected (sequenced + analysis_delay)
  # Intervention: trigger + implementation_delay
  # At IT, we only know about sequences from (IT - analysis_delay) due to lag
  IT <- Inf 
  if (nrow(G1) > 0 && nrow(G1) >= cluster_size) {
    idx <- cluster_size
    t_trigger <- G1$timesequenced[idx] + analysis_delay_days
    IT <- t_trigger + implementation_delay_days
  }
  
  # Filter cluster to those sequenced before (IT - analysis_delay)
  # because at IT we only know about sequences with analysis_delay lag
  G1 <- G1[ G1$timesequenced < (IT - analysis_delay_days) , ]
  
  # -------------------------------------------------------------------------
  # Step 3: Calculate contact network size (both subnetwork assumptions)
  # -------------------------------------------------------------------------
  sum_degrees <- 0
  sum_excess <- 0
  contacts_small <- 0
  contacts_large <- 0
  
  if (nrow(G1) > 0) {
    # Total degree = sum of F (long-term), G (casual), H (one-time) partnerships
    G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
    sum_degrees <- sum(G1$degree)
    n <- nrow(G1)
    
    # Small/Dense subnetwork: cluster members are fully connected internally
    # External contacts = degree minus connections to other cluster members
    # excess = max(degree - (n-1), 0) for each member
    excess <- pmax(G1$degree - (n - 1), 0)
    sum_excess <- sum(excess)
    contacts_small <- n + sum_excess
    
    # Large/Sparse subnetwork: minimal internal connections
    # Only transmission links connect cluster members
    contacts_large <- sum_degrees - (n - 2)
  }
  
  # -------------------------------------------------------------------------
  # Step 4: Calculate intervention outcomes (PIA and PUTA)
  # -------------------------------------------------------------------------
  # PIA: Potential Infections Averted - infections that occur after IT
  # PUTA: Person-years Untreated Averted - time between IT and diagnosis
  #       for people infected before IT but not yet diagnosed
  
  pia <- 0 
  puta <- 0 
  
  # Calculate transmissions from and to cluster members
  if (!is.infinite( IT ) )
  {
    # Find all individuals connected to cluster (donors, recipients, and members)
    piapids <-  D$recipient[D$donor %in% G1$pid] |> 
      union(G1$pid) |> 
      union(D$donor[D$recipient %in% G1$pid]) 
    
    G2  <- G[ G$pid %in% piapids , ]
    
    # PIA: count infections occurring after intervention time
    pia <- sum( G2$timeinfected  > IT )
    
    # PUTA: sum of (diagnosis time - IT) for those infected before IT
    #       but diagnosed after IT
    G3   <- G2[ G2$timeinfected <= IT & G2$timediagnosed > IT , ]
    puta <- sum( G3$timediagnosed - IT )
  }
  
  c(pia = pia, puta = puta, interventiontime = IT, nc = nrow(G1), 
    sum_degrees = sum_degrees, sum_excess = sum_excess,
    contacts_small = contacts_small, contacts_large = contacts_large)
}

  
# -----------------------------------------------------------------------------
# Strategy 1 & 2: Distance-size intervention
# -----------------------------------------------------------------------------

#' Distance-size cluster intervention strategy
#'
#' Triggers intervention when cluster_size genetically-linked cases are detected.
#' A cluster is defined by genetic distance threshold (e.g., 0.005 = 0.5% divergence).
#'
#' @param Ds List of transmission data frames, one per simulation
#' @param Gs List of individual data frames, one per simulation
#' @param Gall Combined individual data for all simulations
#' @param distance_threshold Genetic distance threshold for clustering
#' @param cluster_size Number of cases to trigger intervention
#' @param analysis_delay_days Fixed delay (days) for cluster analysis after detection
#' @param implementation_delay_days Fixed delay (days) for intervention implementation
#'
#' @return List with: o (detailed results), propintervened, n_units,
#'         puta_small, puta_large, pia, total_contacts_small, total_contacts_large
distsize_intervention <- function(
    Ds, Gs, Gall,
    distance_threshold, cluster_size, 
    analysis_delay_days, implementation_delay_days
)
{
  # Process each simulation using proc_cluster
  o <- lapply(seq_along(Ds), function(i) {
    tryCatch({
      proc_cluster(
        D = Ds[[i]], G = Gs[[i]],
        distance_threshold = distance_threshold,
        cluster_size = cluster_size,
        analysis_delay_days = analysis_delay_days,
        implementation_delay_days = implementation_delay_days
      )
    }, error = function(e) {
      # Return default values if processing fails
      c(pia = 0, puta = 0, interventiontime = Inf, nc = 0, 
        sum_degrees = 0, sum_excess = 0, contacts_small = 0, contacts_large = 0)
    })
  })
  
  # Combine results into data frame
  odf  <- do.call( rbind, o ) |> as.data.frame()
  
  # Attach simulation identifiers for traceability
  simids_vec <- names(Ds)
  odf$simid <- if (!is.null(simids_vec)) simids_vec else seq_len(nrow(odf))
  odf <- odf[, c("simid", "pia", "puta", "interventiontime", "nc", 
                 "sum_degrees", "sum_excess", "contacts_small", "contacts_large")]

  # -------------------------------------------------------------------------
  # Sanity-check: Recompute nc and contacts deterministically at recorded IT
  # This ensures consistency when intervention time was sampled stochastically
  # -------------------------------------------------------------------------
  if (nrow(odf) > 0) {
    # Coerce potential character columns
    suppressWarnings({
      odf$interventiontime <- as.numeric(odf$interventiontime)
      odf$nc <- as.numeric(odf$nc)
      odf$sum_degrees <- as.numeric(odf$sum_degrees)
      odf$sum_excess <- as.numeric(odf$sum_excess)
      odf$contacts_small <- as.numeric(odf$contacts_small)
      odf$contacts_large <- as.numeric(odf$contacts_large)
    })
    for (i in seq_len(nrow(odf))) {
      IT <- odf$interventiontime[i]
      if (!is.finite(IT)) next
      sim_key <- as.character(odf$simid[i])
      D <- Ds[[sim_key]]; G <- Gs[[sim_key]]
      if (is.null(D) || is.null(G)) next
      lastgeneration <- max(G$generation)
      D1 <- D[D$distance <= distance_threshold, ]
      keeppids <- "0"
      addpids <- D1$recipient[D1$donor %in% keeppids]
      while (length(addpids) > 0) {
        keeppids <- union(addpids, keeppids)
        addpids <- setdiff(D1$recipient[D1$donor %in% keeppids], keeppids)
      }
      G1 <- G[G$pid %in% keeppids, ]
      G1 <- G1[G1$generation != lastgeneration & G1$timesequenced < IT, ]
      if (nrow(G1) == 0) next
      G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
      n <- nrow(G1)
      sum_deg <- sum(G1$degree)
      excess <- pmax(G1$degree - (n - 1), 0)
      sum_exc <- sum(excess)
      odf$nc[i] <- n
      odf$sum_degrees[i] <- sum_deg
      odf$sum_excess[i] <- sum_exc
      odf$contacts_small[i] <- n + sum_exc
      odf$contacts_large[i] <- sum_deg - (n - 2)
    }
  }
  
  # Filter to simulations where a cluster was found (finite IT)
  odf1 <- odf[ !is.infinite( odf$interventiontime ), ]
  
  lastgen <- max( Gall$generation )
  
  # Handle case where no clusters were found
  if (nrow(odf1) == 0) {
    return(list(
      o = data.frame(pia = numeric(0), puta = numeric(0), interventiontime = numeric(0), 
                     nc = numeric(0), sum_degrees = numeric(0), sum_excess = numeric(0),
                     contacts_small = numeric(0), contacts_large = numeric(0)),
      propintervened = 0,
      n_units = 0,
      puta_small = c(0, 0, 0, 0, 0),
      puta_large = c(0, 0, 0, 0, 0),
      pia_small = c(0, 0, 0, 0, 0),
      pia_large = c(0, 0, 0, 0, 0),
      total_contacts_small = 0,
      total_contacts_large = 0
    ))
  }
  
  # -------------------------------------------------------------------------
  # Compute summary statistics
  # -------------------------------------------------------------------------
  
  # PUTA efficiency for small subnetwork assumption
  e_puta_small <- odf1$puta / odf1$contacts_small
  e_puta_small_valid <- e_puta_small[is.finite(e_puta_small) & !is.na(e_puta_small)]
  med_puta_small <- if (length(e_puta_small_valid) > 0) median(e_puta_small_valid) else NA_real_
  q_puta_small <- if (length(e_puta_small_valid) > 0) quantile(e_puta_small_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
  
  # PUTA efficiency for large subnetwork assumption
  e_puta_large <- odf1$puta / odf1$contacts_large
  e_puta_large_valid <- e_puta_large[is.finite(e_puta_large) & !is.na(e_puta_large)]
  med_puta_large <- if (length(e_puta_large_valid) > 0) median(e_puta_large_valid) else NA_real_
  q_puta_large <- if (length(e_puta_large_valid) > 0) quantile(e_puta_large_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
  
  # PIA efficiency for small subnetwork assumption (same structure as PUTA)
  e_pia_small <- odf1$pia / odf1$contacts_small
  e_pia_small_valid <- e_pia_small[is.finite(e_pia_small) & !is.na(e_pia_small)]
  med_pia_small <- if (length(e_pia_small_valid) > 0) median(e_pia_small_valid) else NA_real_
  q_pia_small <- if (length(e_pia_small_valid) > 0) quantile(e_pia_small_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
  
  # PIA efficiency for large subnetwork assumption
  e_pia_large <- odf1$pia / odf1$contacts_large
  e_pia_large_valid <- e_pia_large[is.finite(e_pia_large) & !is.na(e_pia_large)]
  med_pia_large <- if (length(e_pia_large_valid) > 0) median(e_pia_large_valid) else NA_real_
  q_pia_large <- if (length(e_pia_large_valid) > 0) quantile(e_pia_large_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
  
  # Return comprehensive results
  list(
    o = odf1,
    propintervened = sum(odf1$nc) / sum(Gall$generation > 0 & Gall$generation < lastgen),
    n_units = nrow(odf1),
    puta_small = c(sum(odf1$puta),
         if (sum(odf1$contacts_small) > 0) sum(odf1$puta) / sum(odf1$contacts_small) else NA_real_,
         med_puta_small,
         q_puta_small[1],
         q_puta_small[2]),
    puta_large = c(sum(odf1$puta),
         if (sum(odf1$contacts_large) > 0) sum(odf1$puta) / sum(odf1$contacts_large) else NA_real_,
         med_puta_large,
         q_puta_large[1],
         q_puta_large[2]),
    pia_small = c(sum(odf1$pia),
         if (sum(odf1$contacts_small) > 0) sum(odf1$pia) / sum(odf1$contacts_small) else NA_real_,
         med_pia_small,
         q_pia_small[1],
         q_pia_small[2]),
    pia_large = c(sum(odf1$pia),
         if (sum(odf1$contacts_large) > 0) sum(odf1$pia) / sum(odf1$contacts_large) else NA_real_,
         med_pia_large,
         q_pia_large[1],
         q_pia_large[2]),
    total_contacts_small = sum(odf1$contacts_small),
    total_contacts_large = sum(odf1$contacts_large)
  )
}

# -----------------------------------------------------------------------------
# Strategy 3: Growth-rate cluster intervention
# -----------------------------------------------------------------------------

#' Growth-rate based cluster intervention strategy
#'
#' Triggers intervention when cluster_size cases are sequenced within a sliding
#' time window (lookback_window_months). This detects rapidly growing clusters
#' rather than simply large clusters.
#'
#' Key difference from distance-size:
#'   - Distance-size triggers on k-th sequenced case (any time span)
#'   - Growth-rate triggers when k cases are sequenced within W months
#'
#' Both use sequencing times for trigger and cluster membership.
#'
#' @param Ds List of transmission data frames, one per simulation
#' @param Gs List of individual data frames, one per simulation
#' @param Gall Combined individual data for all simulations
#' @param growth_distance_threshold Genetic distance threshold (typically larger than distsize)
#' @param cluster_size Number of sequenced cases in window to trigger
#' @param lookback_window_months Width of sliding window in months
#' @param analysis_delay_days Fixed delay (days) for cluster analysis after detection
#' @param implementation_delay_days Fixed delay (days) for intervention implementation
#'
#' @return List with same structure as distsize_intervention
growthrate_intervention <- function(
  Ds, Gs, Gall,
  growth_distance_threshold, cluster_size,
  lookback_window_months = 3,
  analysis_delay_days, implementation_delay_days
)
{
  # Convert window to days
  lookback_days <- lookback_window_months * 30

  # Inner function to process one simulation
  process_one <- function(D, G) {
    # Order by SEQUENCING time (like distsize) for growth trigger
    G <- G[order(G$timesequenced), ]
    D <- D[order(D$timetransmission), ]

    lastgeneration <- max(G$generation)

    # -------------------------------------------------------------------------
    # Build genetic cluster (same logic as distsize)
    # -------------------------------------------------------------------------
    D1 <- D[D$distance <= growth_distance_threshold, ]
    keeppids <- "0"
    addpids <- D1$recipient[D1$donor %in% keeppids]
    while (length(addpids) > 0) {
      keeppids <- union(addpids, keeppids)
      addpids <- setdiff(D1$recipient[D1$donor %in% keeppids], keeppids)
    }
    Gcluster <- G[G$pid %in% keeppids, ]

    # -------------------------------------------------------------------------
    # Growth trigger: sliding window over SEQUENCING times
    # -------------------------------------------------------------------------
    # Exclude seed (generation 0) and last generation (not yet observable)
    Gtrig <- Gcluster[Gcluster$generation > 0 & Gcluster$generation != lastgeneration, ]
    if (nrow(Gtrig) == 0) {
      return(c(pia = 0, puta = 0, interventiontime = Inf, nc = 0, 
               sum_degrees = 0, sum_excess = 0, contacts_small = 0, contacts_large = 0))
    }

    # Sliding window algorithm to find earliest time when cluster_size
    # cases are SEQUENCED within lookback_days window
    t <- suppressWarnings(as.numeric(Gtrig$timesequenced))
    if (all(!is.finite(t))) {
      return(c(pia = 0, puta = 0, interventiontime = Inf, nc = 0, 
               sum_degrees = 0, sum_excess = 0, contacts_small = 0, contacts_large = 0))
    }
    # Keep only finite sequencing times for windowing
    keep <- is.finite(t)
    t <- t[keep]
    Gtrig <- Gtrig[keep, ]
    n <- length(t)
    
    # Two-pointer sliding window: j is left edge, i is right edge
    # Find earliest time when k cases are sequenced within lookback_days
    j <- 1
    t_detect <- Inf
    trigger_j <- NA  # Store window indices for later
    trigger_i <- NA
    for (i in seq_len(n)) {
      # Advance left pointer while window exceeds lookback_days
      while (j <= i && (t[i] - t[j]) > lookback_days) {
        j <- j + 1
      }
      # Check if window contains cluster_size sequenced cases
      if ((i - j + 1) >= cluster_size) {
        t_detect <- t[i]  # Trigger at the k-th sequencing time in window
        trigger_j <- j    # Store window boundaries
        trigger_i <- i
        break
      }
    }

    # Set intervention time (trigger + implementation delay)
    # Trigger = detection time + analysis delay
    # At IT, we only know about sequences from (IT - analysis_delay) due to lag
    t_trigger <- if (is.finite(t_detect)) t_detect + analysis_delay_days else Inf
    IT <- if (is.finite(t_trigger)) t_trigger + implementation_delay_days else Inf

    # -------------------------------------------------------------------------
    # Define cluster at IT: only cases we KNOW about at intervention time
    # Due to analysis_delay lag, we only know sequences from before (IT - analysis_delay)
    # -------------------------------------------------------------------------
    G1 <- data.frame()
    if (is.finite(IT) && !is.na(trigger_j)) {
      # Get the start time of the trigger window
      window_start_time <- t[trigger_j]
      # Include cluster members sequenced from window start to (IT - analysis_delay)
      # This is the knowledge cutoff at intervention time
      knowledge_cutoff <- IT - analysis_delay_days
      G1 <- Gcluster[Gcluster$generation != lastgeneration & 
                     Gcluster$timesequenced >= window_start_time & 
                     Gcluster$timesequenced < knowledge_cutoff, ]
    }

    # Compute contacts - both subnetwork assumptions
    sum_degrees <- 0
    sum_excess <- 0
    contacts_small <- 0
    contacts_large <- 0
    if (nrow(G1) > 0 && is.finite(IT)) {
      G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
      n1 <- nrow(G1)
      sum_degrees <- sum(G1$degree)
      excess <- pmax(G1$degree - (n1 - 1), 0)
      sum_excess <- sum(excess)
      contacts_small <- n1 + sum_excess
      contacts_large <- sum_degrees - (n1 - 2)
    }

    # Calculate PIA/PUTA (same logic as distsize)
    pia <- 0
    puta <- 0
    if (is.finite(IT) && nrow(G1) > 0) {
      piapids <- D$recipient[D$donor %in% G1$pid] |>
        union(G1$pid) |>
        union(D$donor[D$recipient %in% G1$pid])
      G2 <- G[G$pid %in% piapids, ]
      pia <- sum(G2$timeinfected > IT)
      G3 <- G2[G2$timeinfected <= IT & G2$timediagnosed > IT, ]
      puta <- sum(G3$timediagnosed - IT)
    }

    c(pia = pia, puta = puta, interventiontime = IT, nc = nrow(G1), 
      sum_degrees = sum_degrees, sum_excess = sum_excess,
      contacts_small = contacts_small, contacts_large = contacts_large)
  }

  # Process all simulations
  o <- lapply(seq_along(Ds), function(i) {
    tryCatch({
      process_one(Ds[[i]], Gs[[i]])
    }, error = function(e) {
      c(pia = 0, puta = 0, interventiontime = Inf, nc = 0, 
        sum_degrees = 0, sum_excess = 0, contacts_small = 0, contacts_large = 0)
    })
  })

  odf <- do.call(rbind, o) |> as.data.frame()
  simids_vec <- names(Ds)
  odf$simid <- if (!is.null(simids_vec)) simids_vec else seq_len(nrow(odf))
  odf <- odf[, c("simid", "pia", "puta", "interventiontime", "nc", 
                 "sum_degrees", "sum_excess", "contacts_small", "contacts_large")]

  # Coerce columns to numeric
  suppressWarnings({
    odf$interventiontime <- as.numeric(odf$interventiontime)
    odf$nc <- as.numeric(odf$nc)
    odf$sum_degrees <- as.numeric(odf$sum_degrees)
    odf$sum_excess <- as.numeric(odf$sum_excess)
    odf$contacts_small <- as.numeric(odf$contacts_small)
    odf$contacts_large <- as.numeric(odf$contacts_large)
  })

  # Filter to simulations where growth trigger fired
  odf1 <- odf[!is.infinite(odf$interventiontime), ]
  lastgen <- max(Gall$generation)

  # Handle case where no clusters triggered
  if (nrow(odf1) == 0) {
    return(list(
      o = data.frame(pia = numeric(0), puta = numeric(0), interventiontime = numeric(0),
                     nc = numeric(0), sum_degrees = numeric(0), sum_excess = numeric(0),
                     contacts_small = numeric(0), contacts_large = numeric(0)),
      propintervened = 0,
      n_units = 0,
      puta_small = c(0, 0, 0, 0, 0),
      puta_large = c(0, 0, 0, 0, 0),
      pia_small = c(0, 0, 0, 0, 0),
      pia_large = c(0, 0, 0, 0, 0),
      total_contacts_small = 0,
      total_contacts_large = 0
    ))
  }

  # -------------------------------------------------------------------------
  # Compute summary statistics (same structure as distsize)
  # -------------------------------------------------------------------------
  
  # PUTA efficiency for small subnetwork
  e_puta_small <- odf1$puta / odf1$contacts_small
  e_puta_small_valid <- e_puta_small[is.finite(e_puta_small) & !is.na(e_puta_small)]
  med_puta_small <- if (length(e_puta_small_valid) > 0) median(e_puta_small_valid) else NA_real_
  q_puta_small <- if (length(e_puta_small_valid) > 0) quantile(e_puta_small_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # PUTA efficiency for large subnetwork
  e_puta_large <- odf1$puta / odf1$contacts_large
  e_puta_large_valid <- e_puta_large[is.finite(e_puta_large) & !is.na(e_puta_large)]
  med_puta_large <- if (length(e_puta_large_valid) > 0) median(e_puta_large_valid) else NA_real_
  q_puta_large <- if (length(e_puta_large_valid) > 0) quantile(e_puta_large_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # PIA efficiency for small subnetwork
  e_pia_small <- odf1$pia / odf1$contacts_small
  e_pia_small_valid <- e_pia_small[is.finite(e_pia_small) & !is.na(e_pia_small)]
  med_pia_small <- if (length(e_pia_small_valid) > 0) median(e_pia_small_valid) else NA_real_
  q_pia_small <- if (length(e_pia_small_valid) > 0) quantile(e_pia_small_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # PIA efficiency for large subnetwork
  e_pia_large <- odf1$pia / odf1$contacts_large
  e_pia_large_valid <- e_pia_large[is.finite(e_pia_large) & !is.na(e_pia_large)]
  med_pia_large <- if (length(e_pia_large_valid) > 0) median(e_pia_large_valid) else NA_real_
  q_pia_large <- if (length(e_pia_large_valid) > 0) quantile(e_pia_large_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  list(
    o = odf1,
    propintervened = sum(odf1$nc) / sum(Gall$generation > 0 & Gall$generation < lastgen),
    n_units = nrow(odf1),
    puta_small = c(sum(odf1$puta),
         if (sum(odf1$contacts_small) > 0) sum(odf1$puta) / sum(odf1$contacts_small) else NA_real_,
         med_puta_small,
         q_puta_small[1],
         q_puta_small[2]),
    puta_large = c(sum(odf1$puta),
         if (sum(odf1$contacts_large) > 0) sum(odf1$puta) / sum(odf1$contacts_large) else NA_real_,
         med_puta_large,
         q_puta_large[1],
         q_puta_large[2]),
    pia_small = c(sum(odf1$pia),
         if (sum(odf1$contacts_small) > 0) sum(odf1$pia) / sum(odf1$contacts_small) else NA_real_,
         med_pia_small,
         q_pia_small[1],
         q_pia_small[2]),
    pia_large = c(sum(odf1$pia),
         if (sum(odf1$contacts_large) > 0) sum(odf1$pia) / sum(odf1$contacts_large) else NA_real_,
         med_pia_large,
         q_pia_large[1],
         q_pia_large[2]),
    total_contacts_small = sum(odf1$contacts_small),
    total_contacts_large = sum(odf1$contacts_large)
  )
}

# -----------------------------------------------------------------------------
# Strategy 4: Random allocation intervention (baseline comparator)
# -----------------------------------------------------------------------------

#' Random allocation intervention strategy
#'
#' Randomly selects individuals for intervention as a baseline comparator.
#' This represents a non-targeted approach where resources are allocated
#' without using cluster or risk information.
#'
#' @param Dall Combined transmission data for all simulations
#' @param Gall Combined individual data for all simulations
#' @param random_sample_size Number of individuals to randomly select
#' @param implementation_delay_days Fixed delay (days) for intervention implementation
#'
#' @return List with: o (detailed results), propintervened (NA for random),
#'         n_units, puta, pia, total_contacts
random_intervention <- function(Dall, Gall, random_sample_size, implementation_delay_days)
{
  lastgeneration <- max( Gall$generation )
  
  # Select from generations 2+ (excluding seed and last generation)
  G1 <- Gall[ Gall$generation > 0 & Gall$generation < lastgeneration, ]
    
    if (nrow(G1) == 0) {
      return(list(
        o = data.frame(pia = numeric(0), puta = numeric(0)),
        propintervened = 0,
        n_units = 0,
        puta = c(0, 0, 0, 0, 0),
        pia = c(0, 0, 0, 0, 0),
        total_contacts = 0
      ))
    }
    
    # Randomly sample cases
    sample_size <- min(random_sample_size, nrow(G1))
    G1 <- G1 |> dplyr::slice_sample(n = sample_size)
    
    # Total degree = F + G + H partnerships
    G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
    
    # Make PIDs unique across simulations by appending simid
    G1$pid <- paste(sep='.', G1$pid, G1$simid )
    D <- Dall; D$donor <- paste(sep='.', D$donor, D$simid )
    D$recipient <- paste(sep='.', D$recipient, D$simid )
    G <- Gall; G$pid <- paste(sep='.', G$pid, G$simid)
    
    # Set intervention time at diagnosis + implementation delay
    G1$IT <- G1$timediagnosed + implementation_delay_days
    
    # Process individual intervention outcomes
    proc_indiv <- function(pid, IT, degree) {
      # Find all contacts: recipients (transmitted to), self, donors (transmitted from)
      piapids <- D$recipient[ D$donor == pid] |> 
        union( c( pid, D$donor[D$recipient==pid] )) 
      
      G2 <- G[ G$pid %in% piapids , ]
      pia <- sum( G2$timeinfected > IT )
      
      G3 <- G2[ G2$timeinfected <= IT & G2$timediagnosed > IT, ]
      puta <- sum( G3$timediagnosed - IT )
      
      # Contacts = degree + 1 (self)
      c( pia, puta, degree + 1 )
    }
    
    mapply(proc_indiv, G1$pid, G1$IT, G1$degree ) -> o 
    o <- as.data.frame( t( o ) )
    colnames(o) <- c('pia', 'puta', 'contacts' )
    
    # Compute summary statistics (same structure as cluster-based)
    e_puta_percontact <- o$puta / o$contacts
    e_puta_valid <- e_puta_percontact[is.finite(e_puta_percontact) & !is.na(e_puta_percontact)]
    med_puta <- if (length(e_puta_valid) > 0) median(e_puta_valid) else NA_real_
    q_puta <- if (length(e_puta_valid) > 0) quantile(e_puta_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
    
    e_pia_percontact <- o$pia / o$contacts
    e_pia_valid <- e_pia_percontact[is.finite(e_pia_percontact) & !is.na(e_pia_percontact)]
    med_pia <- if (length(e_pia_valid) > 0) median(e_pia_valid) else NA_real_
    q_pia <- if (length(e_pia_valid) > 0) quantile(e_pia_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
    
    list(
      o = o,
      propintervened = NA,  # Not applicable for random selection
      n_units = nrow(o),
      puta = c(sum(o$puta),
               sum(o$puta) / sum(o$contacts),
               med_puta,
               q_puta[1],
               q_puta[2]),
      pia = c(sum(o$pia),
              sum(o$pia) / sum(o$contacts),
              med_pia,
              q_pia[1],
              q_pia[2]),
      total_contacts = sum(o$contacts)
    )
}

# -----------------------------------------------------------------------------
# Strategy 5: RITA intervention (Recent Infection Testing Algorithm)
# -----------------------------------------------------------------------------

#' RITA-based intervention strategy
#'
#' Targets individuals identified as recently infected using the Recent
#' Infection Testing Algorithm. RITA detects acute infections based on
#' biomarkers that indicate infection within approximately 6 months.
#'
#' RITA simulation:
#'   - Test is positive if (time_diagnosed - time_infected) < random exponential
#'   - Average detection window is rita_window_months (default 6 months)
#'
#' @param Dall Combined transmission data for all simulations
#' @param Gall Combined individual data for all simulations
#' @param rita_window_months Average RITA detection window in months
#' @param implementation_delay_days Fixed delay (days) for intervention implementation
#'
#' @return List with same structure as random_intervention
rita_intervention <- function(Dall, Gall, rita_window_months, implementation_delay_days)
{
  lastgeneration <- max( Gall$generation )
  
  # Simulate RITA test: positive if diagnosed within random window of infection
  # Average window = rita_window_months * 30 days
  Gall$rita <- with(Gall, (timediagnosed - timeinfected) < rexp(nrow(Gall), 1/(rita_window_months*30)))
  
  # Filter to RITA-positive cases in generations 2+ (excluding last)
  G1 <- Gall[Gall$rita & (Gall$generation > 0) & (Gall$generation < lastgeneration), ]
    
    if (nrow(G1) == 0) {
      return(list(
        o = data.frame(pia = numeric(0), puta = numeric(0), contacts = numeric(0)),
        propintervened = 0,
        n_units = 0,
        puta = c(0, 0, 0, 0, 0),
        pia = c(0, 0, 0, 0, 0),
        total_contacts = 0
      ))
    }
    
    # Total degree for RITA-positive individuals
    G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
    
    # Make PIDs unique across simulations
    G1$pid <- paste(sep = '.', G1$pid, G1$simid)
    D <- Dall; D$donor <- paste(sep='.', D$donor, D$simid )
    D$recipient <- paste(sep='.', D$recipient, D$simid )
    G <- Gall; G$pid <- paste(sep='.', G$pid, G$simid)
    
    # Set intervention time at diagnosis + implementation delay
    G1$IT <- G1$timediagnosed + implementation_delay_days
    
    # Process individual intervention outcomes
    proc_indiv <- function(pid, IT, degree) {
      piapids <- D$recipient[D$donor == pid] |> 
        union(c(pid, D$donor[D$recipient == pid]))
      
      G2 <- G[G$pid %in% piapids, ]
      pia <- sum(G2$timeinfected > IT)
      
      G3 <- G2[G2$timeinfected <= IT & G2$timediagnosed > IT, ]
      puta <- sum(G3$timediagnosed - IT)
      
      c(pia, puta, degree + 1)
    }
    
    # Process all RITA-positive cases
    results <- matrix(nrow = nrow(G1), ncol = 3)
    for (i in seq_len(nrow(G1))) {
      results[i, ] <- proc_indiv(G1$pid[i], G1$IT[i], G1$degree[i])
    }
    
    o <- as.data.frame(results)
    colnames(o) <- c('pia', 'puta', 'contacts')
    
    # Compute summary statistics (same structure as cluster-based)
    e_puta_percontact <- o$puta / o$contacts
    e_puta_valid <- e_puta_percontact[is.finite(e_puta_percontact) & !is.na(e_puta_percontact)]
    med_puta <- if (length(e_puta_valid) > 0) median(e_puta_valid) else NA_real_
    q_puta <- if (length(e_puta_valid) > 0) quantile(e_puta_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
    
    e_pia_percontact <- o$pia / o$contacts
    e_pia_valid <- e_pia_percontact[is.finite(e_pia_percontact) & !is.na(e_pia_percontact)]
    med_pia <- if (length(e_pia_valid) > 0) median(e_pia_valid) else NA_real_
    q_pia <- if (length(e_pia_valid) > 0) quantile(e_pia_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
    
    list(
      o = o,
      propintervened = NA,  # Not directly comparable to cluster-based
      n_units = nrow(o),
      puta = c(sum(o$puta),
               sum(o$puta) / sum(o$contacts),
               med_puta,
               q_puta[1],
               q_puta[2]),
      pia = c(sum(o$pia),
              sum(o$pia) / sum(o$contacts),
              med_pia,
              q_pia[1],
              q_pia[2]),
      total_contacts = sum(o$contacts)
    )
}


# -----------------------------------------------------------------------------
# Strategy 6: Network degree intervention
# -----------------------------------------------------------------------------

#' Network degree-based intervention strategy
#'
#' Targets individuals with high network degree (many sexual partners).
#' This approach prioritizes "superspreaders" who may be more likely to
#' transmit infection to many others.
#'
#' Contact calculation uses the simulated contact counts within the specified
#' partner notification window (90 days for 3 months, 180 days for 6 months):
#'   contacts = Fcontacts_XXd + Gcontacts_XXd + Hcontacts_XXd
#'   where XX is 90 or 180 depending on partner_notification_window_months
#'
#' @param Dall Combined transmission data for all simulations
#' @param Gall Combined individual data for all simulations
#' @param network_degree_threshold Minimum total contacts to trigger intervention
#' @param implementation_delay_days Fixed delay (days) for intervention implementation
#' @param partner_notification_window_months Lookback window: 3 or 6 months (default: 6)
#'
#' @return List with same structure as random_intervention
network_intervention <- function(Dall, Gall, network_degree_threshold, implementation_delay_days,
                                  partner_notification_window_months = 6)
{
  lastgeneration <- max(Gall$generation)
  
  # Make PIDs unique across simulations
  D <- Dall
  D$donor <- paste(sep = '.', D$donor, D$simid)
  D$recipient <- paste(sep = '.', D$recipient, D$simid)
  G <- Gall
  G$pid <- paste(sep = '.', G$pid, G$simid)
  
  # Calculate total contacts based on partner notification window
  # Use 90-day columns for 3-month window, 180-day columns for 6-month window
  if (partner_notification_window_months == 3) {
    G$total_contacts <- with(G, Fcontacts_90d + Gcontacts_90d + Hcontacts_90d)
  } else {
    # Default to 6 months (180 days)
    G$total_contacts <- with(G, Fcontacts_180d + Gcontacts_180d + Hcontacts_180d)
  }
  
  # Filter to high-contact individuals in generations 2+ (excluding last)
  G1 <- G[ (G$total_contacts >= network_degree_threshold) & (G$generation > 0) & (G$generation < lastgeneration), ]
    
    if (nrow(G1) == 0) {
      return(list(
        o = data.frame(pia = numeric(0), puta = numeric(0), contacts = numeric(0)),
        propintervened = 0,
        n_units = 0,
        puta = c(0, 0, 0, 0, 0),
        pia = c(0, 0, 0, 0, 0),
        total_contacts = 0
      ))
    }
    
    # Set intervention time at diagnosis + implementation delay
    G1$IT <- G1$timediagnosed + implementation_delay_days
    
    # Process individual intervention outcomes
    proc_indiv <- function(pid, IT, total_contacts) {
      piapids <- D$recipient[D$donor == pid] |> 
        union(c(pid, D$donor[D$recipient == pid]))
      
      G2 <- G[G$pid %in% piapids, ]
      pia <- sum(G2$timeinfected > IT)
      
      G3 <- G2[G2$timeinfected <= IT & G2$timediagnosed > IT, ]
      puta <- sum(G3$timediagnosed - IT)
      
      c(pia, puta, total_contacts + 1)  # +1 for the index case
    }
    
    mapply(proc_indiv, G1$pid, G1$IT, G1$total_contacts) -> o 
    o <- as.data.frame(t(o))
    colnames(o) <- c('pia', 'puta', 'contacts')
    
    # Compute summary statistics (same structure as cluster-based)
    e_puta_percontact <- o$puta / o$contacts
    e_puta_valid <- e_puta_percontact[is.finite(e_puta_percontact) & !is.na(e_puta_percontact)]
    med_puta <- if (length(e_puta_valid) > 0) median(e_puta_valid) else NA_real_
    q_puta <- if (length(e_puta_valid) > 0) quantile(e_puta_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
    
    e_pia_percontact <- o$pia / o$contacts
    e_pia_valid <- e_pia_percontact[is.finite(e_pia_percontact) & !is.na(e_pia_percontact)]
    med_pia <- if (length(e_pia_valid) > 0) median(e_pia_valid) else NA_real_
    q_pia <- if (length(e_pia_valid) > 0) quantile(e_pia_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
    
    list(
      o = o,
      propintervened = NA,  # Not directly comparable to cluster-based
      n_units = nrow(o),
      puta = c(sum(o$puta),
               sum(o$puta) / sum(o$contacts),
               med_puta,
               q_puta[1],
               q_puta[2]),
      pia = c(sum(o$pia),
              sum(o$pia) / sum(o$contacts),
              med_pia,
              q_pia[1],
              q_pia[2]),
      total_contacts = sum(o$contacts)
    )
}

# =============================================================================
# plot_efficiency_distributions
# =============================================================================
#' Plot distributions of PUTA and PIA efficiency across intervention strategies
#' 
#' Creates 4 plots showing PUTA and PIA efficiency distributions under small 
#' and large subnetwork assumptions for all 6 intervention strategies.
#' 
#' @param results The output from run_interventions()
#' @param plot_type Either "density" for smoothed density plots, "violin" for 
#'   violin plots, or "boxplot" for box plots
#' @param title_prefix Optional prefix for plot titles
#' @return A combined ggplot object with 4 panels
#' 
plot_efficiency_distributions <- function(results, 
                                          plot_type = "density",
                                          title_prefix = "") {
  require(ggplot2)
  require(patchwork)
  
  # Strategy names and display labels
  strategy_names <- c("distsize5", "distsize2", "growth", "random", "rita", "network")
  strategy_labels <- c("Size>=5", "Size>=2", "Growth", "Random", "RITA", "Network")
  
  # Extract per-unit data from each strategy's results
  extract_data <- function(details) {
    dfs <- lapply(seq_along(strategy_names), function(i) {
      sname <- strategy_names[i]
      slabel <- strategy_labels[i]
      
      if (is.null(details[[sname]]) || is.null(details[[sname]]$o)) {
        return(NULL)
      }
      
      o <- details[[sname]]$o
      
      # For network intervention, there's only one contacts column
      if (sname == "network") {
        if (!"contacts" %in% names(o)) return(NULL)
        data.frame(
          strategy = slabel,
          puta = o$puta,
          pia = o$pia,
          contacts_small = o$contacts,  # Use same contacts for both assumptions
          contacts_large = o$contacts
        )
      } else {
        # For cluster-based strategies
        if (!all(c("puta", "pia", "contacts_small", "contacts_large") %in% names(o))) return(NULL)
        data.frame(
          strategy = slabel,
          puta = o$puta,
          pia = o$pia,
          contacts_small = o$contacts_small,
          contacts_large = o$contacts_large
        )
      }
    })
    
    do.call(rbind, Filter(Negate(is.null), dfs))
  }
  
  df <- extract_data(results$details)
  
  if (is.null(df) || nrow(df) == 0) {
    stop("No data available for plotting. Run interventions first.")
  }
  
  # Compute efficiencies
  df$puta_eff_small <- df$puta / df$contacts_small
  df$puta_eff_large <- df$puta / df$contacts_large
  df$pia_eff_small <- df$pia / df$contacts_small
  df$pia_eff_large <- df$pia / df$contacts_large
  
  # Remove invalid values
  df <- df[is.finite(df$puta_eff_small) & is.finite(df$puta_eff_large) &
           is.finite(df$pia_eff_small) & is.finite(df$pia_eff_large), ]
  
  # Set factor levels for ordering
  df$strategy <- factor(df$strategy, levels = strategy_labels)
  
  # Define color palette
  strategy_colors <- c(
    "Size>=5" = "#E41A1C",
    "Size>=2" = "#377EB8", 
    "Growth" = "#4DAF4A",
    "Random" = "#984EA3",
    "RITA" = "#FF7F00",
    "Network" = "#A65628"
  )
  
  # Create plots based on plot_type
  if (plot_type == "density") {
    p1 <- ggplot(df, aes(x = puta_eff_small, fill = strategy, color = strategy)) +
      geom_density(alpha = 0.3) +
      labs(x = "PUTA per contact", y = "Density",
           title = paste0(title_prefix, "PUTA Efficiency (Small Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      scale_color_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "bottom", legend.title = element_blank())
    
    p2 <- ggplot(df, aes(x = puta_eff_large, fill = strategy, color = strategy)) +
      geom_density(alpha = 0.3) +
      labs(x = "PUTA per contact", y = "Density",
           title = paste0(title_prefix, "PUTA Efficiency (Large Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      scale_color_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "bottom", legend.title = element_blank())
    
    p3 <- ggplot(df, aes(x = pia_eff_small, fill = strategy, color = strategy)) +
      geom_density(alpha = 0.3) +
      labs(x = "PIA per contact", y = "Density",
           title = paste0(title_prefix, "PIA Efficiency (Small Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      scale_color_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "bottom", legend.title = element_blank())
    
    p4 <- ggplot(df, aes(x = pia_eff_large, fill = strategy, color = strategy)) +
      geom_density(alpha = 0.3) +
      labs(x = "PIA per contact", y = "Density",
           title = paste0(title_prefix, "PIA Efficiency (Large Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      scale_color_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "bottom", legend.title = element_blank())
    
  } else if (plot_type == "violin") {
    p1 <- ggplot(df, aes(x = strategy, y = puta_eff_small, fill = strategy)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
      labs(x = "", y = "PUTA per contact",
           title = paste0(title_prefix, "PUTA Efficiency (Small Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    p2 <- ggplot(df, aes(x = strategy, y = puta_eff_large, fill = strategy)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
      labs(x = "", y = "PUTA per contact",
           title = paste0(title_prefix, "PUTA Efficiency (Large Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    p3 <- ggplot(df, aes(x = strategy, y = pia_eff_small, fill = strategy)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
      labs(x = "", y = "PIA per contact",
           title = paste0(title_prefix, "PIA Efficiency (Small Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    p4 <- ggplot(df, aes(x = strategy, y = pia_eff_large, fill = strategy)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
      labs(x = "", y = "PIA per contact",
           title = paste0(title_prefix, "PIA Efficiency (Large Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
  } else if (plot_type == "boxplot") {
    p1 <- ggplot(df, aes(x = strategy, y = puta_eff_small, fill = strategy)) +
      geom_boxplot(alpha = 0.7) +
      labs(x = "", y = "PUTA per contact",
           title = paste0(title_prefix, "PUTA Efficiency (Small Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    p2 <- ggplot(df, aes(x = strategy, y = puta_eff_large, fill = strategy)) +
      geom_boxplot(alpha = 0.7) +
      labs(x = "", y = "PUTA per contact",
           title = paste0(title_prefix, "PUTA Efficiency (Large Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    p3 <- ggplot(df, aes(x = strategy, y = pia_eff_small, fill = strategy)) +
      geom_boxplot(alpha = 0.7) +
      labs(x = "", y = "PIA per contact",
           title = paste0(title_prefix, "PIA Efficiency (Small Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    p4 <- ggplot(df, aes(x = strategy, y = pia_eff_large, fill = strategy)) +
      geom_boxplot(alpha = 0.7) +
      labs(x = "", y = "PIA per contact",
           title = paste0(title_prefix, "PIA Efficiency (Large Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # Combine plots
  (p1 + p2) / (p3 + p4) + 
    plot_annotation(title = "Intervention Efficiency Distributions")
}
