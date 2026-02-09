# =============================================================================
# Intervention Analysis for HIV Cluster Detection Strategies
# =============================================================================
#
# This script evaluates the effectiveness of different intervention strategies
# for HIV cluster detection and response. It processes simulated epidemic data
# to compare how various cluster prioritization methods perform in terms of:
#   - Infectious Days Averted (IDA)
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
#   - IDA: Infectious days averted by intervention
#   - IDA/contact: Efficiency metric (IDA per contact traced)
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
library(rlang)  # For .data pronoun in ggplot2 aes()

# Suppress R CMD check notes for ggplot2 non-standard evaluation
utils::globalVariables(c("nc", "strategy", "ida_eff_small", "ida_eff_large",
                         "pia_eff_small", "pia_eff_large"))

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
  random_coverage = 0.30,  # 30% of eligible population
  rita_window_months = 6, # average RITA detection window
  lookback_window_months = 6, # growth-rate trigger window
  growth_distance_threshold = 0.01, # separate D for growth-based trigger
  analysis_delay_days = 14, # 2 weeks for cluster analysis (cluster-based only)
  implementation_delay_days = 14, # 2 weeks to implement intervention (all strategies)
  show_table = TRUE,
  output_dir = "intervention-results"  # if non-NULL, save results to this directory
) {
    # -------------------------------------------------------------------------
    # Validate parameters
    # -------------------------------------------------------------------------
    # Partner notification window must be 3 or 6 months (corresponding to 90d or 180d contact windows)
    if (!partner_notification_window_months %in% c(3, 6)) {
      stop("partner_notification_window_months must be 3 or 6 (corresponding to 90-day or 180-day contact windows)")
    }

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

    # Validate that G contains all simulation IDs present in D
    # (G may have additional sims with no transmissions - this is normal)
    d_simids <- unique(Dall$simid)
    g_simids <- unique(Gall$simid)
    missing_in_g <- setdiff(d_simids, g_simids)
    if (length(missing_in_g) > 0) {
      stop("Simulation ID mismatch: G file missing simids present in D: ",
           paste(head(missing_in_g, 10), collapse = ", "))
    }

    # Filter G to only include simulations present in D
    if (length(g_simids) > length(d_simids)) {
      cat(sprintf("  Note: G has %d simulations, D has %d (filtering G to match D)\n",
                  length(g_simids), length(d_simids)))
      Gall <- Gall[Gall$simid %in% d_simids, ]
    }

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
        implementation_delay_days = implementation_delay_days,
        partner_notification_window_months = partner_notification_window_months
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
        implementation_delay_days = implementation_delay_days,
        partner_notification_window_months = partner_notification_window_months
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
        implementation_delay_days = implementation_delay_days,
        partner_notification_window_months = partner_notification_window_months
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
        implementation_delay_days = implementation_delay_days,
        partner_notification_window_months = partner_notification_window_months
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             ida = c(0, 0, 0, 0, 0), pia = c(0, 0, 0, 0, 0), total_contacts = 0)
      }
    )
    cat(" done (", orand$n_units, " units)\n", sep = "")

    # Strategy 5: RITA (Recent Infection Testing Algorithm)
    cat("  [5/6] RITA (window=", rita_window_months, "mo)...", sep = "")
    orita <- tryCatch(
      rita_intervention(
        Dall = Dall, Gall = Gall,
        rita_window_months = rita_window_months,
        implementation_delay_days = implementation_delay_days,
        partner_notification_window_months = partner_notification_window_months
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             ida = c(0, 0, 0, 0, 0), pia = c(0, 0, 0, 0, 0), total_contacts = 0)
      }
    )
    cat(" done (", orita$n_units, " units)\n", sep = "")

    # Strategy 6: Network degree-based intervention
    cat("  [6/7] Network degree (threshold=", network_degree_threshold, ", window=", partner_notification_window_months, "mo)...", sep = "")
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
             ida = c(0, 0, 0, 0, 0), pia = c(0, 0, 0, 0, 0), total_contacts = 0)
      }
    )
    cat(" done (", onet$n_units, " units)\n", sep = "")

    # Strategy 7: RITA + Secondary Contact Tracing
    cat("  [7/7] RITA+Secondary (window=", rita_window_months, "mo, lookback=", partner_notification_window_months, "mo)...", sep = "")
    oritasec <- tryCatch(
      rita_secondary_intervention(
        Dall = Dall, Gall = Gall,
        rita_window_months = rita_window_months,
        implementation_delay_days = implementation_delay_days,
        partner_notification_window_months = partner_notification_window_months,
        analysis_delay_days = 0
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(o = data.frame(pia = numeric(0), ida = numeric(0), interventiontime = numeric(0),
                           nc = numeric(0), sum_degrees = numeric(0), sum_excess = numeric(0),
                           contacts_small = numeric(0), contacts_large = numeric(0)),
             propintervened = 0, n_units = 0,
             puta_small = c(0, 0, 0, 0, 0), puta_large = c(0, 0, 0, 0, 0),
             pia_small = c(0, 0, 0, 0, 0), pia_large = c(0, 0, 0, 0, 0),
             total_contacts_small = 0, total_contacts_large = 0)
      }
    )
    cat(" done (", oritasec$n_units, " units)\n", sep = "")

    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat("All interventions completed in", round(elapsed, 1), "seconds\n\n")

    # -------------------------------------------------------------------------
    # Compile summary table
    # -------------------------------------------------------------------------
    # Cluster-based strategies get 2 rows each (small/large subnetwork)
    # Individual-based strategies get 1 row (subnetwork = "-")
    
    # Helper to build a row for each strategy+subnetwork combination
    # ida_vec and pia_vec both have 5 elements: Total, Mean/contact, Median, Low, High
    build_row <- function(strategy_name, subnetwork, contacts, ida_vec, pia_vec) {
      c(Strategy = strategy_name, Subnetwork = subnetwork,
        Contacts = contacts, 
        Total_PUTA = ida_vec[1], 
        `IDA/contact` = ida_vec[2],
        Median_PUTA = ida_vec[3], 
        Low_PUTA = ida_vec[4], 
        High_PUTA = ida_vec[5],
        Total_PIA = pia_vec[1], 
        `PIA/contact` = pia_vec[2],
        Median_PIA = pia_vec[3],
        Low_PIA = pia_vec[4], 
        High_PIA = pia_vec[5])
    }

    # Build summary table (wrap in tryCatch to handle potential issues)
    odf <- tryCatch({
      # Build rows list
      rows <- list(
        # Size=5: small and large rows
        build_row(paste0('Size=', cluster_size_5, ',D=', distance_threshold), "small",
                  ods5$total_contacts_small, ods5$ida_small, ods5$pia_small),
        build_row(paste0('Size=', cluster_size_5, ',D=', distance_threshold), "large",
                  ods5$total_contacts_large, ods5$ida_large, ods5$pia_large),
        # Size=2: small and large rows
        build_row(paste0('Size=', cluster_size_2, ',D=', distance_threshold), "small",
                  ods2$total_contacts_small, ods2$ida_small, ods2$pia_small),
        build_row(paste0('Size=', cluster_size_2, ',D=', distance_threshold), "large",
                  ods2$total_contacts_large, ods2$ida_large, ods2$pia_large),
        # Growth: small and large rows
        build_row(paste0('Growth,size=', cluster_size_5, ',W=', lookback_window_months, 'mo,D=', growth_distance_threshold), "small",
                  ogrowth$total_contacts_small, ogrowth$ida_small, ogrowth$pia_small),
        build_row(paste0('Growth,size=', cluster_size_5, ',W=', lookback_window_months, 'mo,D=', growth_distance_threshold), "large",
                  ogrowth$total_contacts_large, ogrowth$ida_large, ogrowth$pia_large),
        # Individual-based strategies: single row each (subnetwork = NA)
        build_row(paste0('Random,', round(100*random_coverage), '%'), "-", orand$total_contacts, orand$ida, orand$pia),
        build_row('RITA', "-", orita$total_contacts, orita$ida, orita$pia),
        build_row(paste0('Network,partners>', network_degree_threshold), "-", onet$total_contacts, onet$ida, onet$pia),
        # RITA+Secondary: cluster-based (small and large rows)
        build_row(paste0('RITA+Secondary,W=', rita_window_months, 'mo'), "small",
                  oritasec$total_contacts_small, oritasec$ida_small, oritasec$pia_small),
        build_row(paste0('RITA+Secondary,W=', rita_window_months, 'mo'), "large",
                  oritasec$total_contacts_large, oritasec$ida_large, oritasec$pia_large)
      )

      # Use incremental rbind instead of do.call(rbind, rows)
      # do.call has issues with named vectors
      odf_mat <- rows[[1]]
      for (i in 2:length(rows)) {
        odf_mat <- rbind(odf_mat, rows[[i]])
      }
      odf_temp <- as.data.frame(odf_mat, stringsAsFactors = FALSE)

      # Convert numeric columns from character
      num_cols <- c("Contacts", "Total_PUTA", "IDA/contact", "Median_PUTA", "Low_PUTA", "High_PUTA",
                    "Total_PIA", "PIA/contact", "Median_PIA", "Low_PIA", "High_PIA")
      odf_temp[num_cols] <- lapply(odf_temp[num_cols], as.numeric)

      odf_temp
    }, error = function(e) {
      cat("\nWarning: Could not build summary table. Details files will still be saved.\n")
      cat("Error message:", e$message, "\n")
      NULL
    })
    # Round numeric columns if summary table was successfully created
    odf1 <- if (!is.null(odf)) {
      # Round IDA metrics to 0dp, PIA per-contact metrics to 2dp
      cols_0dp <- c("Contacts", "Total_PUTA", "IDA/contact", "Median_PUTA", "Low_PUTA", "High_PUTA", "Total_PIA")
      cols_2dp <- c("PIA/contact", "Median_PIA", "Low_PIA", "High_PIA")
      odf_rounded <- odf
      odf_rounded[cols_0dp] <- lapply(odf[cols_0dp], function(x) round(x, 0))
      odf_rounded[cols_2dp] <- lapply(odf[cols_2dp], function(x) round(x, 2))
      odf_rounded
    } else {
      NULL
    }

    # Secondary counts table: units and contacts per strategy
    counts_df <- data.frame(
      Strategy = c(
        paste0('Size=', cluster_size_5, ',D=', distance_threshold),
        paste0('Size=', cluster_size_2, ',D=', distance_threshold),
        paste0('Growth,size=', cluster_size_5, ',W=', lookback_window_months, 'mo,D=', growth_distance_threshold),
        paste0('Random,', round(100*random_coverage), '%'),
        'RITA',
        paste0('Network,partners>', network_degree_threshold),
        paste0('RITA+Secondary,W=', rita_window_months, 'mo')
      ),
      Units = c(ods5$n_units, ods2$n_units, ogrowth$n_units, orand$n_units, orita$n_units, onet$n_units, oritasec$n_units),
      Contacts_Small = round(c(ods5$total_contacts_small, ods2$total_contacts_small, ogrowth$total_contacts_small,
                               orand$total_contacts, orita$total_contacts, onet$total_contacts, oritasec$total_contacts_small), 0),
      Contacts_Large = round(c(ods5$total_contacts_large, ods2$total_contacts_large, ogrowth$total_contacts_large,
                               orand$total_contacts, orita$total_contacts, onet$total_contacts, oritasec$total_contacts_large), 0)
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
      
      # Save summary table (if it was successfully created)
      if (!is.null(odf1)) {
        summary_path <- file.path(output_dir, paste0("summary_", timestamp, ".csv"))
        write.csv(odf1, summary_path, row.names = FALSE)
      }
      
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
      # RITA+Secondary
      if (!is.null(oritasec$o) && nrow(oritasec$o) > 0) {
        oritasec$o$strategy <- paste0("RITA+Secondary,W=", rita_window_months, "mo")
        write.csv(oritasec$o, file.path(output_dir, paste0("details_ritasecondary_", timestamp, ".csv")), row.names = FALSE)
      }

      # Also save a combined long-format details file for easy comparison
      all_details <- list()
      if (!is.null(ods5$o) && nrow(ods5$o) > 0) all_details$distsize5 <- ods5$o
      if (!is.null(ods2$o) && nrow(ods2$o) > 0) all_details$distsize2 <- ods2$o
      if (!is.null(ogrowth$o) && nrow(ogrowth$o) > 0) all_details$growth <- ogrowth$o
      if (!is.null(oritasec$o) && nrow(oritasec$o) > 0) all_details$ritasecondary <- oritasec$o
      # Random/RITA/Network have different column structures; combine cluster-based ones only
      if (length(all_details) > 0) {
        # Use dplyr::bind_rows which handles column mismatches gracefully
        if (requireNamespace("dplyr", quietly = TRUE)) {
          combined_clusters <- dplyr::bind_rows(all_details)
        } else {
          # Fallback: use incremental rbind with error handling
          tryCatch({
            combined_clusters <- all_details[[1]]
            if (length(all_details) > 1) {
              for (i in 2:length(all_details)) {
                combined_clusters <- rbind(combined_clusters, all_details[[i]])
              }
            }
          }, error = function(e) {
            cat("Warning: Could not combine cluster details:", e$message, "\n")
            combined_clusters <- NULL
          })
        }
        if (!is.null(combined_clusters)) {
          write.csv(combined_clusters, file.path(output_dir, paste0("details_all_clusters_", timestamp, ".csv")), row.names = FALSE)
        }
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
      ritasec_nc <- if (!is.null(oritasec$o) && nrow(oritasec$o) > 0) oritasec$o$nc else numeric(0)

      if (length(distsize5_nc) > 0 || length(distsize2_nc) > 0 || length(growth_nc) > 0 || length(ritasec_nc) > 0) {
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
          else NULL,
          if (length(ritasec_nc) > 0)
            data.frame(strategy = paste0('RITA+Secondary (W=', rita_window_months, 'mo)'), nc = ritasec_nc)
          else NULL
        )
        
        if (nrow(plot_df) > 0) {
          p <- ggplot(plot_df, aes(x = .data$nc, fill = .data$strategy)) +
            geom_histogram(binwidth = 1, color = 'black', alpha = 0.7) +
            facet_wrap(~ .data$strategy, scales = 'free_y', ncol = 1) +
            labs(
              title = 'Distribution of Cluster Sizes at Intervention Time',
              x = 'Number of cluster members (nc) at intervention',
              y = 'Count'
            ) +
            theme_bw() +
            theme(legend.position = 'none') +
            scale_x_continuous(breaks = seq(0, max(plot_df$nc, na.rm = TRUE), by = 5))
          
          # Save plot to intervention-plots directory
          plot_dir <- here::here("intervention-plots")
          if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
          plot_path <- file.path(plot_dir, paste0("cluster_sizes_", timestamp, ".png"))
          ggsave(plot_path, p, width = 10, height = 8, dpi = 150)
          cat("  - intervention-plots/cluster_sizes_", timestamp, ".png (cluster size distributions)\n", sep = "")
          
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
        network = onet,
        ritasecondary = oritasec
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
#' @return Named vector: pia, ida, interventiontime, nc, sum_degrees,
#'         sum_excess, contacts_small, contacts_large
proc_cluster <- function(
    D, G, 
    distance_threshold,
    cluster_size, 
    analysis_delay_days,
    implementation_delay_days,
    partner_notification_window_months = 6
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
    # Calculate total contacts based on partner notification window
    if (partner_notification_window_months == 3) {
      G1$total_contacts <- with(G1, Fcontacts_90d + Gcontacts_90d + Hcontacts_90d)
    } else {
      G1$total_contacts <- with(G1, Fcontacts_180d + Gcontacts_180d + Hcontacts_180d)
    }
    sum_degrees <- sum(G1$total_contacts)
    n <- nrow(G1)
    
    # Small/Dense subnetwork: cluster members are fully connected internally
    # External contacts = total_contacts minus connections to other cluster members
    # excess = max(total_contacts - (n-1), 0) for each member
    excess <- pmax(G1$total_contacts - (n - 1), 0)
    sum_excess <- sum(excess)
    contacts_small <- n + sum_excess
    
    # Large/Sparse subnetwork: minimal internal connections
    # Only transmission links connect cluster members
    contacts_large <- sum_degrees - (n - 2)
  }
  
  # -------------------------------------------------------------------------
  # Step 4: Calculate intervention outcomes (PIA and IDA)
  # -------------------------------------------------------------------------
  # PIA: Potential Infections Averted - infections that occur after IT
  # IDA: Infectious Days Averted - time between IT and diagnosis
  #       for people infected before IT but not yet diagnosed
  
  pia <- 0 
  ida <- 0 
  
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
    
    # IDA: sum of (diagnosis time - IT) for those infected before IT
    #       but diagnosed after IT
    G3   <- G2[ G2$timeinfected <= IT & G2$timediagnosed > IT , ]
    ida <- sum( G3$timediagnosed - IT )
  }
  
  c(pia = pia, ida = ida, interventiontime = IT, nc = nrow(G1), 
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
#' @param partner_notification_window_months Lookback window: 3 or 6 months (default: 6)
#'
#' @return List with: o (detailed results), propintervened, n_units,
#'         puta_small, puta_large, pia, total_contacts_small, total_contacts_large
distsize_intervention <- function(
    Ds, Gs, Gall,
    distance_threshold, cluster_size, 
    analysis_delay_days, implementation_delay_days,
    partner_notification_window_months = 6
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
        implementation_delay_days = implementation_delay_days,
        partner_notification_window_months = partner_notification_window_months
      )
    }, error = function(e) {
      # Return default values if processing fails
      c(pia = 0, ida = 0, interventiontime = Inf, nc = 0, 
        sum_degrees = 0, sum_excess = 0, contacts_small = 0, contacts_large = 0)
    })
  })
  
  # Combine results into data frame
  odf  <- do.call( rbind, o ) |> as.data.frame()
  
  # Attach simulation identifiers for traceability
  simids_vec <- names(Ds)
  odf$simid <- if (!is.null(simids_vec)) simids_vec else seq_len(nrow(odf))
  odf <- odf[, c("simid", "pia", "ida", "interventiontime", "nc",
                 "sum_degrees", "sum_excess", "contacts_small", "contacts_large")]

  # Coerce columns to numeric (in case they came through as character from rbind)
  if (nrow(odf) > 0) {
    suppressWarnings({
      odf$interventiontime <- as.numeric(odf$interventiontime)
      odf$nc <- as.numeric(odf$nc)
      odf$sum_degrees <- as.numeric(odf$sum_degrees)
      odf$sum_excess <- as.numeric(odf$sum_excess)
      odf$contacts_small <- as.numeric(odf$contacts_small)
      odf$contacts_large <- as.numeric(odf$contacts_large)
    })
  }

  # Filter to simulations where a cluster was found (finite IT)
  odf1 <- odf[ !is.infinite( odf$interventiontime ), ]
  
  lastgen <- max( Gall$generation )
  
  # Handle case where no clusters were found
  if (nrow(odf1) == 0) {
    return(list(
      o = data.frame(pia = numeric(0), ida = numeric(0), interventiontime = numeric(0), 
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
  # Quantiles: 10th and 90th percentiles used uniformly for all metrics

  # IDA efficiency for small subnetwork assumption
  e_puta_small <- odf1$ida / odf1$contacts_small
  e_puta_small_valid <- e_puta_small[is.finite(e_puta_small) & !is.na(e_puta_small)]
  med_ida_small <- if (length(e_puta_small_valid) > 0) median(e_puta_small_valid) else NA_real_
  q_ida_small <- if (length(e_puta_small_valid) > 0) quantile(e_puta_small_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # IDA efficiency for large subnetwork assumption
  e_puta_large <- odf1$ida / odf1$contacts_large
  e_puta_large_valid <- e_puta_large[is.finite(e_puta_large) & !is.na(e_puta_large)]
  med_ida_large <- if (length(e_puta_large_valid) > 0) median(e_puta_large_valid) else NA_real_
  q_ida_large <- if (length(e_puta_large_valid) > 0) quantile(e_puta_large_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # PIA efficiency for small subnetwork assumption
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
    puta_small = c(sum(odf1$ida),
         if (sum(odf1$contacts_small) > 0) sum(odf1$ida) / sum(odf1$contacts_small) else NA_real_,
         med_ida_small,
         q_ida_small[1],
         q_ida_small[2]),
    puta_large = c(sum(odf1$ida),
         if (sum(odf1$contacts_large) > 0) sum(odf1$ida) / sum(odf1$contacts_large) else NA_real_,
         med_ida_large,
         q_ida_large[1],
         q_ida_large[2]),
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
#' @param partner_notification_window_months Lookback window for contacts: 3 or 6 months (default: 6)
#'
#' @return List with same structure as distsize_intervention
growthrate_intervention <- function(
  Ds, Gs, Gall,
  growth_distance_threshold, cluster_size,
  lookback_window_months = 3,
  analysis_delay_days, implementation_delay_days,
  partner_notification_window_months = 6
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
      return(c(pia = 0, ida = 0, interventiontime = Inf, nc = 0, 
               sum_degrees = 0, sum_excess = 0, contacts_small = 0, contacts_large = 0))
    }

    # Sliding window algorithm to find earliest time when cluster_size
    # cases are SEQUENCED within lookback_days window
    t <- suppressWarnings(as.numeric(Gtrig$timesequenced))
    if (all(!is.finite(t))) {
      return(c(pia = 0, ida = 0, interventiontime = Inf, nc = 0, 
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
    trigger_j <- NA  # Store window index for later
    for (i in seq_len(n)) {
      # Advance left pointer while window exceeds lookback_days
      while (j <= i && (t[i] - t[j]) > lookback_days) {
        j <- j + 1
      }
      # Check if window contains cluster_size sequenced cases
      if ((i - j + 1) >= cluster_size) {
        t_detect <- t[i]  # Trigger at the k-th sequencing time in window
        trigger_j <- j    # Store window boundary
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
      # Calculate total contacts based on partner notification window
      if (partner_notification_window_months == 3) {
        G1$total_contacts <- with(G1, Fcontacts_90d + Gcontacts_90d + Hcontacts_90d)
      } else {
        G1$total_contacts <- with(G1, Fcontacts_180d + Gcontacts_180d + Hcontacts_180d)
      }
      n1 <- nrow(G1)
      sum_degrees <- sum(G1$total_contacts)
      excess <- pmax(G1$total_contacts - (n1 - 1), 0)
      sum_excess <- sum(excess)
      contacts_small <- n1 + sum_excess
      contacts_large <- sum_degrees - (n1 - 2)
    }

    # Calculate PIA/IDA (same logic as distsize)
    pia <- 0
    ida <- 0
    if (is.finite(IT) && nrow(G1) > 0) {
      piapids <- D$recipient[D$donor %in% G1$pid] |>
        union(G1$pid) |>
        union(D$donor[D$recipient %in% G1$pid])
      G2 <- G[G$pid %in% piapids, ]
      pia <- sum(G2$timeinfected > IT)
      G3 <- G2[G2$timeinfected <= IT & G2$timediagnosed > IT, ]
      ida <- sum(G3$timediagnosed - IT)
    }

    c(pia = pia, ida = ida, interventiontime = IT, nc = nrow(G1), 
      sum_degrees = sum_degrees, sum_excess = sum_excess,
      contacts_small = contacts_small, contacts_large = contacts_large)
  }

  # Process all simulations
  o <- lapply(seq_along(Ds), function(i) {
    tryCatch({
      process_one(Ds[[i]], Gs[[i]])
    }, error = function(e) {
      c(pia = 0, ida = 0, interventiontime = Inf, nc = 0, 
        sum_degrees = 0, sum_excess = 0, contacts_small = 0, contacts_large = 0)
    })
  })

  odf <- do.call(rbind, o) |> as.data.frame()
  simids_vec <- names(Ds)
  odf$simid <- if (!is.null(simids_vec)) simids_vec else seq_len(nrow(odf))
  odf <- odf[, c("simid", "pia", "ida", "interventiontime", "nc", 
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
      o = data.frame(pia = numeric(0), ida = numeric(0), interventiontime = numeric(0),
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
  # Quantiles: 10th and 90th percentiles used uniformly for all metrics

  # IDA efficiency for small subnetwork
  e_puta_small <- odf1$ida / odf1$contacts_small
  e_puta_small_valid <- e_puta_small[is.finite(e_puta_small) & !is.na(e_puta_small)]
  med_ida_small <- if (length(e_puta_small_valid) > 0) median(e_puta_small_valid) else NA_real_
  q_ida_small <- if (length(e_puta_small_valid) > 0) quantile(e_puta_small_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # IDA efficiency for large subnetwork
  e_puta_large <- odf1$ida / odf1$contacts_large
  e_puta_large_valid <- e_puta_large[is.finite(e_puta_large) & !is.na(e_puta_large)]
  med_ida_large <- if (length(e_puta_large_valid) > 0) median(e_puta_large_valid) else NA_real_
  q_ida_large <- if (length(e_puta_large_valid) > 0) quantile(e_puta_large_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

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
    puta_small = c(sum(odf1$ida),
         if (sum(odf1$contacts_small) > 0) sum(odf1$ida) / sum(odf1$contacts_small) else NA_real_,
         med_ida_small,
         q_ida_small[1],
         q_ida_small[2]),
    puta_large = c(sum(odf1$ida),
         if (sum(odf1$contacts_large) > 0) sum(odf1$ida) / sum(odf1$contacts_large) else NA_real_,
         med_ida_large,
         q_ida_large[1],
         q_ida_large[2]),
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
#' @param partner_notification_window_months Lookback window: 3 or 6 months (default: 6)
#'
#' @return List with: o (detailed results), propintervened (NA for random),
#'         n_units, ida, pia, total_contacts
random_intervention <- function(Dall, Gall, random_sample_size, implementation_delay_days,
                                partner_notification_window_months = 6)
{
  lastgeneration <- max( Gall$generation )
  
  # Select from generations 2+ (excluding seed and last generation)
  G1 <- Gall[ Gall$generation > 0 & Gall$generation < lastgeneration, ]
    
    if (nrow(G1) == 0) {
      return(list(
        o = data.frame(pia = numeric(0), ida = numeric(0)),
        propintervened = 0,
        n_units = 0,
        ida = c(0, 0, 0, 0, 0),
        pia = c(0, 0, 0, 0, 0),
        total_contacts = 0
      ))
    }
    
    # Randomly sample cases
    sample_size <- min(random_sample_size, nrow(G1))
    G1 <- G1 |> dplyr::slice_sample(n = sample_size)
    
    # Calculate total contacts based on partner notification window
    if (partner_notification_window_months == 3) {
      G1$total_contacts <- with(G1, Fcontacts_90d + Gcontacts_90d + Hcontacts_90d)
    } else {
      G1$total_contacts <- with(G1, Fcontacts_180d + Gcontacts_180d + Hcontacts_180d)
    }

    # Save simid before modifying pid (for paired comparison analysis)
    G1_simid <- G1$simid

    # Make PIDs unique across simulations by appending simid
    G1$pid <- paste(sep='.', G1$pid, G1$simid )
    D <- Dall; D$donor <- paste(sep='.', D$donor, D$simid )
    D$recipient <- paste(sep='.', D$recipient, D$simid )
    G <- Gall; G$pid <- paste(sep='.', G$pid, G$simid)

    # Set intervention time at diagnosis + implementation delay
    G1$IT <- G1$timediagnosed + implementation_delay_days

    # Process individual intervention outcomes
    proc_indiv <- function(pid, IT, total_contacts) {
      # Find all contacts: recipients (transmitted to), self, donors (transmitted from)
      piapids <- D$recipient[ D$donor == pid] |>
        union( c( pid, D$donor[D$recipient==pid] ))

      G2 <- G[ G$pid %in% piapids , ]
      pia <- sum( G2$timeinfected > IT )

      G3 <- G2[ G2$timeinfected <= IT & G2$timediagnosed > IT, ]
      ida <- sum( G3$timediagnosed - IT )

      # Contacts = total_contacts + 1 (self)
      c( pia, ida, total_contacts + 1 )
    }

    mapply(proc_indiv, G1$pid, G1$IT, G1$total_contacts ) -> o
    o <- as.data.frame( t( o ) )
    colnames(o) <- c('pia', 'ida', 'contacts' )
    o$simid <- as.character(G1_simid)  # Add simid for paired comparison

    # Compute summary statistics (same structure as cluster-based)
    # Quantiles: 10th and 90th percentiles used uniformly for all metrics
    e_ida_percontact <- o$ida / o$contacts
    e_ida_valid <- e_ida_percontact[is.finite(e_ida_percontact) & !is.na(e_ida_percontact)]
    med_ida <- if (length(e_ida_valid) > 0) median(e_ida_valid) else NA_real_
    q_ida <- if (length(e_ida_valid) > 0) quantile(e_ida_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
    
    e_pia_percontact <- o$pia / o$contacts
    e_pia_valid <- e_pia_percontact[is.finite(e_pia_percontact) & !is.na(e_pia_percontact)]
    med_pia <- if (length(e_pia_valid) > 0) median(e_pia_valid) else NA_real_
    q_pia <- if (length(e_pia_valid) > 0) quantile(e_pia_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
    
    list(
      o = o,
      propintervened = NA,  # Not applicable for random selection
      n_units = nrow(o),
      ida = c(sum(o$ida),
               sum(o$ida) / sum(o$contacts),
               med_ida,
               q_ida[1],
               q_ida[2]),
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
#' @param partner_notification_window_months Lookback window: 3 or 6 months (default: 6)
#'
#' @return List with same structure as random_intervention
rita_intervention <- function(Dall, Gall, rita_window_months, implementation_delay_days,
                              partner_notification_window_months = 6)
{
  lastgeneration <- max( Gall$generation )
  
  # Simulate RITA test: positive if diagnosed within random window of infection
  # Average window = rita_window_months * 30 days
  Gall$rita <- with(Gall, (timediagnosed - timeinfected) < rexp(nrow(Gall), 1/(rita_window_months*30)))
  
  # Filter to RITA-positive cases in generations 2+ (excluding last)
  G1 <- Gall[Gall$rita & (Gall$generation > 0) & (Gall$generation < lastgeneration), ]
    
    if (nrow(G1) == 0) {
      return(list(
        o = data.frame(pia = numeric(0), ida = numeric(0), contacts = numeric(0)),
        propintervened = 0,
        n_units = 0,
        ida = c(0, 0, 0, 0, 0),
        pia = c(0, 0, 0, 0, 0),
        total_contacts = 0
      ))
    }
    
    # Calculate total contacts based on partner notification window
    if (partner_notification_window_months == 3) {
      G1$total_contacts <- with(G1, Fcontacts_90d + Gcontacts_90d + Hcontacts_90d)
    } else {
      G1$total_contacts <- with(G1, Fcontacts_180d + Gcontacts_180d + Hcontacts_180d)
    }

    # Save simid before modifying pid (for paired comparison analysis)
    G1_simid <- G1$simid

    # Make PIDs unique across simulations
    G1$pid <- paste(sep = '.', G1$pid, G1$simid)
    D <- Dall; D$donor <- paste(sep='.', D$donor, D$simid )
    D$recipient <- paste(sep='.', D$recipient, D$simid )
    G <- Gall; G$pid <- paste(sep='.', G$pid, G$simid)

    # Set intervention time at diagnosis + implementation delay
    G1$IT <- G1$timediagnosed + implementation_delay_days

    # Process individual intervention outcomes
    proc_indiv <- function(pid, IT, total_contacts) {
      piapids <- D$recipient[D$donor == pid] |>
        union(c(pid, D$donor[D$recipient == pid]))

      G2 <- G[G$pid %in% piapids, ]
      pia <- sum(G2$timeinfected > IT)

      G3 <- G2[G2$timeinfected <= IT & G2$timediagnosed > IT, ]
      ida <- sum(G3$timediagnosed - IT)

      c(pia, ida, total_contacts + 1)
    }

    # Process all RITA-positive cases
    results <- matrix(nrow = nrow(G1), ncol = 3)
    for (i in seq_len(nrow(G1))) {
      results[i, ] <- proc_indiv(G1$pid[i], G1$IT[i], G1$total_contacts[i])
    }

    o <- as.data.frame(results)
    colnames(o) <- c('pia', 'ida', 'contacts')
    o$simid <- as.character(G1_simid)  # Add simid for paired comparison

    # Compute summary statistics (same structure as cluster-based)
    # Quantiles: 10th and 90th percentiles used uniformly for all metrics
    e_ida_percontact <- o$ida / o$contacts
    e_ida_valid <- e_ida_percontact[is.finite(e_ida_percontact) & !is.na(e_ida_percontact)]
    med_ida <- if (length(e_ida_valid) > 0) median(e_ida_valid) else NA_real_
    q_ida <- if (length(e_ida_valid) > 0) quantile(e_ida_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
    
    e_pia_percontact <- o$pia / o$contacts
    e_pia_valid <- e_pia_percontact[is.finite(e_pia_percontact) & !is.na(e_pia_percontact)]
    med_pia <- if (length(e_pia_valid) > 0) median(e_pia_valid) else NA_real_
    q_pia <- if (length(e_pia_valid) > 0) quantile(e_pia_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
    
    list(
      o = o,
      propintervened = NA,  # Not directly comparable to cluster-based
      n_units = nrow(o),
      ida = c(sum(o$ida),
               sum(o$ida) / sum(o$contacts),
               med_ida,
               q_ida[1],
               q_ida[2]),
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

  # Calculate total contacts based on partner notification window BEFORE modifying pids
  # Use 90-day columns for 3-month window, 180-day columns for 6-month window
  if (partner_notification_window_months == 3) {
    Gall$total_contacts <- with(Gall, Fcontacts_90d + Gcontacts_90d + Hcontacts_90d)
  } else {
    # Default to 6 months (180 days)
    Gall$total_contacts <- with(Gall, Fcontacts_180d + Gcontacts_180d + Hcontacts_180d)
  }

  # Filter to high-contact individuals in generations 2+ (excluding last) BEFORE modifying pids
  G1_filter <- (Gall$total_contacts >= network_degree_threshold) &
               (Gall$generation > 0) &
               (Gall$generation < lastgeneration)

  if (!any(G1_filter)) {
    return(list(
      o = data.frame(pia = numeric(0), ida = numeric(0), contacts = numeric(0)),
      propintervened = 0,
      n_units = 0,
      ida = c(0, 0, 0, 0, 0),
      pia = c(0, 0, 0, 0, 0),
      total_contacts = 0
    ))
  }

  # Save simid for filtered cases (for paired comparison analysis)
  G1_simid <- Gall$simid[G1_filter]

  # Make PIDs unique across simulations
  D <- Dall
  D$donor <- paste(sep = '.', D$donor, D$simid)
  D$recipient <- paste(sep = '.', D$recipient, D$simid)
  G <- Gall
  G$pid <- paste(sep = '.', G$pid, G$simid)

  # Apply filter after pid modification
  G1 <- G[G1_filter, ]

  # Set intervention time at diagnosis + implementation delay
  G1$IT <- G1$timediagnosed + implementation_delay_days

  # Process individual intervention outcomes
  proc_indiv <- function(pid, IT, total_contacts) {
    piapids <- D$recipient[D$donor == pid] |>
      union(c(pid, D$donor[D$recipient == pid]))

    G2 <- G[G$pid %in% piapids, ]
    pia <- sum(G2$timeinfected > IT)

    G3 <- G2[G2$timeinfected <= IT & G2$timediagnosed > IT, ]
    ida <- sum(G3$timediagnosed - IT)

    c(pia, ida, total_contacts + 1)  # +1 for the index case
  }

  mapply(proc_indiv, G1$pid, G1$IT, G1$total_contacts) -> o
  o <- as.data.frame(t(o))
  colnames(o) <- c('pia', 'ida', 'contacts')
  o$simid <- as.character(G1_simid)  # Add simid for paired comparison

    # Compute summary statistics (same structure as cluster-based)
    # Quantiles: 10th and 90th percentiles used uniformly for all metrics
    e_ida_percontact <- o$ida / o$contacts
    e_ida_valid <- e_ida_percontact[is.finite(e_ida_percontact) & !is.na(e_ida_percontact)]
    med_ida <- if (length(e_ida_valid) > 0) median(e_ida_valid) else NA_real_
    q_ida <- if (length(e_ida_valid) > 0) quantile(e_ida_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
    
    e_pia_percontact <- o$pia / o$contacts
    e_pia_valid <- e_pia_percontact[is.finite(e_pia_percontact) & !is.na(e_pia_percontact)]
    med_pia <- if (length(e_pia_valid) > 0) median(e_pia_valid) else NA_real_
    q_pia <- if (length(e_pia_valid) > 0) quantile(e_pia_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
    
    list(
      o = o,
      propintervened = NA,  # Not directly comparable to cluster-based
      n_units = nrow(o),
      ida = c(sum(o$ida),
               sum(o$ida) / sum(o$contacts),
               med_ida,
               q_ida[1],
               q_ida[2]),
      pia = c(sum(o$pia),
              sum(o$pia) / sum(o$contacts),
              med_pia,
              q_pia[1],
              q_pia[2]),
      total_contacts = sum(o$contacts)
    )
}

#' RITA + Secondary Contact Tracing intervention strategy
#'
#' Combines RITA's early detection with 2-degree contact tracing.
#' Identifies RITA-positive individuals (recently infected) and traces their
#' transmission network backward and forward within the lookback window.
#'
#' Network identification (using transmission tree):
#'   - Primary traced: Individuals with transmission links to/from RITA+ index
#'   - Secondary traced: Individuals with transmission links to/from primary traced
#'   - All transmission events must occur within lookback window
#'   - Transmission links used to IDENTIFY individuals, not count contacts
#'
#' Contact notification (using contact counts):
#'   - Each traced individual has unknown contacts (Fcontacts + Gcontacts + Hcontacts)
#'   - We notify all contacts of all traced individuals
#'   - Small subnetwork: Sum all contact counts (assume no overlap)
#'   - Large subnetwork: Max contact count (assume complete overlap)
#'
#' @param Dall Combined transmission data for all simulations
#' @param Gall Combined individual data for all simulations
#' @param rita_window_months Average RITA detection window in months
#' @param implementation_delay_days Fixed delay (days) for intervention implementation
#' @param partner_notification_window_months Lookback window: 3 or 6 months (default: 6)
#' @param analysis_delay_days Delay for network analysis (default: 0, interviews during implementation)
#'
#' @return List with same structure as cluster-based strategies
rita_secondary_intervention <- function(Dall, Gall, rita_window_months,
                                        implementation_delay_days,
                                        partner_notification_window_months = 6,
                                        analysis_delay_days = 0)
{
  lastgeneration <- max(Gall$generation)
  lookback_days <- partner_notification_window_months * 30

  # Simulate RITA test: positive if diagnosed within random window of infection
  Gall$rita <- with(Gall, (timediagnosed - timeinfected) < rexp(nrow(Gall), 1/(rita_window_months*30)))

  # Filter to RITA-positive cases in generations 1+ (excluding last)
  G_rita <- Gall[Gall$rita & (Gall$generation > 0) & (Gall$generation < lastgeneration), ]

  if (nrow(G_rita) == 0) {
    return(list(
      o = data.frame(pia = numeric(0), ida = numeric(0), interventiontime = numeric(0),
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

  # Save simid before modifying pid (for paired comparison analysis)
  G_rita_simid <- G_rita$simid

  # Make PIDs unique across simulations
  G_rita$pid <- paste(sep = '.', G_rita$pid, G_rita$simid)
  D <- Dall
  D$donor <- paste(sep = '.', D$donor, D$simid)
  D$recipient <- paste(sep = '.', D$recipient, D$simid)
  G <- Gall
  G$pid <- paste(sep = '.', G$pid, G$simid)

  # Calculate total contacts based on partner notification window
  if (partner_notification_window_months == 3) {
    G$total_contacts <- with(G, Fcontacts_90d + Gcontacts_90d + Hcontacts_90d)
  } else {
    G$total_contacts <- with(G, Fcontacts_180d + Gcontacts_180d + Hcontacts_180d)
  }

  # Set intervention time
  G_rita$IT <- G_rita$timediagnosed + analysis_delay_days + implementation_delay_days

  # Helper function: Get transmission-linked individuals within lookback window
  # NOTE: This identifies individuals via transmission tree, NOT their contacts
  # Contacts are unknown - we only know contact COUNTS per person
  get_transmission_linked_individuals <- function(pid, window_start, window_end) {
    # People this person transmitted TO within window
    # (These are a subset of this person's contacts)
    forward <- D$recipient[
      D$donor == pid &
      D$timetransmission >= window_start &
      D$timetransmission <= window_end
    ]

    # People this person was infected BY within window
    # (These are a subset of people who had contact with this person)
    backward <- D$donor[
      D$recipient == pid &
      D$timetransmission >= window_start &
      D$timetransmission <= window_end
    ]

    unique(c(forward, backward))
  }

  # Helper function: Check if person is eligible for intervention at IT
  eligible_for_intervention <- function(person_info, IT) {
    if (nrow(person_info) == 0) return(FALSE)
    person_info$generation[1] > 0 &&
      person_info$generation[1] < lastgeneration &&
      person_info$timediagnosed[1] > IT
  }

  # Process each RITA-positive case
  results <- list()

  for (i in seq_len(nrow(G_rita))) {
    index_pid <- G_rita$pid[i]
    index_dx_time <- G_rita$timediagnosed[i]
    IT <- G_rita$IT[i]

    # Define lookback window anchored to index diagnosis
    lookback_start <- index_dx_time - lookback_days
    lookback_end <- index_dx_time

    # Build traced network using transmission tree
    # NOTE: We use transmission links to IDENTIFY individuals, not to count contacts
    network <- c(index_pid)  # Start with RITA+ index

    # Step 1: Identify primary traced individuals (transmission-linked to index)
    primary_traced <- get_transmission_linked_individuals(index_pid, lookback_start, lookback_end)

    for (primary_pid in primary_traced) {
      primary_info <- G[G$pid == primary_pid, ]
      if (eligible_for_intervention(primary_info, IT)) {
        network <- c(network, primary_pid)

        # Step 2: Identify secondary traced individuals (transmission-linked to primary)
        secondary_traced <- get_transmission_linked_individuals(primary_pid, lookback_start, lookback_end)

        for (secondary_pid in secondary_traced) {
          # Skip if already in network
          if (secondary_pid %in% network) next

          secondary_info <- G[G$pid == secondary_pid, ]
          if (eligible_for_intervention(secondary_info, IT)) {
            network <- c(network, secondary_pid)
          }
        }
      }
    }

    network <- unique(network)

    # Calculate total contacts to notify
    # Each traced individual has their own (unknown) contact list
    # We notify ALL contacts of ALL traced individuals
    network_contacts <- numeric(length(network))
    for (j in seq_along(network)) {
      person <- G[G$pid == network[j], ]
      if (nrow(person) > 0) {
        # total_contacts = Fcontacts_XXd + Gcontacts_XXd + Hcontacts_XXd
        network_contacts[j] <- person$total_contacts[1]
      }
    }

    # Small network: sum all contacts (assume fully unique - no overlap)
    # Total notifications = sum of all contact lists
    contacts_small <- sum(network_contacts)

    # Large network: max contacts (assume fully overlapping/connected)
    # Total unique contacts ≈ max individual contact list
    contacts_large <- max(network_contacts, 0)

    # Calculate PIA and IDA for network
    # PIA/IDA scope: network members + their transmission partners
    # (Following cluster-based strategy pattern)
    piapids <- network
    for (pid in network) {
      # Add all transmission partners of network members
      # (These represent potential infections averted through network intervention)
      piapids <- union(piapids, D$recipient[D$donor == pid])  # People they infected
      piapids <- union(piapids, D$donor[D$recipient == pid])  # People who infected them
    }

    G2 <- G[G$pid %in% piapids, ]

    # PIA: Infections that occurred AFTER intervention time
    pia <- sum(G2$timeinfected > IT)

    # IDA: Infectious days averted
    # (People infected before IT but diagnosed after IT)
    G3 <- G2[G2$timeinfected <= IT & G2$timediagnosed > IT, ]
    ida <- sum(G3$timediagnosed - IT)

    results[[i]] <- data.frame(
      simid = as.character(G_rita_simid[i]),  # Add simid for paired comparison
      pia = pia,
      ida = ida,
      interventiontime = IT,
      nc = length(network),
      sum_degrees = NA_real_,  # Not applicable to RITA+Secondary
      sum_excess = NA_real_,   # Not applicable to RITA+Secondary
      contacts_small = contacts_small,
      contacts_large = contacts_large
    )
  }

  # Combine results - use incremental rbind to avoid do.call issues
  if (length(results) == 0) {
    return(list(
      o = NULL,
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

  odf <- results[[1]]
  if (length(results) > 1) {
    for (i in 2:length(results)) {
      odf <- rbind(odf, results[[i]])
    }
  }
  odf <- as.data.frame(odf, stringsAsFactors = FALSE)

  # Compute summary statistics (matching cluster-based strategies)
  # Quantiles: 10th and 90th percentiles used uniformly for all metrics

  # IDA efficiency for small subnetwork
  e_puta_small <- odf$ida / odf$contacts_small
  e_puta_small_valid <- e_puta_small[is.finite(e_puta_small) & !is.na(e_puta_small)]
  med_ida_small <- if (length(e_puta_small_valid) > 0) median(e_puta_small_valid) else NA_real_
  q_ida_small <- if (length(e_puta_small_valid) > 0) quantile(e_puta_small_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # IDA efficiency for large subnetwork
  e_puta_large <- odf$ida / odf$contacts_large
  e_puta_large_valid <- e_puta_large[is.finite(e_puta_large) & !is.na(e_puta_large)]
  med_ida_large <- if (length(e_puta_large_valid) > 0) median(e_puta_large_valid) else NA_real_
  q_ida_large <- if (length(e_puta_large_valid) > 0) quantile(e_puta_large_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # PIA efficiency for small subnetwork
  e_pia_small <- odf$pia / odf$contacts_small
  e_pia_small_valid <- e_pia_small[is.finite(e_pia_small) & !is.na(e_pia_small)]
  med_pia_small <- if (length(e_pia_small_valid) > 0) median(e_pia_small_valid) else NA_real_
  q_pia_small <- if (length(e_pia_small_valid) > 0) quantile(e_pia_small_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # PIA efficiency for large subnetwork
  e_pia_large <- odf$pia / odf$contacts_large
  e_pia_large_valid <- e_pia_large[is.finite(e_pia_large) & !is.na(e_pia_large)]
  med_pia_large <- if (length(e_pia_large_valid) > 0) median(e_pia_large_valid) else NA_real_
  q_pia_large <- if (length(e_pia_large_valid) > 0) quantile(e_pia_large_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # Return results
  list(
    o = odf,
    propintervened = nrow(odf) / nrow(Gall[Gall$generation > 0 & Gall$generation < lastgeneration, ]),
    n_units = nrow(odf),
    puta_small = c(sum(odf$ida),
         if (sum(odf$contacts_small) > 0) sum(odf$ida) / sum(odf$contacts_small) else NA_real_,
         med_ida_small,
         q_ida_small[1],
         q_ida_small[2]),
    puta_large = c(sum(odf$ida),
         if (sum(odf$contacts_large) > 0) sum(odf$ida) / sum(odf$contacts_large) else NA_real_,
         med_ida_large,
         q_ida_large[1],
         q_ida_large[2]),
    pia_small = c(sum(odf$pia),
         if (sum(odf$contacts_small) > 0) sum(odf$pia) / sum(odf$contacts_small) else NA_real_,
         med_pia_small,
         q_pia_small[1],
         q_pia_small[2]),
    pia_large = c(sum(odf$pia),
         if (sum(odf$contacts_large) > 0) sum(odf$pia) / sum(odf$contacts_large) else NA_real_,
         med_pia_large,
         q_pia_large[1],
         q_pia_large[2]),
    total_contacts_small = sum(odf$contacts_small),
    total_contacts_large = sum(odf$contacts_large)
  )
}


# =============================================================================
# Main execution: Run interventions and generate plots
# =============================================================================
# Only runs when this script is executed directly (not when sourced)

if (sys.nframe() == 0) {
  # Source the plotting functions
  source("src/plot_interventions.R")
  
  # Run interventions with 6-month partner notification window
  cat("Running interventions...\n")
  results <- run_intervention_analysis(
    d_file = "src/experiment1-N10000-gens7-D.csv",
    g_file = "src/experiment1-N10000-gens7-G.csv",
    partner_notification_window_months = 6
  )
  
  # Generate and save violin plots
  cat("\nGenerating efficiency distribution plots...\n")
  p_violin <- plot_efficiency_distributions(results)
  plot_dir <- here::here("intervention-plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  ggsave(file.path(plot_dir, "efficiency_distributions_violin.pdf"), p_violin, width = 12, height = 10)
  cat("Saved: intervention-plots/efficiency_distributions_violin.pdf\n")
  
  cat("\nDone!\n")
}
