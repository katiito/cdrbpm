invisible('
experiment1 intervention model 
')

library(glue)
library(ggplot2)
## Path resolution using the 'here' package (required)
if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required to run this script. Install it with install.packages('here').")
}

## Resolve data file paths anchored at the project root via here::here()
## - Absolute paths are used as-is
## - For bare filenames, prefer here('src', file) if it exists, else here(file)
## - For relative paths with directories (e.g. 'data/x.csv'), use here(file)
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


# Main function to run intervention analysis
run_intervention_analysis <- function(
  d_file = 'experiment1-N10000-gens7-D.csv',
  g_file = 'experiment1-N10000-gens7-G.csv',
  seed = NULL,  # Default to random seed
  cluster_size_5 = 5,
  cluster_size_2 = 2,
  distance_threshold = 0.005,
  network_degree_threshold = 4,
  random_sample_size = 30,
  rita_window_months = 6, # average RITA detection window
  lookback_window_months = 6, # growth-rate trigger window
  growth_distance_threshold = 0.1, # separate D for growth-based trigger
  intervention_rate = 1/90, # average 3 months to intervention
  show_table = TRUE,
  show_debug = FALSE,
  output_dir = "intervention-results"  # if non-NULL, save results to this directory
) {
    # Seed
    if (is.null(seed)) {
      seed <- as.numeric(Sys.time())
      cat("Using random seed:", seed, "\n")
    } else {
      cat("Using fixed seed:", seed, "\n")
    }
    set.seed(seed)

    # Load data
    cat("Loading data...\n")
    Dall <- read.csv(resolve_input_path(d_file), stringsAs = FALSE)
    Gall <- read.csv(resolve_input_path(g_file), stringsAs = FALSE)

    # Split by simid
    simids <- unique(Dall$simid)
    Ds <- split(Dall, Dall$simid)[simids]
    Gs <- split(Gall, Gall$simid)[simids]
    cat("  Loaded", nrow(Dall), "D rows,", nrow(Gall), "G rows across", length(simids), "simulations\n")

    # Run interventions
    cat("Running interventions (6 strategies)...\n")
    t_start <- Sys.time()

    cat("  [1/6] Distance-size (size=", cluster_size_5, ", D=", distance_threshold, ")...", sep = "")
    ods5 <- tryCatch(
      distsize_intervention(
        Ds = Ds, Gs = Gs, Gall = Gall,
        distance_threshold = distance_threshold,
        cluster_size = cluster_size_5,
        intervention_rate = intervention_rate
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta_small = c(0, 0, 0, 0, 0), puta_large = c(0, 0, 0, 0, 0),
             pia = c(0, 0, 0), total_contacts_small = 0, total_contacts_large = 0)
      }
    )
    cat(" done (", ods5$n_units, " units)\n", sep = "")

    cat("  [2/6] Distance-size (size=", cluster_size_2, ", D=", distance_threshold, ")...", sep = "")
    ods2 <- tryCatch(
      distsize_intervention(
        Ds = Ds, Gs = Gs, Gall = Gall,
        distance_threshold = distance_threshold,
        cluster_size = cluster_size_2,
        intervention_rate = intervention_rate
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta_small = c(0, 0, 0, 0, 0), puta_large = c(0, 0, 0, 0, 0),
             pia = c(0, 0, 0), total_contacts_small = 0, total_contacts_large = 0)
      }
    )
    cat(" done (", ods2$n_units, " units)\n", sep = "")

    cat("  [3/6] Growth-rate (size=", cluster_size_5, ", W=", lookback_window_months, "mo, D=", growth_distance_threshold, ")...", sep = "")
    ogrowth <- tryCatch(
      growthrate_intervention(
        Ds = Ds, Gs = Gs, Gall = Gall,
        growth_distance_threshold = growth_distance_threshold,
        cluster_size = cluster_size_5,
        lookback_window_months = lookback_window_months,
        intervention_rate = intervention_rate
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta_small = c(0, 0, 0, 0, 0), puta_large = c(0, 0, 0, 0, 0),
             pia = c(0, 0, 0), total_contacts_small = 0, total_contacts_large = 0)
      }
    )
    cat(" done (", ogrowth$n_units, " units)\n", sep = "")

    cat("  [4/6] Random allocation (n=", random_sample_size, ")...", sep = "")
    orand <- tryCatch(
      random_intervention(
        Dall = Dall, Gall = Gall,
        random_sample_size = random_sample_size,
        intervention_rate = intervention_rate
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta = c(0, 0, 0, 0, 0), pia = c(0, 0, 0), total_contacts = 0)
      }
    )
    cat(" done (", orand$n_units, " units)\n", sep = "")

    cat("  [5/6] RITA (window=", rita_window_months, "mo)...", sep = "")
    orita <- tryCatch(
      rita_intervention(
        Dall = Dall, Gall = Gall,
        rita_window_months = rita_window_months,
        intervention_rate = intervention_rate
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta = c(0, 0, 0, 0, 0), pia = c(0, 0, 0), total_contacts = 0)
      }
    )
    cat(" done (", orita$n_units, " units)\n", sep = "")

    cat("  [6/6] Network degree (threshold=", network_degree_threshold, ")...", sep = "")
    onet <- tryCatch(
      network_intervention(
        Dall = Dall, Gall = Gall,
        network_degree_threshold = network_degree_threshold,
        intervention_rate = intervention_rate
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta = c(0, 0, 0, 0, 0), pia = c(0, 0, 0), total_contacts = 0)
      }
    )
    cat(" done (", onet$n_units, " units)\n", sep = "")

    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat("All interventions completed in", round(elapsed, 1), "seconds\n\n")

    # Compile summary table with both subnetwork calculations
    # For cluster-based strategies (ods5, ods2, ogrowth), we have puta_small and puta_large
    # For individual-based strategies (orand, orita, onet), we only have puta (single value)
    
    strategy_names <- c(
      paste0('Size=', cluster_size_5, ',D=', distance_threshold),
      paste0('Size=', cluster_size_2, ',D=', distance_threshold),
      paste0('Growth, size=', cluster_size_5, ',W=', lookback_window_months, 'mo,D=', growth_distance_threshold),
      'Random allocation',
      'RITA',
      paste0('Network, partners>', network_degree_threshold)
    )
    
    # Build summary with both small and large subnetwork columns
    odf <- rbind(
      # Cluster-based strategies have both subnetwork calculations
      c(ods5$total_contacts_small, ods5$total_contacts_large, ods5$puta_small, ods5$puta_large, ods5$pia),
      c(ods2$total_contacts_small, ods2$total_contacts_large, ods2$puta_small, ods2$puta_large, ods2$pia),
      c(ogrowth$total_contacts_small, ogrowth$total_contacts_large, ogrowth$puta_small, ogrowth$puta_large, ogrowth$pia),
      # Individual-based strategies: total_contacts is same for both, puta also same
      c(orand$total_contacts, orand$total_contacts, orand$puta, orand$puta, orand$pia),
      c(orita$total_contacts, orita$total_contacts, orita$puta, orita$puta, orita$pia),
      c(onet$total_contacts, onet$total_contacts, onet$puta, onet$puta, onet$pia)
    ) |> as.data.frame()
    
    colnames(odf) <- c(
      "Contacts_Small", "Contacts_Large",
      "Total_PUTA", "PUTA/contact_small", "Median_small", "Low_small", "High_small",
      "Total_PUTA_lg", "PUTA/contact_large", "Median_large", "Low_large", "High_large",
      "PIA", "PIA_Low", "PIA_High"
    )
    rownames(odf) <- strategy_names
    odf1 <- round(odf, 2)

    counts_df <- data.frame(
      Strategy = strategy_names,
      Units = c(ods5$n_units, ods2$n_units, ogrowth$n_units, orand$n_units, orita$n_units, onet$n_units),
      Contacts_Small = c(ods5$total_contacts_small, ods2$total_contacts_small, ogrowth$total_contacts_small,
                         orand$total_contacts, orita$total_contacts, onet$total_contacts),
      Contacts_Large = c(ods5$total_contacts_large, ods2$total_contacts_large, ogrowth$total_contacts_large,
                         orand$total_contacts, orita$total_contacts, onet$total_contacts)
    )

    # Save outputs to files if output_dir is specified
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        cat("Created output directory:", output_dir, "\n")
      }
      
      # Generate timestamp for unique filenames
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      
      # Save summary table (wide format with Strategy as rownames)
      summary_df <- cbind(Strategy = rownames(odf), odf)
      rownames(summary_df) <- NULL
      summary_path <- file.path(output_dir, paste0("summary_", timestamp, ".csv"))
      write.csv(summary_df, summary_path, row.names = FALSE)
      
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
                      "growth_distance_threshold", "intervention_rate", "n_simulations"),
        value = c(d_file, g_file, as.character(seed), cluster_size_5, cluster_size_2,
                  distance_threshold, network_degree_threshold,
                  random_sample_size, rita_window_months, lookback_window_months,
                  growth_distance_threshold, intervention_rate, length(simids))
      )
      write.csv(params, file.path(output_dir, paste0("parameters_", timestamp, ".csv")), row.names = FALSE)
      
      cat("\nOutputs saved to:", output_dir, "\n")
      cat("  - summary_", timestamp, ".csv (aggregated metrics)\n", sep = "")
      cat("  - counts_", timestamp, ".csv (units and contacts per strategy)\n", sep = "")
      cat("  - details_*_", timestamp, ".csv (per-simulation/unit details)\n", sep = "")
      cat("  - parameters_", timestamp, ".csv (run parameters for reproducibility)\n", sep = "")
    }

    # Display results
    if (show_table) {
      if (require(knitr, quietly = TRUE)) {
        print(knitr::kable(odf1))
        cat("\nSample sizes and totals (for context):\n")
        print(knitr::kable(counts_df))
        if (show_debug) {
          cat("\nDebug: Cluster-size-5 per-simulation details -> simid, nc, contacts_small, contacts_large, interventiontime, pia, puta\n")
          if (!is.null(ods5$o) && nrow(ods5$o) > 0) {
            debug_cols <- c("simid", "nc", "contacts_small", "contacts_large", "interventiontime", "pia", "puta")
            odf1_5_view <- ods5$o[, debug_cols]
            num_cols_all <- intersect(c("nc", "contacts_small", "contacts_large", "interventiontime", "pia", "puta"), colnames(odf1_5_view))
            odf1_5_view[num_cols_all] <- lapply(odf1_5_view[num_cols_all], function(x) suppressWarnings(as.numeric(x)))
            num_cols <- intersect(c("interventiontime", "pia", "puta", "contacts_small", "contacts_large"), colnames(odf1_5_view))
            odf1_5_view[num_cols] <- lapply(odf1_5_view[num_cols], function(x) if (is.numeric(x)) round(x, 2) else x)
            print(knitr::kable(odf1_5_view))
            cat("\nDebug: Degree breakdown at IT for size-5 clusters\n")
            for (i in seq_len(nrow(ods5$o))) {
              row <- ods5$o[i, ]
              sim_key <- as.character(row$simid)
              IT <- as.numeric(row$interventiontime)
              if (!is.finite(IT)) next
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
              tc_small <- n + sum(excess)
              tc_large <- sum_deg - (n - 2)
              cat(sprintf("simid=%s n=%d contacts_small=%.2f contacts_large=%.2f degrees=%s excess=%s\n",
                          sim_key, n, tc_small, tc_large,
                          paste(round(G1$degree, 2), collapse=","),
                          paste(round(excess, 2), collapse=",")))
            }
          } else {
            cat("(No size-5 clusters met criteria; odf1 is empty.)\n")
          }
        }
      } else {
        cat("knitr package not available, showing table as dataframe\n")
        print(odf1)
        cat("\nSample sizes and totals (for context):\n")
        print(counts_df)
        if (show_debug) {
          cat("\nDebug: Cluster-size-5 per-simulation details -> simid, nc, contacts_small, contacts_large, interventiontime, pia, puta\n")
          if (!is.null(ods5$o) && nrow(ods5$o) > 0) {
            debug_cols <- c("simid", "nc", "contacts_small", "contacts_large", "interventiontime", "pia", "puta")
            odf1_5_view <- ods5$o[, debug_cols]
            num_cols_all <- intersect(c("nc", "contacts_small", "contacts_large", "interventiontime", "pia", "puta"), colnames(odf1_5_view))
            odf1_5_view[num_cols_all] <- lapply(odf1_5_view[num_cols_all], function(x) suppressWarnings(as.numeric(x)))
            num_cols <- intersect(c("interventiontime", "pia", "puta", "contacts_small", "contacts_large"), colnames(odf1_5_view))
            odf1_5_view[num_cols] <- lapply(odf1_5_view[num_cols], function(x) if (is.numeric(x)) round(x, 2) else x)
            print(odf1_5_view)
            cat("\nDebug: Degree breakdown at IT for size-5 clusters\n")
            for (i in seq_len(nrow(ods5$o))) {
              row <- ods5$o[i, ]
              sim_key <- as.character(row$simid)
              IT <- as.numeric(row$interventiontime)
              if (!is.finite(IT)) next
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
              tc_small <- n + sum(excess)
              tc_large <- sum_deg - (n - 2)
              cat(sprintf("simid=%s n=%d contacts_small=%.2f contacts_large=%.2f degrees=%s excess=%s\n",
                          sim_key, n, tc_small, tc_large,
                          paste(round(G1$degree, 2), collapse=","),
                          paste(round(excess, 2), collapse=",")))
            }
          } else {
            cat("(No size-5 clusters met criteria; odf1 is empty.)\n")
          }
        }
      }
    }
  }

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
###### ### ### ### ###  INTERVENTION FUNCTIONS ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
###### ### ### ### ###  INTERVENTION FUNCTIONS ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 



# Process individual cluster for intervention analysis
proc_cluster <- function(
    D, G, 
    distance_threshold,
    cluster_size, 
    intervention_rate
) 
{
  # (growth threshold not implemented)
  
  # Sort data by relevant time variables
  G <- G[ order(G$timesequenced), ]
  D <- D[ order(D$timetransmission),]
  
  lastgeneration <- max( G$generation )
  
  # Retain cases with path to patient 0, within distance threshold, excluding last generation
  D1 <- D[ D$distance <= distance_threshold , ]
  keeppids <- "0" 
  addpids  <- D1$recipient[ D1$donor %in% keeppids ]
  
  while( length( addpids ) > 0 ){
    keeppids <- union( addpids, keeppids )
    addpids  <- setdiff( D1$recipient[ D1$donor %in% keeppids ], keeppids )
  }
  
  G1 <- G[ G$pid %in% keeppids , ]
  G1 <- G1[ G1$generation != lastgeneration, ]
  D1 <- D1[ D1$donor %in% G1$pid & D1$recipient %in% G1$pid, ]
  
  # find intervention time, if there is one 
  IT <- Inf 
  if (nrow(G1) > 0 && nrow(G1) >= cluster_size) {
    idx <- cluster_size
    IT  <- G1$timesequenced[idx] + rexp(1, rate = intervention_rate)
  }
  
  # Exclude clustered infections detected after intervention time
  G1 <- G1[ G1$timesequenced < IT , ]
  
  # Calculate total cluster contact network size (both small and large)
  sum_degrees <- 0
  sum_excess <- 0
  contacts_small <- 0
  contacts_large <- 0
  
  if (nrow(G1) > 0) {
    # Get cluster members' total degrees
    G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
    sum_degrees <- sum(G1$degree)
    n <- nrow(G1)
    
    # Compute both subnetwork assumptions
    excess <- pmax(G1$degree - (n - 1), 0)
    sum_excess <- sum(excess)
    contacts_small <- n + sum_excess          # dense connections inside cluster
    contacts_large <- sum_degrees - (n - 2)   # sparse connections inside cluster
  }
  
  # Calculate potential infections averted (PIA) and person-years untreated averted (PUTA)
  
  pia <- 0 
  puta <- 0 
  
  # calculate transmissions from and to D (this is different to original code where final union was omitted)
  if (!is.infinite( IT ) )
  {
    piapids <-  D$recipient[D$donor %in% G1$pid] |> 
      union(G1$pid) |> 
      union(D$donor[D$recipient %in% G1$pid]) 
    
    G2  <- G[ G$pid %in% piapids , ]
    pia <- sum( G2$timeinfected  > IT )
    
    G3   <- G2[ G2$timeinfected <= IT & G2$timediagnosed > IT , ]
    puta <- sum( G3$timediagnosed - IT )
  }
  c(pia = pia, puta = puta, interventiontime = IT, nc = nrow(G1), 
    sum_degrees = sum_degrees, sum_excess = sum_excess,
    contacts_small = contacts_small, contacts_large = contacts_large)
}

  
### ### ### # Distance-size intervention strategy

distsize_intervention <- function(
    Ds, Gs, Gall,
    distance_threshold, cluster_size, intervention_rate
)
{
  
  # Process all simulations
  o <- lapply(seq_along(Ds), function(i) {
    tryCatch({
      proc_cluster(
        D = Ds[[i]], G = Gs[[i]],
        distance_threshold = distance_threshold,
        cluster_size = cluster_size,
        intervention_rate = intervention_rate
      )
    }, error = function(e) {
      # Return default values if processing fails
      c(pia = 0, puta = 0, interventiontime = Inf, nc = 0, 
        sum_degrees = 0, sum_excess = 0, contacts_small = 0, contacts_large = 0)
    })
  })
  
  odf  <- do.call( rbind, o ) |> as.data.frame()
  # Attach simulation identifiers to each row for debugging/traceability
  simids_vec <- names(Ds)
  odf$simid <- if (!is.null(simids_vec)) simids_vec else seq_len(nrow(odf))
  # Reorder columns to bring simid upfront
  odf <- odf[, c("simid", "pia", "puta", "interventiontime", "nc", 
                 "sum_degrees", "sum_excess", "contacts_small", "contacts_large")]

  # Sanity-correct nc and contacts by recomputing deterministically at the recorded IT
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
  odf1 <- odf[ !is.infinite( odf$interventiontime ), ] # exclude sims where no cluster found 
  
  lastgen <- max( Gall$generation )
  
  if (nrow(odf1) == 0) {
    return(list(
      o = data.frame(pia = numeric(0), puta = numeric(0), interventiontime = numeric(0), 
                     nc = numeric(0), sum_degrees = numeric(0), sum_excess = numeric(0),
                     contacts_small = numeric(0), contacts_large = numeric(0)),
      propintervened = 0,
      n_units = 0,
      puta_small = c(0, 0, 0, 0, 0),
      puta_large = c(0, 0, 0, 0, 0),
      pia = c(0, 0, 0),
      total_contacts_small = 0,
      total_contacts_large = 0
    ))
  }
  
  # Compute PUTA efficiency for small subnetwork
  e_puta_small <- odf1$puta / odf1$contacts_small
  e_puta_small_valid <- e_puta_small[is.finite(e_puta_small) & !is.na(e_puta_small)]
  med_puta_small <- if (length(e_puta_small_valid) > 0) median(e_puta_small_valid) else NA_real_
  q_puta_small <- if (length(e_puta_small_valid) > 0) quantile(e_puta_small_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
  
  # Compute PUTA efficiency for large subnetwork
  e_puta_large <- odf1$puta / odf1$contacts_large
  e_puta_large_valid <- e_puta_large[is.finite(e_puta_large) & !is.na(e_puta_large)]
  med_puta_large <- if (length(e_puta_large_valid) > 0) median(e_puta_large_valid) else NA_real_
  q_puta_large <- if (length(e_puta_large_valid) > 0) quantile(e_puta_large_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
  
  sort_pia <- sort(odf1$pia)
  
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
    pia = c(mean(odf1$pia), 
            sort_pia[max(1, ceiling(0.1 * nrow(odf1)))], 
            sort_pia[min(nrow(odf1), floor(0.9 * nrow(odf1)))]),
    total_contacts_small = sum(odf1$contacts_small),
    total_contacts_large = sum(odf1$contacts_large)
  )
}

### ### ### # Growth-rate cluster intervention

growthrate_intervention <- function(
  Ds, Gs, Gall,
  growth_distance_threshold, cluster_size,
  lookback_window_months = 3,
  intervention_rate
)
{
  lookback_days <- lookback_window_months * 30

  process_one <- function(D, G) {
    # Order by infection time for growth trigger logic
    G <- G[order(G$timeinfected), ]
    D <- D[order(D$timetransmission), ]

    lastgeneration <- max(G$generation)

    # Build cluster reachable from seed via edges within growth_distance_threshold
    D1 <- D[D$distance <= growth_distance_threshold, ]
    keeppids <- "0"
    addpids <- D1$recipient[D1$donor %in% keeppids]
    while (length(addpids) > 0) {
      keeppids <- union(addpids, keeppids)
      addpids <- setdiff(D1$recipient[D1$donor %in% keeppids], keeppids)
    }
    Gcluster <- G[G$pid %in% keeppids, ]

    # Use generation > 0 for growth detection; include last generation for trigger decision
    Gtrig <- Gcluster[Gcluster$generation > 0, ]
    if (nrow(Gtrig) == 0) {
      return(c(pia = 0, puta = 0, interventiontime = Inf, nc = 0, 
               sum_degrees = 0, sum_excess = 0, contacts_small = 0, contacts_large = 0))
    }

    # Sliding window over infection times: earliest time t where count in [t-lookback_days, t] >= cluster_size
    t <- suppressWarnings(as.numeric(Gtrig$timeinfected))
    if (all(!is.finite(t))) {
      # No usable infection times; cannot trigger growth
      return(c(pia = 0, puta = 0, interventiontime = Inf, nc = 0, 
               sum_degrees = 0, sum_excess = 0, contacts_small = 0, contacts_large = 0))
    }
    # Keep only finite infection times for windowing
    keep <- is.finite(t)
    t <- t[keep]
    Gtrig <- Gtrig[keep, ]
    n <- length(t)
    j <- 1
    t_detect <- Inf
    for (i in seq_len(n)) {
      while (j <= i && (t[i] - t[j]) > lookback_days) {
        j <- j + 1
      }
      if ((i - j + 1) >= cluster_size) {
        t_detect <- t[i]
        break
      }
    }

    IT <- if (is.finite(t_detect)) t_detect + rexp(1, rate = intervention_rate) else Inf

    # Define cluster members at IT for outcomes: exclude last generation, require timesequenced < IT
    G1 <- Gcluster[Gcluster$generation != lastgeneration & Gcluster$timesequenced < IT, ]

    # Compute contacts within cluster - both subnetwork assumptions
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

    # PIA/PUTA consistent with other cluster strategy
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

  # Coerce numeric columns
  suppressWarnings({
    odf$interventiontime <- as.numeric(odf$interventiontime)
    odf$nc <- as.numeric(odf$nc)
    odf$sum_degrees <- as.numeric(odf$sum_degrees)
    odf$sum_excess <- as.numeric(odf$sum_excess)
    odf$contacts_small <- as.numeric(odf$contacts_small)
    odf$contacts_large <- as.numeric(odf$contacts_large)
  })

  odf1 <- odf[!is.infinite(odf$interventiontime), ]
  lastgen <- max(Gall$generation)

  if (nrow(odf1) == 0) {
    return(list(
      o = data.frame(pia = numeric(0), puta = numeric(0), interventiontime = numeric(0),
                     nc = numeric(0), sum_degrees = numeric(0), sum_excess = numeric(0),
                     contacts_small = numeric(0), contacts_large = numeric(0)),
      propintervened = 0,
      n_units = 0,
      puta_small = c(0, 0, 0, 0, 0),
      puta_large = c(0, 0, 0, 0, 0),
      pia = c(0, 0, 0),
      total_contacts_small = 0,
      total_contacts_large = 0
    ))
  }

  # Compute PUTA efficiency for small subnetwork
  e_puta_small <- odf1$puta / odf1$contacts_small
  e_puta_small_valid <- e_puta_small[is.finite(e_puta_small) & !is.na(e_puta_small)]
  med_puta_small <- if (length(e_puta_small_valid) > 0) median(e_puta_small_valid) else NA_real_
  q_puta_small <- if (length(e_puta_small_valid) > 0) quantile(e_puta_small_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # Compute PUTA efficiency for large subnetwork
  e_puta_large <- odf1$puta / odf1$contacts_large
  e_puta_large_valid <- e_puta_large[is.finite(e_puta_large) & !is.na(e_puta_large)]
  med_puta_large <- if (length(e_puta_large_valid) > 0) median(e_puta_large_valid) else NA_real_
  q_puta_large <- if (length(e_puta_large_valid) > 0) quantile(e_puta_large_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  sort_pia <- sort(odf1$pia)

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
    pia = c(mean(odf1$pia),
            sort_pia[max(1, ceiling(0.1 * nrow(odf1)))],
            sort_pia[min(nrow(odf1), floor(0.9 * nrow(odf1)))]),
    total_contacts_small = sum(odf1$contacts_small),
    total_contacts_large = sum(odf1$contacts_large)
  )
}

  invisible( '
  	select random case from generation after seed and before last generation 
  	sim intervention time 
  	trace one generation + donor 
  ')
  
  random_intervention <- function(Dall, Gall, random_sample_size, intervention_rate)
  {
    lastgeneration <- max( Gall$generation )
    G1 <- Gall[ Gall$generation > 0 & Gall$generation < lastgeneration, ] # gen 2+ 
    
    if (nrow(G1) == 0) {
      return(list(
        o = data.frame(pia = numeric(0), puta = numeric(0)),
        propintervened = 0,
        n_units = 0,
        puta = c(0, 0, 0, 0, 0),
        pia = c(0, 0, 0),
        total_contacts = 0
      ))
    }
    
    # Randomly sample cases
    sample_size <- min(random_sample_size, nrow(G1))
    G1 <- G1 |> dplyr::slice_sample(n = sample_size)
    
    # Compute weighted degree before making PIDs unique
    G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
    
    
    # make pids unique 
    G1$pid <- paste(sep='.', G1$pid, G1$simid )
    D <- Dall; D$donor <- paste(sep='.', D$donor, D$simid )
    D$recipient <- paste(sep='.', D$recipient, D$simid )
    G <- Gall; G$pid <- paste(sep='.', G$pid, G$simid)
    
    # Set intervention time
    # G1$IT <-  G1$timesequenced + ritdist( nrow(G1 )) 
    G1$IT <- G1$timesequenced + rexp(nrow(G1), rate = intervention_rate)
    
    # Process individual intervention
    proc_indiv <- function(pid, IT, degree )
    {
      piapids <- D$recipient[ D$donor == pid] |> 
        union( c( pid, D$donor[D$recipient==pid] )) 
      
      G2 <- G[ G$pid %in% piapids , ]
      pia <- sum( G2$timeinfected > IT )
      
      G3 <- G2[ G2$timeinfected <= IT & G2$timediagnosed > IT, ]
      puta <- sum( G3$timediagnosed - IT )
      
      c( pia, puta,degree + 1 )
    }
    mapply(proc_indiv, G1$pid, G1$IT, G1$degree ) -> o 
    o <- as.data.frame( t( o ) )
    colnames(o) <- c('pia', 'puta', 'contacts' )
    
    
  e_puta_percontact <- o$puta / o$contacts
  med_puta_percontact <- median(e_puta_percontact, na.rm = TRUE)
  q_puta_percontact <- quantile(e_puta_percontact, probs = c(0.1, 0.9), na.rm = TRUE, names = FALSE)
    sort_pia <- sort(o$pia)
    
    
    list(
      o = o,
      propintervened = NA,
      n_units = nrow(o),
  puta = c(sum(o$puta),
       sum(o$puta) / sum(o$contacts),
       med_puta_percontact,
       q_puta_percontact[1],
       q_puta_percontact[2]),
      pia = c(mean(o$pia), 
              sort_pia[max(1, ceiling(0.01 * nrow(G1)))], 
              sort_pia[min(nrow(G1), floor(0.99 * nrow(G1)))]),
      total_contacts = sum(o$contacts)
    )
  }
  
  invisible('
  Intervene on cases with acute infection 
  Simulated RITA test
  6 month average detection window 
  	  ')
  rita_intervention <- function(Dall, Gall, rita_window_months, intervention_rate)
  {
    lastgeneration <- max( Gall$generation )
    
    # RITA test simulation - 6 month average detection window
    Gall$rita <- with(Gall, (timediagnosed - timeinfected) < rexp(nrow(Gall), 1/(rita_window_months*30)))
    
    # Filter generations & only RITA positive cases
    G1 <- Gall[Gall$rita & (Gall$generation > 0) & (Gall$generation < lastgeneration), ]
    
    if (nrow(G1) == 0) {
      return(list(
        o = data.frame(pia = numeric(0), puta = numeric(0), contacts = numeric(0)),
        propintervened = 0,
        n_units = 0,
        puta = c(0, 0, 0, 0, 0),
        pia = c(0, 0, 0),
        total_contacts = 0
      ))
    }
    
    # Compute weighted degree for RITA-positive individuals
    G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
    
    
    
    # make pids unique 
    G1$pid <- paste(sep = '.', G1$pid, G1$simid)
    D <- Dall; D$donor <- paste(sep='.', D$donor, D$simid )
    D$recipient <- paste(sep='.', D$recipient, D$simid )
    G <- Gall; G$pid <- paste(sep='.', G$pid, G$simid)
    
    
    
    G1$IT <- G1$timesequenced + rexp(nrow(G1), rate = intervention_rate)
    
    # Process individual intervention
    proc_indiv <- function(pid, IT, degree) {
      piapids <- D$recipient[D$donor == pid] |> 
        union(c(pid, D$donor[D$recipient == pid]))
      
      G2 <- G[G$pid %in% piapids, ]
      pia <- sum(G2$timeinfected > IT)
      
      G3 <- G2[G2$timeinfected <= IT & G2$timediagnosed > IT, ]
      puta <- sum(G3$timediagnosed - IT)
      
      c(pia, puta, degree + 1)
    }
    
    # Process all cases
    results <- matrix(nrow = nrow(G1), ncol = 3)
  for (i in seq_len(nrow(G1))) {
      results[i, ] <- proc_indiv(G1$pid[i], G1$IT[i], G1$degree[i])
    }
    
    o <- as.data.frame(results)
    colnames(o) <- c('pia', 'puta', 'contacts')
    
    # Need to make the pia and puta calculations consistent (either efficiency or total)
  e_puta_percontact <- o$puta / o$contacts
  med_puta_percontact <- median(e_puta_percontact, na.rm = TRUE)
  q_puta_percontact <- quantile(e_puta_percontact, probs = c(0.1, 0.9), na.rm = TRUE, names = FALSE)
    sort_pia <- sort(o$pia)
    
    list(
      o = o,
      propintervened = NA,
      n_units = nrow(o),
  puta = c(sum(o$puta),
       sum(o$puta) / sum(o$contacts),
       med_puta_percontact,
       q_puta_percontact[1],
       q_puta_percontact[2]),
      pia = c(mean(o$pia), 
              sort_pia[max(1, ceiling(0.01 * nrow(o)))], 
              sort_pia[min(nrow(o), floor(0.99 * nrow(o)))]),
      total_contacts = sum(o$contacts)
    )
  }
  
  
  network_intervention <- function(Dall, Gall, network_degree_threshold, intervention_rate)
  {
    # do not use weighting and instead just use contact numbers if 
    # w <- 7.0 / 2.0 # ratio duration F to G 
    # ww <-  7.0*30.0 # if unit daiy oo contact rate, E partners in 7 months (F duration) 
    
    lastgeneration <- max(Gall$generation)
    
    # Make pids unique across simulations
    D <- Dall
    D$donor <- paste(sep = '.', D$donor, D$simid)
    D$recipient <- paste(sep = '.', D$recipient, D$simid)
    G <- Gall
    G$pid <- paste(sep = '.', G$pid, G$simid)
    
    # compute weighted degree ( exp partners over 7 months )
    # G$degree  <-  with( G, Fdegree + w*Gdegree + ww*Hdegree)
    G$degree <- with(G, Fdegree + Gdegree + Hdegree)
    
    # filter generations & only degree above threshold  
    G1 <- G[ (G$degree>=network_degree_threshold) & (G$generation > 0) & (G$generation < lastgeneration), ] # gen 2+ 
    
    if (nrow(G1) == 0) {
      return(list(
        o = data.frame(pia = numeric(0), puta = numeric(0), contacts = numeric(0)),
        propintervened = 0,
        n_units = 0,
        puta = c(0, 0, 0, 0, 0),
        pia = c(0, 0, 0),
        total_contacts = 0
      ))
    }
    
    
    
    # Set intervention time
    G1$IT <- G1$timesequenced + rexp(nrow(G1), rate = intervention_rate)
    
    # Process individual intervention
    
    proc_indiv <- function(pid, IT, degree )
    {
      piapids <- D$recipient[D$donor == pid] |> 
        union(c(pid, D$donor[D$recipient == pid]))
      
      G2 <- G[G$pid %in% piapids, ]
      pia <- sum(G2$timeinfected > IT)
      
      G3 <- G2[G2$timeinfected <= IT & G2$timediagnosed > IT, ]
      puta <- sum(G3$timediagnosed - IT)
      
      c(pia, puta, degree + 1)
    }
    mapply(proc_indiv, G1$pid, G1$IT, G1$degree ) -> o 
    o <- as.data.frame( t( o ) )
    colnames(o) <- c('pia', 'puta','contacts' )
  e_puta_percontact <- o$puta / o$contacts
  med_puta_percontact <- median(e_puta_percontact, na.rm = TRUE)
  q_puta_percontact <- quantile(e_puta_percontact, probs = c(0.1, 0.9), na.rm = TRUE, names = FALSE)
  sort_pia <- sort(o$pia / o$contacts)
    
    list(
      o = o,
      propintervened = NA,
      n_units = nrow(o),
  puta = c(sum(o$puta),
       sum(o$puta) / sum(o$contacts),
       med_puta_percontact,
       q_puta_percontact[1],
       q_puta_percontact[2]),
      pia = c(mean(o$pia), 
              sort_pia[max(1, ceiling(0.01 * nrow(o)))], 
              sort_pia[min(nrow(o), floor(0.99 * nrow(o)))]),
      total_contacts = sum(o$contacts)
    )
  }
  
  
  
  