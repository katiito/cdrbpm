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
  d_file = 'experiment1-N100000-gens7-D.csv',
  g_file = 'experiment1-N100000-gens7-G.csv',
  seed = NULL,  # Default to random seed
  cluster_size_5 = 5,
  cluster_size_2 = 2,
  subnetwork = "small", # large or small network (sparse or dense connections inside cluster)
  distance_threshold = 0.005,
  network_degree_threshold = 4,
  random_sample_size = 30,
  rita_window_months = 6, # average RITA detection window
  lookback_window_months = 6, # growth-rate trigger window
  growth_distance_threshold = 0.1, # separate D for growth-based trigger
  intervention_rate = 1/90, # average 3 months to intervention
  show_table = TRUE,
  show_debug = FALSE
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
        subnetwork = subnetwork,
        intervention_rate = intervention_rate
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta = c(0, 0, 0, 0, 0), pia = c(0, 0, 0), total_contacts = 0)
      }
    )
    cat(" done (", ods5$n_units, " units)\n", sep = "")

    cat("  [2/6] Distance-size (size=", cluster_size_2, ", D=", distance_threshold, ")...", sep = "")
    ods2 <- tryCatch(
      distsize_intervention(
        Ds = Ds, Gs = Gs, Gall = Gall,
        distance_threshold = distance_threshold,
        cluster_size = cluster_size_2,
        subnetwork = subnetwork,
        intervention_rate = intervention_rate
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta = c(0, 0, 0, 0, 0), pia = c(0, 0, 0), total_contacts = 0)
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
        subnetwork = subnetwork,
        intervention_rate = intervention_rate
      ),
      error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(propintervened = 0, n_units = 0,
             puta = c(0, 0, 0, 0, 0), pia = c(0, 0, 0), total_contacts = 0)
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

    # Compile summary table
    odf <- rbind(
      with(ods5, c(total_contacts, puta, pia)),
      with(ods2, c(total_contacts, puta, pia)),
      with(ogrowth, c(total_contacts, puta, pia)),
      with(orand, c(total_contacts, puta, pia)),
      with(orita, c(total_contacts, puta, pia)),
      with(onet,  c(total_contacts, puta, pia))
    ) |> as.data.frame()
    colnames(odf) <- c("Contacted Total", "Total PUTA", "PUTA/contacted", "Median", "Low", "High", "PIA", "Low", "High")
    rownames(odf) <- c(
      paste0('Size=', cluster_size_5, ',D=', distance_threshold),
      paste0('Size=', cluster_size_2, ',D=', distance_threshold),
      paste0('Growth, size=', cluster_size_5, ',W=', lookback_window_months, 'mo,D=', growth_distance_threshold),
      'Random allocation',
      'RITA',
      paste0('Network, partners>', network_degree_threshold)
    )
    odf1 <- round(odf, 2)

    counts_df <- data.frame(
      Strategy = c(
        paste0('Size=', cluster_size_5, ',D=', distance_threshold),
        paste0('Size=', cluster_size_2, ',D=', distance_threshold),
        paste0('Growth, size=', cluster_size_5, ',W=', lookback_window_months, 'mo,D=', growth_distance_threshold),
        'Random allocation',
        'RITA',
        paste0('Network, partners>', network_degree_threshold)
      ),
      Units = c(ods5$n_units, ods2$n_units, ogrowth$n_units, orand$n_units, orita$n_units, onet$n_units),
      TotalContacts = c(ods5$total_contacts, ods2$total_contacts, ogrowth$total_contacts, orand$total_contacts, orita$total_contacts, onet$total_contacts)
    )

    # Display results
    if (show_table) {
      if (require(knitr, quietly = TRUE)) {
        print(knitr::kable(odf1))
        cat("\nSample sizes and totals (for context):\n")
        print(knitr::kable(counts_df))
        if (show_debug) {
          cat("\nDebug: Cluster-size-5 per-simulation details (odf1) -> simid, nc, total_contacts, interventiontime, pia, puta\n")
          if (!is.null(ods5$o) && nrow(ods5$o) > 0) {
            debug_cols <- c("simid", "nc", "total_contacts", "interventiontime", "pia", "puta")
            odf1_5_view <- ods5$o[, debug_cols]
            num_cols_all <- intersect(c("nc", "total_contacts", "interventiontime", "pia", "puta"), colnames(odf1_5_view))
            odf1_5_view[num_cols_all] <- lapply(odf1_5_view[num_cols_all], function(x) suppressWarnings(as.numeric(x)))
            num_cols <- intersect(c("interventiontime", "pia", "puta", "total_contacts"), colnames(odf1_5_view))
            odf1_5_view[num_cols] <- lapply(odf1_5_view[num_cols], function(x) if (is.numeric(x)) round(x, 2) else x)
            print(knitr::kable(odf1_5_view))
            cat("\nDebug: Degree breakdown at IT for size-5 clusters (verifies why total_contacts equals its value)\n")
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
              excess <- pmax(G1$degree - (n - 1), 0)
              if (identical(subnetwork, "small")) {
                tc <- n + sum(excess)
              } else {
                total_degree <- sum(G1$degree)
                tc <- total_degree - (n - 2)
              }
              cat(sprintf("simid=%s n=%d nc_recorded=%s total_contacts(recalc)=%.2f recorded=%.2f degrees=%s excess=%s\n",
                          sim_key, n, as.character(row$nc), tc, as.numeric(row$total_contacts),
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
          cat("\nDebug: Cluster-size-5 per-simulation details (odf1) -> simid, nc, total_contacts, interventiontime, pia, puta\n")
          if (!is.null(ods5$o) && nrow(ods5$o) > 0) {
            debug_cols <- c("simid", "nc", "total_contacts", "interventiontime", "pia", "puta")
            odf1_5_view <- ods5$o[, debug_cols]
            num_cols_all <- intersect(c("nc", "total_contacts", "interventiontime", "pia", "puta"), colnames(odf1_5_view))
            odf1_5_view[num_cols_all] <- lapply(odf1_5_view[num_cols_all], function(x) suppressWarnings(as.numeric(x)))
            num_cols <- intersect(c("interventiontime", "pia", "puta", "total_contacts"), colnames(odf1_5_view))
            odf1_5_view[num_cols] <- lapply(odf1_5_view[num_cols], function(x) if (is.numeric(x)) round(x, 2) else x)
            print(odf1_5_view)
            cat("\nDebug: Degree breakdown at IT for size-5 clusters (verifies why total_contacts equals its value)\n")
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
              excess <- pmax(G1$degree - (n - 1), 0)
              if (identical(subnetwork, "small")) {
                tc <- n + sum(excess)
              } else {
                total_degree <- sum(G1$degree)
                tc <- total_degree - (n - 2)
              }
              cat(sprintf("simid=%s n=%d nc_recorded=%s total_contacts(recalc)=%.2f recorded=%.2f degrees=%s excess=%s\n",
                          sim_key, n, as.character(row$nc), tc, as.numeric(row$total_contacts),
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
    subnetwork = "small",
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
  
  # Calculate total cluster contact network size
  if (nrow(G1) > 0) {
    # Get cluster members' total degrees
    G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
    total_degree <- sum(G1$degree)
    n <- nrow(G1)
    
    # Estimate internal connections within cluster
    if (subnetwork == "large"){ # sparse connections inside cluster
      total_contacts <- total_degree - (n - 2)
    } else if (subnetwork == "small"){ # dense connections inside cluster
      excess <- pmax(G1$degree - (n - 1), 0)
      total_contacts <- n + sum(excess)
      
    } else {
      stop("Assumption about how structure of contacts in cluster needs to be defined (large or small subnetwork)")
    }
    
  } else {
    total_contacts <- 0
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
  c(pia = pia, puta = puta, interventiontime = IT, nc = nrow(G1), total_contacts = total_contacts)
}

  
### ### ### # Distance-size intervention strategy

distsize_intervention <- function(
    Ds, Gs, Gall,
    distance_threshold, cluster_size, subnetwork = "small", intervention_rate
)
{
  
  # Process all simulations
  o <- lapply(seq_along(Ds), function(i) {
    tryCatch({
      proc_cluster(
        D = Ds[[i]], G = Gs[[i]],
        distance_threshold = distance_threshold,
        cluster_size = cluster_size,
        subnetwork = subnetwork,
        intervention_rate = intervention_rate
      )
    }, error = function(e) {
      # Return default values if processing fails
      c(pia = 0, puta = 0, interventiontime = Inf, nc = 0, total_contacts = 0)
    })
  })
  
  odf  <- do.call( rbind, o ) |> as.data.frame()
  # Attach simulation identifiers to each row for debugging/traceability
  simids_vec <- names(Ds)
  odf$simid <- if (!is.null(simids_vec)) simids_vec else seq_len(nrow(odf))
  # Reorder columns to bring simid upfront
  odf <- odf[, c("simid", "pia", "puta", "interventiontime", "nc", "total_contacts")]

  # Sanity-correct nc and total_contacts by recomputing deterministically at the recorded IT
  # This avoids any accidental scoping bugs and ensures recorded values match the recomputation used in debug
  if (nrow(odf) > 0) {
    # Coerce potential character columns
    suppressWarnings({
      odf$interventiontime <- as.numeric(odf$interventiontime)
      odf$nc <- as.numeric(odf$nc)
      odf$total_contacts <- as.numeric(odf$total_contacts)
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
      if (identical(subnetwork, "small")) {
        excess <- pmax(G1$degree - (n - 1), 0)
        tc <- n + sum(excess)
      } else if (identical(subnetwork, "large")) {
        total_degree <- sum(G1$degree)
        tc <- total_degree - (n - 2)
      } else {
        tc <- NA_real_
      }
      odf$nc[i] <- n
      odf$total_contacts[i] <- tc
    }
  }
  odf1 <- odf[ !is.infinite( odf$interventiontime ), ] # exclude sims where no cluster found 
  
  lastgen <- max( Gall$generation )
  
  # Do NOT drop clusters with total_contacts==0; include them for counting
  # Guard downstream per-contact computations instead of filtering them out
  
  if (nrow(odf1) == 0) {
    return(list(
      o = data.frame(pia = numeric(0), puta = numeric(0), interventiontime = numeric(0), 
                     nc = numeric(0), total_contacts = numeric(0)),
      propintervened = 0,
      n_units = 0,
      puta = c(0, 0, 0, 0, 0),
      pia = c(0, 0, 0),
      total_contacts = 0
    ))
  }
  
  e_puta_percontact <- odf1$puta / odf1$total_contacts
  e_puta_percontact_valid <- e_puta_percontact[is.finite(e_puta_percontact) & !is.na(e_puta_percontact)]
  med_puta_percontact <- if (length(e_puta_percontact_valid) > 0) median(e_puta_percontact_valid) else NA_real_
  q_puta_percontact <- if (length(e_puta_percontact_valid) > 0) quantile(e_puta_percontact_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
  sort_pia <- sort(odf1$pia)
  
  # output IQ ranges and pia/puta as per contacted individuals in addition to totals
  list(
    o = odf1,
    propintervened = sum(odf1$nc) / sum(Gall$generation > 0 & Gall$generation < lastgen),
    n_units = nrow(odf1),
  puta = c(sum(odf1$puta),
       if (sum(odf1$total_contacts) > 0) sum(odf1$puta) / sum(odf1$total_contacts) else NA_real_,
       med_puta_percontact,
       q_puta_percontact[1],
       q_puta_percontact[2]),
    pia = c(mean(odf1$pia), 
            sort_pia[max(1, ceiling(0.1 * nrow(odf1)))], 
            sort_pia[min(nrow(odf1), floor(0.9 * nrow(odf1)))]),
    total_contacts = sum(odf1$total_contacts)
  )
}

### ### ### # Growth-rate cluster intervention

growthrate_intervention <- function(
  Ds, Gs, Gall,
  growth_distance_threshold, cluster_size,
  lookback_window_months = 3,
  subnetwork = "small",
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
      return(c(pia = 0, puta = 0, interventiontime = Inf, nc = 0, total_contacts = 0))
    }

    # Sliding window over infection times: earliest time t where count in [t-lookback_days, t] >= cluster_size
    t <- suppressWarnings(as.numeric(Gtrig$timeinfected))
    if (all(!is.finite(t))) {
      # No usable infection times; cannot trigger growth
      return(c(pia = 0, puta = 0, interventiontime = Inf, nc = 0, total_contacts = 0))
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

    # Compute contacts within cluster
    if (nrow(G1) > 0 && is.finite(IT)) {
      G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
      n1 <- nrow(G1)
      if (subnetwork == "large") {
        total_degree <- sum(G1$degree)
        total_contacts <- total_degree - (n1 - 2)
      } else if (subnetwork == "small") {
        excess <- pmax(G1$degree - (n1 - 1), 0)
        total_contacts <- n1 + sum(excess)
      } else {
        stop("subnetwork must be 'small' or 'large'")
      }
    } else {
      total_contacts <- 0
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

    c(pia = pia, puta = puta, interventiontime = IT, nc = nrow(G1), total_contacts = total_contacts)
  }

  # Process all simulations
  o <- lapply(seq_along(Ds), function(i) {
    tryCatch({
      process_one(Ds[[i]], Gs[[i]])
    }, error = function(e) {
      c(pia = 0, puta = 0, interventiontime = Inf, nc = 0, total_contacts = 0)
    })
  })

  odf <- do.call(rbind, o) |> as.data.frame()
  simids_vec <- names(Ds)
  odf$simid <- if (!is.null(simids_vec)) simids_vec else seq_len(nrow(odf))
  odf <- odf[, c("simid", "pia", "puta", "interventiontime", "nc", "total_contacts")]

  # Coerce numeric columns
  suppressWarnings({
    odf$interventiontime <- as.numeric(odf$interventiontime)
    odf$nc <- as.numeric(odf$nc)
    odf$total_contacts <- as.numeric(odf$total_contacts)
  })

  odf1 <- odf[!is.infinite(odf$interventiontime), ]
  lastgen <- max(Gall$generation)

  if (nrow(odf1) == 0) {
    return(list(
      o = data.frame(pia = numeric(0), puta = numeric(0), interventiontime = numeric(0),
                     nc = numeric(0), total_contacts = numeric(0)),
      propintervened = 0,
      n_units = 0,
      puta = c(0, 0, 0, 0, 0),
      pia = c(0, 0, 0),
      total_contacts = 0
    ))
  }

  e_puta_percontact <- odf1$puta / odf1$total_contacts
  e_puta_percontact_valid <- e_puta_percontact[is.finite(e_puta_percontact) & !is.na(e_puta_percontact)]
  med_puta_percontact <- if (length(e_puta_percontact_valid) > 0) median(e_puta_percontact_valid) else NA_real_
  q_puta_percontact <- if (length(e_puta_percontact_valid) > 0) quantile(e_puta_percontact_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)
  sort_pia <- sort(odf1$pia)

  list(
    o = odf1,
    propintervened = sum(odf1$nc) / sum(Gall$generation > 0 & Gall$generation < lastgen),
    n_units = nrow(odf1),
    puta = c(sum(odf1$puta),
             if (sum(odf1$total_contacts) > 0) sum(odf1$puta) / sum(odf1$total_contacts) else NA_real_,
             med_puta_percontact,
             q_puta_percontact[1],
             q_puta_percontact[2]),
    pia = c(mean(odf1$pia),
            sort_pia[max(1, ceiling(0.1 * nrow(odf1)))],
            sort_pia[min(nrow(odf1), floor(0.9 * nrow(odf1)))]),
    total_contacts = sum(odf1$total_contacts)
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
  
  
  
  