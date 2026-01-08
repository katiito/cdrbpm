# Test and debug utilities for intervention analysis
# Separate from main intervention0.R for cleaner code organization

source("src/intervention0.R")

# =============================================================================
# Test: Verify simid alignment between D and G data
# =============================================================================
test_simid_alignment <- function(d_file = "src/experiment1-N10000-gens7-D.csv",
                                  g_file = "src/experiment1-N10000-gens7-G.csv") {
  Dall <- read.csv(d_file)
  Gall <- read.csv(g_file)
  
  cat("=== Simid Alignment Test ===\n")
  cat("Unique simids in D:", length(unique(Dall$simid)), "\n")
  cat("Unique simids in G:", length(unique(Gall$simid)), "\n")
  
  # Correct alignment using character keys
  simids <- as.character(unique(Dall$simid))
  Ds <- split(Dall, Dall$simid)[simids]
  Gs <- split(Gall, Gall$simid)[simids]
  
  cat("After alignment:\n")
  cat("  Length of Ds:", length(Ds), "\n")
  cat("  Length of Gs:", length(Gs), "\n")
  cat("  Names match:", identical(names(Ds), names(Gs)), "\n")
  
  # Verify a few random entries
  test_idx <- sample(seq_along(Ds), min(5, length(Ds)))
  all_match <- TRUE
  for (i in test_idx) {
    d_simid <- unique(Ds[[i]]$simid)
    g_simid <- unique(Gs[[i]]$simid)
    if (!identical(d_simid, g_simid)) {
      cat("  MISMATCH at position", i, ": D has", d_simid, ", G has", g_simid, "\n")
      all_match <- FALSE
    }
  }
  if (all_match) cat("  Random sample verification: PASSED\n")
  
  invisible(list(Ds = Ds, Gs = Gs, aligned = identical(names(Ds), names(Gs))))
}


# =============================================================================
# Test: Verify cluster size >= threshold after filtering
# =============================================================================
test_cluster_size_constraint <- function(d_file = "src/experiment1-N10000-gens7-D.csv",
                                          g_file = "src/experiment1-N10000-gens7-G.csv",
                                          cluster_size = 5,
                                          distance_threshold = 0.005,
                                          seed = 123) {
  set.seed(seed)
  
  Dall <- read.csv(d_file)
  Gall <- read.csv(g_file)
  simids <- as.character(unique(Dall$simid))
  Ds <- split(Dall, Dall$simid)[simids]
  Gs <- split(Gall, Gall$simid)[simids]
  
  od <- distsize_intervention(Ds, Gs, Gall, 
                               distance_threshold = distance_threshold, 
                               cluster_size = cluster_size, 
                               intervention_rate = 1/90)
  
  cat("=== Cluster Size Constraint Test ===\n")
  cat("cluster_size threshold:", cluster_size, "\n")
  cat("n_units:", od$n_units, "\n")
  
  if (od$n_units > 0) {
    below_threshold <- sum(od$o$nc < cluster_size)
    cat("Clusters with nc < threshold:", below_threshold, "\n")
    
    if (below_threshold > 0) {
      cat("FAIL: Found clusters below threshold\n")
      cat("  Problem simids:", paste(od$o$simid[od$o$nc < cluster_size], collapse = ", "), "\n")
      return(FALSE)
    } else {
      cat("PASS: All clusters have nc >= threshold\n")
      return(TRUE)
    }
  } else {
    cat("No clusters found (cannot verify constraint)\n")
    return(NA)
  }
}


# =============================================================================
# Debug: Trace through proc_cluster for a specific simulation
# =============================================================================
debug_proc_cluster <- function(simid, 
                                d_file = "src/experiment1-N10000-gens7-D.csv",
                                g_file = "src/experiment1-N10000-gens7-G.csv",
                                distance_threshold = 0.005,
                                cluster_size = 5) {
  Dall <- read.csv(d_file)
  Gall <- read.csv(g_file)
  simids <- as.character(unique(Dall$simid))
  Ds <- split(Dall, Dall$simid)[simids]
  Gs <- split(Gall, Gall$simid)[simids]
  
  simid_char <- as.character(simid)
  if (!(simid_char %in% names(Ds))) {
    cat("Simid", simid, "not found in D data\n")
    return(NULL)
  }
  
  D <- Ds[[simid_char]]
  G <- Gs[[simid_char]]
  
  cat("=== Debug proc_cluster for simid", simid, "===\n")
  
  # Sort by sequencing time
  G <- G[order(G$timesequenced), ]
  D <- D[order(D$timetransmission), ]
  
  lastgen <- max(G$generation)
  cat("Last generation:", lastgen, "\n")
  cat("Total rows in G:", nrow(G), "\n")
  cat("Total rows in D:", nrow(D), "\n")
  
  # Distance filter
  D1 <- D[D$distance <= distance_threshold, ]
  cat("\nEdges within distance", distance_threshold, ":", nrow(D1), "\n")
  
  # Build cluster from seed
  keeppids <- "0"
  addpids <- D1$recipient[D1$donor %in% keeppids]
  iterations <- 0
  while (length(addpids) > 0) {
    iterations <- iterations + 1
    keeppids <- union(addpids, keeppids)
    addpids <- setdiff(D1$recipient[D1$donor %in% keeppids], keeppids)
  }
  cat("PIDs reachable from 0:", length(keeppids), "(", iterations, "iterations)\n")
  cat("  PIDs:", paste(head(keeppids, 10), collapse = ", "), 
      if(length(keeppids) > 10) "..." else "", "\n")
  
  # Filter G to cluster members
  G1_all <- G[G$pid %in% keeppids, ]
  cat("\nCluster members (all gens):", nrow(G1_all), "\n")
  
  G1 <- G1_all[G1_all$generation != lastgen, ]
  cat("Cluster members (excl last gen):", nrow(G1), "\n")
  
  # Trigger check
  cat("\nTrigger check:\n")
  cat("  cluster_size threshold:", cluster_size, "\n")
  cat("  nrow(G1) >= cluster_size:", nrow(G1) >= cluster_size, "\n")
  
  if (nrow(G1) >= cluster_size) {
    cat("  timesequenced[", cluster_size, "]:", G1$timesequenced[cluster_size], "\n")
    set.seed(999)  # Fixed seed for reproducibility
    IT <- G1$timesequenced[cluster_size] + rexp(1, rate = 1/90)
    cat("  IT (with delay):", IT, "\n")
    
    G1_final <- G1[G1$timesequenced < IT, ]
    cat("  Members after IT filter:", nrow(G1_final), "\n")
  } else {
    cat("  Cluster too small, IT = Inf\n")
  }
  
  # Run actual proc_cluster
  cat("\n--- Actual proc_cluster result ---\n")
  set.seed(999)
  result <- proc_cluster(D, G, distance_threshold, cluster_size, intervention_rate = 1/90)
  print(result)
  
  invisible(result)
}


# =============================================================================
# Debug: Show degree breakdown for triggered clusters
# =============================================================================
debug_degree_breakdown <- function(d_file = "src/experiment1-N10000-gens7-D.csv",
                                    g_file = "src/experiment1-N10000-gens7-G.csv",
                                    distance_threshold = 0.005,
                                    cluster_size = 5,
                                    seed = 123,
                                    max_clusters = 10) {
  set.seed(seed)
  
  Dall <- read.csv(d_file)
  Gall <- read.csv(g_file)
  simids <- as.character(unique(Dall$simid))
  Ds <- split(Dall, Dall$simid)[simids]
  Gs <- split(Gall, Gall$simid)[simids]
  
  od <- distsize_intervention(Ds, Gs, Gall, 
                               distance_threshold = distance_threshold, 
                               cluster_size = cluster_size, 
                               intervention_rate = 1/90)
  
  cat("=== Degree Breakdown for Size-", cluster_size, " Clusters ===\n", sep = "")
  cat("Total clusters:", od$n_units, "\n\n")
  
  if (od$n_units == 0) {
    cat("No clusters to analyze\n")
    return(invisible(NULL))
  }
  
  n_show <- min(max_clusters, od$n_units)
  
  for (i in seq_len(n_show)) {
    row <- od$o[i, ]
    sim_key <- as.character(row$simid)
    IT <- as.numeric(row$interventiontime)
    
    D <- Ds[[sim_key]]
    G <- Gs[[sim_key]]
    
    lastgen <- max(G$generation)
    D1 <- D[D$distance <= distance_threshold, ]
    keeppids <- "0"
    addpids <- D1$recipient[D1$donor %in% keeppids]
    while (length(addpids) > 0) {
      keeppids <- union(addpids, keeppids)
      addpids <- setdiff(D1$recipient[D1$donor %in% keeppids], keeppids)
    }
    
    G1 <- G[G$pid %in% keeppids, ]
    G1 <- G1[G1$generation != lastgen & G1$timesequenced < IT, ]
    G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
    
    n <- nrow(G1)
    sum_deg <- sum(G1$degree)
    excess <- pmax(G1$degree - (n - 1), 0)
    tc_small <- n + sum(excess)
    tc_large <- sum_deg - (n - 2)
    
    cat(sprintf("simid=%s n=%d contacts_small=%.2f contacts_large=%.2f\n",
                sim_key, n, tc_small, tc_large))
    cat(sprintf("  degrees: %s\n", paste(round(G1$degree, 2), collapse = ", ")))
    cat(sprintf("  excess:  %s\n\n", paste(round(excess, 2), collapse = ", ")))
  }
  
  if (od$n_units > max_clusters) {
    cat("... (", od$n_units - max_clusters, " more clusters not shown)\n", sep = "")
  }
  
  invisible(od)
}


# =============================================================================
# Run all tests
# =============================================================================
run_all_tests <- function(d_file = "src/experiment1-N10000-gens7-D.csv",
                           g_file = "src/experiment1-N10000-gens7-G.csv",
                           seed = 123) {
  cat("========================================\n")
  cat("Running intervention analysis tests\n")
  cat("========================================\n\n")
  
  # Test 1: Simid alignment
  test_simid_alignment(d_file, g_file)
  cat("\n")
  
  # Test 2: Cluster size constraint for size=5

  test_cluster_size_constraint(d_file, g_file, cluster_size = 5, seed = seed)
  cat("\n")
  
  # Test 3: Cluster size constraint for size=2
  test_cluster_size_constraint(d_file, g_file, cluster_size = 2, seed = seed)
  cat("\n")
  
  cat("========================================\n")
  cat("All tests completed\n")
  cat("========================================\n")
}


# =============================================================================
# Example usage (uncomment to run)
# =============================================================================
# run_all_tests()
# debug_proc_cluster(simid = 1608)
# debug_degree_breakdown(max_clusters = 5)
