# =============================================================================
# plot_interventions.R
# =============================================================================
# Plotting functions for intervention analysis results
#
# Usage:
#   source("src/intervention0.R")
#   source("src/plot_interventions.R")
#   results <- run_interventions(Dall, Gall)
#   plot_efficiency_distributions(results)
# =============================================================================


# =============================================================================
# plot_efficiency_distributions
# =============================================================================
#' Plot distributions of PUTA and PIA efficiency across intervention strategies
#' 
#' Creates 4 plots showing PUTA and PIA efficiency distributions under small 
#' and large subnetwork assumptions for all 6 intervention strategies.
#' Uses violin plots with pseudo-log scale to handle outliers.
#' Strategies are displayed as rows (horizontal orientation).
#' 
#' @param results The output from run_interventions()
#' @param title_prefix Optional prefix for plot titles
#' @return A combined ggplot object with 4 panels
#' 
plot_efficiency_distributions <- function(results, 
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
      
      # For individual-based interventions (network, random, rita), there's only one contacts column
      if (sname %in% c("network", "random", "rita")) {
        if (!"contacts" %in% names(o)) return(NULL)
        data.frame(
          strategy = slabel,
          puta = o$puta,
          pia = o$pia,
          contacts_small = o$contacts,  # Use same contacts for both assumptions
          contacts_large = o$contacts
        )
      } else {
        # For cluster-based strategies (distsize5, distsize2, growth)
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
  
  # Set factor levels for ordering (reversed so top-to-bottom matches desired order)
  df$strategy <- factor(df$strategy, levels = rev(strategy_labels))
  
  # Define color palette
  strategy_colors <- c(
    "Size>=5" = "#E41A1C",
    "Size>=2" = "#377EB8", 
    "Growth" = "#4DAF4A",
    "Random" = "#984EA3",
    "RITA" = "#FF7F00",
    "Network" = "#A65628"
  )
  
  # Pseudo-log transformation using asinh for smoother behavior near zero
  # asinh(x) ≈ x for small x, and ≈ sign(x)*log(2|x|) for large |x|
  # Scale factor adjusts the transition point - higher = more expansion near zero
  
  # PUTA transformation (lower scale to show more detail at higher values)
  puta_scale <- 0.5
  puta_trans <- scales::trans_new(
    name = "puta_pseudo_log",
    transform = function(x) asinh(x * puta_scale) / puta_scale,
    inverse = function(x) sinh(x * puta_scale) / puta_scale
  )
  
  # PIA transformation (higher scale to expand 0-1 range where most density is)
  pia_scale <- 10  # Higher value = more space allocated to values near 0
  pia_trans <- scales::trans_new(
    name = "pia_pseudo_log",
    transform = function(x) asinh(x * pia_scale) / pia_scale,
    inverse = function(x) sinh(x * pia_scale) / pia_scale
  )
  
  # Common theme for horizontal violin plots
  common_theme <- theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      plot.title = element_text(size = 11, face = "bold")
    )
  
  # PUTA plots with pseudo-log scale (horizontal orientation)
  # Extended tick marks to show full range
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
                 color = "white", size = 3, shape = 21, fill = "black", stroke = 1.5) +
    scale_y_continuous(
      trans = puta_trans,
      breaks = c(0, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000),
      labels = c("0", "0.5", "1", "2", "5", "10", "20", "50", "100", "200", "500", "1k", "2k", "5k")
    ) +
    coord_flip() +
    labs(y = "PUTA per contact", x = "",
         title = paste0(title_prefix, "PUTA (Small Subnetwork)")) +
    scale_fill_manual(values = strategy_colors) +
    common_theme
  
  p2 <- ggplot(df, aes(x = strategy, y = puta_eff_large, fill = strategy)) +
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
                 color = "white", size = 3, shape = 21, fill = "black", stroke = 1.5) +
    scale_y_continuous(
      trans = puta_trans,
      breaks = c(0, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000),
      labels = c("0", "0.5", "1", "2", "5", "10", "20", "50", "100", "200", "500", "1k", "2k", "5k")
    ) +
    coord_flip() +
    labs(y = "PUTA per contact", x = "",
         title = paste0(title_prefix, "PUTA (Large Subnetwork)")) +
    scale_fill_manual(values = strategy_colors) +
    common_theme
  
  # PIA plots with pseudo-log scale - truncated to show bulk density
  p3 <- ggplot(df, aes(x = strategy, y = pia_eff_small, fill = strategy)) +
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
                 color = "white", size = 3, shape = 21, fill = "black", stroke = 1.5) +
    scale_y_continuous(
      trans = pia_trans,
      breaks = c(0, 0.01, 0.02, 0.05, 0.1, 0.2),
      labels = c("0", "0.01", "0.02", "0.05", "0.1", "0.2"),
      limits = c(0, 0.2)
    ) +
    coord_flip() +
    labs(y = "PIA per contact", x = "",
         title = paste0(title_prefix, "PIA (Small Subnetwork)")) +
    scale_fill_manual(values = strategy_colors) +
    common_theme
  
  p4 <- ggplot(df, aes(x = strategy, y = pia_eff_large, fill = strategy)) +
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
                 color = "white", size = 3, shape = 21, fill = "black", stroke = 1.5) +
    scale_y_continuous(
      trans = pia_trans,
      breaks = c(0, 0.01, 0.02, 0.05, 0.1, 0.2),
      labels = c("0", "0.01", "0.02", "0.05", "0.1", "0.2"),
      limits = c(0, 0.2)
    ) +
    coord_flip() +
    labs(y = "PIA per contact", x = "",
         title = paste0(title_prefix, "PIA (Large Subnetwork)")) +
    scale_fill_manual(values = strategy_colors) +
    common_theme
  
  # Combine plots: 2x2 grid
  (p1 + p2) / (p3 + p4) + 
    plot_annotation(title = "Intervention Efficiency Distributions")
}


# =============================================================================
# plot_mechanism_analysis
# =============================================================================
#' Plot mechanistic analysis of intervention strategies
#' 
#' Creates a multi-panel figure explaining WHY different strategies perform
#' as they do:
#' (i) Why RITA captures donors so effectively (timing)
#' (ii) How growth clusters identify higher transmitters
#' (iii) How delays are distributed by component
#' 
#' @param D Distance data frame (from experiment CSV)
#' @param G Generation data frame (from experiment CSV)
#' @param distance_threshold_distsize Distance threshold for size clusters (default 0.005, matching intervention0.R)
#' @param distance_threshold_growth Distance threshold for growth clusters (default 0.01, matching intervention0.R)
#' @param n_sims Number of simulations to analyze (default 100)
#' @return A combined ggplot object with mechanism analysis panels
#' 
plot_mechanism_analysis <- function(D, G, 
                                    distance_threshold_distsize = 0.005,
                                    distance_threshold_growth = 0.01,
                                    n_sims = 100) {
  require(ggplot2)
  require(patchwork)
  require(dplyr)
  
  # Color palette
  strategy_colors <- c(
    "Size>=5" = "#E41A1C",
    "Size>=2" = "#377EB8", 
    "Growth" = "#4DAF4A",
    "Random" = "#984EA3",
    "RITA" = "#FF7F00",
    "Network" = "#A65628"
  )
  
  # Consistent ordering for all plots
  strategy_order <- c("RITA", "Growth", "Size>=5", "Random", "Network")
  
  # Parameters - matching those used in intervention0.R analysis
  # (distance thresholds come from function parameters)
  implementation_delay <- 14        # implementation_delay_days
  analysis_delay <- 14              # analysis_delay_days  
  lookback_days <- 6 * 30           # lookback_window_months = 6
  rita_window_days <- 6 * 30        # rita_window_months = 6 (average window for RITA)
  network_degree_threshold <- 4     # network_degree_threshold
  
  # Use all simulations (no filtering/sorting to avoid bias)
  simids <- unique(G$simid)
  if (!is.null(n_sims) && n_sims < length(simids)) {
    simids <- simids[1:n_sims]
  }
  
  # Storage for results
  donor_analysis <- data.frame()
  transmission_analysis <- data.frame()
  delay_results <- data.frame()
  timing_analysis <- data.frame()
  survivorship_results <- data.frame()  # For survivorship bias analysis
  
  cat(sprintf("Analyzing %d simulations...\n", length(simids)))
  
  progress_interval <- max(1, floor(length(simids) / 10))  # Report every 10%
  
  for (sim_idx in seq_along(simids)) {
    simid <- simids[sim_idx]
    
    # Progress update
    if (sim_idx %% progress_interval == 0 || sim_idx == length(simids)) {
      cat(sprintf("  Progress: %d/%d (%.0f%%)\n", sim_idx, length(simids), 100 * sim_idx / length(simids)))
    }
    
    G1 <- G[G$simid == simid, ]
    D1 <- D[D$simid == simid, ]
    lastgen <- max(G1$generation)
    
    # Skip generation 0 and last generation for analysis
    G1_eligible <- G1[G1$generation > 0 & G1$generation < lastgen, ]
    if (nrow(G1_eligible) == 0) next
    
    # =====================================================================
    # RITA strategy: recently infected individuals
    # Matching intervention0.R: RITA positive if (timediagnosed - timeinfected) < rexp(1, 1/rita_window_days)
    # =====================================================================
    G1_eligible_sorted <- G1_eligible[order(G1_eligible$timediagnosed), ]
    
    for (k in 1:nrow(G1_eligible_sorted)) {
      person <- G1_eligible_sorted[k, ]
      pid <- person$pid
      time_dx <- person$timediagnosed
      time_infected <- person$timeinfected
      IT <- time_dx + implementation_delay
      
      # RITA: simulate exponential test window (matching intervention0.R)
      rita_threshold <- rexp(1, 1/rita_window_days)
      if ((time_dx - time_infected) <= rita_threshold) {
        # Use D1 for donor/transmission info
        # Exclude transmissions to recipients in the final generation
        transmissions_from_pid <- D1[D1$donor == pid, ]
        # Filter out recipients in final generation
        recipient_gens <- sapply(transmissions_from_pid$recipient, function(r) {
          g <- G1$generation[G1$pid == r]
          if (length(g) == 0) NA else g[1]
        })
        transmissions_from_pid <- transmissions_from_pid[!is.na(recipient_gens) & recipient_gens < lastgen, ]
        
        is_donor <- nrow(transmissions_from_pid) > 0
        total_trans <- nrow(transmissions_from_pid)
        pia <- sum(transmissions_from_pid$timetransmission > IT)
        
        # Find this person's donor and check if donor is undiagnosed at IT
        donor_row <- D1[D1$recipient == pid, ]
        donor_undiagnosed <- FALSE
        if (nrow(donor_row) > 0) {
          donor_pid <- donor_row$donor[1]
          donor_info <- G1[G1$pid == donor_pid, ]
          if (nrow(donor_info) > 0) {
            donor_undiagnosed <- donor_info$timediagnosed[1] > IT
          }
        }
        
        donor_analysis <- rbind(donor_analysis, data.frame(
          strategy = "RITA", is_donor = is_donor, pia = pia,
          donor_undiagnosed = donor_undiagnosed
        ))
        
        transmission_analysis <- rbind(transmission_analysis, data.frame(
          strategy = "RITA", n_transmitted = total_trans
        ))
        
        if (total_trans > 0) {
          frac_remaining <- pia / total_trans
          timing_analysis <- rbind(timing_analysis, data.frame(
            strategy = "RITA", 
            time_since_infection = IT - time_infected,
            frac_remaining = frac_remaining
          ))
        }
      }
    }
    
    # =====================================================================
    # Random strategy: ALL diagnosed individuals (population baseline)
    # =====================================================================
    if (nrow(G1_eligible_sorted) >= 1) {
      # Use ALL eligible individuals to represent population baseline
      random_sample <- G1_eligible_sorted
      for (k in 1:nrow(random_sample)) {
        person <- random_sample[k, ]
        pid <- person$pid
        time_dx <- person$timediagnosed
        time_infected <- person$timeinfected
        IT <- time_dx + implementation_delay
        
        # Use D1 for donor/transmission info
        # Exclude transmissions to recipients in the final generation
        transmissions_from_pid <- D1[D1$donor == pid, ]
        recipient_gens <- sapply(transmissions_from_pid$recipient, function(r) {
          g <- G1$generation[G1$pid == r]
          if (length(g) == 0) NA else g[1]
        })
        transmissions_from_pid <- transmissions_from_pid[!is.na(recipient_gens) & recipient_gens < lastgen, ]
        
        is_donor <- nrow(transmissions_from_pid) > 0
        total_trans <- nrow(transmissions_from_pid)
        pia <- sum(transmissions_from_pid$timetransmission > IT)
        
        # Find this person's donor and check if donor is undiagnosed at IT
        donor_row <- D1[D1$recipient == pid, ]
        donor_undiagnosed <- FALSE
        if (nrow(donor_row) > 0) {
          donor_pid <- donor_row$donor[1]
          donor_info <- G1[G1$pid == donor_pid, ]
          if (nrow(donor_info) > 0) {
            donor_undiagnosed <- donor_info$timediagnosed[1] > IT
          }
        }
        
        donor_analysis <- rbind(donor_analysis, data.frame(
          strategy = "Random", is_donor = is_donor, pia = pia,
          donor_undiagnosed = donor_undiagnosed
        ))
        
        transmission_analysis <- rbind(transmission_analysis, data.frame(
          strategy = "Random", n_transmitted = total_trans
        ))
        
        if (total_trans > 0) {
          frac_remaining <- pia / total_trans
          timing_analysis <- rbind(timing_analysis, data.frame(
            strategy = "Random",
            time_since_infection = IT - time_infected,
            frac_remaining = frac_remaining
          ))
        }
      }
    }
    
    # =====================================================================
    # Network (high degree) strategy - top decile of contacts (>=12)
    # =====================================================================
    G1_eligible_sorted$total_contacts <- G1_eligible_sorted$Fcontacts_90d + 
                                          G1_eligible_sorted$Gcontacts_90d + 
                                          G1_eligible_sorted$Hcontacts_90d
    high_degree <- G1_eligible_sorted[G1_eligible_sorted$total_contacts >= 8, ]
    
    if (nrow(high_degree) > 0) {
      for (k in 1:nrow(high_degree)) {
        person <- high_degree[k, ]
        pid <- person$pid
        time_dx <- person$timediagnosed
        time_infected <- person$timeinfected
        IT <- time_dx + implementation_delay
        
        # Use D1 for donor/transmission info
        # Exclude transmissions to recipients in the final generation
        transmissions_from_pid <- D1[D1$donor == pid, ]
        recipient_gens <- sapply(transmissions_from_pid$recipient, function(r) {
          g <- G1$generation[G1$pid == r]
          if (length(g) == 0) NA else g[1]
        })
        transmissions_from_pid <- transmissions_from_pid[!is.na(recipient_gens) & recipient_gens < lastgen, ]
        
        is_donor <- nrow(transmissions_from_pid) > 0
        total_trans <- nrow(transmissions_from_pid)
        pia <- sum(transmissions_from_pid$timetransmission > IT)
        
        # Find this person's donor and check if donor is undiagnosed at IT
        donor_row <- D1[D1$recipient == pid, ]
        donor_undiagnosed <- FALSE
        if (nrow(donor_row) > 0) {
          donor_pid <- donor_row$donor[1]
          donor_info <- G1[G1$pid == donor_pid, ]
          if (nrow(donor_info) > 0) {
            donor_undiagnosed <- donor_info$timediagnosed[1] > IT
          }
        }
        
        donor_analysis <- rbind(donor_analysis, data.frame(
          strategy = "Network", is_donor = is_donor, pia = pia,
          donor_undiagnosed = donor_undiagnosed
        ))
        
        transmission_analysis <- rbind(transmission_analysis, data.frame(
          strategy = "Network", n_transmitted = total_trans
        ))
        
        if (total_trans > 0) {
          frac_remaining <- pia / total_trans
          timing_analysis <- rbind(timing_analysis, data.frame(
            strategy = "Network",
            time_since_infection = IT - time_infected,
            frac_remaining = frac_remaining
          ))
        }
      }
    }
    
    # =====================================================================
    # Growth cluster strategy
    # =====================================================================
    D1_filtered <- D1[D1$distance <= distance_threshold_growth, ]
    keeppids <- "0"
    addpids <- D1_filtered$recipient[D1_filtered$donor %in% keeppids]
    while (length(addpids) > 0) {
      keeppids <- union(addpids, keeppids)
      addpids <- setdiff(D1_filtered$recipient[D1_filtered$donor %in% keeppids], keeppids)
    }
    
    Gcluster <- G1[G1$pid %in% keeppids, ]
    Gtrig <- Gcluster[Gcluster$generation > 0 & Gcluster$generation != lastgen, ]
    Gtrig <- Gtrig[order(Gtrig$timesequenced), ]
    
    # Find growth trigger
    t <- Gtrig$timesequenced
    n <- length(t)
    if (n >= 5) {
      j <- 1
      trigger_i <- NA
      window_start_j <- NA
      
      for (ii in 1:n) {
        while (j <= ii && (t[ii] - t[j]) > lookback_days) {
          j <- j + 1
        }
        if ((ii - j + 1) >= 5) {
          trigger_i <- ii
          window_start_j <- j
          break
        }
      }
      
      if (!is.na(trigger_i)) {
        window_start_time <- t[window_start_j]
        IT <- t[trigger_i] + analysis_delay + implementation_delay
        knowledge_cutoff <- IT - analysis_delay
        
        # Get cluster members (matching code logic)
        cluster_at_IT <- Gcluster[Gcluster$generation != lastgen & 
                                  Gcluster$timesequenced >= window_start_time &
                                  Gcluster$timesequenced < knowledge_cutoff, ]
        
        for (k in 1:nrow(cluster_at_IT)) {
          person <- cluster_at_IT[k, ]
          pid <- person$pid
          time_dx <- person$timediagnosed
          time_seq <- person$timesequenced
          time_infected <- person$timeinfected
          
          # Use D1 for transmission info
          # Exclude transmissions to recipients in the final generation
          transmissions_from_pid <- D1[D1$donor == pid, ]
          recipient_gens <- sapply(transmissions_from_pid$recipient, function(r) {
            g <- G1$generation[G1$pid == r]
            if (length(g) == 0) NA else g[1]
          })
          transmissions_from_pid <- transmissions_from_pid[!is.na(recipient_gens) & recipient_gens < lastgen, ]
          
          total_trans <- nrow(transmissions_from_pid)
          pia <- sum(transmissions_from_pid$timetransmission > IT)
          
          # Check if this person is a donor (transmits to anyone not in final gen)
          is_donor <- total_trans > 0
          
          # Check if this person's donor is still undiagnosed at IT
          donor_undiagnosed <- FALSE
          donor_row <- D1[D1$recipient == pid, ]
          if (nrow(donor_row) > 0) {
            donor_pid <- donor_row$donor[1]
            donor_info <- G1[G1$pid == donor_pid, ]
            if (nrow(donor_info) > 0) {
              donor_undiagnosed <- donor_info$timediagnosed[1] > IT
            }
          }
          
          donor_analysis <- rbind(donor_analysis, data.frame(
            strategy = "Growth", is_donor = is_donor, pia = pia,
            donor_undiagnosed = donor_undiagnosed
          ))
          
          transmission_analysis <- rbind(transmission_analysis, data.frame(
            strategy = "Growth", n_transmitted = total_trans
          ))
          
          if (total_trans > 0) {
            frac_remaining <- pia / total_trans
            timing_analysis <- rbind(timing_analysis, data.frame(
              strategy = "Growth",
              time_since_infection = IT - time_infected,
              frac_remaining = frac_remaining
            ))
          }
          
          # Delay components for growth cluster intervention
          # Timeline: Infection → Diagnosis → Sequencing → [Wait for cluster] → Trigger → Analysis → Implementation → IT
          delay_results <- rbind(delay_results, data.frame(
            component = c("Dx to Sequencing\n(sequencing delay)", 
                          "Sequencing to Trigger\n(cluster accumulation)", 
                          "Trigger to Analysis\n(analysis delay)",
                          "Analysis to Intervention\n(implementation delay)"),
            delay = c(
              time_seq - time_dx,                    # Time from diagnosis to sequencing complete
              t[trigger_i] - time_seq,              # Time waiting for 5th case in cluster (varies by member)
              analysis_delay,                        # Fixed 14-day analysis delay after trigger
              implementation_delay                   # Fixed 14-day implementation delay after analysis
            )
          ))
        }
        
        # =====================================================================
        # Survivorship bias analysis: offspring of growth cluster members
        # =====================================================================
        for (k in 1:nrow(cluster_at_IT)) {
          member_pid <- cluster_at_IT$pid[k]
          offspring_pids <- D1$recipient[D1$donor == member_pid]
          for (offspring_pid in offspring_pids) {
            offspring_info <- G1[G1$pid == offspring_pid, ]
            if (nrow(offspring_info) == 0 || offspring_info$generation >= lastgen) next
            
            # Count transmissions from this offspring
            trans_from_offspring <- D1[D1$donor == offspring_pid, ]
            recipient_gens <- sapply(trans_from_offspring$recipient, function(r) {
              g <- G1$generation[G1$pid == r]
              if (length(g) == 0) NA else g[1]
            })
            trans_from_offspring <- trans_from_offspring[!is.na(recipient_gens) & recipient_gens < lastgen, ]
            
            survivorship_results <- rbind(survivorship_results, data.frame(
              group = "Growth cluster offspring",
              generation = offspring_info$generation[1],
              n_transmitted = nrow(trans_from_offspring)
            ))
          }
        }
      }
    }
    
    # =====================================================================
    # Population baseline for survivorship analysis (all eligible individuals)
    # =====================================================================
    for (i in 1:nrow(G1_eligible)) {
      survivorship_results <- rbind(survivorship_results, data.frame(
        group = "Population",
        generation = G1_eligible$generation[i],
        n_transmitted = {
          trans <- D1[D1$donor == G1_eligible$pid[i], ]
          if (nrow(trans) == 0) {
            0
          } else {
            recipient_gens <- sapply(trans$recipient, function(r) {
              g <- G1$generation[G1$pid == r]
              if (length(g) == 0) NA else g[1]
            })
            sum(!is.na(recipient_gens) & recipient_gens < lastgen)
          }
        }
      ))
    }
    
    # =====================================================================
    # Size>=5 cluster strategy
    # =====================================================================
    D1_filtered_size <- D1[D1$distance <= distance_threshold_distsize, ]
    keeppids_size <- "0"
    addpids_size <- D1_filtered_size$recipient[D1_filtered_size$donor %in% keeppids_size]
    while (length(addpids_size) > 0) {
      keeppids_size <- union(addpids_size, keeppids_size)
      addpids_size <- setdiff(D1_filtered_size$recipient[D1_filtered_size$donor %in% keeppids_size], keeppids_size)
    }
    
    Gcluster_size <- G1[G1$pid %in% keeppids_size, ]
    Gcluster_size <- Gcluster_size[Gcluster_size$generation > 0 & Gcluster_size$generation != lastgen, ]
    Gcluster_size <- Gcluster_size[order(Gcluster_size$timesequenced), ]
    
    # Find when cluster reaches size 5
    if (nrow(Gcluster_size) >= 5) {
      time_trigger <- Gcluster_size$timesequenced[5]
      IT <- time_trigger + analysis_delay + implementation_delay
      knowledge_cutoff <- IT - analysis_delay
      
      cluster_at_IT <- Gcluster_size[Gcluster_size$timesequenced < knowledge_cutoff, ]
      
      for (k in 1:nrow(cluster_at_IT)) {
        person <- cluster_at_IT[k, ]
        pid <- person$pid
        
        # Use D1 for transmission info
        # Exclude transmissions to recipients in the final generation
        transmissions_from_pid <- D1[D1$donor == pid, ]
        recipient_gens <- sapply(transmissions_from_pid$recipient, function(r) {
          g <- G1$generation[G1$pid == r]
          if (length(g) == 0) NA else g[1]
        })
        transmissions_from_pid <- transmissions_from_pid[!is.na(recipient_gens) & recipient_gens < lastgen, ]
        
        total_trans <- nrow(transmissions_from_pid)
        pia <- sum(transmissions_from_pid$timetransmission > IT)
        
        # Check if this person is a donor (transmits to anyone not in final gen)
        is_donor <- total_trans > 0
        
        # Check if this person's donor is still undiagnosed at IT
        donor_undiagnosed <- FALSE
        donor_row <- D1[D1$recipient == pid, ]
        if (nrow(donor_row) > 0) {
          donor_pid <- donor_row$donor[1]
          donor_info <- G1[G1$pid == donor_pid, ]
          if (nrow(donor_info) > 0) {
            donor_undiagnosed <- donor_info$timediagnosed[1] > IT
          }
        }
        
        donor_analysis <- rbind(donor_analysis, data.frame(
          strategy = "Size>=5", is_donor = is_donor, pia = pia,
          donor_undiagnosed = donor_undiagnosed
        ))
        
        transmission_analysis <- rbind(transmission_analysis, data.frame(
          strategy = "Size>=5", n_transmitted = total_trans
        ))
        
        if (total_trans > 0) {
          frac_remaining <- pia / total_trans
          timing_analysis <- rbind(timing_analysis, data.frame(
            strategy = "Size>=5",
            time_since_infection = IT - person$timeinfected,
            frac_remaining = frac_remaining
          ))
        }
      }
    }
  }
  
  cat("Creating plots...\n")
  
  # =========================================================================
  # Panel A: RITA's undiagnosed donor advantage
  # Shows what % of targets have an undiagnosed donor at intervention time
  # =========================================================================
  
  if (nrow(donor_analysis) > 0 && "donor_undiagnosed" %in% names(donor_analysis)) {
    # Calculate percentage of targets whose donor is undiagnosed by strategy
    donor_pct <- donor_analysis %>%
      group_by(strategy) %>%
      summarise(
        n_total = n(),
        n_undiag_donor = sum(donor_undiagnosed, na.rm = TRUE),
        pct_undiag_donor = mean(donor_undiagnosed, na.rm = TRUE) * 100,
        .groups = "drop"
      )
    
    donor_pct$strategy <- factor(donor_pct$strategy, levels = strategy_order)
    
    p_donor <- ggplot(donor_pct, aes(x = strategy, y = pct_undiag_donor, fill = strategy)) +
      geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
      geom_text(aes(label = sprintf("%.0f%%", pct_undiag_donor)), vjust = -0.5, size = 3.5) +
      scale_fill_manual(values = strategy_colors) +
      labs(
        x = "", 
        y = "% with undiagnosed donor",
        title = "A. Targeting active transmission chains",
        subtitle = "% of targets whose donor is still undiagnosed"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        axis.text.x = element_text(angle = 15, hjust = 1)
      ) +
      coord_cartesian(ylim = c(0, 100))
  } else {
    p_donor <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No donor data available") +
      theme_void() + labs(title = "A. Targeting active transmission chains")
  }
  
  # =========================================================================
  # Panel B: Timing advantage (RITA intervenes early)
  # =========================================================================
  
  if (nrow(timing_analysis) > 0) {
    timing_summary <- timing_analysis %>%
      group_by(strategy) %>%
      summarise(
        mean_time_since_inf = mean(time_since_infection),
        mean_frac_remaining = mean(frac_remaining) * 100,
        se_frac = sd(frac_remaining) / sqrt(n()) * 100,
        n = n(),
        .groups = "drop"
      )
    
    timing_summary$strategy <- factor(timing_summary$strategy, 
                                      levels = strategy_order)
    
    p_timing <- ggplot(timing_summary, aes(x = strategy, y = mean_frac_remaining, fill = strategy)) +
      geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
      geom_errorbar(aes(ymin = pmax(0, mean_frac_remaining - se_frac), 
                        ymax = pmin(100, mean_frac_remaining + se_frac)),
                    width = 0.2, color = "gray30") +
      geom_text(aes(label = sprintf("%.0f%%", mean_frac_remaining)), vjust = -0.5, size = 3.5) +
      scale_fill_manual(values = strategy_colors) +
      labs(
        x = "", 
        y = "% of transmissions preventable",
        title = "B. RITA intervenes early in transmission",
        subtitle = "Among transmitters: % of transmissions averted by intervention"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        axis.text.x = element_text(angle = 15, hjust = 1)
      ) +
      coord_cartesian(ylim = c(0, 100))
  } else {
    p_timing <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No timing data available") +
      theme_void() + labs(title = "B. RITA intervenes early in transmission career")
  }
  
  # =========================================================================
  # Panel C: Growth clusters identify higher transmitters
  # =========================================================================
  
  if (nrow(transmission_analysis) > 0) {
    trans_summary <- transmission_analysis %>%
      group_by(strategy) %>%
      summarise(
        mean_trans = mean(n_transmitted),
        se_trans = sd(n_transmitted) / sqrt(n()),
        pct_transmitters = mean(n_transmitted > 0) * 100,
        n = n(),
        .groups = "drop"
      )
    
    trans_summary$strategy <- factor(trans_summary$strategy, 
                                     levels = strategy_order)
    
    p_trans <- ggplot(trans_summary, aes(x = strategy, y = mean_trans, fill = strategy)) +
      geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
      geom_errorbar(aes(ymin = mean_trans - se_trans, ymax = mean_trans + se_trans),
                    width = 0.2, color = "gray30") +
      geom_text(aes(label = sprintf("%.2f", mean_trans)), vjust = -0.5, size = 3.5) +
      scale_fill_manual(values = strategy_colors) +
      labs(
        x = "", 
        y = "Mean transmissions per target",
        title = "C. Growth clusters identify higher transmitters",
        subtitle = "Total offspring per intervention target"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        axis.text.x = element_text(angle = 15, hjust = 1)
      ) +
      coord_cartesian(ylim = c(0, max(trans_summary$mean_trans, na.rm = TRUE) * 1.3))
  } else {
    p_trans <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No transmission data available") +
      theme_void() + labs(title = "C. Growth clusters identify higher transmitters")
  }
  
  # =========================================================================
  # Panel D: Survivorship bias - enrichment fades across generations
  # =========================================================================
  
  if (nrow(survivorship_results) > 0 && "Growth cluster offspring" %in% survivorship_results$group) {
    surv_summary <- survivorship_results %>%
      filter(generation >= 1 & generation <= 4) %>%
      group_by(group, generation) %>%
      summarise(
        mean_trans = mean(n_transmitted),
        se_trans = sd(n_transmitted) / sqrt(n()),
        n = n(),
        .groups = "drop"
      )
    
    # Calculate enrichment ratio - rename columns for easier access
    ratio_data <- surv_summary %>%
      select(group, generation, mean_trans) %>%
      tidyr::pivot_wider(names_from = group, values_from = mean_trans)
    
    # Handle column names with spaces
    offspring_col <- names(ratio_data)[grepl("Growth", names(ratio_data))]
    if (length(offspring_col) > 0) {
      ratio_data$offspring_trans <- ratio_data[[offspring_col]]
      ratio_data$ratio <- ratio_data$offspring_trans / ratio_data$Population
    } else {
      ratio_data$offspring_trans <- NA
      ratio_data$ratio <- NA
    }
    
    surv_summary$group <- factor(surv_summary$group, 
                                  levels = c("Population", "Growth cluster offspring"))
    
    p_surv <- ggplot(surv_summary, aes(x = factor(generation), y = mean_trans, fill = group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
      geom_errorbar(aes(ymin = mean_trans - se_trans, ymax = mean_trans + se_trans),
                    position = position_dodge(width = 0.8), width = 0.2, color = "gray30") +
      scale_fill_manual(values = c("Population" = "gray70", "Growth cluster offspring" = "gray30")) +
      # Add ratio annotations
      geom_text(data = ratio_data[!is.na(ratio_data$ratio), ], 
                aes(x = factor(generation), y = offspring_trans + 0.15, 
                    label = sprintf("%.1fx", ratio)),
                inherit.aes = FALSE, size = 3, color = "#2E7D32", fontface = "bold") +
      labs(
        x = "Generation", 
        y = "Mean transmissions",
        fill = NULL,
        title = "D. Survivorship bias in growth clusters",
        subtitle = "Offspring enrichment fades across generations"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40")
      )
  } else {
    p_surv <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No survivorship data available") +
      theme_void() + labs(title = "D. Survivorship bias in growth clusters")
  }
  
  # =========================================================================
  # Panel E: Delay components for growth clusters
  # =========================================================================
  
  if (nrow(delay_results) > 0) {
    delay_summary <- delay_results %>%
      group_by(component) %>%
      summarise(
        mean_delay = mean(delay),
        se_delay = sd(delay) / sqrt(n()),
        .groups = "drop"
      )
    
    delay_summary$component <- factor(delay_summary$component, 
                                      levels = c("Dx to Sequencing\n(sequencing delay)", 
                                                 "Sequencing to Trigger\n(cluster accumulation)", 
                                                 "Trigger to Analysis\n(analysis delay)",
                                                 "Analysis to Intervention\n(implementation delay)"))
  
    p_delay <- ggplot(delay_summary, aes(x = component, y = mean_delay)) +
      geom_bar(stat = "identity", fill = "gray50", color = "black", width = 0.7) +
      geom_errorbar(aes(ymin = mean_delay - se_delay, ymax = mean_delay + se_delay),
                    width = 0.2, color = "black") +
      geom_text(aes(label = sprintf("%.0f days", mean_delay)), vjust = -0.5, size = 3.5) +
      labs(
        x = "", 
        y = "Days",
        title = "E. Growth cluster delay components",
        subtitle = "Timeline: Diagnosis → Sequencing → Cluster trigger → Analysis → Intervention"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        axis.text.x = element_text(size = 8, angle = 15, hjust = 1)
      )
  } else {
    p_delay <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No delay data available") +
      theme_void() + labs(title = "E. Growth cluster delay components")
  }
  
  # =========================================================================
  # Combine all panels (left to right, top to bottom: A B / C D / E)
  # =========================================================================
  
  combined <- (p_donor + p_timing) / (p_trans + p_surv) / (p_delay + plot_spacer()) +
    plot_annotation(
      title = "Mechanism Analysis: Why Strategies Differ in Efficiency",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
  
  return(combined)
}


# =============================================================================
# run_mechanism_analysis
# =============================================================================
#' Convenience function to run the full mechanism analysis
#' 
#' @param D_path Path to distance CSV file
#' @param G_path Path to generation CSV file  
#' @param save_path Optional path to save the figure
#' @param n_sims Number of simulations to analyze (default 100)
#' @param width Figure width in inches (default 12)
#' @param height Figure height in inches (default 10)
#' @return The ggplot object
#' 
run_mechanism_analysis <- function(D_path = "src/experiment1-N10000-gens7-D.csv",
                                   G_path = "src/experiment1-N10000-gens7-G.csv",
                                   save_path = "intervention-plots/mechanism_analysis.png",
                                   n_sims = 100,
                                   width = 12,
                                   height = 12) {
  
  cat("Loading data...\n")
  D <- read.csv(D_path)
  G <- read.csv(G_path)
  
  cat("Generating mechanism analysis plots...\n")
  p <- plot_mechanism_analysis(D, G, n_sims = n_sims)
  
  if (!is.null(save_path)) {
    # Ensure directory exists
    save_dir <- dirname(save_path)
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    cat(sprintf("Saving figure to %s...\n", save_path))
    ggsave(save_path, p, width = width, height = height, dpi = 300)
  }
  
  return(invisible(p))
}
