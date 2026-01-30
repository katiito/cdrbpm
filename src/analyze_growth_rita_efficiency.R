# Analyze relationship between RITA cases in growth clusters and IDA per contact
# This script examines whether:
# 1. Number of RITA+ cases in a growth cluster affects IDA/contact efficiency
# 2. Timing of first RITA+ case relative to growth trigger affects efficiency

library(dplyr)
library(ggplot2)
library(tidyr)

# Source the intervention code to get RITA definitions
source('src/intervention0.R')

# Load the simulation data
cat('Loading simulation data...\n')
Dall <- read.csv('src/experiment1-N10000-gens7-D.csv')
Gall <- read.csv('src/experiment1-N10000-gens7-G.csv')

# Load the most recent growth cluster results
growth_files <- list.files('intervention-results/',
                           pattern = 'details_growth_.*\\.csv',
                           full.names = TRUE)
latest_growth_file <- sort(growth_files, decreasing = TRUE)[1]
cat('Loading growth cluster results from:', latest_growth_file, '\n')
growth_results <- read.csv(latest_growth_file)

# Parameters (matching intervention0.R defaults)
growth_distance_threshold <- 0.01
cluster_size <- 5
lookback_window_months <- 6
lookback_days <- lookback_window_months * 30
rita_window_months <- 6
analysis_delay_days <- 14
implementation_delay_days <- 14

cat('Parameters:\n')
cat('  Growth distance threshold:', growth_distance_threshold, '\n')
cat('  Cluster size:', cluster_size, '\n')
cat('  Lookback window:', lookback_window_months, 'months\n')
cat('  RITA window:', rita_window_months, 'months\n')
cat('  Total delays:', analysis_delay_days + implementation_delay_days, 'days\n\n')

# Set seed for reproducibility
set.seed(123)

# Simulate RITA status for all individuals
cat('Simulating RITA status...\n')
Gall$rita <- with(Gall, (timediagnosed - timeinfected) < rexp(nrow(Gall), 1/(rita_window_months*30)))

# Split data by simulation
simids <- as.character(unique(Dall$simid))
Ds <- split(Dall, Dall$simid)[simids]
Gs <- split(Gall, Gall$simid)[simids]

# Function to find growth clusters and analyze RITA within them
analyze_growth_clusters <- function() {
  results_list <- list()

  for (sim_key in simids) {
    D <- Ds[[sim_key]]
    G <- Gs[[sim_key]]

    # Order by sequencing time
    G <- G[order(G$timesequenced), ]
    lastgen <- max(G$generation)

    # Build cluster based on genetic distance
    D1 <- D[D$distance <= growth_distance_threshold, ]
    keeppids <- '0'
    addpids <- D1$recipient[D1$donor %in% keeppids]
    while (length(addpids) > 0) {
      keeppids <- union(addpids, keeppids)
      addpids <- setdiff(D1$recipient[D1$donor %in% keeppids], keeppids)
    }
    Gcluster <- G[G$pid %in% keeppids, ]

    # Get trigger candidates (exclude seed and last gen)
    Gtrig <- Gcluster[Gcluster$generation > 0 & Gcluster$generation != lastgen, ]
    if (nrow(Gtrig) < cluster_size) next

    # Find trigger time using sliding window
    t <- as.numeric(Gtrig$timesequenced)
    n <- length(t)
    j <- 1
    t_detect <- Inf
    trigger_idx <- NA

    for (i in seq_len(n)) {
      while (j <= i && (t[i] - t[j]) > lookback_days) {
        j <- j + 1
      }
      if ((i - j + 1) >= cluster_size) {
        t_detect <- t[i]
        trigger_idx <- i
        break
      }
    }

    if (!is.finite(t_detect)) next

    # Calculate intervention time (with delays)
    analysis_delay <- analysis_delay_days
    implementation_delay <- implementation_delay_days
    t_trigger <- t_detect + analysis_delay
    IT <- t_trigger + implementation_delay

    # Get the growth cluster result for this simulation
    growth_row <- growth_results[growth_results$simid == sim_key, ]
    if (nrow(growth_row) == 0) next

    # Calculate IDA per contact for both network assumptions
    ida <- growth_row$puta  # IDA (was called PUTA)
    ida_per_contact_small <- ida / growth_row$contacts_small
    ida_per_contact_large <- ida / growth_row$contacts_large

    # Identify cluster members who were sequenced before IT
    # (These are the ones we "know about" at intervention time)
    Gcluster_known <- Gcluster[Gcluster$generation != lastgen &
                                Gcluster$timesequenced < IT, ]

    # Also get the full cluster (all generations except last, all time)
    Gcluster_full <- Gcluster[Gcluster$generation != lastgen, ]

    # Count RITA+ cases in known cluster (at IT)
    rita_in_known <- Gcluster_known[Gcluster_known$rita, ]
    n_rita_known <- nrow(rita_in_known)
    prop_rita_known <- n_rita_known / nrow(Gcluster_known)

    # Count RITA+ cases in full cluster (eventual)
    rita_in_full <- Gcluster_full[Gcluster_full$rita, ]
    n_rita_full <- nrow(rita_in_full)
    prop_rita_full <- n_rita_full / nrow(Gcluster_full)

    # Find timing of first RITA+ case
    if (n_rita_full > 0) {
      # Earliest RITA diagnosis in the cluster
      first_rita_time <- min(rita_in_full$timediagnosed)
      first_rita_detection <- first_rita_time - analysis_delay_days  # When we'd know about it
      time_to_growth_trigger <- t_detect - first_rita_detection  # Days RITA was earlier

      # Find earliest RITA case known at IT
      rita_known_at_IT <- rita_in_full[rita_in_full$timediagnosed < (IT - analysis_delay_days), ]
      if (nrow(rita_known_at_IT) > 0) {
        earliest_rita_known_at_IT <- min(rita_known_at_IT$timediagnosed)
        time_advantage_at_IT <- IT - earliest_rita_known_at_IT
      } else {
        earliest_rita_known_at_IT <- NA
        time_advantage_at_IT <- NA
      }
    } else {
      first_rita_time <- NA
      first_rita_detection <- NA
      time_to_growth_trigger <- NA
      earliest_rita_known_at_IT <- NA
      time_advantage_at_IT <- NA
    }

    # Store results
    results_list[[length(results_list) + 1]] <- data.frame(
      simid = sim_key,
      # Growth cluster info
      IT = IT,
      t_detect = t_detect,
      t_trigger = t_trigger,
      nc_known = nrow(Gcluster_known),  # Cluster size at IT
      nc_full = nrow(Gcluster_full),     # Full cluster size
      # IDA efficiency metrics
      ida = ida,
      pia = growth_row$pia,
      ida_per_contact_small = ida_per_contact_small,
      ida_per_contact_large = ida_per_contact_large,
      contacts_small = growth_row$contacts_small,
      contacts_large = growth_row$contacts_large,
      # RITA metrics at IT (what we know)
      n_rita_known = n_rita_known,
      prop_rita_known = prop_rita_known,
      # RITA metrics eventual (full cluster)
      n_rita_full = n_rita_full,
      prop_rita_full = prop_rita_full,
      # Timing metrics
      first_rita_time = first_rita_time,
      first_rita_detection = first_rita_detection,
      time_to_growth_trigger = time_to_growth_trigger,  # Days RITA faster than growth
      earliest_rita_known_at_IT = earliest_rita_known_at_IT,
      time_advantage_at_IT = time_advantage_at_IT,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results_list)
}

# Run the analysis
cat('Analyzing growth clusters for RITA content and timing...\n')
results <- analyze_growth_clusters()

# Save results
output_file <- 'intervention-results/growth_rita_efficiency_analysis.csv'
write.csv(results, output_file, row.names = FALSE)
cat('Results saved to:', output_file, '\n\n')

# Print summary statistics
cat('=== SUMMARY STATISTICS ===\n\n')
cat('Total growth clusters analyzed:', nrow(results), '\n\n')

cat('RITA presence in growth clusters (at intervention time):\n')
cat('  Clusters with ≥1 RITA+ case:', sum(results$n_rita_known > 0), '/', nrow(results),
    sprintf('(%.1f%%)\n', 100 * mean(results$n_rita_known > 0)))
cat('  Mean RITA+ cases per cluster:', round(mean(results$n_rita_known), 2), '\n')
cat('  Range of RITA+ cases:', min(results$n_rita_known), '-', max(results$n_rita_known), '\n')
cat('  Mean proportion RITA+:', sprintf('%.1f%%\n', 100 * mean(results$prop_rita_known)), '\n\n')

cat('RITA presence in growth clusters (eventual/full cluster):\n')
cat('  Clusters with ≥1 RITA+ case:', sum(results$n_rita_full > 0), '/', nrow(results),
    sprintf('(%.1f%%)\n', 100 * mean(results$n_rita_full > 0)))
cat('  Mean RITA+ cases per cluster:', round(mean(results$n_rita_full), 2), '\n')
cat('  Range of RITA+ cases:', min(results$n_rita_full), '-', max(results$n_rita_full), '\n')
cat('  Mean proportion RITA+:', sprintf('%.1f%%\n\n', 100 * mean(results$prop_rita_full)))

cat('Timing advantage of RITA:\n')
valid_timing <- results[!is.na(results$time_to_growth_trigger), ]
cat('  Clusters with RITA timing data:', nrow(valid_timing), '\n')
cat('  Mean days RITA faster than growth:',
    round(mean(valid_timing$time_to_growth_trigger), 1), '\n')
cat('  Range:', round(min(valid_timing$time_to_growth_trigger), 1), 'to',
    round(max(valid_timing$time_to_growth_trigger), 1), 'days\n')
cat('  Clusters where RITA is faster:',
    sum(valid_timing$time_to_growth_trigger > 0), '/', nrow(valid_timing), '\n\n')

cat('IDA per contact efficiency (small network assumption):\n')
cat('  Mean:', round(mean(results$ida_per_contact_small, na.rm = TRUE), 2), 'days/contact\n')
cat('  Median:', round(median(results$ida_per_contact_small, na.rm = TRUE), 2), 'days/contact\n')
cat('  Range:', round(min(results$ida_per_contact_small, na.rm = TRUE), 2), 'to',
    round(max(results$ida_per_contact_small, na.rm = TRUE), 2), 'days/contact\n\n')

cat('IDA per contact efficiency (large network assumption):\n')
cat('  Mean:', round(mean(results$ida_per_contact_large, na.rm = TRUE), 2), 'days/contact\n')
cat('  Median:', round(median(results$ida_per_contact_large, na.rm = TRUE), 2), 'days/contact\n')
cat('  Range:', round(min(results$ida_per_contact_large, na.rm = TRUE), 2), 'to',
    round(max(results$ida_per_contact_large, na.rm = TRUE), 2), 'days/contact\n\n')

# Correlation analysis
cat('=== CORRELATION ANALYSIS ===\n\n')

# Number of RITA cases vs IDA per contact
cor_n_rita_small <- cor(results$n_rita_known, results$ida_per_contact_small,
                        use = "complete.obs")
cor_n_rita_large <- cor(results$n_rita_known, results$ida_per_contact_large,
                        use = "complete.obs")

cat('Number of RITA+ cases (at IT) vs IDA/contact:\n')
cat('  Small network: r =', round(cor_n_rita_small, 3), '\n')
cat('  Large network: r =', round(cor_n_rita_large, 3), '\n\n')

# Proportion RITA vs IDA per contact
cor_prop_small <- cor(results$prop_rita_known, results$ida_per_contact_small,
                      use = "complete.obs")
cor_prop_large <- cor(results$prop_rita_known, results$ida_per_contact_large,
                      use = "complete.obs")

cat('Proportion RITA+ (at IT) vs IDA/contact:\n')
cat('  Small network: r =', round(cor_prop_small, 3), '\n')
cat('  Large network: r =', round(cor_prop_large, 3), '\n\n')

# Timing advantage vs IDA per contact
valid_timing_subset <- results[!is.na(results$time_to_growth_trigger), ]
if (nrow(valid_timing_subset) > 2) {
  cor_timing_small <- cor(valid_timing_subset$time_to_growth_trigger,
                          valid_timing_subset$ida_per_contact_small,
                          use = "complete.obs")
  cor_timing_large <- cor(valid_timing_subset$time_to_growth_trigger,
                          valid_timing_subset$ida_per_contact_large,
                          use = "complete.obs")

  cat('RITA timing advantage (days faster) vs IDA/contact:\n')
  cat('  Small network: r =', round(cor_timing_small, 3), '\n')
  cat('  Large network: r =', round(cor_timing_large, 3), '\n\n')
}

# Create visualizations
cat('Creating visualizations...\n')

# Plot 1: Number of RITA cases vs IDA per contact
p1 <- ggplot(results, aes(x = n_rita_known, y = ida_per_contact_small)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "IDA per Contact vs Number of RITA+ Cases in Growth Cluster",
    subtitle = "Small network assumption (at intervention time)",
    x = "Number of RITA+ cases in cluster (at IT)",
    y = "IDA per contact (days/contact)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave('intervention-results/rita_count_vs_efficiency.png', p1, width = 8, height = 6)

# Plot 2: Proportion RITA vs IDA per contact
p2 <- ggplot(results, aes(x = prop_rita_known, y = ida_per_contact_small)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "IDA per Contact vs Proportion of RITA+ Cases in Growth Cluster",
    subtitle = "Small network assumption (at intervention time)",
    x = "Proportion of cluster that is RITA+ (at IT)",
    y = "IDA per contact (days/contact)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave('intervention-results/rita_proportion_vs_efficiency.png', p2, width = 8, height = 6)

# Plot 3: RITA timing advantage vs IDA per contact
if (nrow(valid_timing_subset) > 2) {
  p3 <- ggplot(valid_timing_subset, aes(x = time_to_growth_trigger, y = ida_per_contact_small)) +
    geom_point(alpha = 0.6, size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
    labs(
      title = "IDA per Contact vs RITA Timing Advantage",
      subtitle = "Small network assumption (positive = RITA detected earlier than growth)",
      x = "Days RITA detected earlier than growth trigger",
      y = "IDA per contact (days/contact)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

  ggsave('intervention-results/rita_timing_vs_efficiency.png', p3, width = 8, height = 6)
}

# Plot 4: Comparison across network assumptions
results_long <- results %>%
  select(simid, n_rita_known, ida_per_contact_small, ida_per_contact_large) %>%
  pivot_longer(cols = c(ida_per_contact_small, ida_per_contact_large),
               names_to = "network_assumption",
               values_to = "ida_per_contact") %>%
  mutate(network_assumption = ifelse(network_assumption == "ida_per_contact_small",
                                     "Small network", "Large network"))

p4 <- ggplot(results_long, aes(x = n_rita_known, y = ida_per_contact, color = network_assumption)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "IDA per Contact vs RITA+ Cases: Network Assumption Comparison",
    x = "Number of RITA+ cases in cluster (at IT)",
    y = "IDA per contact (days/contact)",
    color = "Network\nAssumption"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave('intervention-results/network_assumption_comparison.png', p4, width = 10, height = 6)

cat('\nPlots saved to intervention-results/\n')
cat('  - rita_count_vs_efficiency.png\n')
cat('  - rita_proportion_vs_efficiency.png\n')
cat('  - rita_timing_vs_efficiency.png\n')
cat('  - network_assumption_comparison.png\n\n')

cat('Analysis complete!\n')
