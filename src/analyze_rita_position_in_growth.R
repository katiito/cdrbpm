# Analyze RITA position within growth cluster trigger window
# Does having RITA+ cases in positions 1-5 of the triggering window affect IDA/contact?

library(dplyr)
library(ggplot2)
library(tidyr)

# Source the intervention code
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

# Parameters
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
cat('  RITA window:', rita_window_months, 'months\n\n')

# Set seed for reproducibility
set.seed(123)

# Simulate RITA status
cat('Simulating RITA status...\n')
Gall$rita <- with(Gall, (timediagnosed - timeinfected) < rexp(nrow(Gall), 1/(rita_window_months*30)))

# Split data by simulation
simids <- as.character(unique(Dall$simid))
Ds <- split(Dall, Dall$simid)[simids]
Gs <- split(Gall, Gall$simid)[simids]

# Function to analyze RITA position in trigger window
analyze_rita_positions <- function() {
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

    # Get trigger candidates
    Gtrig <- Gcluster[Gcluster$generation > 0 & Gcluster$generation != lastgen, ]
    if (nrow(Gtrig) < cluster_size) next

    # Find trigger time using sliding window
    t <- as.numeric(Gtrig$timesequenced)
    n <- length(t)
    j <- 1
    t_detect <- Inf
    trigger_start_idx <- NA
    trigger_end_idx <- NA

    for (i in seq_len(n)) {
      while (j <= i && (t[i] - t[j]) > lookback_days) {
        j <- j + 1
      }
      if ((i - j + 1) >= cluster_size) {
        t_detect <- t[i]
        trigger_start_idx <- j
        trigger_end_idx <- i
        break
      }
    }

    if (!is.finite(t_detect)) next

    # Calculate intervention time
    analysis_delay <- analysis_delay_days
    implementation_delay <- implementation_delay_days
    t_trigger <- t_detect + analysis_delay
    IT <- t_trigger + implementation_delay

    # Get the growth cluster result
    growth_row <- growth_results[growth_results$simid == sim_key, ]
    if (nrow(growth_row) == 0) next

    # Calculate IDA per contact
    ida <- growth_row$puta
    ida_per_contact_small <- ida / growth_row$contacts_small
    ida_per_contact_large <- ida / growth_row$contacts_large

    # Get the 5 cases in the trigger window
    trigger_window_cases <- Gtrig[trigger_start_idx:trigger_end_idx, ]

    # Verify we have exactly 5 cases (or at least cluster_size)
    if (nrow(trigger_window_cases) < cluster_size) {
      cat('Warning: sim', sim_key, 'has', nrow(trigger_window_cases), 'trigger cases\n')
      next
    }

    # Take first 5 cases in the trigger window
    trigger_window_cases <- trigger_window_cases[1:5, ]

    # For each position (1-5), record if it's RITA+
    rita_positions <- trigger_window_cases$rita

    # Count RITA+ by position
    n_rita_pos_1to2 <- sum(rita_positions[1:2])  # Early positions
    n_rita_pos_3to5 <- sum(rita_positions[3:5])  # Later positions
    n_rita_total <- sum(rita_positions)

    # Check if specific positions are RITA+
    has_rita_pos1 <- rita_positions[1]
    has_rita_pos2 <- rita_positions[2]
    has_rita_pos3 <- rita_positions[3]
    has_rita_pos4 <- rita_positions[4]
    has_rita_pos5 <- rita_positions[5]

    # Store detailed position info
    positions_detail <- paste(which(rita_positions), collapse = ',')

    # Get timing of earliest RITA case in trigger window
    if (n_rita_total > 0) {
      earliest_rita_pos <- min(which(rita_positions))
      earliest_rita_time <- trigger_window_cases$timediagnosed[earliest_rita_pos]
    } else {
      earliest_rita_pos <- NA
      earliest_rita_time <- NA
    }

    # Store results
    results_list[[length(results_list) + 1]] <- data.frame(
      simid = sim_key,
      IT = IT,
      ida = ida,
      pia = growth_row$pia,
      ida_per_contact_small = ida_per_contact_small,
      ida_per_contact_large = ida_per_contact_large,
      contacts_small = growth_row$contacts_small,
      contacts_large = growth_row$contacts_large,
      # Total RITA counts
      n_rita_total = n_rita_total,
      n_rita_pos_1to2 = n_rita_pos_1to2,
      n_rita_pos_3to5 = n_rita_pos_3to5,
      # Individual positions
      rita_pos1 = has_rita_pos1,
      rita_pos2 = has_rita_pos2,
      rita_pos3 = has_rita_pos3,
      rita_pos4 = has_rita_pos4,
      rita_pos5 = has_rita_pos5,
      # Position details
      rita_positions_list = positions_detail,
      earliest_rita_position = earliest_rita_pos,
      earliest_rita_time = earliest_rita_time,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results_list)
}

# Run the analysis
cat('Analyzing RITA positions in growth trigger windows...\n')
results <- analyze_rita_positions()

# Save results
output_file <- 'intervention-results/rita_position_analysis.csv'
write.csv(results, output_file, row.names = FALSE)
cat('Results saved to:', output_file, '\n\n')

# Print summary statistics
cat('=== SUMMARY STATISTICS ===\n\n')
cat('Total growth clusters analyzed:', nrow(results), '\n\n')

cat('RITA positions in trigger window (positions 1-5):\n')
cat('  Clusters with RITA+ in position 1:', sum(results$rita_pos1), '/', nrow(results),
    sprintf('(%.1f%%)\n', 100 * mean(results$rita_pos1)))
cat('  Clusters with RITA+ in position 2:', sum(results$rita_pos2), '/', nrow(results),
    sprintf('(%.1f%%)\n', 100 * mean(results$rita_pos2)))
cat('  Clusters with RITA+ in position 3:', sum(results$rita_pos3), '/', nrow(results),
    sprintf('(%.1f%%)\n', 100 * mean(results$rita_pos3)))
cat('  Clusters with RITA+ in position 4:', sum(results$rita_pos4), '/', nrow(results),
    sprintf('(%.1f%%)\n', 100 * mean(results$rita_pos4)))
cat('  Clusters with RITA+ in position 5:', sum(results$rita_pos5), '/', nrow(results),
    sprintf('(%.1f%%)\n\n', 100 * mean(results$rita_pos5)))

cat('RITA+ in early vs late positions:\n')
cat('  Mean RITA+ in positions 1-2 (early):', round(mean(results$n_rita_pos_1to2), 2), '\n')
cat('  Mean RITA+ in positions 3-5 (later):', round(mean(results$n_rita_pos_3to5), 2), '\n\n')

cat('Earliest RITA position distribution:\n')
for (pos in 1:5) {
  n_earliest_at_pos <- sum(results$earliest_rita_position == pos, na.rm = TRUE)
  cat('  Earliest RITA at position', pos, ':', n_earliest_at_pos, 'clusters\n')
}
cat('\n')

# Compare IDA/contact by earliest RITA position
cat('IDA per contact by earliest RITA position (small network):\n')
for (pos in 1:5) {
  subset_data <- results[results$earliest_rita_position == pos, ]
  if (nrow(subset_data) > 0) {
    mean_ida <- mean(subset_data$ida_per_contact_small, na.rm = TRUE)
    median_ida <- median(subset_data$ida_per_contact_small, na.rm = TRUE)
    cat('  Position', pos, '(n=', nrow(subset_data), '):',
        'mean =', round(mean_ida, 1), 'days/contact,',
        'median =', round(median_ida, 1), 'days/contact\n')
  }
}
cat('\n')

# Compare IDA/contact: has RITA in position 1 vs doesn't
cat('IDA per contact by whether RITA+ in position 1 (small network):\n')
cat('  With RITA+ in pos 1 (n=', sum(results$rita_pos1), '):',
    'mean =', round(mean(results$ida_per_contact_small[results$rita_pos1], na.rm = TRUE), 1),
    'days/contact\n')
cat('  Without RITA+ in pos 1 (n=', sum(!results$rita_pos1), '):',
    'mean =', round(mean(results$ida_per_contact_small[!results$rita_pos1], na.rm = TRUE), 1),
    'days/contact\n\n')

# Compare early vs late RITA positions
cat('IDA per contact by early vs late RITA positions (small network):\n')
results$has_early_rita <- results$n_rita_pos_1to2 > 0
results$has_late_rita <- results$n_rita_pos_3to5 > 0
cat('  With RITA+ in positions 1-2 (n=', sum(results$has_early_rita), '):',
    'mean =', round(mean(results$ida_per_contact_small[results$has_early_rita], na.rm = TRUE), 1),
    'days/contact\n')
cat('  Without RITA+ in positions 1-2 (n=', sum(!results$has_early_rita), '):',
    'mean =', round(mean(results$ida_per_contact_small[!results$has_early_rita], na.rm = TRUE), 1),
    'days/contact\n\n')

# Correlation analysis
cat('=== CORRELATION ANALYSIS ===\n\n')

# Early positions (1-2) vs IDA per contact
cor_early_small <- cor(results$n_rita_pos_1to2, results$ida_per_contact_small, use = "complete.obs")
cor_early_large <- cor(results$n_rita_pos_1to2, results$ida_per_contact_large, use = "complete.obs")
cat('RITA+ in early positions (1-2) vs IDA/contact:\n')
cat('  Small network: r =', round(cor_early_small, 3), '\n')
cat('  Large network: r =', round(cor_early_large, 3), '\n\n')

# Late positions (3-5) vs IDA per contact
cor_late_small <- cor(results$n_rita_pos_3to5, results$ida_per_contact_small, use = "complete.obs")
cor_late_large <- cor(results$n_rita_pos_3to5, results$ida_per_contact_large, use = "complete.obs")
cat('RITA+ in later positions (3-5) vs IDA/contact:\n')
cat('  Small network: r =', round(cor_late_small, 3), '\n')
cat('  Large network: r =', round(cor_late_large, 3), '\n\n')

# Earliest RITA position vs IDA per contact
cor_earliest_small <- cor(results$earliest_rita_position, results$ida_per_contact_small,
                          use = "complete.obs")
cor_earliest_large <- cor(results$earliest_rita_position, results$ida_per_contact_large,
                          use = "complete.obs")
cat('Earliest RITA position (1=early, 5=late) vs IDA/contact:\n')
cat('  Small network: r =', round(cor_earliest_small, 3), '\n')
cat('  Large network: r =', round(cor_earliest_large, 3), '\n\n')

# Create visualizations
cat('Creating visualizations...\n')

# Plot 1: IDA per contact by earliest RITA position
p1 <- ggplot(results[!is.na(results$earliest_rita_position), ],
             aes(x = factor(earliest_rita_position), y = ida_per_contact_small)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red") +
  labs(
    title = "IDA per Contact by Position of Earliest RITA+ Case",
    subtitle = "Red diamond = mean; Small network assumption",
    x = "Position of earliest RITA+ case in trigger window (1=first sequenced)",
    y = "IDA per contact (days/contact)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave('intervention-results/rita_earliest_position_boxplot.png', p1, width = 10, height = 6)

# Plot 2: Number of RITA+ in early (1-2) vs late (3-5) positions
p2 <- ggplot(results, aes(x = n_rita_pos_1to2, y = ida_per_contact_small)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "IDA per Contact vs RITA+ Cases in Early Positions (1-2)",
    subtitle = "Small network assumption",
    x = "Number of RITA+ cases in positions 1-2 (of 5 trigger cases)",
    y = "IDA per contact (days/contact)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave('intervention-results/rita_early_positions_scatter.png', p2, width = 8, height = 6)

# Plot 3: Number of RITA+ in late positions
p3 <- ggplot(results, aes(x = n_rita_pos_3to5, y = ida_per_contact_small)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "IDA per Contact vs RITA+ Cases in Later Positions (3-5)",
    subtitle = "Small network assumption",
    x = "Number of RITA+ cases in positions 3-5 (of 5 trigger cases)",
    y = "IDA per contact (days/contact)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave('intervention-results/rita_late_positions_scatter.png', p3, width = 8, height = 6)

# Plot 4: Comparison by whether position 1 is RITA+
results$pos1_label <- ifelse(results$rita_pos1, "RITA+ in position 1", "No RITA+ in position 1")
p4 <- ggplot(results, aes(x = pos1_label, y = ida_per_contact_small, fill = pos1_label)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red") +
  scale_fill_manual(values = c("lightcoral", "lightgreen")) +
  labs(
    title = "IDA per Contact: First Case RITA+ vs Not",
    subtitle = "Red diamond = mean; Small network assumption",
    x = "",
    y = "IDA per contact (days/contact)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")

ggsave('intervention-results/rita_position1_comparison.png', p4, width = 8, height = 6)

# Plot 5: Heatmap of RITA positions
# Create a matrix showing which positions are RITA+ for each cluster
position_matrix <- data.frame(
  cluster_id = 1:nrow(results),
  simid = results$simid,
  ida_per_contact = results$ida_per_contact_small,
  pos1 = results$rita_pos1,
  pos2 = results$rita_pos2,
  pos3 = results$rita_pos3,
  pos4 = results$rita_pos4,
  pos5 = results$rita_pos5
)

# Order by IDA per contact
position_matrix <- position_matrix[order(position_matrix$ida_per_contact), ]
position_matrix$cluster_order <- 1:nrow(position_matrix)

# Reshape for heatmap
position_long <- position_matrix %>%
  select(cluster_order, ida_per_contact, pos1, pos2, pos3, pos4, pos5) %>%
  pivot_longer(cols = starts_with("pos"),
               names_to = "position",
               values_to = "is_rita") %>%
  mutate(position = as.numeric(gsub("pos", "", position)))

p5 <- ggplot(position_long, aes(x = factor(position), y = factor(cluster_order), fill = is_rita)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_manual(values = c("white", "darkblue"),
                    labels = c("Not RITA+", "RITA+"),
                    name = "Status") +
  labs(
    title = "RITA+ Position Patterns Across Growth Clusters",
    subtitle = "Clusters ordered by IDA per contact (bottom = lowest, top = highest)",
    x = "Position in trigger window (1 = first sequenced)",
    y = "Cluster (ordered by IDA/contact efficiency)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave('intervention-results/rita_position_heatmap.png', p5, width = 8, height = 10)

# Plot 6: Scatter plot with earliest RITA position colored
p6 <- ggplot(results[!is.na(results$earliest_rita_position), ],
             aes(x = n_rita_total, y = ida_per_contact_small,
                 color = factor(earliest_rita_position))) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_viridis_d(option = "plasma", name = "Earliest\nRITA+\nPosition") +
  labs(
    title = "Total RITA+ Cases vs IDA per Contact",
    subtitle = "Colored by position of earliest RITA+ case (small network)",
    x = "Total number of RITA+ cases in trigger window (of 5)",
    y = "IDA per contact (days/contact)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave('intervention-results/rita_total_vs_efficiency_colored.png', p6, width = 10, height = 6)

cat('\nPlots saved to intervention-results/\n')
cat('  - rita_earliest_position_boxplot.png\n')
cat('  - rita_early_positions_scatter.png\n')
cat('  - rita_late_positions_scatter.png\n')
cat('  - rita_position1_comparison.png\n')
cat('  - rita_position_heatmap.png\n')
cat('  - rita_total_vs_efficiency_colored.png\n\n')

cat('Analysis complete!\n')
