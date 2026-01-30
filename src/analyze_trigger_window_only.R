# Compare RITA vs Growth delays for TRIGGER WINDOW members only
# Focus only on the 5 cases that meet the growth cluster criterion (5-in-6-months)
# This is a fair comparison of the same population under two strategies

library(dplyr)
library(ggplot2)
library(tidyr)

source('src/intervention0.R')

cat('Loading simulation data...\n')
Dall <- read.csv('src/experiment1-N10000-gens7-D.csv')
Gall <- read.csv('src/experiment1-N10000-gens7-G.csv')

# Parameters
growth_distance_threshold <- 0.01
cluster_size <- 5
lookback_window_months <- 6
lookback_days <- lookback_window_months * 30
rita_window_months <- 6
analysis_delay_days <- 14
implementation_delay_days <- 14

set.seed(123)

# Simulate RITA status
cat('Simulating RITA status...\n')
Gall$rita <- with(Gall, (timediagnosed - timeinfected) < rexp(nrow(Gall), 1/(rita_window_months*30)))

# Split data
simids <- as.character(unique(Dall$simid))
Ds <- split(Dall, Dall$simid)[simids]
Gs <- split(Gall, Gall$simid)[simids]

cat('Analyzing trigger window members only...\n')

results_list <- list()

for (sim_key in simids) {
  D <- Ds[[sim_key]]
  G <- Gs[[sim_key]]

  G <- G[order(G$timesequenced), ]
  lastgen <- max(G$generation)

  # Build cluster
  D1 <- D[D$distance <= growth_distance_threshold, ]
  keeppids <- '0'
  addpids <- D1$recipient[D1$donor %in% keeppids]
  while (length(addpids) > 0) {
    keeppids <- union(addpids, keeppids)
    addpids <- setdiff(D1$recipient[D1$donor %in% keeppids], keeppids)
  }
  Gcluster <- G[G$pid %in% keeppids, ]

  # Find growth trigger
  Gtrig <- Gcluster[Gcluster$generation > 0 & Gcluster$generation != lastgen, ]
  if (nrow(Gtrig) < cluster_size) next

  t <- as.numeric(Gtrig$timesequenced)
  n <- length(t)
  j <- 1
  t_detect <- Inf
  trigger_idx <- NA
  window_start_idx <- NA

  for (i in seq_len(n)) {
    while (j <= i && (t[i] - t[j]) > lookback_days) {
      j <- j + 1
    }
    if ((i - j + 1) >= cluster_size) {
      t_detect <- t[i]
      trigger_idx <- i
      window_start_idx <- j
      break
    }
  }

  if (!is.finite(t_detect)) next

  # Growth intervention time
  IT_growth <- t_detect + analysis_delay_days + implementation_delay_days

  # FOCUS ON TRIGGER WINDOW ONLY: the 5 cases that triggered
  trigger_window <- Gtrig[window_start_idx:trigger_idx, ]

  # For each member of the trigger window
  for (k in 1:nrow(trigger_window)) {
    member <- trigger_window[k, ]

    # Position in trigger window (1 = first sequenced, 5 = last/trigger)
    position <- k

    # RITA intervention time (if they are RITA+)
    IT_rita <- member$timediagnosed + implementation_delay_days

    # Calculate delays
    # 1. From INFECTION
    delay_from_infection_growth <- IT_growth - member$timeinfected
    delay_from_infection_rita <- IT_rita - member$timeinfected
    diff_from_infection <- delay_from_infection_growth - delay_from_infection_rita

    # 2. From DIAGNOSIS (modifiable)
    delay_from_diagnosis_growth <- IT_growth - member$timediagnosed
    delay_from_diagnosis_rita <- IT_rita - member$timediagnosed
    diff_from_diagnosis <- delay_from_diagnosis_growth - delay_from_diagnosis_rita

    # 3. From SEQUENCING (for growth)
    delay_from_sequencing_growth <- IT_growth - member$timesequenced

    # 4. How much earlier would RITA have intervened?
    if (member$rita) {
      time_saved_by_rita <- IT_growth - IT_rita
    } else {
      time_saved_by_rita <- NA  # Not RITA+, so RITA wouldn't catch them
    }

    results_list[[length(results_list) + 1]] <- data.frame(
      simid = sim_key,
      pid = member$pid,
      position = position,
      is_rita = member$rita,
      # Absolute times
      time_infected = member$timeinfected,
      time_diagnosed = member$timediagnosed,
      time_sequenced = member$timesequenced,
      IT_growth = IT_growth,
      IT_rita = IT_rita,
      # Delays from different origins
      delay_from_infection_growth = delay_from_infection_growth,
      delay_from_infection_rita = delay_from_infection_rita,
      diff_from_infection = diff_from_infection,
      delay_from_diagnosis_growth = delay_from_diagnosis_growth,
      delay_from_diagnosis_rita = delay_from_diagnosis_rita,
      diff_from_diagnosis = diff_from_diagnosis,
      delay_from_sequencing_growth = delay_from_sequencing_growth,
      time_saved_by_rita = time_saved_by_rita,
      stringsAsFactors = FALSE
    )
  }
}

results <- do.call(rbind, results_list)

# Save results
output_file <- 'intervention-results/trigger_window_comparison.csv'
write.csv(results, output_file, row.names = FALSE)
cat('Results saved to:', output_file, '\n\n')

# Summary statistics
cat('=== TRIGGER WINDOW ONLY: RITA vs GROWTH COMPARISON ===\n\n')

cat('Sample: Only the 5 cases in each trigger window (5-in-6-months)\n')
cat('Total trigger window members analyzed:', nrow(results), '\n')
cat('RITA+ members:', sum(results$is_rita), '/', nrow(results),
    sprintf('(%.1f%%)\n\n', 100 * mean(results$is_rita)))

cat('1. DELAYS FROM DIAGNOSIS (modifiable by strategies):\n')
cat('   All trigger window members:\n')
cat('     Growth delay: mean =', round(mean(results$delay_from_diagnosis_growth), 1), 'days\n')
cat('     RITA delay (if RITA+): mean =', round(mean(results$delay_from_diagnosis_rita[results$is_rita]), 1), 'days\n\n')

cat('   RITA+ members only:\n')
rita_members <- results[results$is_rita, ]
cat('     Growth delay: mean =', round(mean(rita_members$delay_from_diagnosis_growth), 1), 'days (',
    round(mean(rita_members$delay_from_diagnosis_growth) / 365, 2), 'years)\n')
cat('     RITA delay: mean =', round(mean(rita_members$delay_from_diagnosis_rita), 1), 'days (',
    round(mean(rita_members$delay_from_diagnosis_rita) / 365, 2), 'years)\n')
cat('     Time saved by RITA: mean =', round(mean(rita_members$time_saved_by_rita, na.rm = TRUE), 1), 'days (',
    round(mean(rita_members$time_saved_by_rita, na.rm = TRUE) / 365, 2), 'years)\n\n')

cat('2. DELAYS FROM SEQUENCING (for growth strategy):\n')
cat('   Mean delay from sequencing to intervention:', round(mean(results$delay_from_sequencing_growth), 1), 'days\n\n')

cat('3. BY POSITION IN TRIGGER WINDOW:\n')
for (pos in 1:5) {
  pos_data <- results[results$position == pos, ]
  n_rita <- sum(pos_data$is_rita)
  cat('   Position', pos, '(n=', nrow(pos_data), ', RITA+=', n_rita, '):\n')
  cat('     Delay from diagnosis (growth): mean =', round(mean(pos_data$delay_from_diagnosis_growth), 1), 'days\n')
  if (n_rita > 0) {
    cat('     Time saved by RITA (RITA+ only): mean =',
        round(mean(pos_data$time_saved_by_rita, na.rm = TRUE), 1), 'days\n')
  }
}

cat('\n4. CLUSTER ACCUMULATION TIME (trigger window only):\n')
cluster_stats <- results %>%
  group_by(simid) %>%
  summarise(
    first_diagnosed = min(time_diagnosed),
    first_sequenced = min(time_sequenced),
    trigger_time = first(IT_growth) - analysis_delay_days - implementation_delay_days,
    diagnosis_to_trigger = trigger_time - first_diagnosed,
    sequencing_to_trigger = trigger_time - first_sequenced,
    .groups = 'drop'
  )

cat('   From first diagnosis in window to trigger: mean =',
    round(mean(cluster_stats$diagnosis_to_trigger), 1), 'days (',
    round(mean(cluster_stats$diagnosis_to_trigger) / 365, 2), 'years)\n')
cat('   From first sequenced in window to trigger: mean =',
    round(mean(cluster_stats$sequencing_to_trigger), 1), 'days (',
    round(mean(cluster_stats$sequencing_to_trigger) / 365, 2), 'years)\n\n')

# Create visualizations
cat('Creating visualizations...\n')

# Plot 1: Time saved by RITA (for RITA+ cases in trigger window)
p1 <- ggplot(rita_members, aes(x = time_saved_by_rita / 365)) +
  geom_histogram(fill = "orange", alpha = 0.7, bins = 30) +
  geom_vline(aes(xintercept = mean(time_saved_by_rita, na.rm = TRUE) / 365),
             color = "red", linetype = "dashed", size = 1) +
  annotate("text",
           x = mean(rita_members$time_saved_by_rita, na.rm = TRUE) / 365,
           y = Inf,
           label = sprintf("Mean: %.2f years", mean(rita_members$time_saved_by_rita, na.rm = TRUE) / 365),
           vjust = 2, hjust = -0.1, color = "red", fontface = "bold") +
  labs(
    title = "Time Saved by RITA vs Growth (Trigger Window Members Only)",
    subtitle = "For RITA+ cases in the 5-member trigger window",
    x = "Years earlier RITA would intervene compared to growth cluster",
    y = "Count"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave('intervention-results/time_saved_trigger_window.png', p1, width = 10, height = 6)

# Plot 2: Delays by position in trigger window
position_summary <- results %>%
  group_by(position, is_rita) %>%
  summarise(
    mean_delay = mean(delay_from_diagnosis_growth),
    se_delay = sd(delay_from_diagnosis_growth) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )

p2 <- ggplot(position_summary, aes(x = factor(position), y = mean_delay, fill = is_rita)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_delay - se_delay, ymax = mean_delay + se_delay),
                position = position_dodge(0.8), width = 0.3) +
  scale_fill_manual(values = c("gray70", "orange"),
                    labels = c("Non-RITA", "RITA+"),
                    name = "Status") +
  labs(
    title = "Growth Cluster Delay by Position in Trigger Window",
    subtitle = "Time from diagnosis to growth intervention (trigger window members only)",
    x = "Position in trigger window (1 = first sequenced, 5 = triggers cluster)",
    y = "Days from diagnosis to intervention"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave('intervention-results/delay_by_position_trigger_window.png', p2, width = 10, height = 6)

# Plot 3: Comparison of delays from diagnosis
comparison_data <- rita_members %>%
  select(delay_from_diagnosis_growth, delay_from_diagnosis_rita) %>%
  pivot_longer(everything(), names_to = "strategy", values_to = "delay") %>%
  mutate(strategy = ifelse(strategy == "delay_from_diagnosis_growth",
                           "Growth Cluster", "RITA"))

p3 <- ggplot(comparison_data, aes(x = strategy, y = delay / 365, fill = strategy)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red") +
  scale_fill_manual(values = c("steelblue", "orange")) +
  labs(
    title = "RITA vs Growth: Diagnosis to Intervention Delay",
    subtitle = "RITA+ cases in trigger window only (modifiable delays)",
    x = "",
    y = "Years from diagnosis to intervention"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")

ggsave('intervention-results/rita_vs_growth_trigger_window.png', p3, width = 8, height = 6)

cat('\nPlots saved to intervention-results/\n')
cat('  - time_saved_trigger_window.png\n')
cat('  - delay_by_position_trigger_window.png\n')
cat('  - rita_vs_growth_trigger_window.png\n\n')

cat('Analysis complete!\n')
