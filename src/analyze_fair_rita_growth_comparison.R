# Fair comparison of RITA vs Growth detection delays
# Examines whether the apparent ~4 year advantage of RITA is due to:
# 1. Actual strategy differences (diagnosis → intervention delays)
# 2. Population differences (infection → diagnosis delays for RITA vs non-RITA cases)

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

cat('Analyzing delays from different time origins...\n')

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

  # Growth intervention time
  IT_growth <- t_detect + analysis_delay_days + implementation_delay_days

  # For each cluster member, compare RITA vs Growth delays
  Gcluster_members <- Gcluster[Gcluster$generation != lastgen, ]

  for (k in 1:nrow(Gcluster_members)) {
    member <- Gcluster_members[k, ]

    # Check if this member is RITA+
    is_rita <- member$rita

    # RITA intervention time (if they are RITA+)
    IT_rita <- member$timediagnosed + implementation_delay_days

    # Calculate delays from different origins
    # 1. From INFECTION (includes non-modifiable infection→diagnosis delay)
    delay_from_infection_growth <- IT_growth - member$timeinfected
    delay_from_infection_rita <- IT_rita - member$timeinfected

    # 2. From DIAGNOSIS (only modifiable delays)
    delay_from_diagnosis_growth <- IT_growth - member$timediagnosed
    delay_from_diagnosis_rita <- IT_rita - member$timediagnosed

    # 3. From SEQUENCING (for growth members who are sequenced)
    if (!is.na(member$timesequenced) && member$timesequenced < IT_growth) {
      delay_from_sequencing_growth <- IT_growth - member$timesequenced
    } else {
      delay_from_sequencing_growth <- NA
    }

    # 4. Infection to diagnosis delay (non-modifiable)
    infection_to_diagnosis <- member$timediagnosed - member$timeinfected

    results_list[[length(results_list) + 1]] <- data.frame(
      simid = sim_key,
      pid = member$pid,
      is_rita = is_rita,
      # Absolute times
      time_infected = member$timeinfected,
      time_diagnosed = member$timediagnosed,
      time_sequenced = member$timesequenced,
      IT_growth = IT_growth,
      IT_rita = IT_rita,
      # Non-modifiable delay
      infection_to_diagnosis = infection_to_diagnosis,
      # Total delays from infection (includes non-modifiable component)
      delay_from_infection_growth = delay_from_infection_growth,
      delay_from_infection_rita = delay_from_infection_rita,
      diff_from_infection = delay_from_infection_growth - delay_from_infection_rita,
      # Modifiable delays from diagnosis
      delay_from_diagnosis_growth = delay_from_diagnosis_growth,
      delay_from_diagnosis_rita = delay_from_diagnosis_rita,
      diff_from_diagnosis = delay_from_diagnosis_growth - delay_from_diagnosis_rita,
      # From sequencing (for growth)
      delay_from_sequencing_growth = delay_from_sequencing_growth,
      stringsAsFactors = FALSE
    )
  }
}

results <- do.call(rbind, results_list)

# Save results
output_file <- 'intervention-results/fair_rita_growth_comparison.csv'
write.csv(results, output_file, row.names = FALSE)
cat('Results saved to:', output_file, '\n\n')

# Summary statistics
cat('=== SUMMARY: RITA vs GROWTH DELAY COMPARISON ===\n\n')

cat('1. NON-MODIFIABLE DELAY (Infection → Diagnosis):\n')
cat('   RITA+ cases: mean =', round(mean(results$infection_to_diagnosis[results$is_rita]), 1), 'days\n')
cat('   Non-RITA cases: mean =', round(mean(results$infection_to_diagnosis[!results$is_rita]), 1), 'days\n')
cat('   Difference:', round(mean(results$infection_to_diagnosis[!results$is_rita]) -
                             mean(results$infection_to_diagnosis[results$is_rita]), 1), 'days\n\n')

cat('2. DELAYS FROM INFECTION (includes non-modifiable component):\n')
rita_cases <- results[results$is_rita, ]
cat('   For RITA+ cases:\n')
cat('     Growth delay: mean =', round(mean(rita_cases$delay_from_infection_growth), 1), 'days\n')
cat('     RITA delay: mean =', round(mean(rita_cases$delay_from_infection_rita), 1), 'days\n')
cat('     Difference (Growth - RITA): mean =', round(mean(rita_cases$diff_from_infection), 1), 'days\n')
cat('     Difference in years:', round(mean(rita_cases$diff_from_infection) / 365, 2), 'years\n\n')

cat('3. DELAYS FROM DIAGNOSIS (only modifiable by strategies):\n')
cat('   For RITA+ cases:\n')
cat('     Growth delay: mean =', round(mean(rita_cases$delay_from_diagnosis_growth), 1), 'days\n')
cat('     RITA delay: mean =', round(mean(rita_cases$delay_from_diagnosis_rita), 1), 'days\n')
cat('     Difference (Growth - RITA): mean =', round(mean(rita_cases$diff_from_diagnosis), 1), 'days\n')
cat('     Difference in years:', round(mean(rita_cases$diff_from_diagnosis) / 365, 2), 'years\n\n')

cat('4. GROWTH CLUSTER DELAY FROM SEQUENCING:\n')
valid_seq <- results[!is.na(results$delay_from_sequencing_growth), ]
cat('   Mean delay from sequencing to IT:', round(mean(valid_seq$delay_from_sequencing_growth), 1), 'days\n\n')

# Create visualizations
cat('Creating visualizations...\n')

# Plot 1: Comparison of delays from infection vs diagnosis
plot_data <- results %>%
  filter(is_rita) %>%
  select(simid, pid, diff_from_infection, diff_from_diagnosis) %>%
  pivot_longer(cols = c(diff_from_infection, diff_from_diagnosis),
               names_to = "comparison_type",
               values_to = "delay_difference") %>%
  mutate(comparison_type = ifelse(comparison_type == "diff_from_infection",
                                  "From Infection\n(includes non-modifiable delay)",
                                  "From Diagnosis\n(modifiable by strategies)"))

p1 <- ggplot(plot_data, aes(x = comparison_type, y = delay_difference / 365)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "darkblue") +
  labs(
    title = "Growth vs RITA Delay: Does the Comparison Include Non-Modifiable Delays?",
    subtitle = "For RITA+ cases in growth clusters (positive = Growth slower than RITA)",
    x = "",
    y = "Delay difference: Growth - RITA (years)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave('intervention-results/fair_comparison_delays.png', p1, width = 10, height = 6)

# Plot 2: Infection-to-diagnosis delay by RITA status
p2 <- ggplot(results, aes(x = is_rita, y = infection_to_diagnosis / 365, fill = is_rita)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red") +
  scale_fill_manual(values = c("gray70", "orange"),
                    labels = c("Non-RITA", "RITA+"),
                    name = "Status") +
  scale_x_discrete(labels = c("Non-RITA\n(not recently infected)", "RITA+\n(recently infected)")) +
  labs(
    title = "Infection-to-Diagnosis Delay: RITA+ vs Non-RITA Cases",
    subtitle = "This delay is NOT modifiable by intervention strategies",
    x = "",
    y = "Time from infection to diagnosis (years)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")

ggsave('intervention-results/infection_to_diagnosis_by_rita.png', p2, width = 8, height = 6)

# Plot 3: Breakdown of total delay components
delay_components <- rita_cases %>%
  select(infection_to_diagnosis, delay_from_diagnosis_growth) %>%
  summarise(
    infection_to_diagnosis = mean(infection_to_diagnosis),
    diagnosis_to_IT = mean(delay_from_diagnosis_growth)
  ) %>%
  pivot_longer(everything(), names_to = "component", values_to = "days") %>%
  mutate(
    component = factor(component,
                      levels = c("infection_to_diagnosis", "diagnosis_to_IT"),
                      labels = c("Infection → Diagnosis\n(non-modifiable)",
                                "Diagnosis → Intervention\n(growth strategy)")),
    modifiable = c(FALSE, TRUE)
  )

p3 <- ggplot(delay_components, aes(x = component, y = days, fill = modifiable)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.0f days\n(%.1f years)", days, days/365)),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("gray50", "steelblue"),
                    labels = c("Non-modifiable", "Modifiable"),
                    name = "Strategy Control") +
  labs(
    title = "Growth Cluster Delay Components for RITA+ Cases",
    subtitle = "Total time from infection to intervention",
    x = "",
    y = "Days"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave('intervention-results/delay_components_breakdown.png', p3, width = 8, height = 6)

cat('\nPlots saved to intervention-results/\n')
cat('  - fair_comparison_delays.png\n')
cat('  - infection_to_diagnosis_by_rita.png\n')
cat('  - delay_components_breakdown.png\n\n')

cat('CONCLUSION:\n')
cat('When comparing RITA vs Growth from INFECTION, the difference includes non-modifiable\n')
cat('population differences (RITA+ have shorter infection-to-diagnosis by definition).\n')
cat('A fairer comparison measures from DIAGNOSIS, showing only strategy-modifiable delays.\n\n')

cat('Analysis complete!\n')
