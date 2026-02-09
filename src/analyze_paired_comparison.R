# Paired comparison of strategies vs random baseline
# For each simulation ID, compare strategy outcomes to random strategy
# This controls for simulation-specific network effects

library(dplyr)
library(ggplot2)
library(tidyr)
library(here)

# Source required scripts to get load_cached_results() function
source(here::here("src", "run_analysis.R"))

cat('Loading cached intervention results...\n')
cached_results <- load_cached_results()

# Extract details and calculate efficiency metrics for large subnetwork
cat('\nCalculating efficiency metrics...\n')
results_list <- list()

for (strategy_name in names(cached_results$details)) {
  data <- cached_results$details[[strategy_name]]$o

  cat('Processing strategy:', strategy_name, 'with', nrow(data), 'rows\n')
  cat('  Columns:', paste(names(data)[1:min(5, ncol(data))], collapse = ', '), '...\n')

  # Handle different column names: cluster-based strategies use "puta", individual-based use "ida"
  ida_col <- if ("puta" %in% names(data)) "puta" else "ida"

  # For individual-based strategies, contacts_large doesn't exist, just "contacts"
  contacts_col_large <- if ("contacts_large" %in% names(data)) "contacts_large" else "contacts"

  # puta/ida = IDA (Infectious Days Averted), pia = PIA (Possible Infections Averted)
  # Filter out rows with 0 contacts (efficiency per contact is undefined)
  data_with_efficiency <- data %>%
    filter(.data[[contacts_col_large]] > 0) %>%
    mutate(
      intervention = strategy_name,
      ida_per_contact = .data[[ida_col]] / .data[[contacts_col_large]],
      pia_per_contact = pia / .data[[contacts_col_large]]
    ) %>%
    select(simid, intervention, ida_per_contact, pia_per_contact)

  results_list[[strategy_name]] <- data_with_efficiency
}

results_large <- bind_rows(results_list)

cat('Loaded strategies:\n')
print(table(results_large$intervention))
cat('\n')

# IMPORTANT: Aggregate to simulation level first
# Each strategy may have multiple observations per simid (e.g., multiple clusters, multiple random samples)
# We need to average within each simulation before doing paired comparison
cat('Aggregating to simulation level...\n')
results_by_sim <- results_large %>%
  group_by(simid, intervention) %>%
  summarise(
    ida_per_contact = mean(ida_per_contact, na.rm = TRUE),
    pia_per_contact = mean(pia_per_contact, na.rm = TRUE),
    n_obs = n(),  # Track how many observations were averaged
    .groups = 'drop'
  )

cat('After aggregation:\n')
print(results_by_sim %>% group_by(intervention) %>% summarise(n_sims = n(), mean_obs_per_sim = mean(n_obs)))
cat('\n')

# Get random strategy results (now one observation per simid)
random_results <- results_by_sim %>%
  filter(intervention == "random") %>%
  select(simid,
         ida_random = ida_per_contact,
         pia_random = pia_per_contact)

# For each non-random strategy, calculate paired differences
paired_results_list <- list()

strategies <- unique(results_by_sim$intervention)
strategies <- strategies[strategies != "random"]

for (strat in strategies) {
  strat_data <- results_by_sim %>%
    filter(intervention == strat) %>%
    select(simid, intervention, ida_per_contact, pia_per_contact)

  # Join with random results (only keeps simids present in both)
  paired <- inner_join(strat_data, random_results, by = "simid")

  # Calculate differences (strategy - random)
  paired <- paired %>%
    mutate(
      ida_diff = ida_per_contact - ida_random,
      pia_diff = pia_per_contact - pia_random,
      # Only calculate ratios/percentages when baseline is non-zero
      ida_ratio = ifelse(ida_random > 0, ida_per_contact / ida_random, NA),
      pia_ratio = ifelse(pia_random > 0, pia_per_contact / pia_random, NA),
      ida_pct_change = ifelse(ida_random > 0, (ida_per_contact - ida_random) / ida_random * 100, NA),
      pia_pct_change = ifelse(pia_random > 0, (pia_per_contact - pia_random) / pia_random * 100, NA)
    )

  paired_results_list[[strat]] <- paired
}

# Combine all strategies
paired_results <- bind_rows(paired_results_list)

# Save results
output_file <- 'intervention-results/paired_comparison_results.csv'
write.csv(paired_results, output_file, row.names = FALSE)
cat('Paired comparison results saved to:', output_file, '\n\n')

# Summary statistics
cat('=== PAIRED COMPARISON SUMMARY (vs Random Baseline) ===\n\n')

# Report baseline zero counts
cat('Random baseline characteristics:\n')
cat('  IDA = 0:', sum(paired_results$ida_random == 0), '/', nrow(paired_results),
    sprintf('(%.1f%%)\n', 100 * mean(paired_results$ida_random == 0)))
cat('  PIA = 0:', sum(paired_results$pia_random == 0), '/', nrow(paired_results),
    sprintf('(%.1f%%)\n', 100 * mean(paired_results$pia_random == 0)))
cat('  Note: Ratios and percent changes are only calculated for non-zero baseline cases\n\n')

cat('Strategy-specific sample sizes (simulations with both random and strategy data):\n')
sample_sizes <- paired_results %>%
  group_by(intervention) %>%
  summarise(n_paired = n(), .groups = 'drop')
print(sample_sizes)
cat('\n')

summary_stats <- paired_results %>%
  group_by(intervention) %>%
  summarise(
    n = n(),
    # Absolute differences (strategy - random) - uses all cases
    mean_ida_diff = mean(ida_diff),
    median_ida_diff = median(ida_diff),
    mean_pia_diff = mean(pia_diff),
    median_pia_diff = median(pia_diff),
    # Ratios (strategy / random) - only non-zero baseline
    n_nonzero_ida = sum(!is.na(ida_ratio)),
    mean_ida_ratio = mean(ida_ratio, na.rm = TRUE),
    median_ida_ratio = median(ida_ratio, na.rm = TRUE),
    n_nonzero_pia = sum(!is.na(pia_ratio)),
    mean_pia_ratio = mean(pia_ratio, na.rm = TRUE),
    median_pia_ratio = median(pia_ratio, na.rm = TRUE),
    # Percent improvements - only non-zero baseline
    mean_ida_pct = mean(ida_pct_change, na.rm = TRUE),
    median_ida_pct = median(ida_pct_change, na.rm = TRUE),
    mean_pia_pct = mean(pia_pct_change, na.rm = TRUE),
    median_pia_pct = median(pia_pct_change, na.rm = TRUE),
    .groups = 'drop'
  )

cat('Mean absolute differences (positive = strategy better than random):\n')
cat('Uses ALL paired cases\n')
print(summary_stats %>% select(intervention, n, mean_ida_diff, mean_pia_diff))
cat('\n')

cat('Mean ratios (>1 = strategy better than random):\n')
cat('Only uses cases where random baseline > 0\n')
print(summary_stats %>% select(intervention, n, n_nonzero_ida, mean_ida_ratio, n_nonzero_pia, mean_pia_ratio))
cat('\n')

cat('Mean percent improvement over random:\n')
cat('Only uses cases where random baseline > 0\n')
print(summary_stats %>% select(intervention, n, n_nonzero_ida, mean_ida_pct, n_nonzero_pia, mean_pia_pct))
cat('\n')

# Create visualizations
cat('Creating paired comparison plots...\n')

# Prepare data for plotting (using differences)
plot_data_ida <- paired_results %>%
  select(simid, intervention, ida_diff) %>%
  mutate(metric = "Infectious Days Averted")

plot_data_pia <- paired_results %>%
  select(simid, intervention, pia_diff) %>%
  rename(ida_diff = pia_diff) %>%
  mutate(metric = "Possible Infections Averted")

plot_data <- bind_rows(plot_data_ida, plot_data_pia) %>%
  rename(value = ida_diff)

# Define intervention order and labels (matching existing code)
# Map internal names to display labels
intervention_name_map <- c(
  "growth" = "Growth",
  "rita" = "RITA",
  "ritasecondary" = "RITA+Secondary",
  "distsize5" = "Size>=5",
  "distsize2" = "Size>=2",
  "network" = "Network"
)

# Order for plotting
intervention_order <- c("rita", "ritasecondary", "growth", "distsize5", "network")
intervention_labels <- intervention_name_map[intervention_order]

# Colors (matching existing plots)
intervention_colors <- c(
  "rita" = "#E41A1C",
  "ritasecondary" = "#F781BF",
  "growth" = "#377EB8",
  "distsize5" = "#4DAF4A",
  "network" = "#984EA3"
)

plot_data <- plot_data %>%
  filter(intervention %in% intervention_order) %>%
  mutate(
    intervention = factor(intervention, levels = intervention_order),
    intervention_label = intervention_name_map[as.character(intervention)]
  )

# Violin plot (similar to efficiency_distributions)
p_paired <- ggplot(plot_data, aes(x = intervention_label, y = value, fill = intervention)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  scale_fill_manual(values = intervention_colors) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(
    title = "Paired Comparison: Strategy Performance Relative to Random Baseline",
    subtitle = "Difference per contact (Strategy - Random) | Large subnetwork | Matched simulations",
    x = "",
    y = "Difference from random baseline\n(positive = better than random)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

# Save plot
ggsave('intervention-results/paired_comparison_distributions.png',
       p_paired, width = 10, height = 8, dpi = 300)

cat('\nPlot saved to: intervention-results/paired_comparison_distributions.png\n')

# Additional plot: Percent improvement
plot_data_pct_ida <- paired_results %>%
  select(simid, intervention, ida_pct_change) %>%
  mutate(metric = "Infectious Days Averted")

plot_data_pct_pia <- paired_results %>%
  select(simid, intervention, pia_pct_change) %>%
  rename(ida_pct_change = pia_pct_change) %>%
  mutate(metric = "Possible Infections Averted")

plot_data_pct <- bind_rows(plot_data_pct_ida, plot_data_pct_pia) %>%
  rename(value = ida_pct_change) %>%
  filter(intervention %in% intervention_order) %>%
  mutate(
    intervention = factor(intervention, levels = intervention_order),
    intervention_label = intervention_name_map[as.character(intervention)]
  )

p_pct <- ggplot(plot_data_pct, aes(x = intervention_label, y = value, fill = intervention)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  scale_fill_manual(values = intervention_colors) +
  coord_cartesian(ylim = c(-1000, 1000)) +  # Focus on main distribution, extreme outliers excluded from view
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(
    title = "Paired Comparison: Percent Improvement Over Random Baseline",
    subtitle = "Percent change relative to random | Large subnetwork | Matched simulations\nBlack diamond = mean, Boxplot line = median | Y-axis limited to [-1000%, 1000%]",
    x = "",
    y = "Percent improvement over random\n(positive = better than random)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

ggsave('intervention-results/paired_comparison_percent.png',
       p_pct, width = 10, height = 8, dpi = 300)

cat('Plot saved to: intervention-results/paired_comparison_percent.png\n')

cat('\nAnalysis complete!\n')
cat('\nInterpretation:\n')
cat('- Positive values = strategy performs better than random baseline\n')
cat('- Negative values = strategy performs worse than random baseline\n')
cat('- Dashed red line at 0 = no difference from random\n')
cat('- Black diamond = mean difference\n')
