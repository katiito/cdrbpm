# Check the 97.5th percentile values to inform axis limits

source("src/intervention0.R")
source("src/plot_interventions.R")
source("src/run_analysis.R")

# Load cached results
results <- load_cached_results()

# Extract data (same as in plot function)
strategy_names <- c("distsize5", "distsize2", "growth", "random", "rita", "network")
strategy_labels <- c("Size>=5", "Size>=2", "Growth", "Random", "RITA", "Network")

extract_data <- function(details) {
  dfs <- lapply(seq_along(strategy_names), function(i) {
    sname <- strategy_names[i]
    slabel <- strategy_labels[i]

    if (is.null(details[[sname]]) || is.null(details[[sname]]$o)) {
      return(NULL)
    }

    o <- details[[sname]]$o

    if (sname %in% c("network", "random", "rita")) {
      if (!"contacts" %in% names(o)) return(NULL)
      data.frame(
        strategy = slabel,
        puta = o$puta,
        pia = o$pia,
        contacts_small = o$contacts,
        contacts_large = o$contacts
      )
    } else {
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

# Compute efficiencies
df$puta_eff_small <- df$puta / df$contacts_small
df$puta_eff_large <- df$puta / df$contacts_large
df$pia_eff_small <- df$pia / df$contacts_small
df$pia_eff_large <- df$pia / df$contacts_large

# Remove invalid values
df <- df[is.finite(df$puta_eff_small) & is.finite(df$puta_eff_large) &
         is.finite(df$pia_eff_small) & is.finite(df$pia_eff_large), ]

cat("97.5th Percentile Values:\n")
cat("=========================\n\n")

cat("PUTA Small Subnetwork:\n")
cat(sprintf("  Max: %.2f\n", max(df$puta_eff_small)))
cat(sprintf("  97.5th percentile: %.2f\n\n", quantile(df$puta_eff_small, 0.975)))

cat("PUTA Large Subnetwork:\n")
cat(sprintf("  Max: %.2f\n", max(df$puta_eff_large)))
cat(sprintf("  97.5th percentile: %.2f\n\n", quantile(df$puta_eff_large, 0.975)))

cat("PIA Small Subnetwork:\n")
cat(sprintf("  Max: %.2f\n", max(df$pia_eff_small)))
cat(sprintf("  97.5th percentile: %.2f\n\n", quantile(df$pia_eff_small, 0.975)))

cat("PIA Large Subnetwork:\n")
cat(sprintf("  Max: %.2f\n", max(df$pia_eff_large)))
cat(sprintf("  97.5th percentile: %.2f\n\n", quantile(df$pia_eff_large, 0.975)))
