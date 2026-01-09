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
#' 
#' @param results The output from run_interventions()
#' @param plot_type Either "density" for smoothed density plots, "violin" for 
#'   violin plots, or "boxplot" for box plots
#' @param title_prefix Optional prefix for plot titles
#' @return A combined ggplot object with 4 panels
#' 
plot_efficiency_distributions <- function(results, 
                                          plot_type = "density",
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
  
  # Set factor levels for ordering
  df$strategy <- factor(df$strategy, levels = strategy_labels)
  
  # Define color palette
  strategy_colors <- c(
    "Size>=5" = "#E41A1C",
    "Size>=2" = "#377EB8", 
    "Growth" = "#4DAF4A",
    "Random" = "#984EA3",
    "RITA" = "#FF7F00",
    "Network" = "#A65628"
  )
  
  # Create plots based on plot_type
  if (plot_type == "density") {
    p1 <- ggplot(df, aes(x = puta_eff_small, fill = strategy, color = strategy)) +
      geom_density(alpha = 0.3) +
      labs(x = "PUTA per contact", y = "Density",
           title = paste0(title_prefix, "PUTA Efficiency (Small Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      scale_color_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "bottom", legend.title = element_blank())
    
    p2 <- ggplot(df, aes(x = puta_eff_large, fill = strategy, color = strategy)) +
      geom_density(alpha = 0.3) +
      labs(x = "PUTA per contact", y = "Density",
           title = paste0(title_prefix, "PUTA Efficiency (Large Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      scale_color_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "bottom", legend.title = element_blank())
    
    p3 <- ggplot(df, aes(x = pia_eff_small, fill = strategy, color = strategy)) +
      geom_density(alpha = 0.3) +
      labs(x = "PIA per contact", y = "Density",
           title = paste0(title_prefix, "PIA Efficiency (Small Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      scale_color_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "bottom", legend.title = element_blank())
    
    p4 <- ggplot(df, aes(x = pia_eff_large, fill = strategy, color = strategy)) +
      geom_density(alpha = 0.3) +
      labs(x = "PIA per contact", y = "Density",
           title = paste0(title_prefix, "PIA Efficiency (Large Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      scale_color_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "bottom", legend.title = element_blank())
    
  } else if (plot_type == "violin") {
    p1 <- ggplot(df, aes(x = strategy, y = puta_eff_small, fill = strategy)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
      labs(x = "", y = "PUTA per contact",
           title = paste0(title_prefix, "PUTA Efficiency (Small Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    p2 <- ggplot(df, aes(x = strategy, y = puta_eff_large, fill = strategy)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
      labs(x = "", y = "PUTA per contact",
           title = paste0(title_prefix, "PUTA Efficiency (Large Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    p3 <- ggplot(df, aes(x = strategy, y = pia_eff_small, fill = strategy)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
      labs(x = "", y = "PIA per contact",
           title = paste0(title_prefix, "PIA Efficiency (Small Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    p4 <- ggplot(df, aes(x = strategy, y = pia_eff_large, fill = strategy)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
      labs(x = "", y = "PIA per contact",
           title = paste0(title_prefix, "PIA Efficiency (Large Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
  } else if (plot_type == "boxplot") {
    p1 <- ggplot(df, aes(x = strategy, y = puta_eff_small, fill = strategy)) +
      geom_boxplot(alpha = 0.7) +
      labs(x = "", y = "PUTA per contact",
           title = paste0(title_prefix, "PUTA Efficiency (Small Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    p2 <- ggplot(df, aes(x = strategy, y = puta_eff_large, fill = strategy)) +
      geom_boxplot(alpha = 0.7) +
      labs(x = "", y = "PUTA per contact",
           title = paste0(title_prefix, "PUTA Efficiency (Large Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    p3 <- ggplot(df, aes(x = strategy, y = pia_eff_small, fill = strategy)) +
      geom_boxplot(alpha = 0.7) +
      labs(x = "", y = "PIA per contact",
           title = paste0(title_prefix, "PIA Efficiency (Small Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    p4 <- ggplot(df, aes(x = strategy, y = pia_eff_large, fill = strategy)) +
      geom_boxplot(alpha = 0.7) +
      labs(x = "", y = "PIA per contact",
           title = paste0(title_prefix, "PIA Efficiency (Large Subnetwork)")) +
      scale_fill_manual(values = strategy_colors) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # Combine plots
  (p1 + p2) / (p3 + p4) + 
    plot_annotation(title = "Intervention Efficiency Distributions")
}
