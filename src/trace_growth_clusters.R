# Trace growth-rate cluster formation to verify trigger vs intervention time

source('src/intervention0.R')

# Load data
Dall <- read.csv('src/experiment1-N10000-gens7-D.csv')
Gall <- read.csv('src/experiment1-N10000-gens7-G.csv')

simids <- as.character(unique(Dall$simid))
Ds <- split(Dall, Dall$simid)[simids]
Gs <- split(Gall, Gall$simid)[simids]

# Parameters
growth_distance_threshold <- 0.1
cluster_size <- 5
lookback_days <- 6 * 30  # 6 months
intervention_rate <- 1/90

set.seed(123)

# Find a few growth-triggered clusters and trace them
cat('=== Tracing Growth-Rate Cluster Formation ===\n\n')

n_examples <- 0
max_examples <- 5

for (sim_key in simids) {
  if (n_examples >= max_examples) break
  
  D <- Ds[[sim_key]]
  G <- Gs[[sim_key]]
  
  # Build cluster
  G <- G[order(G$timesequenced), ]
  lastgen <- max(G$generation)
  
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
  
  # Check sliding window
  t <- as.numeric(Gtrig$timesequenced)
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
  
  if (!is.finite(t_detect)) next
  
  # Found one!
  n_examples <- n_examples + 1
  IT <- t_detect + rexp(1, rate = intervention_rate)
  
  # Count at different times
  n_at_trigger <- sum(Gcluster$generation != lastgen & Gcluster$timesequenced <= t_detect)
  n_at_IT <- sum(Gcluster$generation != lastgen & Gcluster$timesequenced < IT)
  n_total_cluster <- sum(Gcluster$generation != lastgen)
  
  cat('--- Simulation', sim_key, '---\n')
  cat('Trigger time (t_detect):', round(t_detect, 1), 'days\n')
  cat('Intervention time (IT):', round(IT, 1), 'days\n')
  cat('Delay (IT - t_detect):', round(IT - t_detect, 1), 'days\n')
  cat('Cases sequenced at trigger:', n_at_trigger, '\n')
  cat('Cases sequenced at IT (nc):', n_at_IT, '\n')
  cat('Total cluster size (all time):', n_total_cluster, '\n')
  
  # Show the sequencing times in the trigger window
  window_cases <- Gtrig[j:i, ]
  cat('Trigger window cases (', nrow(window_cases), ' cases in ', 
      round(t[i] - t[j], 1), ' days):\n', sep='')
  for (k in 1:nrow(window_cases)) {
    cat('  pid', window_cases$pid[k], ': sequenced day', round(window_cases$timesequenced[k], 1), '\n')
  }
  cat('\n')
}
