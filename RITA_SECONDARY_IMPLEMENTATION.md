# RITA+Secondary Strategy Implementation Plan

## Date: 2026-01-27

## Strategy Overview

**RITA+Secondary** combines RITA's early detection advantage with 2-degree contact tracing to intervene on the transmission network around recently infected individuals.

### Conceptual Framework: Transmission Tree vs Contact Networks

**Key distinction:**

1. **Transmission tree** (what we observe in simulation data):
   - Who infected whom (donor → recipient relationships)
   - Timing of transmission events
   - This is a SUBSET of actual contact relationships (most contacts don't result in transmission)

2. **Contact networks** (what we don't observe, only counts):
   - For each individual: total number of contacts (Fcontacts + Gcontacts + Hcontacts)
   - We DON'T know the identity of these contacts
   - We DON'T know if contacts overlap between individuals

**Strategy implementation:**

- **Use transmission tree to IDENTIFY network members** (who to intervene on)
- **Use contact counts to ESTIMATE notification burden** (how many contacts to notify)
- **Small/large assumptions to MODEL uncertainty** in contact overlap

**Example:**
```
Person A has 10 contacts (unknown identities)
Person B has 8 contacts (unknown identities)
A infected B (observed in transmission tree)

We intervene on both A and B.
We notify all contacts of A (10 people) + all contacts of B (8 people).

Small network assumption: 10 + 8 = 18 total notifications (no overlap)
Large network assumption: max(10, 8) = 10 total notifications (complete overlap)
Truth: Somewhere between 10 and 18 (partial overlap)
```

### Algorithm

1. **Identify RITA-positive individuals** (recently infected, diagnosed within ~6 months)
2. **Trace network within lookback window** (same as `partner_notification_window_months`):
   - **Primary traced individuals**: People identified through transmission tree as having transmitted to/from the index within the lookback window
   - **Secondary traced individuals**: People identified through transmission tree as having transmitted to/from the primary traced individuals within the lookback window
3. **Intervene on traced network**:
   - Each traced individual has an unknown set of contacts (Fcontacts + Gcontacts + Hcontacts)
   - We notify all contacts of all traced individuals
4. **Calculate total contacts notified**:
   - Small network: Assume all contacts are unique across traced individuals (sum all)
   - Large network: Assume contacts are fully overlapping (max of traced individuals)

### Key Distinction: Transmission Tree vs Contact Networks

**What we OBSERVE:**
- Transmission tree: We see who infected whom (D data)
- This allows us to identify individuals in the network

**What we DON'T KNOW:**
- Identity of contacts: We only know the COUNT of contacts per person
- Contact overlap: Whether Person A's contacts include Person B's contacts

**Strategy logic:**
- Use transmission tree to identify network members
- Use contact counts to estimate notification burden
- Small/large assumptions model uncertainty in contact overlap

### Lookback Window Logic

```
Index diagnosed at time T_dx
Lookback window: [T_dx - lookback_days, T_dx]
Intervention time (IT): T_dx + implementation_delay_days

Individual included in traced network if their transmission link to/from
a traced person occurred within lookback window.

Note: Transmission links are a SUBSET of contact relationships.
We use transmission links to identify network members, then notify
their entire contact lists.
```

### Contact Calculation

Following existing cluster-based strategy pattern:
- **Small network (contacts_small)**: Sum of all contacts (fully unique)
- **Large network (contacts_large)**: Max individual contacts (fully connected)

---

## Implementation Files to Modify

### 1. `src/intervention0.R`

**Add new function after `network_intervention()`** (around line 1450):

```r
#' RITA + Secondary Contact Tracing intervention strategy
#'
#' Combines RITA's early detection with 2-degree contact tracing.
#' Identifies RITA-positive individuals (recently infected) and traces their
#' transmission network backward and forward within the lookback window.
#'
#' Network identification (using transmission tree):
#'   - Primary traced: Individuals with transmission links to/from RITA+ index
#'   - Secondary traced: Individuals with transmission links to/from primary traced
#'   - All transmission events must occur within lookback window
#'   - Transmission links used to IDENTIFY individuals, not count contacts
#'
#' Contact notification (using contact counts):
#'   - Each traced individual has unknown contacts (Fcontacts + Gcontacts + Hcontacts)
#'   - We notify all contacts of all traced individuals
#'   - Small subnetwork: Sum all contact counts (assume no overlap)
#'   - Large subnetwork: Max contact count (assume complete overlap)
#'
#' @param Dall Combined transmission data for all simulations
#' @param Gall Combined individual data for all simulations
#' @param rita_window_months Average RITA detection window in months
#' @param implementation_delay_days Fixed delay (days) for intervention implementation
#' @param partner_notification_window_months Lookback window: 3 or 6 months (default: 6)
#' @param analysis_delay_days Delay for network analysis (default: 0, interviews during implementation)
#'
#' @return List with same structure as cluster-based strategies
rita_secondary_intervention <- function(Dall, Gall, rita_window_months,
                                        implementation_delay_days,
                                        partner_notification_window_months = 6,
                                        analysis_delay_days = 0)
{
  lastgeneration <- max(Gall$generation)
  lookback_days <- partner_notification_window_months * 30

  # Simulate RITA test: positive if diagnosed within random window of infection
  Gall$rita <- with(Gall, (timediagnosed - timeinfected) < rexp(nrow(Gall), 1/(rita_window_months*30)))

  # Filter to RITA-positive cases in generations 1+ (excluding last)
  G_rita <- Gall[Gall$rita & (Gall$generation > 0) & (Gall$generation < lastgeneration), ]

  if (nrow(G_rita) == 0) {
    return(list(
      o = data.frame(pia = numeric(0), puta = numeric(0),
                     contacts_small = numeric(0), contacts_large = numeric(0)),
      propintervened = 0,
      n_units = 0,
      puta_small = c(0, 0, 0, 0, 0),
      puta_large = c(0, 0, 0, 0, 0),
      pia_small = c(0, 0, 0, 0, 0),
      pia_large = c(0, 0, 0, 0, 0),
      total_contacts_small = 0,
      total_contacts_large = 0
    ))
  }

  # Make PIDs unique across simulations
  G_rita$pid <- paste(sep = '.', G_rita$pid, G_rita$simid)
  D <- Dall
  D$donor <- paste(sep = '.', D$donor, D$simid)
  D$recipient <- paste(sep = '.', D$recipient, D$simid)
  G <- Gall
  G$pid <- paste(sep = '.', G$pid, G$simid)

  # Calculate total contacts based on partner notification window
  if (partner_notification_window_months == 3) {
    G$total_contacts <- with(G, Fcontacts_90d + Gcontacts_90d + Hcontacts_90d)
  } else {
    G$total_contacts <- with(G, Fcontacts_180d + Gcontacts_180d + Hcontacts_180d)
  }

  # Set intervention time
  G_rita$IT <- G_rita$timediagnosed + analysis_delay_days + implementation_delay_days

  # Helper function: Get transmission-linked individuals within lookback window
  # NOTE: This identifies individuals via transmission tree, NOT their contacts
  # Contacts are unknown - we only know contact COUNTS per person
  get_transmission_linked_individuals <- function(pid, window_start, window_end) {
    # People this person transmitted TO within window
    # (These are a subset of this person's contacts)
    forward <- D$recipient[
      D$donor == pid &
      D$timetransmission >= window_start &
      D$timetransmission <= window_end
    ]

    # People this person was infected BY within window
    # (These are a subset of people who had contact with this person)
    backward <- D$donor[
      D$recipient == pid &
      D$timetransmission >= window_start &
      D$timetransmission <= window_end
    ]

    unique(c(forward, backward))
  }

  # Helper function: Check if person is eligible for intervention at IT
  eligible_for_intervention <- function(person_info, IT) {
    if (nrow(person_info) == 0) return(FALSE)
    person_info$generation[1] > 0 &&
      person_info$generation[1] < lastgeneration &&
      person_info$timediagnosed[1] > IT
  }

  # Process each RITA-positive case
  results <- list()

  for (i in seq_len(nrow(G_rita))) {
    index_pid <- G_rita$pid[i]
    index_dx_time <- G_rita$timediagnosed[i]
    IT <- G_rita$IT[i]

    # Define lookback window anchored to index diagnosis
    lookback_start <- index_dx_time - lookback_days
    lookback_end <- index_dx_time

    # Build traced network using transmission tree
    # NOTE: We use transmission links to IDENTIFY individuals, not to count contacts
    network <- c(index_pid)  # Start with RITA+ index

    # Step 1: Identify primary traced individuals (transmission-linked to index)
    primary_traced <- get_transmission_linked_individuals(index_pid, lookback_start, lookback_end)

    for (primary_pid in primary_traced) {
      primary_info <- G[G$pid == primary_pid, ]
      if (eligible_for_intervention(primary_info, IT)) {
        network <- c(network, primary_pid)

        # Step 2: Identify secondary traced individuals (transmission-linked to primary)
        secondary_traced <- get_transmission_linked_individuals(primary_pid, lookback_start, lookback_end)

        for (secondary_pid in secondary_traced) {
          # Skip if already in network
          if (secondary_pid %in% network) next

          secondary_info <- G[G$pid == secondary_pid, ]
          if (eligible_for_intervention(secondary_info, IT)) {
            network <- c(network, secondary_pid)
          }
        }
      }
    }

    network <- unique(network)

    # Calculate total contacts to notify
    # Each traced individual has their own (unknown) contact list
    # We notify ALL contacts of ALL traced individuals
    network_contacts <- numeric(length(network))
    for (j in seq_along(network)) {
      person <- G[G$pid == network[j], ]
      if (nrow(person) > 0) {
        # total_contacts = Fcontacts_XXd + Gcontacts_XXd + Hcontacts_XXd
        network_contacts[j] <- person$total_contacts[1]
      }
    }

    # Small network: sum all contacts (assume fully unique - no overlap)
    # Total notifications = sum of all contact lists
    contacts_small <- sum(network_contacts)

    # Large network: max contacts (assume fully overlapping/connected)
    # Total unique contacts ≈ max individual contact list
    contacts_large <- max(network_contacts, 0)

    # Calculate PIA and IDA for network
    # PIA/IDA scope: network members + their transmission partners
    # (Following cluster-based strategy pattern)
    piapids <- network
    for (pid in network) {
      # Add all transmission partners of network members
      # (These represent potential infections averted through network intervention)
      piapids <- union(piapids, D$recipient[D$donor == pid])  # People they infected
      piapids <- union(piapids, D$donor[D$recipient == pid])  # People who infected them
    }

    G2 <- G[G$pid %in% piapids, ]

    # PIA: Infections that occurred AFTER intervention time
    pia <- sum(G2$timeinfected > IT)

    # IDA: Person-years of untreated infection averted
    # (People infected before IT but diagnosed after IT)
    G3 <- G2[G2$timeinfected <= IT & G2$timediagnosed > IT, ]
    puta <- sum(G3$timediagnosed - IT)

    results[[i]] <- data.frame(
      pia = pia,
      puta = puta,
      contacts_small = contacts_small,
      contacts_large = contacts_large
    )
  }

  # Combine results
  odf <- do.call(rbind, results)

  # Compute summary statistics (matching cluster-based strategies)
  # Quantiles: 10th and 90th percentiles used uniformly for all metrics

  # IDA efficiency for small subnetwork
  e_puta_small <- odf$puta / odf$contacts_small
  e_puta_small_valid <- e_puta_small[is.finite(e_puta_small) & !is.na(e_puta_small)]
  med_puta_small <- if (length(e_puta_small_valid) > 0) median(e_puta_small_valid) else NA_real_
  q_puta_small <- if (length(e_puta_small_valid) > 0) quantile(e_puta_small_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # IDA efficiency for large subnetwork
  e_puta_large <- odf$puta / odf$contacts_large
  e_puta_large_valid <- e_puta_large[is.finite(e_puta_large) & !is.na(e_puta_large)]
  med_puta_large <- if (length(e_puta_large_valid) > 0) median(e_puta_large_valid) else NA_real_
  q_puta_large <- if (length(e_puta_large_valid) > 0) quantile(e_puta_large_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # PIA efficiency for small subnetwork
  e_pia_small <- odf$pia / odf$contacts_small
  e_pia_small_valid <- e_pia_small[is.finite(e_pia_small) & !is.na(e_pia_small)]
  med_pia_small <- if (length(e_pia_small_valid) > 0) median(e_pia_small_valid) else NA_real_
  q_pia_small <- if (length(e_pia_small_valid) > 0) quantile(e_pia_small_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # PIA efficiency for large subnetwork
  e_pia_large <- odf$pia / odf$contacts_large
  e_pia_large_valid <- e_pia_large[is.finite(e_pia_large) & !is.na(e_pia_large)]
  med_pia_large <- if (length(e_pia_large_valid) > 0) median(e_pia_large_valid) else NA_real_
  q_pia_large <- if (length(e_pia_large_valid) > 0) quantile(e_pia_large_valid, probs = c(0.1, 0.9), names = FALSE) else c(NA_real_, NA_real_)

  # Return results
  list(
    o = odf,
    propintervened = nrow(odf) / nrow(Gall[Gall$generation > 0 & Gall$generation < lastgeneration, ]),
    n_units = nrow(odf),
    puta_small = c(sum(odf$puta), med_puta_small, q_puta_small),
    puta_large = c(sum(odf$puta), med_puta_large, q_puta_large),
    pia_small = c(sum(odf$pia), med_pia_small, q_pia_small),
    pia_large = c(sum(odf$pia), med_pia_large, q_pia_large),
    total_contacts_small = sum(odf$contacts_small),
    total_contacts_large = sum(odf$contacts_large)
  )
}
```

**Call in `run_intervention_analysis()`** (around line 450, after network intervention):

```r
# RITA+Secondary intervention
cat("  Running RITA+Secondary intervention...\n")
oritasec <- tryCatch({
  rita_secondary_intervention(
    Dall, Gall,
    rita_window_months = rita_window_months,
    implementation_delay_days = implementation_delay_days,
    partner_notification_window_months = partner_notification_window_months,
    analysis_delay_days = 0  # No extra delay for contact tracing
  )
}, error = function(e) {
  cat("  ERROR in RITA+Secondary intervention:", e$message, "\n")
  list(o = data.frame(), propintervened = 0, n_units = 0,
       puta_small = rep(0, 5), puta_large = rep(0, 5),
       pia_small = rep(0, 5), pia_large = rep(0, 5),
       total_contacts_small = 0, total_contacts_large = 0)
})
```

---

### 2. `src/intervention0.R` - Update Summary Table

Add RITA+Secondary row to summary output (around line 460):

```r
odf <- data.frame(
  strategy = c("distsize5", "distsize2", "growth", "random", "rita", "network", "ritasecondary"),
  # ... rest of columns
)
```

---

### 3. `src/intervention0.R` - Update File Saving

Save RITA+Secondary details (around line 490):

```r
if (!is.null(output_dir) && nrow(oritasec$o) > 0) {
  write.csv(oritasec$o, file.path(output_dir, paste0("details_ritasecondary_", timestamp, ".csv")), row.names = FALSE)
}
```

---

### 4. `src/plot_interventions.R` - Add Strategy

Update strategy lists (around line 39):

```r
strategy_names <- c("distsize5", "distsize2", "growth", "random", "rita", "network", "ritasecondary")
strategy_labels <- c("Size>=5", "Size>=2", "Growth", "Random", "RITA", "Network", "RITA+Secondary")
```

Add color for RITA+Secondary (around line 100):

```r
strategy_colors <- c(
  "Size>=5" = "#E41A1C",
  "Size>=2" = "#377EB8",
  "Growth" = "#4DAF4A",
  "Random" = "#984EA3",
  "RITA" = "#FF7F00",
  "Network" = "#A65628",
  "RITA+Secondary" = "#F781BF"  # Pink/magenta
)
```

---

### 5. `src/run_analysis.R` - Update Loading

Add RITA+Secondary to cached results loading (around line 79):

```r
strategy_files <- list(
  distsize5 = find_most_recent("^details_distsize5_"),
  distsize2 = find_most_recent("^details_distsize2_"),
  growth = find_most_recent("^details_growth_"),
  random = find_most_recent("^details_random_"),
  rita = find_most_recent("^details_rita_"),
  network = find_most_recent("^details_network_"),
  ritasecondary = find_most_recent("^details_ritasecondary_")
)
```

---

## Testing Plan

### 1. Run Fresh Analysis
```r
source("src/run_analysis.R")
run_full_analysis(use_cached = FALSE, n_sims = 100)
```

Expected output:
- New `details_ritasecondary_TIMESTAMP.csv` file
- RITA+Secondary row in summary table
- 7 strategies in efficiency plots

### 2. Verify Plots
```r
run_full_analysis(use_cached = TRUE, n_sims = NULL)
```

Expected:
- 7 strategies in violin plots
- RITA+Secondary in mechanism analysis (if applicable)

### 3. Compare Results

RITA+Secondary should show:
- **Higher mean transmissions per target** than RITA alone (Panel C)
- **Lower efficiency** than RITA alone (more contacts)
- **Competitive total impact** (catches network around early infections)

---

## Expected Results

Based on the strategy design:

| Metric | RITA | RITA+Secondary | Expectation |
|--------|------|----------------|-------------|
| Mean offspring (Panel C) | 0.29 | ~0.8-1.2 | Higher (includes donor + siblings) |
| % Undiagnosed donor (Panel A) | 84% | ~70-80% | Slightly lower (diluted by secondary contacts) |
| % Preventable (Panel B) | 38% | ~20-30% | Lower (mixed early + later cases) |
| IDA efficiency | High | Medium | Lower (more contacts per unit) |
| Total IDA | Medium | High | Higher (more people intervened on) |

---

## Files to Create/Modify Summary

1. ✅ `RITA_SECONDARY_IMPLEMENTATION.md` (this file)
2. ⏳ `src/intervention0.R` - Add function + call
3. ⏳ `src/plot_interventions.R` - Add strategy name/color
4. ⏳ `src/run_analysis.R` - Add to cache loading

Total new lines of code: ~200 lines
