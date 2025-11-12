# Intervention model – functional changes (new → old)

1. Overall workflow

- Single entry point: New code wraps everything in run_intervention_analysis(...) (loads data, runs all strategies, prints a table). Old code runs top-level script blocks and prints a different table at the end.

- Error handling: New code uses tryCatch around each strategy and returns zeroed summaries on failure; old version would error/stop.

2. Randomness & seeds

- Seed control: New accepts seed (default = current system time; printed) and calls set.seed(seed). Old has no seed control beyond using rexp(...) via ritdist.

3. Inputs & parameterization

- File I/O: New loads inside the main function via d_file, g_file args. Old loads at top level.

- Renamed/normalized parameters:
    - thdist → distance_threshold
    - thsize → cluster_size
    - thdegree → network_degree_threshold
    - ritdist (function) → intervention_rate (scalar rate passed to rexp)
    - New options: subnetwork = "small"|"large" (dense vs. sparse within-cluster contacts), random_sample_size, rita_window_months.

4. Cluster definition & intervention time (distance–size strategy)

Same concept (grow from index “0” along edges ≤ threshold, drop last generation), but:

Intervention time:

New: if nrow(G1) ≥ cluster_size, set IT = G1$timesequenced[cluster_size] + rexp(1, intervention_rate).

Old: loop until size reaches thsize, then timesequenced[i] + ritdist().
→ Functionally equivalent trigger; new version is simpler and explicit.

Post-IT filtering: Both exclude cases with timesequenced ≥ IT.

5. Who counts toward PIA/PUTA (key behavioral change)

New: For cluster interventions, piapids includes:

recipients of cluster donors ∪ all cluster members ∪ donors of cluster recipients (one “generation” outward in both directions).

Old: piapids included only recipients of cluster donors ∪ cluster members (the donor-of-recipient link was omitted).
→ New version expands the affected set; PIA/PUTA can be larger.

6. Contact modeling (new behavior)

New computes contacts for strategies that need them.

In proc_cluster, each member’s degree = Fdegree + Gdegree + Hdegree.

Dense (“small”) subnetwork: total_contacts = nrow(G1) + sum(pmax(degree - (nrow(G1)-1), 0)).

Sparse (“large”) subnetwork: total_contacts = total_degree - (nrow(G1)-2).

Old had no contact accounting; most outputs were per cluster member (/ nc) rather than per contact.

7. Strategy-specific changes

A.  Distance–size (distsize_intervention)

Output scale:

New: reports totals and per-contact efficiency; returns quantiles.

Old: reported means/variances per member (/nc) and computed Wald CIs.

Intervened proportion: Both compute sum(nc) / sum(generation>0 & < last).

Zero-contact guard: New drops rows with total_contacts == 0 and returns all-zeros if nothing remains; old had no such guard.

B. Random selection (random_intervention)

Sampling:

New: samples up to random_sample_size (default 30) from eligible cases (gens 1..last-1).

Old: used all eligible cases.

Contacts & degree:

New: computes degree = F+G+H, sets contacts = degree + 1 per individual, and aggregates total PUTA, PUTA/contact, and quantiles; PIA quantiles at 1%/99%.

Old: no contacts; reported mean/variance of PIA/PUTA per person.

Prop intervened: Both mark as not meaningful (NA).

C. RITA-based (rita_intervention)

Window:

New: rita_window_months parameter (default 6) → rexp(1/(months*30)).

Old: fixed 6 months.

Contacts & outputs:

New: adds degree-based contacts, reports total PUTA, per-contact PUTA, and quantiles; returns sum of contacts.

Old: per-person mean/variance.

Prop intervened:

New: NA.

Old: nrow(G1) / nrow(Gall).

D. Network-degree targeting (network_intervention)

Degree weighting:

New: unweighted degree = F + G + H.

Old: weighted (F + (7/2)*G + (7*30)*H) to approximate partners over 7 months.
→ New materially changes who qualifies at a given threshold and contact counts.

Threshold name: network_degree_threshold (new) vs thdegree (old).

Contacts & outputs:

New: per-contact PUTA efficiency + totals; PIA list mixes mean of raw pia with per-contact quantiles (inconsistency noted in a code comment).

Old: per-person mean/variance and propintervened = nrow(G1)/nrow(Gall).

8. Aggregation & summary table

New table columns:
"Contacted Total", "Total PUTA", "PUTA/contacted", "Low", "High", "PIA", "Low", "High"
Rows: size/threshold labels plus ‘Random’, ‘RITA’, ‘Network, partners>threshold’.
Values are rounded and printed via knitr::kable if available (else print).

Old table columns:
"Proportion intervened", "Mean", "Variance", "lb", "ub" — a CI around mean PUTA (per member), not contacts-based.

9. Edge cases & defaults

Empty selections: New functions consistently return empty frames plus zeroed summaries (and total_contacts = 0) if nothing qualifies; old functions returned means/variances that could be undefined for empty sets.

Quantiles vs. variances: New reports empirical quantiles (10%/90% for PUTA/contact; 10%/90% or 1%/99% for PIA depending on strategy) instead of normal-approximation CIs.

10. Backward-incompatible behavior (pay attention)

Metrics now emphasize per-contact efficiency and totals, not per-member means/variances.

PIA/PUTA population widened in distance–size strategy by including donors of cluster recipients (one extra edge) → larger effects possible.

Network targeting uses unweighted degrees (F, G, H equally), changing who is targeted at the same numeric threshold.

Random strategy samples a fixed number (default 30) rather than the full eligible set.

RITA “proportion intervened” is no longer reported (now NA).

Results table shape & semantics changed (downstream code that reads old columns will break).