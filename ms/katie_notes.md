# Katie's notes on results so far and possible directions

## Key changes in the code after forking and some notes on the results

1. Maximum number of contacts available to be traced calculated over 3mo and 6mo period (as sum of F, G, and H type contacts)
2. This maximum number of contacts across the 3mo/6mo lookback window for ecah individual becomes the denominator in the total contacts traced for an index case. 
3. Because we track the number of contacts as a proxy on the amount of effort (to obtain an efficiency metric), we must make assumptions about the identity of the contacts traced. We take two exteme assumptions: that they form the smallest or the largest subnetwork possible. The variation produced across these assumptions is relatively small compared to the overall variability between simulations, with the means of the efficiency metrics only stightly shifted. 
3. Inclusion of a 'growth cluster' strategy that extends the 5 person cluster strategy to 5 person clusters that are diagnosed within a specific time interval
4. An additional strategy that combines RITA positive tracing with also tracing their contacts. Note that this not increasing the generational depth of the searching but increasing the number of transmission branches. The thinking was that the RITA strategy was very good at finding the donor compared to the other strategies (due to the infection happening within the look back window of contact tracing). Much of the infectious days averted indeed came from identification of an undiagnosed donor. Finding the donor meant a chance for 'backward contact tracing' - that is, tracing a known transmitter (ie a donor) was associated with identifying more individuals that that person infected (through the friendship paradox). However, the infectious days averted per contact DID NOT INCREASE for the infectious days averted (diganosing individuals faster) because the number of infections found quicker scaled with the number of contacts. It is likely that the IDA per contact would increase if there was heterogeneous transmission risk (so donor had higher number of transmission relative to contacts compared with rest of population). However the PIA per contact did increase as the number of infections averted increased across the new chains of transmissions that would't have been captured by RITA alone. This is the benefit of secondary contact tracing (averted future infections rather than already infected).


## TO do items

0. Check on the delays plotted for growth cluster.
1. Work out a RITA + growth cluster metric (all cluster members have at least one RITA - should you act quickly when you get one?)
1. Plot a relative effectiveness of interventions to reduce the reduce the variability by the different simulations
2. Implement a secondary contacts only (on top of random) as a fair assessment of contact tracing 
3. Implement a. differential transmission risk for individuals and b. transmission risk assortatitively, c. contact rate assortatitivy


