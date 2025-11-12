invisible('
experiment1 intervention model 
')

library( glue )
library( ggplot2 )

    
    

# Main function to run intervention analysis
run_intervention_analysis <- function(
    d_file = 'experiment1-D.csv',
    g_file = 'experiment1-G.csv',
    seed = NULL,  # Default to random seed
    cluster_size_5 = 5,
    cluster_size_2 = 2,
    subnetworksize = "small", # large of small network (sparse or dense connections inside cluster)
    distance_threshold = 0.005,
    network_degree_threshold = 4,
    random_sample_size = 30,
    rita_window_months = 6,
    intervention_rate = 1/90,
    show_table = TRUE) 
{
 
  
  
  # Handle seed setting
  if (is.null(seed)) {
    # Use system time for random seed
    seed <- as.numeric(Sys.time())
    cat("Using random seed:", seed, "\n")
  } else {
    cat("Using fixed seed:", seed, "\n")
  }
  # set my seed
  set.seed(seed)

  # Load data
  cat("Loading data...\n")
  Dall = read.csv( d_file , stringsAs=FALSE)
  Gall = read.csv( g_file, stringsAs=FALSE)
  
  # Split data by simulation ID
  simids = unique( Dall$simid )
  Ds = split( Dall, Dall$simid )[simids]
  Gs = split( Gall, Gall$simid )[ simids]
  
  # Run all interventions with error handling
  cat("Running interventions...\n")
  
  ods5 <- tryCatch(distsize_intervention(thsize = cluster_size_5, subnetwork = subnetworksize), error = function(e) {
    cat("Error in ods5:", e$message, "\n")
    list(propintervened = 0, puta = c(0, 0, 0, 0), pia = c(0, 0, 0), total_contacts = 0)
  })
  
  ods2 <- tryCatch(distsize_intervention(thsize = cluster_size_2, subnetwork = subnetworksize), error = function(e) {
    cat("Error in ods2:", e$message, "\n")
    list(propintervened = 0, puta = c(0, 0, 0, 0), pia = c(0, 0, 0), total_contacts = 0)
  })
  
  orand <- tryCatch(random_intervention(), error = function(e) {
    cat("Error in orand:", e$message, "\n")
    list(propintervened = 0, puta = c(0, 0, 0, 0), pia = c(0, 0, 0), total_contacts = 0)
  })
  
  orita <- tryCatch(rita_intervention(), error = function(e) {
    cat("Error in orita:", e$message, "\n")
    list(propintervened = 0, puta = c(0, 0, 0, 0), pia = c(0, 0, 0), total_contacts = 0)
  })
  
  onet <- tryCatch(network_intervention(), error = function(e) {
    cat("Error in onet:", e$message, "\n")
    list(propintervened = 0, puta = c(0, 0, 0, 0), pia = c(0, 0, 0), total_contacts = 0)
  })
  
  # Compile results
  odf <- rbind(
    with(ods5, c(total_contacts, puta, pia)),
    with(ods2, c(total_contacts, puta, pia)),
    with(orand, c(total_contacts, puta, pia)),
    with(orita, c(total_contacts, puta, pia)),
    with(onet, c(total_contacts, puta, pia))
  ) |> as.data.frame()
  
  colnames(odf) <- c("Contacted Total", "Total PUTA", "PUTA/contacted", "Low", "High", "PIA", "Low", "High")
  
  rownames(odf) <- c(
    paste0('Size=', cluster_size_5, ',D=', distance_threshold),
    paste0('Size=', cluster_size_2, ',D=', distance_threshold),
    'Random allocation',
    'RITA',
    paste0('Network, partners>', network_degree_threshold)
  )
  
  odf1 <- round(odf, 2)
  
  # Display results
  if (show_table) {
    if (require(knitr, quietly = TRUE)) {
      print(knitr::kable(odf1))
    } else {
      cat("knitr package not available, showing table as dataframe\n")
      print(odf1)
    }
  }
}
 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
###### ### ### ### ###  INTERVENTION FUNCTIONS ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


  # Process individual cluster for intervention analysis
  proc_cluster <- function(D, G, thdist = distance_threshold, thsize = cluster_size_5, 
                           thgrowth = NA, subnetwork = "small", ritdist = function() rexp(1, rate = intervention_rate)) 
  {
    
    # Validate input
  	stopifnot( is.na( thgrowth )) # growth threshold not implemented
  	
    # Sort data by relevant time variables
    G = G[ order(G$timesequenced), ]
  	D = D[ order(D$timetransmission),]
  
  	lastgeneration <- max( G$generation )
  
  
  	# Retain cases with path to patient 0, within distance threshold, excluding last generation
  	D1 <- D[ D$distance <= thdist , ]
  	keeppids = "0" 
  	addpids = D1$recipient[ D1$donor %in% keeppids ]
  	
  	while( length( addpids ) > 0 ){
  		keeppids <- union( addpids, keeppids )
  		addpids = setdiff( D1$recipient[ D1$donor %in% keeppids ], keeppids )
  	}
  	
  	G1 <- G[ G$pid %in% keeppids , ]
  	G1 <- G1[ G1$generation != lastgeneration, ]
  	D1 <- D1[ D1$donor %in% G1$pid & D1$recipient %in% G1$pid, ]
  	
  	# find intervention time, if there is one 
  	csize <- 0 
  	IT <- Inf 
  	
  	if (nrow(G1) > 0) {
  	  for (i in 1:nrow(G1)) {
  	    csize <- csize + 1
  	    if (csize >= thsize) {
  	      IT <- G1$timesequenced[i] + ritdist()
  	      break
  	    }
  	  }
  	}
  	
  	# Exclude clustered infections detected after intervention time
  	G1 <- G1[ G1$timesequenced < IT , ]
  	
  	# Calculate total cluster contact network size
  	if (nrow(G1) > 0) {
  	  # Get cluster members' total degrees
  	  G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
  	  total_degree <- sum(G1$degree)
  	  
  	  # Estimate internal connections within cluster
  	  
  	  if(subnetwork == "large"){ # sparse connections inside cluster
  	    total_contacts <- total_degree - (nrow(G1) - 2)
  	  }else if(subnetwork == "small"){# dense connections inside cluster
  	    total_contacts <- nrow(G1) + sum(pmax(G1$degree - (nrow(G1)-1), 0))
  	  }else{
  	    error("Assumption about how structure of contacts in cluster needs to be defined (large or small subnetwork)")
  	  }
  	  
  	} else {
  	  total_contacts <- 0
  	}
  	
  	# Calculate potential infections averted (PIA) and person-years under treatment averted (PUTA)
  	
  	pia = 0 
  	puta = 0 
  	piapids <- c() 
  	
  	# calculate transmissions from and to D (this is different to original code where final union was omitted)
  	if (length(IT) > 0 && !is.infinite( IT ) )
  	{
  		piapids <-  D$recipient[D$donor %in% G1$pid] |> 
  		      union(G1$pid) |> 
  		      union(D$donor[D$recipient %in% G1$pid]) 
  		
  		
  		G2 <- G[ G$pid %in% piapids , ]
  		pia <- sum( G2$timeinfected  > IT )
  
  		G3 <- G2[ G2$timeinfected <= IT & G2$timediagnosed > IT , ]
  		puta <- sum( G3$timediagnosed - IT )
  	}
  	c(pia = pia, puta = puta, interventiontime = IT, nc = nrow(G1), total_contacts = total_contacts)
  }
  
  ### ### ### # Distance-size intervention strategy
  
  distsize_intervention <- function( thsize = cluster_size_5, subnetwork = "small")
  {

    # Process all simulations
    o <- lapply(1:length(Ds), function(i) {
      tryCatch({
        proc_cluster(Ds[[i]], Gs[[i]], thsize = thsize, subnetwork = subnetwork)
      }, error = function(e) {
        # Return default values if processing fails
        c(pia = 0, puta = 0, interventiontime = Inf, nc = 0, total_contacts = 0)
      })
    })
    
  	odf = do.call( rbind, o ) |> as.data.frame()
  	odf1 <- odf[ !is.infinite( odf$interventiontime ), ] # exclude sims where no cluster found 
  
  	
  	lastgen <- max( Gall$generation )
  
  	# Filter out cases where total_contacts is 0 to avoid division by zero
  	odf1 <- odf1[odf1$total_contacts > 0, ]
  	
  	if (nrow(odf1) == 0) {
  	  return(list(
  	    o = data.frame(pia = numeric(0), puta = numeric(0), interventiontime = numeric(0), 
  	                   nc = numeric(0), total_contacts = numeric(0)),
  	    propintervened = 0,
  	    puta = c(0, 0, 0, 0),
  	    pia = c(0, 0, 0),
  	    total_contacts = 0
  	  ))
  	}
  	
  	sort_puta_percontact <- sort(odf1$puta / odf1$total_contacts)
  	sort_pia <- sort(odf1$pia)
  	
  	# output IQ ranges and pia/puta as per contacted individuals in addition to totals
  	list(
  	  o = odf1,
  	  propintervened = sum(odf1$nc) / sum(Gall$generation > 0 & Gall$generation < lastgen),
  	  puta = c(sum(odf1$puta),
  	           sum(odf1$puta) / sum(odf1$total_contacts), 
  	           sort_puta_percontact[max(1, ceiling(0.1 * nrow(odf1)))], 
  	           sort_puta_percontact[min(nrow(odf1), floor(0.9 * nrow(odf1)))]),
  	  pia = c(mean(odf1$pia), 
  	          sort_pia[max(1, ceiling(0.1 * nrow(odf1)))], 
  	          sort_pia[min(nrow(odf1), floor(0.9 * nrow(odf1)))]),
  	  total_contacts = sum(odf1$total_contacts)
  	)
  }
  
  invisible( '
  	select random case from generation after seed and before last generation 
  	sim intervention time 
  	trace one generation + donor 
  ')
  random_intervention <- function(
  				# ritdist = function(n) rexp(n,rate=1/90) 
  )
  {
  	lastgeneration <- max( Gall$generation )
  	G1 <- Gall[ Gall$generation > 0 & Gall$generation < lastgeneration, ] # gen 2+ 
  	
  	if (nrow(G1) == 0) {
  	  return(list(
  	    o = data.frame(pia = numeric(0), puta = numeric(0)),
  	    propintervened = 0,
  	    puta = c(0, 0, 0, 0),
  	    pia = c(0, 0, 0),
  	    total_contacts = 0
  	  ))
  	}
  	
  	# Randomly sample cases
  	sample_size <- min(random_sample_size, nrow(G1))
  	G1 <- G1 |> dplyr::slice_sample(n = sample_size)
  	
  	# Compute weighted degree before making PIDs unique
  	G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
  	
  	
  	# make pids unique 
  	G1$pid <- paste(sep='.', G1$pid, G1$simid )
  	D <- Dall; D$donor <- paste(sep='.', D$donor, D$simid )
  	D$recipient <- paste(sep='.', D$recipient, D$simid )
  	G <- Gall; G$pid <- paste(sep='.', G$pid, G$simid)
  	
  	# Set intervention time
  	# G1$IT <-  G1$timesequenced + ritdist( nrow(G1 )) 
  	G1$IT <- G1$timesequenced + rexp(nrow(G1), rate = intervention_rate)
  	
  	# Process individual intervention
  	proc_indiv <- function(pid, IT, degree )
  	{
  		piapids <- D$recipient[ D$donor == pid] |> 
  		        union( c( pid, D$donor[D$recipient==pid] ))	
  		
  		G2 <- G[ G$pid %in% piapids , ]
  		pia <- sum( G2$timeinfected > IT )
  		
  		G3 <- G2[ G2$timeinfected <= IT & G2$timediagnosed > IT, ]
  		puta <- sum( G3$timediagnosed - IT )
  		
  		c( pia, puta,degree + 1 )
  	}
  	mapply(proc_indiv, G1$pid, G1$IT, G1$degree ) -> o 
  	o <- as.data.frame( t( o ) )
  	colnames(o) <- c('pia', 'puta', 'contacts' )
  
  	
  	sort_puta_percontact <- sort(o$puta / o$contacts)
  	sort_pia <- sort(o$pia)
  	
  	
  	list(
  	  o = o,
  	  propintervened = NA,
  	  puta = c(sum(o$puta),
  	           sum(o$puta) / sum(o$contacts), 
  	           sort_puta_percontact[max(1, ceiling(0.1 * nrow(G1)))], 
  	           sort_puta_percontact[min(nrow(G1), floor(0.9 * nrow(G1)))]),
  	  pia = c(mean(o$pia), 
  	          sort_pia[max(1, ceiling(0.01 * nrow(G1)))], 
  	          sort_pia[min(nrow(G1), floor(0.99 * nrow(G1)))]),
  	  total_contacts = sum(o$contacts)
  	)
  }
  
  invisible('
  Intervene on cases with acute infection 
  Simulated RITA test
  6 month average detection window 
  	  ')
  rita_intervention <- function(
  				# ritdist = function(n) rexp(n,rate=1/rita_window_months) 
  			       )
  {
  	lastgeneration <- max( Gall$generation )
  	
  	# RITA test simulation - 6 month average detection window
  	Gall$rita <- with(Gall, (timediagnosed - timeinfected) < rexp(nrow(Gall), 1/(rita_window_months*30)))
  	
  	# Filter generations & only RITA positive cases
  	G1 <- Gall[Gall$rita & (Gall$generation > 0) & (Gall$generation < lastgeneration), ]
  	
  	if (nrow(G1) == 0) {
  	  return(list(
  	    o = data.frame(pia = numeric(0), puta = numeric(0), contacts = numeric(0)),
  	    propintervened = 0,
  	    puta = c(0, 0, 0, 0),
  	    pia = c(0, 0, 0),
  	    total_contacts = 0
  	  ))
  	}
  	
  	# Compute weighted degree for RITA-positive individuals
  	G1$degree <- with(G1, Fdegree + Gdegree + Hdegree)
  	
  	
  
  	# make pids unique 
  	G1$pid <- paste(sep = '.', G1$pid, G1$simid)
  	D <- Dall; D$donor <- paste(sep='.', D$donor, D$simid )
  	D$recipient <- paste(sep='.', D$recipient, D$simid )
  	G <- Gall; G$pid <- paste(sep='.', G$pid, G$simid)
  	
  	
  	
  	G1$IT <- G1$timesequenced + rexp(nrow(G1), rate = intervention_rate)
  	
  	# Process individual intervention
  	proc_indiv <- function(pid, IT, degree) {
  	  piapids <- D$recipient[D$donor == pid] |> 
  	    union(c(pid, D$donor[D$recipient == pid]))
  	  
  	  G2 <- G[G$pid %in% piapids, ]
  	  pia <- sum(G2$timeinfected > IT)
  	  
  	  G3 <- G2[G2$timeinfected <= IT & G2$timediagnosed > IT, ]
  	  puta <- sum(G3$timediagnosed - IT)
  	  
  	  c(pia, puta, degree + 1)
  	}
  	
  	# Process all cases
  	results <- matrix(nrow = nrow(G1), ncol = 3)
  	for (i in 1:nrow(G1)) {
  	  results[i, ] <- proc_indiv(G1$pid[i], G1$IT[i], G1$degree[i])
  	}
  	
  	o <- as.data.frame(results)
  	colnames(o) <- c('pia', 'puta', 'contacts')
  	
  	# Need to make the pia and puta calculations consistent (either efficiency or total)
  	sort_puta_percontact <- sort(o$puta / o$contacts)
  	sort_pia <- sort(o$pia)
  	
  	list(
  	  o = o,
  	  propintervened = NA,
  	  puta = c(sum(o$puta),
  	           sum(o$puta) / sum(o$contacts), 
  	           sort_puta_percontact[max(1, ceiling(0.1 * nrow(o)))], 
  	           sort_puta_percontact[min(nrow(o), floor(0.9 * nrow(o)))]),
  	  pia = c(mean(o$pia), 
  	          sort_pia[max(1, ceiling(0.01 * nrow(o)))], 
  	          sort_pia[min(nrow(o), floor(0.99 * nrow(o)))]),
  	  total_contacts = sum(o$contacts)
  	)
  }
  
  
  network_intervention <- function( 
    # network_degree_threshold = 30, ritdist = function(n) rexp(n,rate=1/90) 
  )
  {
  	# do not use weighting and instead just use contact numbers if 
    # w <- 7.0 / 2.0 # ratio duration F to G 
  	# ww <-  7.0*30.0 # if unit daiy oo contact rate, E partners in 7 months (F duration) 
  
    lastgeneration <- max(Gall$generation)
    
    # Make pids unique across simulations
    D <- Dall
    D$donor <- paste(sep = '.', D$donor, D$simid)
    D$recipient <- paste(sep = '.', D$recipient, D$simid)
    G <- Gall
    G$pid <- paste(sep = '.', G$pid, G$simid)
  	
  	# compute weighted degree ( exp partners over 7 months )
  	# G$degree  <-  with( G, Fdegree + w*Gdegree + ww*Hdegree)
  	G$degree <- with(G, Fdegree + Gdegree + Hdegree)
  
  	# filter generations & only degree above threshold  
  	G1 <- G[ (G$degree>=network_degree_threshold) & (G$generation > 0) & (G$generation < lastgeneration), ] # gen 2+ 
  	
  	if (nrow(G1) == 0) {
  	  return(list(
  	    o = data.frame(pia = numeric(0), puta = numeric(0), contacts = numeric(0)),
  	    propintervened = 0,
  	    puta = c(0, 0, 0, 0),
  	    pia = c(0, 0, 0),
  	    total_contacts = 0
  	  ))
  	}
  	
  	
  	# intervention time 
  	# G1$IT <-  G1$timesequenced + ritdist( nrow(G1 )) 
  	# Set intervention time
  	G1$IT <- G1$timesequenced + rexp(nrow(G1), rate = intervention_rate)
  	
  	# Process individual intervention
  	
  	proc_indiv <- function(pid, IT, degree )
  	{
  	  piapids <- D$recipient[D$donor == pid] |> 
  	    union(c(pid, D$donor[D$recipient == pid]))
  	  
  	  G2 <- G[G$pid %in% piapids, ]
  	  pia <- sum(G2$timeinfected > IT)
  	  
  	  G3 <- G2[G2$timeinfected <= IT & G2$timediagnosed > IT, ]
  	  puta <- sum(G3$timediagnosed - IT)
  	  
  	  c(pia, puta, degree + 1)
  	}
  	mapply(proc_indiv, G1$pid, G1$IT, G1$degree ) -> o 
  	o <- as.data.frame( t( o ) )
  	colnames(o) <- c('pia', 'puta','contacts' )
  	sort_puta_percontact <- sort(o$puta / o$contacts)
  	sort_pia <- sort(o$pia / o$contacts)
  
  	list(
  	  o = o,
  	  propintervened = NA,
  	  puta = c(sum(o$puta),
  	           sum(o$puta) / sum(o$contacts), 
  	           sort_puta_percontact[max(1, ceiling(0.1 * nrow(o)))], 
  	           sort_puta_percontact[min(nrow(o), floor(0.9 * nrow(o)))]),
  	  pia = c(mean(o$pia), 
  	          sort_pia[max(1, ceiling(0.01 * nrow(o)))], 
  	          sort_pia[min(nrow(o), floor(0.99 * nrow(o)))]),
  	  total_contacts = sum(o$contacts)
  	)
  }
  
  
  
  