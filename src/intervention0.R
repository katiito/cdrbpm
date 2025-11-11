invisible('
experiment1 intervention model 
')

library( glue )
library( ggplot2 )



  Dall = read.csv( 'experiment1-D.csv' , stringsAs=FALSE)
  Gall = read.csv( 'experiment1-G.csv' , stringsAs=FALSE)
  
  simids = unique( Dall$simid )
  Ds = split( Dall, Dall$simid )[simids]
  Gs = split( Gall, Gall$simid )[ simids]
  
  proc_cluster <- function( D, G, thdist = 0.005, thsize = 5, thgrowth=NA
   			  , ritdist = function() rexp(1,rate=1/90) )
  {
  
  	stopifnot( is.na( thgrowth )) # not impl 
  	G = G[ order(G$timesequenced), ]
  	D = D[ order(D$timetransmission),]
  
  	lastgeneration <- max( G$generation )
  
  
  	# retain cases with a path to pid 0, links < thdist ; exclude last generation 
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
  	for (i in 1:nrow(G1)){
  		csize <- csize + 1
  		if ( csize >= thsize ){
  			IT <- G1$timesequenced[i] + ritdist() 
  			break
  		}
  	}
  	# NOTE excluding here clustered infections detected after intervention time
  	# because intervention would not act directly on these; they may count as infections averted etc;
  	G1 <- G1[ G1$timesequenced < IT , ]
  	
  	# potential infections averted; should count clustered _and_ nonclustered cases; should follow only one generation 
  	pia = 0 
  	# potential undiagnosed time averted (cumulative across all cases); should follow only one generation 
  	puta = 0 
  	piapids <- c() 
  	if (!is.infinite( IT ) )
  	{
  		piapids <-  D$recipient[D$donor %in% G1$pid] |> union(G1$pid) 
  		G2 <- G[ G$pid %in% piapids , ]
  		pia <- sum( G2$timeinfected  > IT )
  
  		G3 <- G2[ G2$timeinfected <= IT & G2$timediagnosed > IT , ]
  		puta <- sum( G3$timediagnosed - IT )
  	}
  	c(  pia = pia 
  	  , puta = puta 
  	  , interventiontime = IT 
  	  , nc = nrow( G1 )
  	)
  }
  
  distsize_intervention <- function(thdist = 0.005, thsize = 5, thgrowth=NA
   			  , ritdist = function() rexp(1,rate=1/90) )
  {
  	o = lapply( 1:length( Ds ), function(i) proc_cluster( Ds[[i]], Gs[[i]]
  							     , thdist
  							     , thsize
  							     , thgrowth
  							     , ritdist))
  	odf = do.call( rbind, o ) |> as.data.frame()
  	odf1 <- odf[ !is.infinite( odf$interventiontime ), ] # exclude sims where no cluster found 
  
  	lastgen <- max( Gall$generation )
  
  	list(
  		o = odf1 
  		, propintervened = sum( odf1$nc ) / sum(Gall$generation>0 & Gall$generation<lastgen)# TODO remove 1st and last gen from denominator 
  		, puta = c(mean( odf1$puta / odf1$nc ), var( odf1$puta / odf1$nc  ) / nrow(odf1) )
  		, pia =  c( mean( odf1$pia / odf1$nc ), var( odf1$pia / odf1$nc ) / nrow(odf1 ))
  	     )
  }
  
  invisible( '
  	select random case from generation after seed and before last generation 
  	sim intervention time 
  	trace one generation + donor 
  ')
  random_intervention <- function(
  				ritdist = function(n) rexp(n,rate=1/90) 
  )
  {
  	lastgeneration <- max( Gall$generation )
  	G1 <- Gall[ Gall$generation > 0 & Gall$generation < lastgeneration, ] # gen 2+ 
  	# make pids unique 
  	G1$pid <- paste(sep='.', G1$pid, G1$simid )
  	D <- Dall; D$donor <- paste(sep='.', D$donor, D$simid )
  	D$recipient <- paste(sep='.', D$recipient, D$simid )
  	G <- Gall; G$pid <- paste(sep='.', G$pid, G$simid)
  	#intervention time 
  	G1$IT <-  G1$timesequenced + ritdist( nrow(G1 )) 
  	
  	proc_indiv <- function(pid, IT )
  	{
  		piapids <- D$recipient[ D$donor == pid] |> union( c( pid, D$donor[D$recipient==pid] ))	
  		G2 <- G[ G$pid %in% piapids , ]
  		pia <- sum( G2$timeinfected > IT )
  		G3 <- G2[ G2$timeinfected <= IT & G2$timediagnosed > IT, ]
  		puta <- sum( G3$timediagnosed - IT )
  		c( pia, puta )
  	}
  	mapply(proc_indiv, G1$pid, G1$IT ) -> o 
  	o <- as.data.frame( t( o ) )
  	colnames(o) <- c('pia', 'puta' )
  	o
  
  	list(
  		o = o 
  		, propintervened = NA # output is not meaningful for this intervention 
  		, puta = c( mean( o$puta ), var(o$puta)/nrow(o) )
  		, pia =  c( mean( o$pia ), var(o$puta)/nrow(o))
  	     )
  }
  
  invisible('
  Intervene on cases with acute infection 
  Simulated RITA test
  6 month average detection window 
  	  ')
  rita_intervention <- function(
  				ritdist = function(n) rexp(n,rate=1/90) 
  			       )
  {
  	lastgeneration <- max( Gall$generation )
  
  	# make pids unique 
  	D <- Dall; D$donor <- paste(sep='.', D$donor, D$simid )
  	D$recipient <- paste(sep='.', D$recipient, D$simid )
  	G <- Gall; G$pid <- paste(sep='.', G$pid, G$simid)
  	
  	# rita test, 6 month average detection window 
  	G$rita <- with( G, (timediagnosed - timeinfected) < rexp( nrow(G), 1/(6*30)) )
  
  	# filter generations & only rita+ 
  	G1 <- G[ G$rita & (G$generation > 0) & (G$generation < lastgeneration), ] # gen 2+ 
  	
  	# intervention time 
  	G1$IT <-  G1$timesequenced + ritdist( nrow(G1 )) 
  	
  	proc_indiv <- function(pid, IT )
  	{
  		piapids <- D$recipient[ D$donor == pid] |> union( c( pid, D$donor[D$recipient==pid] ))	
  		G2 <- G[ G$pid %in% piapids , ]
  		pia <- sum( G2$timeinfected > IT )
  		G3 <- G2[ G2$timeinfected <= IT & G2$timediagnosed > IT, ]
  		puta <- sum( G3$timediagnosed - IT )
  		c( pia, puta )
  	}
  	mapply(proc_indiv, G1$pid, G1$IT ) -> o 
  	o <- as.data.frame( t( o ) )
  	colnames(o) <- c('pia', 'puta' )
  	o
  
  	list(
  		o = o 
  		, propintervened = nrow(G1) / nrow(Gall) 
  		, puta = c( mean( o$puta ), var(o$puta)/nrow(o) )
  		, pia =  c( mean( o$pia ), var(o$puta)/nrow(o))
  	     )
  }
  
  
  network_intervention <- function( thdegree = 30.
  	, ritdist = function(n) rexp(n,rate=1/90) 
  )
  {
  	w <- 7.0 / 2.0 # ratio duration F to G 
  	ww <-  7.0*30.0 # if unit daiy oo contact rate, E partners in 7 months (F duration) 
  
  	lastgeneration <- max( Gall$generation )
  
  	# make pids unique 
  	D <- Dall; D$donor <- paste(sep='.', D$donor, D$simid )
  	D$recipient <- paste(sep='.', D$recipient, D$simid )
  	G <- Gall; G$pid <- paste(sep='.', G$pid, G$simid)
  	
  	# compute weighted degree ( exp partners over 7 months )
  	G$degree  <-  with( G, Fdegree + w*Gdegree + ww*Hdegree)
  
  	# filter generations & only degree above threshold  
  	G1 <- G[ (G$degree>=thdegree) & (G$generation > 0) & (G$generation < lastgeneration), ] # gen 2+ 
  	
  	# intervention time 
  	G1$IT <-  G1$timesequenced + ritdist( nrow(G1 )) 
  	
  	proc_indiv <- function(pid, IT )
  	{
  		piapids <- D$recipient[ D$donor == pid] |> union( c( pid, D$donor[D$recipient==pid] ))	
  		G2 <- G[ G$pid %in% piapids , ]
  		pia <- sum( G2$timeinfected > IT )
  		G3 <- G2[ G2$timeinfected <= IT & G2$timediagnosed > IT, ]
  		puta <- sum( G3$timediagnosed - IT )
  		c( pia, puta )
  	}
  	mapply(proc_indiv, G1$pid, G1$IT ) -> o 
  	o <- as.data.frame( t( o ) )
  	colnames(o) <- c('pia', 'puta' )
  	o
  
  	list(
  		o = o 
  		, propintervened = nrow(G1) / nrow(Gall) 
  		, puta = c( mean( o$puta ), var(o$puta)/nrow(o) )
  		, pia =  c( mean( o$pia ), var(o$puta)/nrow(o))
  	     )
  }
  
  
  
  # Run all interventions with error handling
  cat("Running interventions...\n")
  
  ods5 <- tryCatch(distsize_intervention(thsize = cluster_size_5), error = function(e) {
    cat("Error in ods5:", e$message, "\n")
    list(propintervened = 0, puta = c(0, 0, 0, 0), pia = c(0, 0, 0), total_contacts = 0)
  })
  
  ods2 <- tryCatch(distsize_intervention(thsize = cluster_size_2), error = function(e) {
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
  
  rbind(
  with( ods5, c( propintervened, puta ))
  , with( ods2, c( propintervened, puta ))
  , with( orand, c( propintervened, puta ))
  , with( orita, c( propintervened, puta ))
  , with( onet, c( propintervened, puta ))
  ) |> as.data.frame() -> odf 
  colnames( odf ) <- c( "Proportion intervened", "Mean", "Variance" )
  odf$lb <- with( odf, Mean - 1.96*sqrt(Variance) )
  odf$ub <- with( odf, Mean + 1.96*sqrt(Variance) )
  odf
  rownames(odf) <- c(
  'Size=5,D=0.005'
  , 'Size=2,D=0.005'
  , 'Random allocation'
  , 'RITA'
  , 'Network, partners>30'
  )
  odf
  odf1 <- round( odf, 2 )
  knitr::kable(odf1)
