rates.within.bins <-
function(pP, begin=max(pP$st), end=min(pP$en), bin.length = 1){
	#Takes any paleoPhylo object and computes per-lineage origination and extinction rates within the set of bins
	#specified by the user.  Lineages persisting to the specified end are not considered to have gone extinct. Events
	#taking place exactly on bin boundaries are taken to occur in the later bin.
	
	#check it's a paleoPhylo object
  if (class(pP) !="paleoPhylo") stop("Object passed to prune.to.date is not of class paleoPhylo!")
  
  #check there's an interval within which rates will be calculated
  if (end >= begin) stop("end must be later than begin!")
  if (bin.length <= 0) stop("bin.length must be positive!")
  
  cb <-createBifurcate(pP) #Can't be sure the paleoPhylo object passed will be bifurcating, so make it so
  
  df <- with(cb, data.frame(nm, pn, st, en, xx, label, grp))
    #Convert to data frame
  bin.start <- seq(from = begin, to = end+bin.length, by = -bin.length)
  rwb <- data.frame(bin.start)
  rwb$mu <- rwb$lambda <- rwb$extinctions <- rwb$originations <- rwb$branch.length <- rwb$N.at.start <- NA
	
  #Identify first lineage in createBifurcate object, whose origination should not be counted
  given <- cb$nm[1]

  for (i in 1:length(rwb$bin.start))
    {
    bs <- bin.start[i]
      #Identify lineages that do something at or after the start of the bin, but before the end of the bin
      #(Note that this means lineages surviving to the specified end time are not viewed as having gone extinct)
    whoDoesSomething <- cb$nm[cb$en <= bs & cb$en > (bs-bin.length)]
      #Work out which of these disappear because they originate new lineages
    newLineages <- cb$nm[cb$st <= bs & cb$st > (bs-bin.length)]
    whoOriginates <- unique(cb$pn[is.element(cb$nm, newLineages)])
      #The remainder have died out
    whoDies <- setdiff(whoDoesSomething, whoOriginates)
    whoOriginates <- setdiff(whoOriginates, given)
      #Don't count origin of first lineage in createBifurcate object
    
    rwb$originations[i] <- length(whoOriginates)
    rwb$extinctions[i] <- length(whoDies)

	  #Identify and sum branch lengths that fall within the bin
	cb.copy <- cb
	cb.copy$st <- pmin(cb.copy$st, bs)
	cb.copy$en <- pmax(cb.copy$en, bs-bin.length)
	bl <- pmax(0, cb.copy$st - cb.copy$en)
	rwb$branch.length[i] <- sum(bl)

    #Find diversity at start of bins
    ##should it not be length(cb$nm[cb$st>=bs & cb$en<bs])? ANDY SAYS NO!
    rwb$N.at.start[i] <- length(prune.to.date(cb, bs, let.speciate=FALSE, let.die=FALSE)$phylo.tree$tip.label)
    }
  
  #Compute per-lineage rates from the totals computed in the loop
  rwb$lambda <- rwb$originations/rwb$branch.length
  rwb$mu <- rwb$extinctions/rwb$branch.length	

  return(rwb)
  }

