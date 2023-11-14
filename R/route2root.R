route2root <- function(pP, focLin) 
  {
  if (class(pP) != "paleoPhylo") stop("pP is not of class 'paleoPhylo'")
  
  x <- which(pP$st==max(pP$st))
  rtDur <- pP$st[x]-pP$en[x]
  rt <- pP$nm[x]
  
  allAnc <- immAnc <- focLin
  tm2rt  <- pP$en[pP$nm==focLin]
  
  while(immAnc != rt)
    {
    i <- which(pP$nm == immAnc)
    tm2rt <- c(tm2rt, pP$st[i])
    immAnc <- pP$pn[i]
    allAnc <- c(allAnc,immAnc)
    }
    
  if(length(allAnc)!=length(tm2rt)) stop ("Number of lineages does not match number of durations.")
  tm2rt <- diff(c(tm2rt, pP$st[x]))
  out <- list(path=rev(allAnc), duration=rev(tm2rt), nNode=length(allAnc))
  return(out)
  }
