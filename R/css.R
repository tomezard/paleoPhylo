css <- function(pP, lastDate=min(pP$en))
  {
  if(class(pP)!="paleoPhylo") stop("object is not of class 'paleoPhylo'")
  
  ln <- length(pP$st)
  spec <- sapply(1:ln, function(i) length(which(pP$pn==pP$nm[i]))>0)
  dur <- with(pP, st-en)
  extRtotR <- numeric(length(spec))
  notExtC <- 1 - (pP$en==0 | spec==1)
  notExtant <- 1 - (pP$en==lastDate)
  
  dat <- with(pP, data.frame(nm=nm,st=st,en=en))
  for (i in 1:ln) 
	{
    tmp <- dat[dat$en<=dat$st[i] & dat$st>=dat$en[i], ]
    extRtotR[i] <- (length(tmp$nm[tmp$en >= pP$en[i]])/ln) / (1-mean(notExtant))
 	}

  css <- list(nm=pP$nm, css=extRtotR*dur, Duration=dur, notExtant=notExtant, notExtC=notExtC)
  return(css)
  }