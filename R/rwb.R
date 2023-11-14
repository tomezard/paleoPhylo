rwb <- function (pP, bL = 1, st = max(pP$st), en = min(pP$en)) 
  {
  if (class(pP) != "paleoPhylo") 
       stop("Object is not of class paleoPhylo")
  if (en >= st) 
        stop("End must be later than start.")

  bpp <- with(pP, data.frame(nm, pn, st, en))
   
  if(length(bL) > 1) {bins <- bL}
  if(length(bL) == 1) {bins <- seq(st, en, -bL)}
  nb <- length(bins)
  if (nb <= 1) 
        stop("Bin Length must be larger than 1.")
  #check there's an interval within which rates will be calculated
  
  # Lineages persisting to the specified end are not considered to have gone extinct.
  #Events taking place exactly on bin boundaries are taken to occur in the later bin.

  rwb <- data.frame(binStart = bins, N = NA, branchLength = NA, orig = NA, extn = NA, lambda = NA, mu = NA)
  for (k in 1:(nb - 1)) 
     {
     if(k!=1)
       {
       rm(orig)
       rm(extn)
       rm(nn)
       rm(brln)
       }
     inBin <- bpp[bpp$st > bins[k + 1] & bpp$en <= bins[k], ]
     nn <- sum(inBin$st > bins[k])
     orig <- sum(inBin$st <= bins[k])
     extn <- sum(inBin$en > bins[k + 1])
     inBin$st[inBin$st > bins[k]] <- bins[k]
     inBin$en[inBin$en < bins[k + 1]] <- bins[k + 1]
     brln <- sum(inBin$st - inBin$en)
     rwb[k, 2:7] <- c(nn, brln, orig, extn, orig/brln, extn/brln)
     }
  rwb <- rwb[-nb, ]
  rwb <- rwb[rev(order(rwb$binStart)), ]
  rwb$logN <- log(rwb$N)
  rwb$growth <- c(NA, rwb$N[-(nb-1)]/rwb$N[-1])
  return(rwb)
  }



