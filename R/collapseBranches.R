collapseBranches <- function(pP)
  {
  if (class(pP) != "paleoPhylo") stop(" object is not of class 'paleoPhylo'")
  
  df <- with(pP, data.frame(nm, pn, st, en, label, grp))
  df$nm <- as.character(df$nm)
  df$pn <- as.character(df$pn)
  df$label <- as.character(df$label)
  
  nn <- length(df$nm)
  cnts <- ancs <- logical(nn)
  cnts[1] <- NA		#root
  for(k in 2:nn)
    {
    cnts[k] <- df$st[k]==df$en[which(df$nm==df$pn[k])]
    ancs[k] <- sum(df$pn==df$pn[k], na.rm=TRUE)==1
    }
  whrs <- cnts & ancs
  	#whrs contains a vector saying whether the focal species is a continuation of the previous one
  	#and can therefore be collapsed into it
  
  n2c <- which(whrs)
  lc <- length(n2c)
  if(lc>0)
    {
    for(k in 1:lc)
      {
      dec <- as.character(df$nm[n2c[k]])
      whrDec <- n2c[k]
      anc <- as.character(df$pn[whrDec])
      whrAnc <- which(df$nm==anc)

      df$en[whrAnc] <- df$en[whrDec]
      df$nm[whrAnc] <- paste(anc, dec, sep="_")
      df$label[whrAnc] <- paste(df$label[whrAnc], df$label[whrDec], sep="_")
      df$pn[which(df$pn==anc | df$pn==dec)] <- paste(anc, dec, sep="_")
      }
    dfc <- df[-n2c,] 
    }
  if(lc==0) dfc <- df
  
  cpP <- with(dfc, as.paleoPhylo(nm, pn, st, en, label=label, grp=grp))
  return(cpP)
  }