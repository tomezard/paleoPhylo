stratUnc <- function (uSR = NULL, pP, lwdLin = 1, ltyLin=1, clCd=1:length(unique(pP$grp)), style=NULL)
  {
  if (class(pP) != "paleoPhylo") stop(" object is not of class 'paleoPhylo'")

  if(is.null(style))
    {styles <- list(certain = c(lwdLin, "black", 1), lazarus = c(lwdLin, "grey60", 1), 
      extension = c(lwdLin, "black", 2), CI = c(lwdLin * 0.4, "black", 1), point=c(lwdLin * 4, "black", 1))}
  if(!is.null(style))
    {
    nu <- length(unique(style))
    styles <- vector("list", nu)
    for(k in 1:nu) styles[[k]] <- c(lwdLin, clCd[k], ltyLin)
    }

  ##define some quantities to reduce the amount of characters required later on.
  nLin <- length(pP$nm)
  nLwdLin <- length(lwdLin)
  nLtyLin <- length(ltyLin)
  nUnqGrp <- length(unique(pP$grp))
  
  dts <- vector("list", nLin)
  for (i in 1:nLin) {dts[[i]] <- c(pP$st[i], pP$en[i]) }
  usr <- list(id = as.list(pP$nm), dates = dts, types = as.list(rep(1, nLin)), styles = styles)

  if((sum(is.na(pP$grp))==0) && is.null(uSR))
    {
    for (i in 1:nLin) {usr$types[[i]] <- as.numeric(pP$grp)[i]} 
    if(nLwdLin < nUnqGrp) lwdLin <- rep(lwdLin, ceiling(nUnqGrp/nLwdLin))
    if(nLtyLin < nUnqGrp) ltyLin <- rep(ltyLin, ceiling(nUnqGrp/nLtyLin))
    usr <- list(id=usr$id,dates=usr$dates,types=usr$types, styles=NULL)
    for (i in 1:nUnqGrp) usr$styles[[i]] <- c(lwdLin[i], clCd[i], ltyLin[i])
    }
    
  if (!is.null(uSR) & !is.data.frame(uSR)) 
    {
    reqTypes <- unique(unlist(uSR$types))
    nRqTypes <- length(reqTypes)
    nDates   <- length(uSR$dates)
    
    if (sum(sapply(1:nRqTypes, function(i) {is.null(uSR$styles[[reqTypes[i]]]) })) != 0)
        stop("Not all types have a defined style")
        cat("\n")
    badd <- which(sapply(1:nDates, function(i) {length(uSR$dates[[i]]) != (length(uSR$types[[i]]) + 1)}))
        if (length(badd) > 0) {cat(paste("Mismatch at lineage ", badd, sep = ""))}
        cat("\n")
    if (sum(sapply(1:nDates, function(i) { length(uSR$dates[[i]]) != (length(uSR$types[[i]]) +  1)})) != 0)
      {
      ##which is the problem?	
      stop("The number of dates and types do not match. N(dates)!=(N(types)+1) for all lineages.")
      }
    usr <- uSR
    }
    
  if(is.data.frame(uSR))
    {
    for (k in 1:length(uSR$id))
      {
      flag <- FALSE
      focID <- which(usr$id==uSR$id[k])
      oldDates <- as.numeric(usr$dates[[focID]])
      splitLoc <- which(oldDates<=uSR$st[k] & oldDates>=uSR$en[k])[1]
      if (is.na(splitLoc)) {splitLoc <- which(oldDates<=uSR$st[k])  ;  flag<-TRUE}
      newDates <- rev(sort(c(oldDates[-splitLoc],uSR$st[k],uSR$en[k])))
      if (flag) newDates <- c(newDates, oldDates[length(oldDates)])
      
      oldTypes <- newTypes <- usr$types[[focID]]
      if(splitLoc==1) {newTypes <- c(uSR$type[k],oldTypes)}
      if(splitLoc==(length(newDates)-1)) {newTypes <- c(oldTypes,uSR$type[k])}
      if(splitLoc!=(length(newDates)-1) & splitLoc!=1) 
        {
        if(uSR$st[k]!=uSR$en[k])  {newTypes <- c(oldTypes[1:(splitLoc-1)],uSR$type[k],oldTypes[(splitLoc+1):length(oldTypes)])}
        if(uSR$st[k]==uSR$en[k])  {newTypes <- c(oldTypes[1:(splitLoc-1)],uSR$type[k],oldTypes[(splitLoc-1):(length(oldTypes)-1)])}
        }
      if((oldDates[1]==uSR$st[k] & oldDates[2]==uSR$en[k] & length(oldDates)==2)==FALSE)
        {
        usr$types[[focID]] <- as.numeric(newTypes)
        usr$dates[[focID]] <- as.numeric(newDates)
        }
      if((oldDates[1]==uSR$st[k] & oldDates[2]==uSR$en[k] & length(oldDates)==2)==TRUE)  {usr$types[[focID]] <- as.numeric(uSR$type[k])}
      }
    }
  return(usr)
  }