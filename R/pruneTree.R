pruneTree <- function(pP, focLin=NULL, focDate=NULL,
  keepTips=TRUE, keepFocLin=TRUE,
  letSpeciate=FALSE, letDie=FALSE, pruneDead=FALSE,
  outPhylo=FALSE, collapseBranches=FALSE)
  {
  #pP <- p93; focLin <- "90" ; focDate <- 51; keepTips=FALSE; keepFocLin=TRUE
  #pruneDead <- TRUE; letSpeciate <- TRUE; letDie <- TRUE;outPhylo=FALSE
  #test <- pruneTree1(p93, "90", 51, keepTips=FALSE)$paleoPhylo
  #with(test, tapply(nm, pn, length))
  
  if (class(pP) != "paleoPhylo") stop(" object is not of class 'paleoPhylo'")
  nn <- length(pP$nm)
  p2r <- vector("list", nn)
  for(n in 1:nn) p2r[[n]] <- route2root(pP, pP$nm[n])$path
  prT  <- with(pP, data.frame(nm, pn, st, en, label, grp))
  #
  
  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
  #*#*#*if pruning around a focal lineage
  if(!is.null(focLin))
    {
    kT <- ifelse(keepTips, 1, 0)
    kp   <- function(n) length(intersect(focLin, p2r[[n]]))==kT
    	#could be amended to include all descendants of this lineage, too
    	#that way can prune2end of lineage
    keep <- sapply(1:nn, kp)

    if(keepFocLin) keep[which(pP$nm==focLin)] <- TRUE
    prT  <- prT[keep,]

    if(keepTips)
      {
      #remove the ancestor of the oldest lineage as not in the pruned tree
      whr <- which(prT$st==max(prT$st))
      prT$pn[whr] <- NA 
      }

     if(pruneDead)
        {
        if(letDie) extant <- which(prT$en<focDate) else extant <- which(prT$en<=focDate)
        nExt             <- length(extant)
        fromExtantTips   <- vector("list", nExt)  
        for(n in 1:nExt)
          {
          whr <- which(pP$nm==prT$nm[extant[n]])
          if(length(whr)>0)
            {
            fromExtantTips[[n]] <- p2r[[whr]]
            }
          }
          
        unqFromExtant    <- unique(unlist(fromExtantTips))
        nn <- length(unqFromExtant)
        whrs <- numeric(nn)  
        for(n in 1:nn)
          {
          if(length(intersect(prT$nm, unqFromExtant[n]))>0)
            {
            whrs[n] <- which(prT$nm==unqFromExtant[n])
            }
          }  
        whrs <- whrs[whrs!=0]
        prT <- prT[whrs,]
        if(letDie) prT$en[prT$en<focDate] <- focDate else prT$en[prT$en<=focDate] <- focDate
        }

    }  
  prT <- prT[rev(order(prT$st)),]

  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
  #*#*#*if pruning around a focal date
  if(!is.null(focDate) & is.null(focLin))
    {
    if(!keepTips)
      {
      focSt <- prT$st[which(pP$nm==focLin)]
      #if(focSt<focDate)
       #    stop(paste("the prune2date (", focDate, ") is before the start focDate (", focSt, ") of focLin (", as.character(focLin), ")", sep=""))

      if(letSpeciate) prT <- prT[prT$st>=focDate,] else prT <- prT[prT$st>focDate,]
          #remove species that started too late
       
      if(pruneDead)
        {
        if(letDie) extant <- which(prT$en<focDate) else extant <- which(prT$en<=focDate)
        nExt              <- length(extant)
        fromExtantTips    <- vector("list", nExt)  
        for(n in 1:nExt)
          {
          whr <- which(pP$nm==prT$nm[extant[n]])
          if(length(whr)>0)
            {
            fromExtantTips[[n]] <- p2r[[whr]]
            }
          }
          
        unqFromExtant    <- unique(unlist(fromExtantTips))
        nn <- length(unqFromExtant)
        whrs <- numeric(nn)  
        for(n in 1:nn)
          {
          if(length(intersect(prT$nm, unqFromExtant[n]))>0)
            {
            whrs[n] <- which(prT$nm==unqFromExtant[n])
            }
          }  
        whrs <- whrs[whrs!=0]
        prT <- prT[whrs,]
        if(letDie) prT$en[prT$en<focDate] <- focDate else prT$en[prT$en<=focDate] <- focDate
        }
     if(!pruneDead) if(letDie) prT$en[prT$en<focDate] <- focDate else prT$en[prT$en<=focDate] <- focDate
      }
    }
    if(!is.null(focLin) & !is.null(focDate))
      {
      if(!keepTips)
        {stop("prune2focDate currently coded only for the 'root-end' of the tree when focal lineage specified.")}
      if(keepTips)
        {
        prT$st[which(prT$nm==focLin)] <- focDate
        keep <- prT$st<=focDate
        prT <- prT[keep,]
        }
      } 
    


  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
  #*#*#*output stuff 
  prTpP <- with(prT, as.paleoPhylo(nm, pn, st, en, label=label, grp=grp))
  if(collapseBranches) prTpP <- collapseBranches(prTpP)
  if(outPhylo) prTpa <- reorder(buildApe(createBifurcate(prTpP)))  else prTpa <- NULL
  
  if(!outPhylo) return(prTpP) else return(list(paleoPhylo=prTpP, phylo=prTpa))	
  }