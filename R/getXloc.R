getXloc <- function (pP) 
  {
  	if (class(pP) != "paleoPhylo") stop("object is not of class 'paleoPhylo'")
     {
     dat <- data.frame(nm = pP$nm, pn = pP$pn, st = pP$st, en = pP$en)
     #dat <- dat[rev(order(dat$st, dat$en)),]
     pos <- 0.5
     time <- minTime <- min(-pP$st)
     ids <- pP$nm[1]
     tmFrm <- round(-dat$st, 5)
     nn <- length(pP$nm)
     dex <- numeric(nn)
     nAncs <- vector("list", nn)
     for (i in 1:nn) nAncs[[i]] <- route2root(pP, pP$nm[i])$path
     for (j in 1:nn) dex[j] <- sum(sapply(1:nn, function(k) pP$nm[j] %in% nAncs[[k]]))
     
     while (time != max(tmFrm))
       {
       time <- tmFrm[(round(tmFrm, 8) > round(time, 8)) == TRUE][1]
       event <- dat[which(-round(dat$st, 5) == round(time, 5)), ]
       parents <- unique(as.character(event$pn))
       inIds <- sapply(1:length(parents), function(i) sum(intersect(ids, parents) == as.list(parents)[[i]]) > 0)
       parents <- c(parents[inIds == TRUE], parents[inIds == FALSE])
       for (i in 1:length(parents))
         {
         singleEvent <- event[as.character(event$pn) == parents[i], ]
         if (length(singleEvent[, 1]) > 1) evnt <- 1
         if (length(singleEvent[, 1]) == 1) evnt <- 0
         focInd <- as.character(singleEvent$nm)
         pntLoc <- which(ids == parents[i])
         if (time != minTime)
           {
           fcPrnt <- as.character(unique(singleEvent$pn))
           fcSist <- as.character(unique(singleEvent$nm))
           if (length(pntLoc) == 0 & length(pos) > 1) 
             stop(paste("There is a lack of congruence in the tree around", focInd, "\nHave a look at the ancestor", 
               fcPrnt, "but also the sister species", fcSist[1], "and", fcSist[2], "\n"))
           }
         parentPos <- pos[pntLoc]
         ids <- c(ids, focInd)
         sortPos <- sort(pos)
         if (length(pos) > 1)
           {
           whr <- which(sortPos == parentPos)
           lowPos <- (sortPos[whr - 1] + sortPos[whr])/2
           hghPos <- (sortPos[whr + 1] + sortPos[whr])/2
           if (evnt == 0) 
             {
             pta <- dex[which(pP$nm == parents[i])]
             ida <- dex[which(pP$nm == focInd)]
                    
             if (pta > ida & parentPos >= 0.5) newPos <- lowPos else newPos <- hghPos
             if (pta > ida & parentPos < 0.5)  newPos <- hghPos else newPos <- lowPos
                    
             if (length(newPos) == 0 & parentPos == max(pos)) newPos <- extendrange(pos)[2]
             if (length(newPos) == 0 & parentPos == min(pos)) newPos <- extendrange(pos)[1]
             }
          if (evnt == 1)
            {
             if (parentPos == min(pos)) lowPos <- extendrange(pos)[1]
             if (parentPos == max(pos)) hghPos <- extendrange(pos)[2]
             newPos <- c(lowPos, hghPos)
             if (length(focInd) > 2) 
               {for (j in 3:length(focInd)) newPos <- c(newPos, (newPos[j - 2] + parentPos)/2)}
             }
           pos <- c(pos, newPos)
           }
              
         if (length(pos) == 1)  ifelse(evnt == 0, pos <- seq(0, 1, 1), pos <- c(0.5, 0, 1))
              
         lnsp <- length(sortPos)
         lnfi <- length(focInd)
         if (evnt == 1 & lnfi > 4 & lnsp == 1)  stop("Cannot yet handle polytomies of order > 4 immediately post root.")
         if (evnt == 1 & lnfi == 3 & lnsp == 1) pos <- c(0.5, 0, 1, 0.25)
         if (evnt == 1 & lnfi == 4 & lnsp == 1) pos <- c(0.5, 0, 1, 0.25, 0.75)
         if (length(pos) != length(ids)) 
            stop(paste("There is not a position for all individuals.\n The problem is something to with", focInd))
         pos <- as.numeric(as.factor(pos))/length(pos)
         }
       }
     locDat <- data.frame(ids, pos)
     locDat$xx <- pos #as.numeric(as.factor(locDat$pos))/length(locDat$pos)
     #print(length(locDat$xx))
     #print(length(unique(locDat$xx)))
     if (length(locDat$xx) != length(unique(locDat$xx))) stop("non-unique locations")
     cmbDat <- merge(dat, locDat, by.x = c("nm"), by.y = c("ids"))
     }
  return(cmbDat)
  }
  