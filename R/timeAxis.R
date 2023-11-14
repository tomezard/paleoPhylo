timeAxis <-
function (tmScl, whatTime, addTimeLine = "classic", l2r = FALSE, nmLim = 2, cexText = 0.5, srtText=0, cexTime = NULL, cexLab = 0.5, whSpc=0.01, lwdLin=1, dumpLast=FALSE, sz=0.2) 
  {
  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
  #*#*#some initial conditions and initialize the plot
  if(l2r) xLm <- extendrange(-tmScl$MA, f=whSpc) else xLm <- range(0,1)
  if(l2r) yLm <- range(0,1) else yLm <- extendrange(-tmScl$MA, f=whSpc)
  prntName <- c(1, 0+(abs(diff(tmScl$MA))>nmLim))
  if(is.null(cexTime)) cexTime <- cexText[length(cexText)]
   
  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
  #*#*#draw time options, CLASSIC
  if(length(grep("c",addTimeLine))>0)
	{
	##an arbitrary equation to divide space between timescale and phylogeny; fix at 1/3 for now
	  if(l2r) fig.mat <- matrix(c(0,0,1,1,0,sz,sz,1), nrow=2) else fig.mat <- matrix(c(0,sz,sz,1,0,0,1,1), nrow=2)
	  split.screen(fig.mat)
	  screen(1)
	  par(mar = rep(0.2, 4))
	  plot(0:1,0:1,type='n',xlab="",ylab="",xlim=xLm,ylim=yLm,axes=FALSE)

      xv <- seq(0, 0.9, length.out = (length(whatTime) + 1))
      for (i in 1:length(whatTime))
	    {
	    ii <- which(colnames(tmScl) == whatTime[i])
 		st <- sort(tapply(tmScl$MA, as.character(tmScl[, ii]), max))
 		st <- st[-length(st)]
		en <- sort(tapply(tmScl$MA, as.character(tmScl[, ii]), max))
 		en <- en[-1]
 		lb <- names(en) 			   		
 		
 		if(min(c(st,en))!=min(tmScl$MA) | max(c(st,en))!=max(tmScl$MA))
 		{
		en <- c(min(st),en)
		st <- c(min(tmScl$MA),st)
		lb <- unique(names(c(st,en)))
		if(length(lb)!=length(en)) lb <- lb[lb!=""]
		}
	 
	   dL <-length(st)
	   if(dumpLast) {st <- st[-dL] ; en <- en[-dL] ; lb <- lb[-dL]}
	   if(length(cexText)<length(whatTime)) cexText <- rep(cexText, length(whatTime)/length(cexText))
	   if(length(srtText)<length(whatTime)) srtText <- rep(srtText, length(whatTime)/length(srtText))

       vis <- (abs(st - en) > nmLim)
       par(srt = 90 - (l2r*90) + srtText[i])
       if (l2r == FALSE)
         {
       	  rect(xv[i], -st, xv[i + 1], -en)
       	  text((xv[i] + xv[i + 1])/2, y=-c((st + en)/2), lb, col=vis, cex=cexText[i], adj=c(0.5,0.5))
       	  segments(xv[i], -st, xv[i + 1], -st, lwd=lwdLin)
       	  }
       if (l2r == TRUE)
         {
         rect(-st, xv[i], -en, xv[i + 1])
         text(-c((st + en)/2), y=(xv[i] + xv[i + 1])/2, lb, col=vis, cex=cexText[i], adj=c(0.5,0.5))
         segments(-st, xv[i], -st, xv[i + 1], lwd=lwdLin)
         }
       }

    tmAxs <- unique(c(st,en))
    vis <- c(1, 0 + (abs(diff(0-tmAxs)) > nmLim))
    par(srt = 90 - (l2r*90) + srtText[length(srtText)])
    if (l2r) text(-tmAxs, 0.99, abs(tmAxs), cex=cexTime, col=vis) else text(0.99, -tmAxs, abs(tmAxs), cex=cexTime, col=vis)
      
    screen(2, FALSE)
    par(mar = rep(0.2, 4))
    plot(0:1,0:1,type='n',xlab="",ylab="",xlim=xLm,ylim=yLm,axes=FALSE)
    thk <- rep(1,length(tmAxs))
    ##if more than one time box is to be drawn, split based on previous two
    if(length(whatTime)>1)
      {
      i1 <-  which(colnames(tmScl) == whatTime[i])
      i2 <-  which(colnames(tmScl) == whatTime[i-1])
      th1 <- tmScl$MA[(!duplicated(tmScl[,i1][-1]))]
      th2 <- tmScl$MA[(!duplicated(tmScl[,i2][-1]))]
      thk <- duplicated(sort(c(th1,th2))) + 1
      thk <- thk[-(which(thk==2)-1)]
      }
    if (l2r) arrows(-tmAxs, 0.05, -tmAxs, yLm[2], col = "slategray4", lwd = thk, length=0) else arrows(0.05, -tmAxs, yLm[2], -tmAxs, col = "slategray4", lwd = thk, length=0)
  }
  
  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
  #*#*#draw time options, TUBE
    if(length(grep("t",addTimeLine))>0)
      {
      if (length(whatTime) != 2) stop("if addTimeLine=='tube' then whatTime currently must have 2 levels.")
      
      par(mar = rep(0.2, 4))
      plot(0:1,0:1,type='n',xlab="",ylab="",xlim=xLm,ylim=yLm,axes=FALSE)
	  
      col1 <- which(colnames(tmScl) ==  whatTime[1])
      col2 <- which(colnames(tmScl) ==  whatTime[2])
      splits <- paste(tmScl[, col1], tmScl[, col2], sep="  ")
      boundaries <- sort(tapply(tmScl$MA, splits, max))
      if(dumpLast) boundaries <- boundaries[-length(boundaries)]
      st <- boundaries[1:(length(boundaries) - 1)]
      en <- boundaries[2:length(boundaries)]
      
      cls <- c("white", "gray75")
      if(l2r==TRUE) rect(-st, 0, -en, 1, col=cls, border=cls) else rect(0, -st, 1, -en, col=cls, border=cls)
      nms1 <- sapply(1:length(en), function(i) strsplit(names(en),"  ")[[i]][1])
      nms2 <- sapply(1:length(en), function(i) strsplit(names(en),"  ")[[i]][2])
      
      par(srt=0 + 270*l2r)
      if (l2r == FALSE)
        {
        text(0.49, y =-(st + en)/2, nms1, col=rev(cls), cex=cexText, adj=c(1, 0.5))
        text(0.51, y =-(st + en)/2, nms2, col=rev(cls), cex=cexText, adj=c(0, 0.5))
        }
      if (l2r == TRUE)
        {
        text(-(st + en)/2, y=0.51, nms1, col=rev(cls), cex=cexText, adj=c(1, 0.5))
        text(-(st + en)/2, y=0.49, nms2, col=rev(cls), cex=cexText, adj=c(0, 0.5))
        }
    tmLab <- -tmScl$MA[prntName == 1]
    mtext(abs(tmLab), side =2-l2r, at=tmLab, cex=cexTime, line=-1)
    }
  }

