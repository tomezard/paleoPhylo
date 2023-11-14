drawPhylo <- function (pP, uSR = NULL, addTimeLine = "none", tmScl, whatTime, l2r = FALSE, nmLim = 2, cexText = 0.5, srtText=0, cexTime = NULL, cexLab = 0.5, lwdLin = 2, hlty = NULL, barLen=0, whSpc=0.05, dumpLast=FALSE, sz=NULL) 
  {
  if (class(pP) != "paleoPhylo") stop(" object is not of class 'paleoPhylo'")
    {
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #*#*#some initial conditions and initialize the plot
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #*#*#stratigraphic uncertainty
    usrFlg <- is.null(uSR)
  	 uSR <- stratUnc(uSR, pP, lwdLin = lwdLin)
  	 if(l2r) xx <- -range(unlist(uSR$dates)) else xx <- range(pP$xx)
  	 if(l2r) yy <- range(pP$xx) else yy <- -range(unlist(uSR$dates))
  	 if(l2r)
  	   {
  	  xLm <- extendrange(xx, f=whSpc)
  	  xLm <- c(-max(pP$st),xLm[2])
  	  yLm <- yy
  	  }
  	if(!l2r)
  	  {
  	  yLm <- extendrange(yy, f=whSpc)
  	  yLm <- c(-max(pP$st),yLm[2])
  	  xLm <- xx
  	  }

  	  
	#*#*#the time scale used
   if ((length(grep("t", addTimeLine)) > 0 | length(grep("c", addTimeLine)) > 0)) 
     {
     	rngSt <- which(tmScl$MA <= range(unlist(uSR$dates))[1])
     	startPoint <- rngSt[length(rngSt)]
     	endPoint <- which(tmScl$MA >= range(unlist(uSR$dates))[2])[1]
     	tmScl <- tmScl[startPoint:endPoint, ]
     	if (max(pP$st) < 0) tmScl$MA <- 0 - tmScl$MA
     	prntName <- c(1, 0 + (abs(diff(tmScl$MA)) > nmLim))
     }
	  
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #*#*#draw time options
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #op <- par(no.readonly=TRUE)
    par(mar=rep(0, 4), srt=srtText + 90*(!l2r))
	if(is.null(cexTime)) cexTime <- cexText[length(cexText)]
	if(length(grep("c",addTimeLine))>0)
	{
	  ##an arbitrary equation to divide space between timescale and phylogeny
	  if(is.null(sz)) sz <- exp(-(length(pP$nm) + 45)/50) + 0.08

            if (l2r) fig.mat <- matrix(c(0, 0, 1, 1, 0, sz, sz, 1), nrow = 2) else fig.mat <- matrix(c(0, sz, sz, 1, 0, 0, 1, 1), nrow = 2)
            split.screen(fig.mat)
            screen(1)
            plot(xx, yy, type = "n", xlab = "", ylab = "", xlim = xLm, ylim = yLm, axes = FALSE)
            xv <- seq(0.05, 0.7, length.out = (length(whatTime) + 1))
            for (i in 1:length(whatTime))
                {
                ii <- which(colnames(tmScl) == whatTime[i])
                st <- sort(tapply(tmScl$MA, as.character(tmScl[, ii]), max))
                st <- st[-length(st)]
                en <- sort(tapply(tmScl$MA, as.character(tmScl[, ii]), max))
                en <- en[-1]
                lb <- names(en)
                if (min(c(st, en)) != min(tmScl$MA) | max(c(st, en)) != max(tmScl$MA))
                  {
                  en <- c(min(st), en)
                  st <- c(min(tmScl$MA), st)
                  lb <- unique(names(c(st, en)))
                  if (length(lb) != length(en)) lb <- lb[lb != ""]
                  }
                dL <- length(st)
                if (dumpLast) {
                  st <- st[-dL]
                  en <- en[-dL]
                  lb <- lb[-dL]
                }
                if (length(cexText) < length(whatTime)) 
                  cexText <- rep(cexText, length(whatTime)/length(cexText))
                if (length(srtText) < length(whatTime)) 
                  srtText <- rep(srtText, length(whatTime)/length(srtText))
                vis <- (abs(st - en) > nmLim)
                
                if (!l2r) {
                  rect(xv[i], -st, xv[i + 1], -en, lwd=lwdLin)
                  text((xv[i] + xv[i + 1])/2, y = -c((st + en)/2), lb,
                     col = vis, cex = cexText[i], adj = c(0.5, 0.5))
                  segments(xv[i], -st, xv[i + 1], -st, lwd = lwdLin)
                }
                if (l2r)
                  {
                  rect(-st, xv[i], -en, xv[i + 1], lwd=lwdLin)
                  text(-c((st + en)/2), y = (xv[i] + xv[i + 1])/2, lb,
                     col = vis, cex = cexText[i], adj = c(0.5, 0.5))
                  segments(-st, xv[i], -st, xv[i + 1], lwd = lwdLin)
                }
            }
            tmAxs <- unique(c(st, en))
            vis <- c(1, 0 + (abs(diff(0 - tmAxs)) > nmLim))
            
            if (l2r) text(-tmAxs, 1, abs(tmAxs), cex = cexTime, col = vis) else text(1, -tmAxs, abs(tmAxs), cex = cexTime, col = vis, adj=1)
                
            screen(2, FALSE)
            plot(xx, yy, type = "n", xlab = "", ylab = "", ylim = yLm, axes = FALSE)
            thk <- rep(1, length(tmAxs))
            ##if more than one time box is to be drawn, split based on previous two
            if (length(whatTime) > 1) {
                i1 <- which(colnames(tmScl) == whatTime[i])
                i2 <- which(colnames(tmScl) == whatTime[i - 1])
                th1 <- tmScl$MA[(!duplicated(tmScl[, i1][-1]))]
                th2 <- tmScl$MA[(!duplicated(tmScl[, i2][-1]))]
                thk <- duplicated(sort(c(th1, th2))) + 1
                thk <- thk[-(which(thk == 2) - 1)]
            }
            if (l2r) abline(v = -tmAxs, col = "slategray4", lwd = lwdLin) else abline(h = -tmAxs, col = "slategray4", lwd = lwdLin)
        }

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #*#*# straightforward arrow
    
    if(length(grep("ar",addTimeLine))>0)
      {
      plot(xx,yy,type='n',xlab="",ylab="",xlim=xLm,ylim=yLm,axes=FALSE)
      if(l2r==TRUE)  arrows(xx[2],(yLm[1]+yy[1])/2,xx[1],(yLm[1]+yy[1])/2,lwd=lwdLin)
      if(l2r==FALSE) arrows((xLm[1]+xx[1])/2,yy[2],(xLm[1]+xx[1])/2,yy[1],lwd=lwdLin)
      mtext("Time", 2-l2r, line=-1, cex=cexText)
      }
      
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #*#*# funny tube background bit
    if(length(grep("t",addTimeLine))>0)
      {
      lwt <- length(whatTime) 
      if (lwt>2) stop("if addTimeLine=='tube' then whatTime currently must have 2 or fewer levels.")
      plot(xx,yy,type='n',xlab="",ylab="",xlim=xLm,ylim=yLm,axes=FALSE)
      
      col1 <- which(colnames(tmScl) ==  whatTime[1])
      col2 <- which(colnames(tmScl) ==  whatTime[2])
      splits <- paste(tmScl[, col1], tmScl[, col2], sep="  ")
      boundaries <- sort(tapply(tmScl$MA, splits, max))
      if(dumpLast) boundaries <- boundaries[-length(boundaries)]
      st <- boundaries[1:(length(boundaries) - 1)]
      en <- boundaries[2:length(boundaries)]
      
      cls <- c("white", "gray75")
      if(l2r) rect(-st, 0, -en, 1, col=cls, border=cls) else rect(0, -st, 1, -en, col=cls, border=cls)
      nms1 <- sapply(1:length(en), function(i) strsplit(names(en),"  ")[[i]][1])
      nms2 <- sapply(1:length(en), function(i) strsplit(names(en),"  ")[[i]][2])
      
      if(lwt==2)
        {
        if (!l2r)
          {
          text(0.49, y =-(st + en)/2, nms1, col=rev(cls), cex=cexText, adj=c(1, 0.5))
          text(0.51, y =-(st + en)/2, nms2, col=rev(cls), cex=cexText, adj=c(0, 0.5))
          }
        if (l2r)
          {
          text(-(st + en)/2, y=0.51, nms1, col=rev(cls), cex=cexText, adj=c(1, 0.5))
          text(-(st + en)/2, y=0.49, nms2, col=rev(cls), cex=cexText, adj=c(0, 0.5))
          }
        }
      if(lwt==1 & !l2r) {text(0.5, y =-(st + en)/2, nms1, col=rev(cls), cex=cexText)}
      if(lwt==1 & l2r) {text(-(st + en)/2, y=0.5, nms1, col=rev(cls), cex=cexText)}
      
      tmLab <- -tmScl$MA[prntName == 1]
      if(dumpLast) tmLab <- tmLab[-length(tmLab)]
      mtext(abs(tmLab), side =2-l2r, at=tmLab, cex=cexTime, line=-2, las=0+(!l2r))
      #axis(1, at=tmLab, label=abs(tmLab), outer=TRUE, line=-whSpc*100, col="white",  col.ticks="white", cex=cexTime, las=0+(!l2r))
    }
       
    if(length(grep("n",addTimeLine))>0) plot(xx, yy, type='n', xlab="", ylab="",xlim=xLm, ylim=yLm, axes=FALSE)
    
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #*#*#draw the phylogeny itself
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #*#*#begin by sorting out stratigraphic uncertainty & drawing lineages independently
        for (n in 1:length(pP$nm)) {
            mm <- sapply(uSR$types, length)
            for (m in 1:mm[n]) {
                lwdI <- as.numeric(uSR$styles[[uSR$types[[n]][m]]][1])
                clrI <- uSR$styles[[uSR$types[[n]][m]]][2]
                ltyI <- as.numeric(uSR$styles[[uSR$types[[n]][m]]][3])
                if (l2r == FALSE) {
                  arrows(pP$xx[n], -uSR$dates[[n]][m], pP$xx[n], 
                    -uSR$dates[[n]][m + 1], col = clrI, lwd = lwdI, 
                    lty = ltyI, length = (uSR$types[[n]][m] == 
                      4) * barLen, code = 3, angle = 90)
                }
                if (l2r == TRUE) {
                  arrows(-uSR$dates[[n]][m], pP$xx[n], -uSR$dates[[n]][m + 
                    1], pP$xx[n], col = clrI, lwd = lwdI, lty = ltyI, 
                    length = (uSR$types[[n]][m] == 4) * barLen, 
                    code = 3, angle = 90)
                }
            }
        }
        par(srt = 90 * (!l2r))
        x0s <- y0s <- x1s <- cls <- hlts <- hlds <- c()
        prnts <- as.character(unique(pP$pn[!is.na(pP$pn)]))
	
	#*#*#get the formats for links between lineages
      for (i in 1:length(prnts)) {
            offsprng <- which(pP$pn == prnts[i])
            for (j in 1:length(offsprng)) {
            #make this part more computationally efficient
               x0s <- c(x0s, pP$xx[which(pP$nm == prnts[i])])
               x1s <- c(x1s, pP$xx[offsprng[j]])
               y0s <- c(y0s, pP$st[offsprng[j]])
               hclO <- uSR$styles[[uSR$types[[offsprng[j]]][1]]][2]
               #hclO <- uSR$styles[[uSR$types[[which(pP$nm == prnts[i])]][1]]][2]
               cls <- c(cls, hclO)
               #hltO <- uSR$styles[[uSR$types[[which(pP$nm == prnts[i])]][1]]][3]
               hltO <- uSR$styles[[uSR$types[[offsprng[j]]][1]]][3]
               hlts <- c(hlts, hltO)
               #hldO <- uSR$styles[[uSR$types[[which(pP$nm == prnts[i])]][1]]][1]
               hldO <- uSR$styles[[uSR$types[[offsprng[j]]][1]]][1]
               hlds <- c(hlds, hldO)
            }
        }
        clCd <- sapply(1:length(pP$nm), function(i) uSR$styles[[min(uSR$types[[i]])]][2])
        if (!is.null(hlty)) hlts <- rep(hlty, length(pP$nm))
        #*#*#link up lineages by joining ancestors and add tip and node labels
        if (!l2r) {
            segments(x0s, -y0s, x1s, -y0s, col = "grey70", lwd = .5, 
                lty = as.numeric(hlts))
            text(x = pP$xx, y = -pP$en, pP$label, adj = 0, font = 3, 
                cex = cexLab, col = clCd)
        }
        if (l2r) 
          {
          	segments(-y0s, x0s, -y0s, x1s, col = cls, lwd = as.numeric(hlds), lty = as.numeric(hlts))
          	text(x = -pP$en, y = pP$xx, pP$label, adj = 0 + (min(pP$st < 0)), font = 3, cex = cexLab, col = clCd)
           }
        	#*#*#add point occurences if applicable
        	if (!usrFlg)
          	{
           rdPtXX <- which(sapply(1:length(uSR$dates), function(i) which(diff(uSR$dates[[i]]) == 0)) > 0)
           rdPt <- unlist(sapply(1:length(uSR$dates), function(j) uSR$dates[[j]][which(diff(uSR$dates[[j]]) == 0)]))
           ifelse(l2r, xx <- -rdPt, xx <- pP$xx[rdPtXX])
           ifelse(l2r, yy <- pP$xx[rdPtXX], yy <- -rdPt)
           if (length(xx) != 0) arrows(xx, yy, xx, yy, length = 0, lwd = uSR$styles$point[1])
           }
        if (length(grep("c", addTimeLine)) > 0) close.screen(all = TRUE)
        #par(op)
    }
}
