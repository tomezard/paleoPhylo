`createBifurcate` <- function (pP) 
	{
	mergeSplit <- function (dt, spL) 	{dt <- dt[dt$cd != spL$cd[1],]  ;  dt <- dt[dt$cd != spL$cd[3],]  ;  dt <- rbind(dt,spL)}

	splitLineage <- function (prnt, off, nxtCd, ltCd)
		{
		ifelse (ltCd==1 ,  
			splName<-paste(substr(prnt$nm,1,nchar(as.character(prnt$nm))),LETTERS[ltCd],sep=""),  
			splName<-paste(substr(prnt$nm,1,nchar(as.character(prnt$nm))-1),LETTERS[ltCd],sep=""))
		sp1 <- data.frame(nm=as.character(prnt$nm),pn=as.character(prnt$pn),st=prnt$st,en=off$st,cd=prnt$cd,label=prnt$label)
		sp2 <- data.frame(nm=splName,pn=sp1$nm,st=off$st,en=prnt$en,cd=nxtCd,label=prnt$label)
		sp3 <- data.frame(nm=off$nm,pn=sp1$nm,st=off$st,en=off$en,cd=off$cd,label=off$label)
		spL <- data.frame(rbind(sp1,sp2,sp3))
		return(spL)
		}
	
	if(class(pP)!="paleoPhylo") stop("object is not of class 'paleoPhylo'")
		{
		dat <- with(pP, data.frame(nm=as.character(nm),pn=as.character(pn),st=st,en=en,cd=1:length(nm),label=label))
		Nms <- data.frame(nms=as.character(pP$nm),st=pP$st)  ;  nms <- as.character(Nms[rev(order(pP$st)),1])

		for (k in 1:length(unique(dat$nm)))
			{
			prnt <- spL <- dat[which(dat$nm == nms[k]),]
			off <- dat[c(which(as.character(prnt$nm)==as.character(dat$pn))),]
			off <- off[rev(order(off$st)),]
	
			if (length(off[,1]) != 0)
				{
				evnts <- unique(off$st)
				for (i in 1:length(evnts))
					{
					off1 <- off[off$st == evnts[i],]
					if (sum(off$st == evnts[i]) == 2) 
						{for (j in 1:length(off1$nm)) 
							{
							dat$pn <- as.character(dat$pn)
							dat$pn[which(as.character(dat$nm) == as.character(off1$nm[j]))] <- as.character(prnt$nm)}
							}
					if (sum(off$st == evnts[i]) == 1) 
						{
						off1$pn <- prnt$nm
						spL <- splitLineage(prnt,off1,max(dat$cd)+1,i)
						prnt <- spL[2,]
						dat <- mergeSplit (dat, spL)
						}
					}
				}
			}
		dat <- dat[rev(order(dat$st)),]
		if (sum(with(dat[dat$st < max(dat$st),],tapply(st,paste(st,pn),length))!=2) !=0) {cat("splitting of lineages to form bifurcating tree not exhaustive")}
		YY <- with(dat,as.paleoPhylo(nm,pn,st,en,label=label))
		return(YY)
		}	
	}