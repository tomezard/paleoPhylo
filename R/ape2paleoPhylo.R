ape2paleoPhylo <- function(phy,retainNodeLabels=TRUE, nC=0, getloc=TRUE)
	{
	if(class(phy)!="phylo") stop("object is not of class 'phylo'")
		{
		if(is.null(phy$edge.length)) stop("No edge lengths in tree, so cannot calculate start and end points.")
			{
			nm <- phy$edge[,2]
			cd <- sapply(1:dim(phy$edge)[1], function(i) ifelse(!is.na(phy$tip.label[phy$edge[i,2]]), phy$tip.label[phy$edge[i,2]], phy$edge[i,2]))
			cd[nchar(cd)<nC] <- ""
			pn <- phy$edge[,1]
			en <- rep(NA,length(nm))
			for (k in 1:length(nm))
				{
				xx <- x<- nm[k]
				while (x != length(phy$tip.label) + 1) {phy$edge[which(phy$edge[, 2] == x),1] ; x <- phy$edge[which(phy$edge[, 2] == x),1] ; xx <- c(xx,x)}
				xxEdge <- matrix(c(xx[2:length(xx)],xx[1:(length(xx)-1)]),ncol=2)
				xxDist <- 0
				if(sum(is.na(xxEdge))==0)
					{en[k] <- sum(sapply(1:dim(xxEdge)[1], function(i) xxDist <- xxDist + phy$edge.length[which(phy$edge[,1]==xxEdge[i,1] & phy$edge[,2]==xxEdge[i,2])]))}
				if(sum(is.na(xxEdge))>0) {en[k] <- 0}
				}
			st <- en - phy$edge.length
			root <- as.character(length(phy$tip.label)+1)
			if(getloc) pP <- as.paleoPhylo(c(root,nm),c(NA,pn),abs(c(-.0001,st)-max(en)),abs(c(0,en)-max(en)),label=c(NA,cd))
			if(!getloc) pP <- as.paleoPhylo(c(root,nm),c(NA,pn),abs(c(-.0001,st)-max(en)),abs(c(0,en)-max(en)),xx=runif(length(c(root,nm))), label=c(NA,cd))
			return(pP)
			}
		}	
	}