prune.to.date <-
function(pP, date, prune.dead.branches=TRUE, let.speciate=FALSE, let.die=FALSE){
	#Is passed a paleoPhylo object and a cutoff date
	#If prune.dead.branches, returns the phylogeny of just the species extant at the cutoff
	#If !prune.dead.branches, returns phylogeny of all species that started by the cutoff
	#The phylogeny is returned in two formats: paleoPhylo and phylo.
	#If let.speciate = TRUE, speciations at exactly the cutoff date are taken to have happened.
	#If let.die = TRUE, extinctions at exactly the cutoff date are taken to have happened.
	#let.die is only relevant if prune.dead.branches=TRUE

	#check it's a paleoPhylo object
	if (class(pP) !="paleoPhylo")
		stop("Object passed to prune.to.date is not of class paleoPhylo!")
	
	dF <- with(pP, data.frame(nm, pn, st, en, xx, label, grp))
	
	#Get rid of lineages that started too late
	if (let.speciate) pruned <- dF[dF$st>=date,]  else  pruned <- dF[dF$st>date,]

	if (length(pruned$nm) == 0){
		warning(paste("No lineages at time", date))
		return(NULL)
	}
	
	if (length(pruned$nm) == 1){
		warning(paste("Only one lineage at time", date))
		return(NULL)
	}

	if (let.die) labels.to.keep <- pruned$label[pruned$en < date] else labels.to.keep <- pruned$label[pruned$en <= date]
	#The line above is needed because label might not be unique; some instances of it might meet the criteria for being
	#dropped, but there might be one extant at the time of interest.
	
	lis <- with(pruned, as.paleoPhylo(nm, pn, st, en, xx, label, grp))
	cb <-createBifurcate(lis) #Create structures that will let extinct lineages be identified unambiguously
	cb$en[cb$en < date] <- 0 #Pad out all terminal edges to the present day
	tmp <- buildApe(cb, label=TRUE) #Produce corresponding phylo object
	tmpTree <- write.tree(tmp)
	pruned.ape <- read.tree(, text= tmpTree)
	
	#Identify lineages to be pruned (if there are any) and prune them from the tree
	extinct.by.cutoff <- labels.to.drop <- NULL
	pruned.extant.ape <- pruned.ape #Initiate phylo object for selected taxa, prior to any pruning
	if (prune.dead.branches){
		extinct.by.cutoff <- cb$nm[cb$en > date]
		if (let.die) extinct.by.cutoff <- c(extinct.by.cutoff, cb$nm[cb$en == date & !is.element(cb$nm, dF$pn)])
		labels.to.drop <- cb$label[cb$nm %in% extinct.by.cutoff]
		labels.to.drop <- setdiff(labels.to.drop, labels.to.keep) #labels.to.keep overrides labels.to.drop, because labels might not be unique
		to.drop <- (intersect(pruned.ape$tip.label, labels.to.drop))
		if (length(to.drop)>0) pruned.extant.ape <- drop.tip(pruned.extant.ape,to.drop)
	}

	#Now go back to paleoPhylo to trim terminal edges to make tree ultrametric and ending at cutoff
	app <- ape2paleoPhylo(pruned.extant.ape, retainNodeLabels=TRUE, nC=0)
	app$en[app$en < date] <- date
	
	#Make ultrametric phylo tree
	pruned.extant.ape2 <-reorder(buildApe(app))
		
	#Set up and fill return structure
	to.return <- list(app, pruned.extant.ape2)
	names(to.return)<- c("paleoPhylo.tree","phylo.tree")
	
	return(to.return)
}
