\name{paleoPhylo-package}
\alias{paleoPhylo-package}
\alias{paleoPhylo}
\keyword{package}

\docType{package}
\title{Biodiversity Analyses in a Paleontological and Phylogenetic Context}
\description{
paleoPhylo contains a suite of functions to analyze how speciation, extinction & character change contribute to biodiversity as well as tools to display these relationships graphically.  The functions are designed to inform realistic models of clade evolution; examples are drawn from data on planktonic foraminifera or generated randomly.  The package is still in development.  It will be added to in due course - the last edit was in 2013, but I have now put it onto Github!
}
\details{
\tabular{ll}{
Package: \tab paleoPhylo\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-03-11\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Thomas H.G. Ezard and Andy Purvis.

Maintainer: Thomas H.G. Ezard \email{t.ezard@soton.ac.uk}
}

\references{
It would be really nice if you'd be kind enough to cite the manuscript! 

Ezard, T.H.G. & Purvis, A. (2009). paleoPhylo: free software to draw paleobiological phylogenies. Paleobiology 35, 460-464. doi: 10.1666/0094-8373-35.3.460.
}

\keyword{ package }
\seealso{Click on the index link below to see the diverse range of things that \code{paleoPhylo} does.}
\examples{
##plot a lineage phylogeny of planktonic foraminifera
data(pearson93)
data(berggren95)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))
drawPhylo(p93, tmScl=berggren95, whatTime=c("epoch","zone"), addTimeLine="c", 
  cexText=c(1, .5), l2r=TRUE, dumpLast=TRUE, whSpc=.1, nmLim=.8)

##convert the lineage phylogeny with budding and bifurcating speciation 
##	events to class "phylo" format with only bifurcating speciations.
apeTree <- buildApe(createBifurcate(p93))
apeTree

##fig. 2 in Ezard & Purvis (2009); incorporating stratigraphic uncertainty
eg1 <- data.frame(nm=paste("id",1:8,sep=""),
	pn=c("id6","id7","id7","id5","id8",NA,"id8","id6"),
	st=c(3.1,2,2,1,3,3.2,3,3.1),en=c(0,1.1,0.8,0,0,3.1,2,3),
	label=c("E","B","A","C","D","","",""))
pP <- as.paleoPhylo(eg1$nm,eg1$pn,eg1$st,eg1$en,label=eg1$label)
unc <- data.frame(id=c("id1","id1","id3","id4","id4","id2","id2"),
	st=c(3.1,1,2,1.1,0.1,2,1.7),
	en=c(2,0,1.1,0.9,0,1.1,1.7),
	type=c(2,2,3,4,4,2,5))
uSR <- stratUnc(unc,pP)
drawPhylo(pP,uSR,cexLab=1.5, addTimeLine="n", cexText=1.5, hlty=3, l2r=TRUE)

}
