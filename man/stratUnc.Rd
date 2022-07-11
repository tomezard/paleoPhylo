\name{stratUnc}
\alias{stratUnc}

\title{Include stratigraphic uncertainty in paleoPhylo plots.}
\description{
Incorporate different types of stratigraphic and mathematical uncertainty in paleoPhylo plots.  By default, these are envisaged as excellent resolution, range extensions, (dis)appearances of a lineage within part of its overall range, single point occurrences and confidence intervals around first and last appearance dates.  Full flexibility is possible however.}
\usage{
stratUnc(uSR=NULL, pP, lwdLin=1, ltyLin=1, clCd=1:length(unique(pP$grp)), style=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{uSR}{A list (see below).  Defaults to NULL.  If uSR is NULL, then a default list is created for pP which assumes that there is no stratigraphic uncertainty.}
\item{pP}{An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}). }
\item{lwdLin}{Base line thickness (is recycled).}
\item{clCd}{Defaults to base colors per group (see \code{\link{as.paleoPhylo}}).}
\item{ltyLin}{defaults to 1; base line type (is recycled).}
\item{style}{A vector of default style types.}
}
\details{
If a non-null data-frame is passed to \code{stratUnc}, this function checks its suitability.  Specifically, are all types and dates congruent and are all types and styles congruent?  The data-frame should contain the unique ID code for the individuals along with start and end dates for the uncertainty.  

The default settings are intended to reflect various types, but are fully editable by the end-user.  In the details below, lineage letter codes relate to the example below.  These types are excellent resolution (lineage D), range extensions (also referred to as ghost lineages, lineage A), (dis)appearances of a lineage within part of its overall range (lineage E) and single point occurrences (lineage B).  Furthermore, confidence intervals around first and last appearance dates (e.g. Strauss and Sadler 1989, Marshall 1994, Solow 1996, Wagner 2000, Wang and Marshall 2004) can be included (lineage C in fig. 2).  Note that calculation of confidence intervals is not yet implemented.

Any number of uncertainties can be appended within a lineage, although point occurrences can, by definition, only occur once.  
}

\value{
A list containing:
\item{ids}{Another list of unique ID codes per lineage.}
\item{dates}{A list of vectors per lineage for the timings of divisions between the different uncertainty types.}
\item{types}{A list of vectors per lineage for the types.  Has length == length(dates)[[x]]-1 for lineage x.}
\item{styles}{The style types, which are vectors per type containg default line thickness (can be adapted in \code{\link{drawPhylo}}, color and line type.  See \code{\link{par}}.}
}

\references{
Ezard, T.H.G. & Purvis, A. (2009). paleoPhylo: free software to draw paleobiological phylogenies. Paleobiology 35, 460-464. doi: 10.1666/0094-8373-35.3.460.

To calculate confidence intervals on lineage duration, see (for example):

Marshall, C. R. 1994. Confidence intervals on stratigraphic ranges: partial relaxation of the assumption of randomly distributed fossil horizons. Paleobiology 20:459-469.

Solow, A. R. 1996. Tests and confidence intervals for a common upper endpoint in fossil taxa. Paleobiology 22:406-410.

Strauss, D., and P. M. Sadler. 1989. Classical confidence intervals and Bayesian probability estimates for ends of local taxon ranges. . Mathematical Geology 21:411-427.

Wagner, P. J. 2000. Likelihood tests of hypothesized durations: determining and accommodating biasing factors. Paleobiology 26:431-449.

Wang, S. C., and C. R. Marshall. 2004. Improved confidence intervals for estimating the position of a mass extinction boundary. Paleobiology 30:5-18.


}
\author{Thomas H.G. Ezard \email{t.ezard@soton.ac.uk}}

\seealso{\code{\link{createBifurcate}}, \code{\link{drawPhylo}}, \code{\link{as.paleoPhylo}}}

\examples{
op <- par(no.readonly=TRUE)
eg1 <- data.frame(nm=paste("id",1:8,sep=""),
	pn=c("id6","id7","id7","id5","id8",NA,"id8","id6"),
	st=c(3.1,2,2,1,3,3.2,3,3.1),
	en=c(0,1.1,0.8,0,0,3.1,2,3),
	label=c("E","B","A","C","D","","",""))
pP <- with(eg1, as.paleoPhylo(nm, pn, st, en,label=label))

par(mfrow=c(1,2))
##example 1; no uncertainty
drawPhylo(pP, 
  cexLab=1.5, addTimeLine="n", cexText=1.5, hlty=3, l2r=TRUE)

##example 2; hypothetical uncertainty
##this is figure 2 in Ezard & Purvis (2009).
unc <- data.frame(id=c("id1","id1","id3","id4","id4","id2","id2"),
	st=c(3.1,1,2,1.1,0.1,2,1.7),
	en=c(2,0,1.1,0.9,0,1.1,1.7),
	type=c(2,2,3,4,4,2,5))
uSR <- stratUnc(unc,pP)
drawPhylo(pP, uSR, 
  cexLab=1.5, addTimeLine="n", cexText=1.5, hlty=3, l2r=TRUE)

par(op)
}

