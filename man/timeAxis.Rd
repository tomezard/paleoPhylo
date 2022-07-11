\name{timeAxis}
\alias{timeAxis}

\title{
Plot a time axis to append to, e.g., diversity curve plots.
}

\description{
Works in the same way as the \code{addTimeLine} part of \code{\link{drawPhylo}}, only ....
}
\usage{
timeAxis(tmScl, whatTime, addTimeLine = "classic", l2r = FALSE, nmLim = 2, cexText = 0.5,
  srtText = 0, cexTime = NULL, cexLab = 0.5, whSpc = 0.01, lwdLin = 1, dumpLast = FALSE, sz=0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tmScl}{The time scale for the time axis.}
  \item{whatTime}{Which columns of the time scale should be plotted}
  \item{addTimeLine}{Two options (or only their first letter) are accepted: 
  	"classic" (default, plots a set of rectangular boxes with labels on the left-hand side of the cladogram);
  	"tube" (plots grey and white rectangles underneath the cladogram, as London tube maps plot travelcard zones).}
  \item{l2r}{Should the time axis run from left to right or top to bottom?  Defaults to FALSE, i.e. top to bottom.}
 \item{nmLim}{A number of million years that serves as the limit for when the name of the period/epoch/stage/etc... should be printed in the 'classic' time axis.  Defaults to 2.}
  \item{cexText}{Change the text size of time labels.  See \code{\link{par}}. Can be a vector for different text sizes for each time classification (see example 2a below).}
  \item{srtText}{Change the text direction of time labels.  See \code{\link{par}}. Can be a vector for different rotations for each time classification. It probably makes sense to stick to either 0 (parallel to phylogeny) or 90 (perpendicular to phylogeny).}
  \item{cexTime}{Change the text size of time numbers.  Defaults to the size of \code{cexText} used for the narrowest time intervals unless otherwise stated. See \code{\link{par}}.}
  \item{cexLab}{Change the text size of species lineage names.  See \code{\link{par}}.}
 \item{whSpc}{Controls the amount of white space around the phlyogram.  Increasing values of \code{whSpc} will yield more white space, which may be useful when formatting smaller phylograms or if the tip labels are particularly long. It works using the \code{f} in \code{\link{extendrange}}.  Note that \code{mar()} in \code{\link{extendrange}} is not supported.}
  \item{lwdLin}{Change the line thickness.  See \code{\link{par}}.}
  \item{dumpLast}{If \code{addTimeLine=="classic" or "tube"} then should the oldest time box be drawn? Removing it might be desired if, eg., an ancestor is shown from an earlier time period but only just which would mean "messy" unclosed time boxes.}
  \item{sz}{Controls the proportion of the plot given over to a 'classic' time axis. Defaults to 0.2.}

}

\section{Warning }{need to use \code{close.screen()} after the plot if \code{addTimeLine="classic"} as \code{\link{split.screen}} is used to partition up the device.}

\references{Ezard, T.H.G. & Purvis, A. (2009). paleoPhylo: free software to draw paleobiological phylogenies. Paleobiology 35, 460-464. doi: 10.1666/0094-8373-35.3.460.}

\author{Thomas H.G. Ezard \email{t.ezard@soton.ac.uk}}
\seealso{\code{\link{drawPhylo}}, \code{\link{split.screen}}\code{\link{pearson93}}}


\examples{
##draw a diversity curve for pearson93 phylogeny
data(pearson93)
op <- par(no.readonly=TRUE)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))

data(berggren95)
tmScl <- berggren95
tmScl <- tmScl[tmScl$MA>=min(p93$en),]
tmScl <- tmScl[tmScl$MA<=max(p93$st),]

dvrstyTime <- seq(max(p93$st),min(p93$en), -0.2)
abndncTime <- numeric(length(dvrstyTime))
for (k in 1:length(dvrstyTime))
	{abndncTime[k] <- sum(p93$st>dvrstyTime[k] & p93$en<=dvrstyTime[k])}


timeAxis(tmScl, c("epoch", "zone"), l2r=TRUE, dumpLast=FALSE, 
  cexText=c(1, 0.5), srtText=c(0,90, 0), nmLim=0.3, whSpc=0.12)
lines(-dvrstyTime ,abndncTime/max(abndncTime), lwd=2)

#currently the y-axis scales from 0 to 1, hence the use of proportions here for the diversity curve
abdnLab <- round(seq(0, max(abndncTime), length.out=5),0)
axis(2, at= seq(0,1,.25), labels=abdnLab, line=-3)
mtext("Abundance",2, line=-1, outer=FALSE, cex=1.5)
close.screen(2,all=TRUE)
##the 2 refers to the second panel, which is active and draws the (here) diversity curve

par(op)
}


