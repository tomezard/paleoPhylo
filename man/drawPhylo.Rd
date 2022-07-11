\name{drawPhylo}
\alias{drawPhylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a lineage phylogeny.}
\description{Generates a cladogram that can distinguish between budding and bifurcating speciations.  Hypothetical relationships and stratigraphic uncertainty can be incorporated using \code{\link{stratUnc}}.  Further details on the distinction can be found in \code{\link{as.paleoPhylo}}.  Additional arguments enable clipping of time axes in variable formats and from different geological scales to the main plot.}
\usage{
drawPhylo(pP, 
	uSR = NULL,
	addTimeLine = "none", tmScl, whatTime,
	l2r=FALSE, nmLim=2, cexText = 0.5, srtText=0, cexTime = NULL, cexLab=0.5,
	lwdLin = 2, hlty=NULL, barLen=0, whSpc=0.05, dumpLast=FALSE, sz=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pP}{An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}). }
  \item{uSR}{A list generated through \code{\link{stratUnc}} to incorporate stratigraphic uncertainty.  Defaults to NULL. }
  \item{addTimeLine}{Add a time axis to the plot.  Four options (or only their first letter) are accepted: 
  	"none" (plots no time axis); 
  	"classic" (default, plots a set of rectangular boxes with labels on the left-hand side of the cladogram);
  	"tube" (plots grey and white rectangles underneath the cladogram, as London tube maps plot travelcard zones);
  	"arrow" (plots an arrow in the direction of the time axis)}
  \item{tmScl}{A time scale for the time axis.  The "Gradstein 2004" (see \code{\link{gradstein04}}) and "Berggren 1995" (see \code{\link{berggren95}}) Geologic Time Scales are included, but this argument can take *any* user-read in data frame (see details below).  Required if addTimeLine is either classic or tube.}
  \item{whatTime}{A vector of column names to be plotted in the time axis on the plot.  Defaults to 'epoch' and 'age', as used in \code{\link{berggren95}}.}
  \item{l2r}{Should the time axis run from left to right or top to bottom?  Defaults to FALSE, i.e. top to bottom.}
 \item{nmLim}{A number of million years that serves as the limit for when the name of the period/epoch/stage/etc... should be printed in the 'classic' time axis.  Defaults to 2.}
  \item{cexText}{Change the text size of time labels.  See \code{\link{par}}. Can be a vector for different text sizes for each time classification (see example 2a below).}
  \item{srtText}{Change the text direction of time labels in the classic time axis.  See \code{\link{par}}.}
  \item{cexTime}{Change the text size of time numbers.  Defaults to the size of \code{cexText} used for the narrowest time intervals unless otherwise stated. See \code{\link{par}}.}
  \item{cexLab}{Change the text size of species lineage names.  See \code{\link{par}}.}
  \item{lwdLin}{Change the line thickness.  See \code{\link{par}}.}
  \item{hlty}{Change the horizontal line type if not otherwise specified using \code{\link{stratUnc}}.  Takes values for lty in \code{\link{par}}.}
 \item{barLen}{The length of the horizontal bars if confidence intervals are being drawn using \code{\link{stratUnc}}.}
 \item{whSpc}{Controls the amount of white space around the phlyogram.  Increasing values of \code{whSpc} will yield more white space, which may be useful when formatting smaller phylograms or if the tip labels are particularly long. It works using the \code{f} in \code{\link{extendrange}}.  Note that \code{mar()} in \code{\link{extendrange}} is not supported.}
  \item{dumpLast}{If \code{addTimeLine=="classic" or "tube"} then should the oldest time box be drawn? Removing it might be desired if, eg., an ancestor is shown from an earlier time period but only just which would mean "messy" unclosed time boxes.}
 \item{sz}{Controls the proportion of the plot given over to a 'classic' time axis. }
}
\details{
In essence, all this does is generate a blank window and put some lines on top.  As such, its use is not restricted to 'paleoPhylo' objects.  It can be used in conjunction with fully bifurcating trees (\code{\link{ape2paleoPhylo}} can convert an object of class 'phylo' to class 'paleoPhylo' using  \code{\link{as.paleoPhylo}}) to add time axes in multiple forms to plots, and in conjunction with the pdf graphics device (see \code{\link{pdf}}) to output high-quality figures in pdf format (other graphics formats are available; see links from the pdf page).

Currently, any user read-in time scale MUST have a column called "MA", which contains age estimates for each period/epoch/stage/etc....  See \code{\link{gradstein04}} or \code{\link{berggren95}} for ways of *ensuring* your data frame will read in correctly.

\code{addTimeLines=="classic"} uses \code{\link{split.screen}}.  The left-hand portion containing the axis has two fixed widths depending on the number of species lineages to be drawn.  In due course, its width will become fully dynamic.
}
\value{
The function opens a graphics window and plots the lineage phylogeny (see \code{\link{as.paleoPhylo}} for details).  It does not return a list, since all the information is contained within input data.}
\references{Ezard, T.H.G. & Purvis, A. (2009). paleoPhylo: free software to draw paleobiological phylogenies. Paleobiology 35, 460-464. doi: 10.1666/0094-8373-35.3.460.}
\author{Thomas H.G. Ezard \email{t.ezard@soton.ac.uk}}
\seealso{\code{\link{as.paleoPhylo}}, \code{\link{gradstein04}}, \code{\link{pearson93}}, \code{\link{berggren95}}, \code{\link{stratUnc}}}
\examples{
data(pearson93)
data(berggren95)
op <- par(no.readonly=TRUE)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))

##example 1: without a time axis.
drawPhylo (p93)

##example 2: with a "classic" time axis.
drawPhylo (p93,
  tmScl=berggren95, addTimeLine="classic", whatTime=c("epoch","zone"),
  l2r=TRUE, cexText=c(1, 0.5), dumpLast=TRUE,  nmLim=0.3)
	
##example 3: a fully bifucating version with different text sizes for time
p93cb <- createBifurcate(p93)
drawPhylo (p93cb, 
  tmScl=berggren95, addTimeLine="cl", whatTime=c("epoch","zone"),
  dumpLast=TRUE, cexText=c(1, 0.5), nmLim=0.3, srtText=0)
	
#example 4: with species arranged on the x-axis according ...
## .... to some random trait data from a uniform distribution.
## .... THIS DELIBERATELY LOOKS A TOTAL MESS
##don't include latin names by leaving the "label" argument free.
p93rnd <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End,
	 xx=runif(length(Code),-1,1)))
drawPhylo (p93rnd, tmScl=berggren95, addTimeLine="classic", whatTime=c("epoch"),
	l2r=TRUE, cexTime=0.8, dumpLast=TRUE, cexText=1, nmLim=0.3)
	
##THERE ARE MORE EXAMPLES IN  ape2paleoPhylo
par(op)

}

