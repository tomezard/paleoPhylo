\name{rwb}
\Rdversion{1.1}
\alias{rwb}
\alias{rates.within.bins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate per-lineage origination and extinction rates within a set of bins
}
\description{
Takes any object of class \code{paleoPhylo} and computes per-lineage origination and extinction rates within the set of bins specified by the user.  Lineages persisting to the specified end are not considered to have gone extinct. Events taking place exactly on bin boundaries are taken to occur in the later bin.}
\usage{
rwb(pP, bL=1, st=max(pP$st), en=min(pP$en))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pP}{An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}).}
  \item{bL}{The length of bin; defaults to 1.}
  \item{st}{Start date for binning; defaults to start of the \code{\link{as.paleoPhylo}} object.}
  \item{en}{End date for binning;  defaults to end of the \code{\link{as.paleoPhylo}} object.}
}
\details{
Can be used to calculate time-series of rates for correlation analysis, etc. The function identifies lineages that do something at or after the start of the bin, but before the end of the bin, meaning that lineages surviving to the specified end time are not viewed as having gone extinct)

rates.within.bins does the same thing, with slightly different variable names. This function is intended to be a faster replacement, but both are currently (version 1.0-80) still in the package. The big change is in the usage order; in rates.within.bins() it is pP, begin, end, bin.length)
}

\value{
rwb is a data.frame, containing the following columns:
  \item{binStart}{The start date of the bin.}
  \item{nStart}{Number of extant lineages present at the start of the bin.}
  \item{branchLength}{The branch length of the phylogeny within the bin.}
  \item{orig}{The number of originations in the bin.}
  \item{extn}{The number of extinctions in the bin.}
  \item{lambda}{Speciation rate across the bin.}
  \item{mu}{Extinction rate across the bin.}
  }

\references{Ezard, T.H.G. & Purvis, A. (2009). paleoPhylo: free software to draw paleobiological phylogenies. Paleobiology 35, 460-464. doi: 10.1666/0094-8373-35.3.460.}

\author{Thomas H.G. Ezard and Andy Purvis \email{t.ezard@soton.ac.uk}}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{rates.within.bins}}, \code{\link{as.paleoPhylo}}, \code{\link{pearson93}}}

\examples{
data(pearson93)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))

p93RWB <- rwb(p93, bL=5)
with(p93RWB, plot(-binStart, lambda, col=4, lwd=2, type='l',
  xlab="Start of bin (ma)", ylab="Rate", axes=FALSE))
with(p93RWB, lines(-binStart, mu, col=2, lwd=2))
axis(1, at=seq(-65,-25,5), label=abs(seq(-65,-25,5)))
axis(2, at=seq(0,0.25,.05))

}