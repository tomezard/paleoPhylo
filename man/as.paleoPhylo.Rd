\name{as.paleoPhylo}
\alias{as.paleoPhylo}
\alias{getXloc}

\title{The 'paleoPhylo' class.}
\description{
A class to contain core information on a lineage phylogeny.
}
\usage{
as.paleoPhylo(nm, pn, st, en, xx=NA, label=nm, grp=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nm}{A unique identifying ID code for each species lineage.}
  \item{pn}{A unique identifying ID code for the ancestor of each species lineage.  Only the root is currently permitted "NA".}
  \item{st}{The date of first appearance for each species lineage.}
  \item{en}{The date of last appearance for each species lineage.}
  \item{xx}{An optional vector of some morphological trait, for use in plotting (see \code{\link{drawPhylo}}).  Defaults to NA, which is overwritten using \code{getXloc}.}
  \item{label}{An optional vector of labels, for use in annotating plots (see \code{\link{drawPhylo}}).  Defaults to the values given for nm.}
  \item{grp}{An optional vector of groups, should different (morpho)groups wish to be denoted differently (see \code{\link{as.paleoPhylo}}).}
}
\details{
Provided to ease use for various data formats.  Other functions in the package work with this class to do the interesting and/or useful bits.  For example: convert to a fully bifurcating tree (see \code{\link{createBifurcate}}) of class 'phylo' (see \code{as.phylo} or \code{write.tree} in  the \code{ape} package) or plot lineage cladograms (see \code{getXloc} and \code{\link{drawPhylo}}).	
 
\code{getXloc} walks iteratively through the lineage phylogeny from root to tip.  It allocates positions along the x-axis for each species lineage such that when plotted (see \code{\link{drawPhylo}}), they are evenly spaced.  If \code{xx} is not assigned, then \code{getXloc} is called.  Alternatively, any continuous morphological trait (a mixture of positive and negative values can be accomodated) can be used. To ensure that equalities are resolved, the dates are rounded to 8 decimal places in \code{getXloc}.
}

\value{
An object of class 'paleoPhylo', which is a list containing the five vectors defined above.}

\references{
Pearson, P. N. 1993 A lineage phylogeny for the Paleogene planktonic foraminifera. Micropaleontology 39, 193-222.
}

\author{Thomas H.G. Ezard \email{t.ezard@soton.ac.uk}}

\seealso{\code{\link{pearson93}}, \code{\link{createBifurcate}}, \code{\link{drawPhylo}}, \code{\link{ape2paleoPhylo}}}

\examples{
data(pearson93)
##not-defining xx leads to evenly-spaced positions along the x-axis in drawPhylo
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))
str(p93)

}

