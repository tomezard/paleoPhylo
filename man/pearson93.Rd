\name{pearson93}
\alias{pearson93}
\docType{data}
\title{An example lineage phylogeny.}
\description{Pearson (1993) published a lineage phylogeny for the Paleogene planktonic framinifera, which grouped 800 morphospecies into six principal clades and 134 lineages.  This data frame contains one of these clades (Cancellate spinose) and is included primarily as example data. A lineage phylogeny differs from a *normal* phylogeny in that speciation events can do not necessarily result in the termination of the ancestral lineage.}
\usage{data(pearson93)}
\format{
  A data frame containing details on 33 species lineages.
  \describe{
    \item{\code{Name}}{Names of the 33 species lineages within the focal clade.}
    \item{\code{Code}}{Unique ID codes for the 33 species.}
    \item{\code{Start}}{Date of first apperance of each species lineage.}
    \item{\code{End}}{Date of last apperance of each species lineage.}
    \item{\code{Parent}}{Name of ancestor for each species lineage.}
  }
}

\source{
Pearson, P. N. 1993 A lineage phylogeny for the Paleogene planktonic foraminifera. Micropaleontology 39, 193-222.

Pearson, P. N., R. K. Olsson, B. T. Huber, C. Hemleben, and W. A. Berggren. 2006. Atlas of Eocene Planktonic Foraminifera. Cushman Foundation for Foraminiferal Research.
}

\seealso{\code{\link{as.paleoPhylo}}, \code{\link{drawPhylo}}, \code{\link{berggren95}}}
\examples{
data(pearson93)
str(pearson93)

}
\keyword{datasets}
