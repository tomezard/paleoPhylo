\name{prune.to.date}
\alias{prune.to.date}

\title{Returns the phylogeny of the clade to some cut-off point.}

\description{
If passed an object of class "paleoPhylo" and a particular cutoff date, then \code{prune.to.date} returns the phylogeny of the clade to that point, either excluding (the default) or including the extinct lineages within the clade. The phylogeny is returned in two formats: paleoPhylo and phylo.
}

\usage{
prune.to.date(pP, date, 
  prune.dead.branches=TRUE, let.speciate = FALSE, let.die = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pP}{An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}).}
  \item{date}{A particular time in the past, to which the phylogeny is cut to.}
  \item{prune.dead.branches}{If TRUE, prune the whole history to \code{date} rather than only the then-extant phylogeny, i.e. do not prune the extinct lineages out. Takes precedence over \code{let.die}.}
  \item{let.speciate}{If TRUE, speciations at exactly the cutoff date are taken to have happened.}
  \item{let.die}{IF TRUE, extinctions at exactly the cutoff date are taken to have happened.}
  }
\details{
Enables an event-based protocol to be implemented, through iterative application along a phylogeny.
}

\value{
  \item{paleoPhylo.tree}{The pruned tree as an object of class "paleoPhylo".}
  \item{phylo.tree}{The pruned tree as an object of class "phylo".}
}

\references{Ezard, T.H.G. & Purvis, A. (2009). paleoPhylo: free software to draw paleobiological phylogenies. Paleobiology 35, 460-464. doi: 10.1666/0094-8373-35.3.460.}

\author{Andy Purvis \email{a.purvis@imperial.ac.uk}}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{as.paleoPhylo}}, \code{\link{pearson93}} }
\examples{
data(pearson93)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))

par(mfrow=c(2,4), mar=rep(1,4), cex=0.8)
eg1 <- prune.to.date(p93, 54.3)
plot(eg1$phylo.tree, main = "Default")
eg2 <- prune.to.date(p93, 54.3, prune.dead.branches=FALSE)
plot(eg2$phylo.tree, main = "Default, keeping the dead")
eg3 <- prune.to.date(p93, 54.3, let.speciate=TRUE)
plot(eg3$phylo.tree, main="Let speciate")
eg4 <- prune.to.date(p93, 54.3, prune.dead.branches=FALSE, let.speciate=TRUE) 
plot(eg4$phylo.tree, main="Let speciate, keeping the dead")
eg5 <- prune.to.date(p93, 54.3, let.die=TRUE)
plot(eg5$phylo.tree, main = "Let die")
eg6 <- prune.to.date(p93, 54.3, prune.dead.branches=FALSE, let.die=TRUE)
plot(eg6$phylo.tree, main = "Let die, keeping the dead")
eg7 <- prune.to.date(p93, 54.3, let.speciate = TRUE, let.die=TRUE)
plot(eg7$phylo.tree, main = "Live and let die")
eg8 <- prune.to.date(p93, 54.3, prune.dead.branches=FALSE, let.speciate = TRUE, let.die=TRUE)
plot(eg8$phylo.tree, main = "Live and let die, keeping the dead")

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), cex=1)
  ##RESTORE DEFAULT SETTINGS
}
