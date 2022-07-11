\name{createBifurcate}
\alias{createBifurcate}
\alias{buildApe}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a fully bifurcating object of class 'phylo'.}
\description{createBifurcate and buildApe take a 'paleoPhylo' object and convert it to its fully bifurcating equivalent (i.e. all ancestors terminate if they speciate). \code{createBifurcate} converts the \code{\link{as.paleoPhylo}} object into one with only bifurcating speciations; \code{buildApe} does the conversion to ape.}
\usage{
createBifurcate(pP)
buildApe(pP, label=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pP}{An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}).}
  \item{label}{Should the tip and node labels be labels (\code{TRUE}) or ID codes (\code{FALSE}).}
}
\details{
\code{createBifurcate} works by *walking* through the lineage phylogeny until a speciation event occurs.  If the speciation event is a budding-type (where the ancestor lineage persists simultaneously with the descendant), then \code{createBifurcate} splits the ancestor such that the speciation is a bifurcating-type (new, but , unique IDs are generated).  This process is repeated iteratively until extant (or most recently deceased) tips have been dealt with.  The goal is to enable easy conversion to class 'phylo' for use with \code{\link[ape]{read.tree}} formats.

\code{createBifurcate} does not work for phylogenies with polytomies since the fully-birfurcating version would require resolution of said polytomy.

\code{buildApe} converts this into an object of class 'phylo' for use in \code{ape} (e.g. \code{\link[ape]{read.tree}} or \code{\link[ape]{plot.phylo}}), i.e. assigning of node labels and branch lengths, edge matrix and root definitions. The root is usually the first row of the data frame and needs to be unique. Phylogenies containing more than 500,000 species will not \code{buildApe} properly.

\code{createBifurcate} and \code{buildApe} are kept separate to enable less manipulation before plotting fully bifurcating trees with arguments in \code{\link{drawPhylo}}.  
}
\value{
An object of class 'paleoPhylo' (see \code{\link{as.paleoPhylo}}).}

\references{Ezard, T.H.G. & Purvis, A. (2009). paleoPhylo: free software to draw paleobiological phylogenies. Paleobiology 35, 460-464. doi: 10.1666/0094-8373-35.3.460.}

\author{Thomas H.G. Ezard \email{t.ezard@soton.ac.uk}}

\seealso{\code{\link{as.paleoPhylo}}}
\examples{
data(pearson93)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))
apeTree <- buildApe(createBifurcate(p93))
apeTree

##using ape, do some simple illustratory stuff
##plot a fully bifurcating phylogeny using plot.phylo
plot(apeTree)
##write apeTree using the Newick (aka New Hampshire) format.
write.tree(apeTree)

}


