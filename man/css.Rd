\name{css}
\alias{css}

\title{Calculate Corrected Survivorship Scores.}
\description{Currently, calculates Pearson's (1992) corrected survivorship scores to account for variation in variable extinction rates.
}
\usage{
css(pP, lastDate=min(pP$en))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pP}{An object of class \code{paleoPhylo}.}
  \item{lastDate}{The date of the end of the experiment, for use in censoring in survivorship analysis. Defaults to the minumum age in the \code{paleoPhylo} object implying that EVERYTHING DIED.}
}
\details{
Pearson (1992) criticized survivorship analysis when conducted on the fossil record on the grounds that taxa arising at different periods must have experienced different extinction probabilities  and therefore cannot be considered homogeneous, implying that they are inappropriate for survival analysis
}

\value{
An list containing:
 \item{nm}{A unique identifying ID code for each species lineage, ripped from \code{\link{as.paleoPhylo}}.}
 \item{css}{The corrected survivorship score (Pearson 1992).}
 \item{Duration}{The raw duration}
 \item{notExtant}{Censorship if \code{nm} persists until \code{lastDate}.}
 \item{notExtC}{Censorship if \code{nm} persists until \code{lastDate} or speciates.}
}

\references{
Pearson, P. (1992). Survivorship analysis of fossil taxa when real-time extinction rates vary: the paleogene planktonic foraminifera. Paleobiology 18 (2), 115-131. 
}

\author{Thomas H.G. Ezard \email{t.ezard@soton.ac.uk}}

\seealso{\code{\link{pearson93}}, \code{\link{as.paleoPhylo}}}

\examples{
data(pearson93)
p93 <- with(pearson93, as.paleoPhylo(Code, Parent, Start, End, label=Name))
p93css <- css(p93)
with(p93css, plot(css, Duration))

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
##NOT RUN
##To do a survivorship analysis, use 
##library(survival)
##m1 <- summary(survfit(Surv(css,notExtant) ~ 1, data=p93css))
##plot(m1)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
}

