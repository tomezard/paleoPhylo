\name{berggren95}
\alias{berggren95}
\alias{bp05}
\docType{data}
\title{W.A. Berggren's 1995 geologic timescale}
\description{W.A. Berggren's 1995 geologic timescale for the Cenozoic, giving tropical to subtropical planktonic foraminiferal zonation.

Additionally, the Berggren and Pearson revision \code{bp05} for the Paleogene in the same format with the 1995 timescale for more recent times.}
\usage{
data(berggren95)
data(bp05)
}
\format{
 A data frame containing the following columns: 
  \describe{
    \item{\code{epoch}}{Epoch name.}
    \item{\code{zone}}{Tropical to subtropical planktonic foraminiferal zone code.}
    \item{\code{MA}}{Estimated boundary age for the associated interval.}
    \item{\code{epch}}{A handy shortened (in an \code{\link{nchar}} sense) form of epoch (which I use often).}
 }
}

\source{
Berggren, W.A., Kent, D.V., Swisher III, C.C. & Aubry, M.-P., 1995.	A Revised Cenozoic Geochronology and Chronostratigraphy, in Berggren, W.A., Kent, D.V., Swisher III, C.C., Aubry, M.-P. & Hardenbol, J., (eds), Geochronology, Time Scales and Global Stratigraphic Correlation: SEPM (Society for Sedimentary Geology), Special Publication No. 54, p. 129-212.

Berggren, W.A. & Pearson, P.N., 2005. A revised tropical to subtropical Paleogene planktonic foraminiferal zonation. Journal of Foraminiferal Research, 35, 279-298
}
 
\seealso{\code{\link{gradstein04}}, \code{\link{drawPhylo}}}

\examples{
data(berggren95)
str(berggren95)

data(bp05)
str(bp05)

}
\keyword{datasets}
