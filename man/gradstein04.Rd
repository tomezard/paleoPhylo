\name{gradstein04}
\alias{gradstein04}
\docType{data}
\title{Gradstein et al.'s 2004 Geologic Time Scale.}
\description{A data frame that gives the geologic time scale published in 2004 by Gradstein and co-workers.  Most column names are self-explanatory}
\usage{data(gradstein04)}
\format{
 A data frame containing the following columns: 
  \describe{
    \item{\code{eon}}{Eon name.}
    \item{\code{era}}{Era name.}
    \item{\code{period}}{Period name.}
    \item{\code{epoch}}{Epoch name.}
    \item{\code{stage}}{Stage name.}
    \item{\code{MA}}{Estimated boundary age for the associated interval.}
    \item{\code{error}}{Estimated errors associated with each age estimate.}
    \item{\code{GSSP}}{Binary response denoting whether the age estimate is defined by a basal Global Standard Section and Point}
  }
}
\details{
Full details on derivation of the time scale are given in the references below.  The following excerpt is the ad verbatim legend from the published figure:
"Subdivisions of the global geologic record are formally defined by their lower boundary. Each unit of the Phanerozoic (~542 Ma to Present) and the base of Ediacaran are defined by a basal Global Standard Section and Point (GSSP), whereas Precambrian units are formally subdivided by absolute age (Global Standard Stratigraphic Age, GSSA). Details of each GSSP are posted on the ICS website \url{www.stratigraphy.org}.

Numerical ages of the unit boundaries in the Phanerozoic are subject to revision. Some stages within the Cambrian will be formally named upon international agreement on their GSSP limits. Most sub-Series boundaries (e.g., Middle and Upper Aptian) are not formally defined."

Versions of R older than 2.8 on Windows contain a bug (this data-frame is erroneously titled gst).  Either substitute gradstein04 in examples for gst, or update to the latest release.
}
\source{
Gradstein F.M., Ogg J.G., Smith A.G. 2004 A Geologic Time Scale 2004. Cambridge University Press, Cambridge, UK.
\url{www.cambridge.org/uk/catalogue/catalogue.asp?isbn=0521786738}

Ogg, J.G., Ogg, G. and Gradstein, F.M. 2008 The Concise Geologic Time Scale.  Cambridge University Press, Cambridge, UK. 
\url{www.cambridge.org/catalogue/catalogue.asp?isbn=0521898498}

Gradstein F.M. and Ogg J.G. Geologic Time Scale 2004: Why, how and where next!
\url{www.stratigraphy.org/scale04.pdf}
}

\seealso{\code{\link{drawPhylo}}}

\examples{
data(gradstein04)
str(gradstein04)
}
\keyword{datasets}
