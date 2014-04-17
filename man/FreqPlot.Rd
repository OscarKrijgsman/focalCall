\name{FreqPlot.Rd}
\alias{FreqPlot.Rd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{frequency plot copy number aberrations.}
\description{
Frequency plot of all aberrations in the dataset.  
}
\usage{
  FreqPlot(calls, header)
}
\arguments{
  \item{calls}{ Object of class \code{focalCall}, contains focalCall output.}
  \item{header}{ String with header name of frequencyPlot.}
}
\details{
  FreqPlot is applied to the calls of an object of class \code{cghCall} (as returned by focalCall version 0.0.99 or higher). 
}
\value{
  See \code{plot}.
}
\references{ 
Krijgsman O., Benner C., Meijer G.A., Van de Wiel, M.A., Ylstra, B. (2013), "FocalCall: an R-package to detect genomic focal aberrations.", \emph{Submitted}
}
\author{ Oscar Krijgsman: \email{o.krijgsman@vumc.nl} }
\seealso{ \code{focalCall}, \code{CGHcall} }
\examples{

# generate frequency plot of copy number data
\dontrun{FreqPlot(calls, header)}

}
