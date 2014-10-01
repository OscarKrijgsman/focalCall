\name{FreqPlot}
\alias{FreqPlot}
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
  FreqPlot is applied to the calls of an object of class \code{cghCall} (as returned by focalCall). 
}

\value{
 A single plot is returned with the frequency of aberrations in the complete dataset. The x-axis represents the chromosomes, the y-axis the percentage of samples with a gain or loss at that genomic position. 
}

\references{ 
Krijgsman O., Benner C., Meijer G.A., Van de Wiel, M.A., Ylstra, B. (2013), "FocalCall: an R-package to detect genomic focal aberrations.", \emph{Submitted}
}
\author{ Oscar Krijgsman: \email{o.krijgsman@vumc.nl} }
\seealso{ \code{focalCall}, \code{CGHcall} }
\examples{

# Due to runtime and space limitations example data only contains complete data for chromosome 2
data(BierkensCNA)

# generate frequency plot of copy number data
FreqPlot(calls_focals, "FrequencyPlot_BierkensSamples")

}
