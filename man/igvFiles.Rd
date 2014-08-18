\name{igvFiles.Rd}
\alias{igvFiles.Rd}
\title{Match Ensembl genes to focal aberrations.}
\description{
Generates three tracks that can be loaded into IGV ().  
1) Segmented data per sample 2) frequency plot based on the calls from CGHcall 
3) Frequency plot of focal aberrations as generated with focalCall.
}
\usage{
igvFiles(CGHset)
}
\arguments{
  \item{CGHset}{ Object of class \code{focalCall}, contains focalCall output.}
}
\details{
  igvFocal is applied to a  (as returned by focalCall version 0.0.99 or higher). 
}
\value{
  See \code{plot}.
}
\references{ 
Krijgsman O., Benner C., Meijer G.A., Van de Wiel, M.A., Ylstra, B. (2013), "FocalCall: an R-package to detect genomic focal aberrations.", \emph{Submitted}
}
\author{ Oscar Krijgsman: \email{o.krijgsman@nki.nl} }
\seealso{ \code{focalCall}, \code{CGHcall} }
\examples{

# generate IGV readable files
\dontrun{igvFiles(CGHset)}
}
