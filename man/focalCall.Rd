\name{focalCall}
\alias{focalCall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Detection genomic focal aberrations }
\description{
Detection of genomic focal copy number aberrations from high-resolution arrayCGH data or shallow-Seq copy number data. 
}
\usage{
  focalCall(CGHset, CNVset , focalSize = 3, minFreq=2)
}
\arguments{
  \item{CGHset}{ Object of class \code{cghCall}, contains CGHcall output of tumor data.}
  \item{CNVset}{ Object of class \code{cghCall}, contains CGHcall output of matched normal data or .bed file with CNV locations.}
  \item{focalSize}{ Threshold for the function \code{focalCall}. Sets max size of a focal aberration/ }
  \item{minFreq}{ Threshold for the function \code{focalCall}. Sets minimal number of samples that contain the aberration/ }
}
\details{
  focalCall is applied to the calls of an object of class \code{cghCall} (as returned by CGHcall version 1.2.0 or higher). FocalCall extract all focal aberrations from the tumor copy number data and 
}
\value{
  A \code{list}-object with the following entries:
  \item{focal_List }{\code{data.frame} containing the annotation data of the detected focal aberrations. The first column contains the chromosome information, the second and third its start and end base pair, and the fourth column the number of features constituting the region.}
  \item{focalcall }{Object with focal and CNV data \code{cghCall}.}
}
\references{ 
Krijgsman O., Benner C., Meijer G.A., Van de Wiel, M.A., Ylstra, B. (2013), "FocalCall: an R-package to detect genomic focal aberrations.", \emph{Submitted}
}
\author{ Oscar Krijgsman: \email{o.krijgsman@vumc.nl} }
\seealso{ \code{focalCall}, \code{CGHcall} }
\examples{

# generate object of cghCall class object and input files
# Due to runtime and space limitations example data only contains complete data for chromosome 2
data(BierkensCNA)

# Extract focal aberrations and distinguish somatic from germ-line copy number aberrations (CNA)
focalCall(CGHset, CNVset, focalSize=3, minFreq=2)

}
