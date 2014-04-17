\name{singleSample.Rd}
\alias{singleSample.Rd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Detection focal aberrations with 1 sample only}
\description{
Detection and plotting of focal aberrations and distinguish somatic from germ-line in a single sample. This function, as \code{focalCall} can be used for samples with and without patient matched normal data.
}
\usage{
  singleSample(CGHset, CNVset, focalSize=3, OverlapPerc=0.2)
}
\arguments{
  \item{CGHset}{ Object of class \code{cghCall}, contains CGHcall output of tumor data.}
  \item{CNVset}{ Object of class \code{cghCall}, contains CGHcall output of matched normal data or .bed file with CNV locations.}
  \item{focalSize}{ Threshold for the function \code{singleSample}. Sets max size of a focal aberration/ }
  \item{OverlapPerc}{ Threshold for the function \code{singleSample}. sets minimal overlap of focal aberration with CNVs from in CNdata.normal to be classified as germ-line./}
}

\details{
  singleSample is applied to the calls of an object of class \code{cghCall}.
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

# generate object of cghCall class object and input files
# Due to runtime and space limitations example data only contains complete data for chromosome 2
data(BierkensCNA)

# Extract focal aberrations and distinguish somatic from germ-line copy number aberrations (CNA) in a single sample
singleSample(CGHset[,1], CNVset, focalSize=3)

}
