\name{focalCall}
\alias{focalCall}
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
	Two files are returned - an object of class \code{cghCall} and a text files with all focal aberrations detected in all sample. 
	\itemize{
	\item "focalCall.Rdata" is an object of class \code{cghCall} and contains all DNA copy number information including the focal aberrations.
	\item "focalList_sampleName.txt" is a text file with all focal aberrations detected in the set of sample.  
	}
}


\references{ 
Krijgsman O., Benner C., Meijer G.A., Van de Wiel, M.A., Ylstra, B. (2013), "FocalCall: an R-package to detect genomic focal aberrations.", \emph{Submitted}
}

\author{ Oscar Krijgsman: \email{o.krijgsman@vumc.nl} }
\seealso{ \code{focalCall}, \code{CGHcall} }

\examples{
# Generate object of cghCall class object and input files
data(BierkensCNA)

# Extract focal aberrations and distinguish somatic from germ-line copy number aberrations (CNA)
ExampleRun<-focalCall(CGHset, CNVset, focalSize=3, minFreq=2)

}
