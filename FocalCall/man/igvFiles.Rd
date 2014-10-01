\name{igvFiles}
\alias{igvFiles}
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
 Three IGV readable files are returned - Frequency aberrations, Frequency focal aberrations and all segmented values for all samples.
 \itemize{
    \item "FrequencyPlot.igv" is a IGV readable files and contains the frequency plot of all samples and all aberrations.
    \item "FrequencyPlotfocals.igv" is a IGV readable files and contains the frequency plot of all samples and all focal aberrations.
    \item "Overview_segments.igv" is a IGV readable files and contains the segmentation log2ratio values.
 }
}

\references{ 
Krijgsman O., Benner C., Meijer G.A., Van de Wiel, M.A., Ylstra, B. (2013), "FocalCall: an R-package to detect genomic focal aberrations.", \emph{Submitted}
}

\author{ Oscar Krijgsman: \email{o.krijgsman@nki.nl} }
\seealso{ \code{focalCall}, \code{CGHcall} }
\examples{

# Due to runtime and space limitations example data only contains complete data for chromosome 2
data(BierkensCNA)

# generate IGV readable files
igvFiles(calls_focals)
}
