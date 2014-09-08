\name{Bierkens}
\docType{data}
\alias{Bierkens}
\title{Cervical cancer arrayCGH data}
\description{
  Example dataset with arrayCGH data of cervical pre-cancer lesions. OF NOTE: this set contains only 8 samples from the data and only the complete probeset for chromosome 2. Complete dataset is available from GEO (GSE34575). 
}	
}
\usage{Bierkens}
\format{

  This dataset contain two .Rdata files:
  
   \describe{
   \item{CGHset}{a CGHcall object containing  8 cervical pre-curses cancer lesions. These samples were taken from  Bierkens et al. 2013 and available under GEO accession number GSE34575.}
   \item{CNVset}{a bed file with CNV locations as provided by \url{http://dgv.tcag.ca/.}
  }

  CGHset: A CGHcall object containing the following 8 samples:

  \describe{
   \item{hgCIN.2.10_16}{high-grade Cervical Intra-epithelial Neoplasia sample: hgCIN.2.10_16.}
   \item{hgCIN.2.13_16}{high-grade Cervical Intra-epithelial Neoplasia sample: hgCIN.2.13_16.}
   \item{hgCIN.2.8_31}{high-grade Cervical Intra-epithelial Neoplasia sample: hgCIN.2.8_31.}
   \item{hgCIN.2.16_31}{high-grade Cervical Intra-epithelial Neoplasia sample: hgCIN.2.16_31.}
   \item{CIN3.1.2_16}{Cervical Intra-epithelial Neoplasia sample: CIN3.1.2_16.}
   \item{hgCIN.2.4_16}{high-grade Cervical Intra-epithelial Neoplasia sample: hgCIN.2.4_16.}
   \item{CIN3.1.8_16}{Cervical Intra-epithelial Neoplasia sample: CIN3.1.8_16.}
   \item{hgCIN.2.19_58}{high-grade Cervical Intra-epithelial Neoplasia sample: hgCIN.2.19_58.}
  }
  CNVset: A bed file containing the following columns:

  \describe{
   \item{CNV ID}{Identification ID for copy number variation (CNV)}
   \item{start}{Basepair position that indicates start of CNV}
   \item{end}{Basepair position that indicates end of CNV}
  }
}
\source{
 Bierkens, M., Krijgsman, O., Wilting, S.M., Bosch, L., Jaspers, A., Meijer, G.A., Meijer, C.J., Snijders, P.J., Ylstra, B., Steenbergen, R.D. (2013). Focal aberrations indicate EYA2 and hsa-miR-375 as oncogene and tumor suppressor in cervical carcinogenesis. \emph{Genes Chromosomes Cancer 52},56-68. 
}
\keyword{datasets}