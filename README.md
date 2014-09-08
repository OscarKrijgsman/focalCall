focalCall
=========

FocalCall: An R package for the annotation of focal copy number aberrations

To identify somatic focal copy number aberrations (CNAs) in cancer specimens and distinguish them from germ-line copy number variations (CNVs) we developed a software package named, FocalCall. FocalCall permits user-defined size cut-offs to recognize focal aberrations and builds on established array CGH segmentation and calling algorithms. To differentiate CNAs from CNVs the algorithm can either make use of a matched normal reference signal or, if this is not available, an external CNV track. FocalCall furthermore differentiates between homozygous and hemizygous deletions as well as gains and amplifications and is applicable to high-resolution array and sequencing data. 


FocalCall has been submitted to [Bioconductor](http://bioconductor.org/). 

# Requirements:
CGHcall and associated packages available from BioConductor [Wiel et al. 2011].

Installation of the package follows standard R installation guidelines:
UNIX/LINUX/MAC:
`> R CMD install focalCall_1.0.tar.gz`

# Running FocalCall 
Here we  use focalCall to analyse the dataset previously published by 
Bierkens et al. (2013) and preprocessed using CGHcall by Wiel et al. (2007). The example set used here only contains complete chromosome 2. 
For the other chromosomes only a small portion of the CGH probes are included. The 
complete dataset can be downloaded from the NCBI Gene Expression Omnibus (GEO), 
accession number GSE34575.

	> library(focalCall)
	> data(BierkensCNA)
	> calls_focals<-focalCall(CGHset, CNVset, focalSize = 3, minFreq=2)
	> FreqPlot(calls_focals, header="FrequencyPlot all aberrations")
	> FreqPlotfocal(calls_focals, header="FrequencyPlot all aberrations")
	> igvFiles(calls_focals)

An R-vignette with more details on how to analyse this dataset is available with the package. 

## Contact

We have tried to make the FocalCall readable and its use as easy as possible. If any questions arise regarding the package, or if you want to report any bugs, please do not hesitate and contact:

- [Oscar Krijgsman](mailto:oscarkrijgsman@gmail.com)

This package was developed in the lab of Prof Dr. Bauke Ylstra, [website] (http://www.vumc.nl/afdelingen/microarrays/).


# References:
van de Wiel M.A. et al. (2011) Preprocessing and downstream analysis of microarray DNA copy number profiles. Brief Bioinform., 12:10-21. 

Bierkens M. et al. (2013) Focal aberrations indicate EYA2 and hsa-miR-375 as oncogene and tumor suppressor in cervical carcinogenesis. Genes Chromosomes Cancer, 52:56-68.

Krijgsman O. et al. (2014) Focal chromosomal copy number aberrations in cancer - Needles in the genome haystack. BBA- molecular cell research.
