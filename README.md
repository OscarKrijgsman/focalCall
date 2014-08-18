focalCall
=========

FocalCall: An R package for the annotation of focal copy number aberrations

To identify somatic focal copy number aberrations (CNAs) in cancer specimens and distinguish them from germ-line copy number variations (CNVs) we developed a software package named, FocalCall. FocalCall permits user-defined size cut-offs to recognize focal aberrations and builds on established array CGH segmentation and calling algorithms. To differentiate CNAs from CNVs the algorithm can either make use of a matched normal reference signal or, if this is not available, an external CNV track. FocalCall furthermore differentiates between homozygous and hemizygous deletions as well as gains and amplifications and is applicable to high-resolution array and sequencing data. 

The tools described here will be made available through the Bioconductor website after final revisions are approved. 


The script Example_run.R describes the complete analysis of a dataset available in GEO (GSE34575) and published by Bierkens et al. (2013). 


# Requirements:
CGHcall and associated packages available from BioConductor [Wiel et al. 2011].

Installation of the package follows standard R installation guidelines:
UNIX/LINUX/MAC:
`> R CMD install focalCall_0.0.99.tar.gz`


# References:
van de Wiel M.A. et al. (2011) Preprocessing and downstream analysis of microarray DNA copy number profiles. Brief Bioinform., 12:10-21. 

Bierkens M. et al. (2013) Focal aberrations indicate EYA2 and hsa-miR-375 as oncogene and tumor suppressor in cervical carcinogenesis. Genes Chromosomes Cancer, 52:56-68.

Krijgsman O. et al. (2014) Focal chromosomal copy number aberrations in cancer - Hot needles in the genome stack. BBA- molecular cell research.