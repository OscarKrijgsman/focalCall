# Example script for the detection of focal chromosomal aberrations with FocalCall

# FocalCall: An R package for the annotation of focal copy number aberrations
# Krijgsman et al. 

library(GEOquery)
library(CGHcall)
library(focalCall)

## Download data from ncbi's Gene Expression Omnibus (GEO)
## This data was published by Bierkens et al. in 2013
## FocalCall was used to identify genes on focal aberrations of which multiple were validated in this study.
Bierkens_CGH<-getGEO("GSE34575")



## Pre-processing using CGHcall
## This pipeline was previously described by Wiel et al 2012. 

raw <- make_cghRaw(stanford)
prep <- preprocess(raw, maxmiss = 30, nchrom = 23)
nor <-  normalize(prep,method = "median", smoothOutliers = TRUE)  

seg <-  segmentData(nor, method = "DNAcopy",nperm=2000,undo.splits="sdundo",min.width=5,undo.SD=3,
clen=25, relSDlong=5)
segnorm <- postsegnormalize(seg,inter=c(-1,1))
	
listcalls <- CGHcall(segnorm,nclass=5,robustsig="yes",cellularity=1,ncpus=2)
calls <- ExpandCGHcall(listcalls,segnorm, divide=5, memeff=FALSE)
save(calls, file="calls.Rdata") 

## Analysis with focalCall
focalCall(calls)