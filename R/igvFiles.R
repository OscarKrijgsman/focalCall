# Generates three IGV files
# 1 - CBS file with segmentation values all samples
# 2 - Frequencyplot (all aberrations)
# 3 - Frequencyplot (focals only)

igvFiles<-function(CGHset){

### 1) Segmentation file
	IGV.cbs<-matrix(data=0, ncol=ncol(calls(CGHset))+4, nrow=nrow(calls(CGHset)))
	IGV.cbs[,1]<-paste("chr", chromosomes(CGHset), sep="")
	IGV.cbs[which(IGV.cbs[,1]=="chr23"),1]<-"chrX"
	IGV.cbs[,2]<-bpstart(CGHset)
	IGV.cbs[,3]<-bpend(CGHset)
	IGV.cbs[,4]<-featureNames(CGHset)
	IGV.cbs[,-c(1:4)]<-segmented(CGHset)
	
	# Insert .igv header for formatting
	tmp<-rep(0,ncol(IGV.cbs))
	tmp[1]<-"#type=COPY_NUMBER"
	IGV.cbs<-rbind(tmp, c("Chromosome","Start","End", "Feature", colnames(segmented(CGHset))), IGV.cbs)
	IGV.cbs[1,which(IGV.cbs[1,]==0)]<-""
	write.table(IGV.cbs, file="Overview_segments.igv", row.names=FALSE, quote=FALSE, append=FALSE, sep="\t", col.names=FALSE)
	cat("Generated .SEG file for loading in IGV", "\n")
	
### 2) Frequency plot (all aberrations)
	# Test length probes
	half.length<-round((bpend(CGHset)[1] - bpstart(CGHset)[1])/2)

	IGV.fp<-matrix(data=0, ncol=5, nrow=2*nrow(calls(CGHset)))
	IGV.fp[1:(nrow(IGV.fp)/2),1]<-paste("chr", chromosomes(CGHset), sep="")
	IGV.fp[((nrow(IGV.fp)/2)+1):nrow(IGV.fp),1]<-paste("chr", chromosomes(CGHset), sep="")	
	
	IGV.fp[1:(nrow(IGV.fp)/2),2]<-bpstart(CGHset)
	IGV.fp[((nrow(IGV.fp)/2)+1):nrow(IGV.fp),2]<-bpstart(CGHset)+(half.length+1)
	IGV.fp[1:(nrow(IGV.fp)/2),3]<-bpend(CGHset)-half.length
	IGV.fp[((nrow(IGV.fp)/2)+1):nrow(IGV.fp),3]<-bpend(CGHset)
	
	IGV.fp[1:(nrow(IGV.fp)/2),4]<-featureNames(CGHset)
	IGV.fp[((nrow(IGV.fp)/2)+1):nrow(IGV.fp),4]<-featureNames(CGHset)
	
	## Get frequencies
	Gains<-matrix(data=0, ncol=ncol(calls(CGHset)), nrow=nrow(calls(CGHset)))
	Losses<-matrix(data=0, ncol=ncol(calls(CGHset)), nrow=nrow(calls(CGHset)))
	Gains[which(calls(CGHset)>0)]<- 1
	Losses[which(calls(CGHset)<0)]<- 1
	IGV.fp[1:(nrow(IGV.fp)/2),5]<-round((rowSums(Gains)/ncol(calls(CGHset)))*100)
	IGV.fp[((nrow(IGV.fp)/2)+1):nrow(IGV.fp),5]<-round(((rowSums(Losses)/ncol(calls(CGHset)))*100)*-1)
	
	IGV.fp<-IGV.fp[order(as.numeric(gsub("chr","",IGV.fp[,1])), partial=as.numeric(IGV.fp[,2])),]
	IGV.fp[which(IGV.fp[,1]=="chr23"),1]<-"chrX"
	
	# Insert .igv header for formatting
	trackIDs<-matrix(data=0, ncol=5, nrow=2)
	trackIDs[1,1]<-"#type=OTHER"
	trackIDs[2,1]<-"#track color=150,0,0 altColor=0,0,150 autoScale=off graphType=bar scaleType=linear"
	IGV.fp<-rbind(trackIDs, c("Chromosome","Start","End","Feature", "Frequency"), IGV.fp)
	IGV.fp[1,which(IGV.fp[1,]==0)]<-""
	IGV.fp[2,which(IGV.fp[2,]==0)]<-""
	write.table(IGV.fp, file="FrequencyPlot.igv",row.names=FALSE, quote=FALSE, append=FALSE, col.names=FALSE, sep="\t")
	cat("Generated .IGV file for loading in IGV", "\n")
	
	
### 3) Frequency plot (focal aberrations)
	# Test length probes
	half.length<-round((bpend(CGHset)[1] - bpstart(CGHset)[1])/2)

	IGV.focal<-matrix(data=0, ncol=5, nrow=2*nrow(calls(CGHset)))
	IGV.focal[1:(nrow(IGV.focal)/2),1]<-paste("chr", chromosomes(CGHset), sep="")
	IGV.focal[((nrow(IGV.focal)/2)+1):nrow(IGV.focal),1]<-paste("chr", chromosomes(CGHset), sep="")	
	
	IGV.focal[1:(nrow(IGV.focal)/2),2]<-bpstart(CGHset)
	IGV.focal[((nrow(IGV.focal)/2)+1):nrow(IGV.focal),2]<-bpstart(CGHset)+(half.length+1)
	IGV.focal[1:(nrow(IGV.focal)/2),3]<-bpend(CGHset)-half.length
	IGV.focal[((nrow(IGV.focal)/2)+1):nrow(IGV.focal),3]<-bpend(CGHset)
	
	IGV.focal[1:(nrow(IGV.focal)/2),4]<-featureNames(CGHset)
	IGV.focal[((nrow(IGV.focal)/2)+1):nrow(IGV.focal),4]<-featureNames(CGHset)
	
	## Get frequencies
	Gains<-matrix(data=0, ncol=ncol(assayDataElement(CGHset, 'focal')), nrow=nrow(assayDataElement(CGHset, 'focal')))
	Losses<-matrix(data=0, ncol=ncol(assayDataElement(CGHset, 'focal')), nrow=nrow(assayDataElement(CGHset, 'focal')))
	Gains[which(assayDataElement(CGHset, 'focal')>0)]<- 1
	Losses[which(assayDataElement(CGHset, 'focal')<0)]<- 1
	
	IGV.focal[1:(nrow(IGV.focal)/2),5]<-round((rowSums(Gains)/ncol(calls(CGHset)))*100)
	IGV.focal[((nrow(IGV.focal)/2)+1):nrow(IGV.focal),5]<-round(((rowSums(Losses)/ncol(calls(CGHset)))*100)*-1)
	
	IGV.focal<-IGV.focal[order(as.numeric(gsub("chr","",IGV.focal[,1])), partial=as.numeric(IGV.focal[,2])),]
	IGV.focal[which(IGV.focal[,1]=="chr23"),1]<-"chrX"
	
	# Insert .igv header for formatting
	trackIDs<-matrix(data=0, ncol=5, nrow=2)
	trackIDs[1,1]<-"#type=OTHER"
	trackIDs[2,1]<-"#track color=150,0,0 altColor=0,0,150 autoScale=off graphType=bar scaleType=linear viewLimits=-100:100"
	IGV.focal<-rbind(trackIDs, c("Chromosome","Start","End","Feature", "Frequency_focals"), IGV.focal)
	IGV.focal[1,which(IGV.focal[1,]==0)]<-""
	IGV.focal[2,which(IGV.focal[2,]==0)]<-""

	
	write.table(IGV.focal, file="FrequencyPlotfocals.igv",row.names=FALSE, quote=FALSE, append=FALSE, col.names=FALSE, sep="\t")
	cat("Generated .IGV file for loading in IGV", "\n")
	
}


