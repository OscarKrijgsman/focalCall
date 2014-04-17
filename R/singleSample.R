singleSample<-function(CGHset, CNVset, focalSize=3, OverlapPerc=0.2){


	timeStarted <- proc.time()


##########################################
# Input control: Dimensions, format etc. #
##########################################

# Test whether 1 or two eSets

	if (class(CGHset) == "cghCall"){
		ncol.cgh<-ncol(calls(CGHset))
	} else stop("Input should be a CGHcall object.")
	
	if (class(CNVset) == "cghCall"){
		CNVcalls <- TRUE
	} else if (class(CNVset) == "data.frame"){
		CNVcalls <- FALSE
	}
	
# Check probe names CGHset and CNVset
# In case featureNames don't match provide matching and report percentage matching 
	if (CNVcalls==TRUE){
		if (sum(featureNames(CGHset)==featureNames(CGHset))!=length(featureNames(CGHset))){
			featureOverlap<-intersect(featureNames(CGHset), featureNames(CNVset))
			length_CGH<-length(featureOverlap)/length(featureNames(CGHset))*100
			length_CNV<-length(featureOverlap)/length(featureNames(CNVset))*100
			CGHset<-CGHset[match(featureOverlap, featureNames(CGHset)),]
			CNVset<-CNVset[match(featureOverlap, featureNames(CNVset)),]
			cat(paste("Number of featureNames not equal between CGHcall sets. \n", "Matching between the CGHcall sets performed: \n",length_CGH ,"% of tumor set matched \n", length_CNV, "% of normal set matched\n" , sep=""))
		}
	}


# Check dimensions CNVset (.bed file)

	if (CNVcalls==FALSE){
		class(CNVset) == "data.frame"
	}	
	
# Check resolution after matching CGHset and CNVset
	if (nrow(calls(CGHset))<30000){
		cat(paste("Array resolution too low for calling aberrations smaller than ", focalSize , "MB. \n", sep=""))
	}
	

############
# Analysis #
############

# Add CNV data to CGHset
if (CNVcalls==FALSE){ 
	fData(CGHset)$CNV <- .match_CNV2CGH(CGHset, CNVset)
}

if (CNVcalls==TRUE){ 
	assayDataElement(CGHset, 'matchedNormal')<-calls(CNVset)
	cat("Added matched normal calls to the tumor R-object.", "\n")
} 

# Count segments
	Seg_count<-.getSegments(segmented(CGHset), chromosomes(CGHset))
	cat("Counted number of segments for tumor sample", "\n")


## Focal calls are stored in cghCalls object as 'focal'
	focal.call<-matrix(data=0,ncol=ncol(calls(CGHset)), nrow=nrow(calls(CGHset)))
	for (k in 1:ncol(focal.call)){
		uni.seg<-unique(Seg_count[,k])
		for (i in 1:length(uni.seg)){
			temp<-which(Seg_count[,k]==uni.seg[i])
			if ((bpend(CGHset)[max(temp)]-bpstart(CGHset)[min(temp)])<(1000000*focalSize)){
				focal.call[temp,k]<-calls(CGHset[temp,k])
			}	
		}
	}
	assayDataElement(CGHset, 'focal')<-focal.call
	rm(focal.call)
	cat(paste("Generated matrix with detected aberrations <", focalSize , "MB. \n", sep=""))

#### Identify focal regions
tmp.focal<-matrix(data=0, ncol=ncol(calls(CGHset)), nrow=nrow(calls(CGHset)))
tmp.focal[which(assayDataElement(CGHset, 'focal')!=0)]<-1
focalList<-data.frame(chromosomes(CGHset), bpstart(CGHset), bpend(CGHset),tmp.focal)
focalList<-data.frame(focalList, rep(0,nrow(focalList)))
focalList[,5]<-1:nrow(focalList)

### Only get focals that are recurrent
focalList_short<-focalList[which(focalList[,4]!=0),]
focalList_short[,6]<-.getSegments_increase(focalList_short[,5], focalList_short[,1])
colnames(focalList_short)<-c("Chromo", "BPstart", "BPend", "focal", "index", "focalSegment")

### Generate matrix to export
MaxPeaks<-matrix(data=0, ncol=7, nrow=max(focalList_short[6]))
colnames(MaxPeaks)<-c("Chromo","BPstart","BPend","MB","Start_index","End_index","Type")

uni.focal<-unique(focalList_short[,6])

for (i in 1:max(focalList_short[6])){
	MaxPeaks[i,1]<-focalList_short[min(which(focalList_short[,6]==uni.focal[i])),1]
	MaxPeaks[i,2]<-focalList_short[min(which(focalList_short[,6]==uni.focal[i])),2]
	MaxPeaks[i,3]<-focalList_short[max(which(focalList_short[,6]==uni.focal[i])),3]
	MaxPeaks[i,4]<-round((MaxPeaks[i,3]-MaxPeaks[i,2])/1000000,3)
	MaxPeaks[i,5]<-focalList_short[min(which(focalList_short[,6]==uni.focal[i])),5]
	MaxPeaks[i,6]<-focalList_short[max(which(focalList_short[,6]==uni.focal[i])),5]
	if (CNVcalls==TRUE){
		tmp.calls<-assayDataElement(CGHset, 'matchedNormal')[MaxPeaks[i,5]:MaxPeaks[i,6]]}
	if (CNVcalls==FALSE){		
		tmp.calls<-fData(CGHset)$CNV[MaxPeaks[i,5]:MaxPeaks[i,6]]}
# Adding type of aberation		
	if (median(calls(CGHset)[MaxPeaks[i,5]:MaxPeaks[i,6]])==2){
		MaxPeaks[i,7]<- 2}
	if (median(calls(CGHset)[MaxPeaks[i,5]:MaxPeaks[i,6]])== -2){
		MaxPeaks[i,7]<- -2}
	if (median(calls(CGHset)[MaxPeaks[i,5]:MaxPeaks[i,6]])== 1){
		MaxPeaks[i,7]<-1}
	if (median(calls(CGHset)[MaxPeaks[i,5]:MaxPeaks[i,6]])== -1){
		MaxPeaks[i,7]<- -1}
	if ((length(which(tmp.calls!=0))/length(tmp.calls)) > 0.2){
		MaxPeaks[i,7]<-99}
}


# Plotting
if (CNVcalls==TRUE){
	
	png(file=paste("Tumor_",colnames(calls(CGHset[,1])), ".png", sep=""),height=480, width=2*480 )
		plot(CGHset[,1], dotres=1, ylim=c(-2,4))
	dev.off()
	
	png(file=paste("Normal_",colnames(calls(CGHset[,1])), ".png", sep=""),height=480, width=2*480 )
		plot(CNVset[,1], dotres=1, ylim=c(-2,3))
	dev.off()	

}
if (CNVcalls==FALSE){
	png(file=paste("Tumor_",colnames(calls(CGHset[,1])), ".png", sep=""), height=480, width=2*480)
		plot(CGHset[,1], dotres=1, ylim=c(-2,4))
	dev.off()
}


#### Output (.txt & .RData)
cat("Generating output files...", "\n")

# Focal.call eSet (CGHset)
save(CGHset, file=paste("focalCall",colnames(calls(CGHset[,1])), ".Rdata", sep=""))
# Focallist
write.table(MaxPeaks, file=paste("focalList",colnames(calls(CGHset[,1])), ".txt", sep=""), sep="\t", row.names=FALSE)
############### End analysis section ######	
	
# Tidy up nicely ;-)
	cat("FINISHED!\n")
	timeFinished <- round((proc.time() - timeStarted)[1] / 60)
	cat("Total time:", timeFinished, "minutes\n")	
	
#return(CGHset)
}
#EOF
