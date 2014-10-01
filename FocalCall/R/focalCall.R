focalCall<-function(CGHset, CNVset, focalSize=3, minFreq=2){

	timeStarted <- proc.time()


#################################################
# Input control: Dimensions, format, resolution #
#################################################

# Test whether 1 or two eSets

	if (class(CGHset) == "cghCall"){
		ncol.cgh<-ncol(calls(CGHset))
	} else stop("Input should be a CGHcall object.")
	
	if (class(CNVset) == "cghCall"){
		CNVcalls <- TRUE
	} else if (class(CNVset) == "data.frame"){
		CNVcalls <- FALSE
	}
	
# Check dimensions CNVset (.bed file)

	if (CNVcalls==FALSE){
		class(CNVset) == "data.frame"
	}	
	
# Check dimensions CGHset and CNVset
	if (CNVcalls==TRUE){
		if((ncol(calls(CGHset))==ncol(calls(CNVset)))!=TRUE){
			stop("Number of samples in eSets are not equal! \n")	
		}
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

# Check resolution after matching CGHset and CNVset
	if (nrow(calls(CGHset))<30000){
		cat(paste("Array resolution too low for calling aberrations smaller than ", focalSize , "MB. \n", sep=""))
	}
	
##########################
# Start analysis section #
##########################
	
# Count segments
	Seg_count<-.getSegments(segmented(CGHset), chromosomes(CGHset))
	cat("Counted number of segments for each tumor sample", "\n")

## Focal calls are stored in cghCalls object as 'focal'
	focal.call<-matrix(data=0,ncol=ncol(calls(CGHset)), nrow=nrow(calls(CGHset)))
	for (k in 1:ncol(focal.call)){
		uni.seg<-unique(Seg_count[,k])
		for (i in 1:length(uni.seg)){
			temp<-which(Seg_count[,k]==uni.seg[i])
			if ((bpstart(CGHset)[max(temp)]-bpstart(CGHset)[min(temp)])<1000000*focalSize){
				focal.call[temp,k]<-calls(CGHset[temp,k])
			}	
		}
	}
	assayDataElement(CGHset, 'focal')<-focal.call
	rm(focal.call)
	cat("Generated matrix with detected aberrations <", focalSize, "Mb...", "\n")

############## Add calls normals to cancer set! 

if (CNVcalls==TRUE){ 
	assayDataElement(CGHset, 'matchedNormal')<-calls(CNVset)
	rm(CNVset)
	cat("Added matched normal calls to the tumor R-object.", "\n")
} 

#### Identify focal regions
tmp.focal<-matrix(data=0, ncol=ncol(calls(CGHset)), nrow=nrow(calls(CGHset)))
tmp.focal[which(assayDataElement(CGHset, 'focal')!=0)]<-1
focalList<-data.frame(chromosomes(CGHset), bpstart(CGHset), bpend(CGHset),rowSums(tmp.focal))
focalList<-data.frame(focalList, rep(0,nrow(focalList)))
focalList[,5]<-1:nrow(focalList)

### Only get focals that are recurrent
focalList_short<-focalList[which(focalList[,4]>(minFreq-1)),]
focalList_short[,6]<-.getSegments_increase(focalList_short[,5], focalList_short[,1])
colnames(focalList_short)<-c("Chromo", "BPstart", "BPend", "focal", "index", "focalSegment")

#######################################
# Matrix generation focal aberrations #
#######################################

MaxPeaks<-matrix(data=0, ncol=16, nrow=2)
colnames(MaxPeaks)<-c("index","Chromo","BPstart","BPend","MB","Peak_start","Peak_end","Start_index","End_index","Freq_gain","Freq_loss", "Total_gain","Total_loss","freq_Amps","Freq_HDs","CNV_call")

cat("Start detection of peak regions...", "\n")
# get each region
uni.peak<-unique(focalList_short[,"focalSegment"])
for (i in 1:max(focalList_short[,"focalSegment"])){
	index_temp<-which(focalList_short[,"focalSegment"]==uni.peak[i])
# Set k as one when initiating loop
	ifelse(i==1, k<-1, k<-k+1)
	#cat("Display_",i, k, sep=" ", "\n")
	MaxPeaks<-rbind(MaxPeaks, rep(0,16))
	MaxPeaks[k,"index"]<-focalList_short[index_temp[1],"index"]
	MaxPeaks[k,"Chromo"]<-focalList_short[index_temp[1],"Chromo"]
	MaxPeaks[k,"BPstart"]<-focalList_short[index_temp[1],"BPstart"]
	MaxPeaks[k,"BPend"]<-focalList_short[max(index_temp),"BPend"]
	MaxPeaks[k,"MB"]<-(as.numeric(MaxPeaks[k,"BPend"])-as.numeric(MaxPeaks[k,"BPstart"]))/1000000

	if(length(unique(focalList_short[index_temp,"focal"]))==1){
		cat("Simple_", i, "\n")
		MaxPeaks[k,"Peak_start"]<-MaxPeaks[k,"BPstart"]
		MaxPeaks[k,"Peak_end"]<-MaxPeaks[k,"BPend"]
		MaxPeaks[k,"Start_index"]<-focalList_short[index_temp[1],"index"]
		MaxPeaks[k,"End_index"]<-focalList_short[max(index_temp),"index"]
	}

	if(length(unique(focalList_short[index_temp,4]))!=1){
		cat("Complex_", i, "\n")
		#rm(peak_region)
		peak_region<-.getPeaks(focalList_short[index_temp,])
		if (length(peak_region)==2){
			MaxPeaks[k,"Peak_start"]<-focalList_short[index_temp[peak_region[1]],"BPstart"]
			MaxPeaks[k,"Peak_end"]<-focalList_short[index_temp[peak_region[2]],"BPend"]
			MaxPeaks[k,"Start_index"]<-focalList_short[index_temp[peak_region[1]],"index"]
			MaxPeaks[k,"End_index"]<-focalList_short[index_temp[peak_region[2]],"index"]
		}
		if (length(peak_region)>2){
			for (j in 1:nrow(peak_region)){
				if (j==1){
					MaxPeaks[k,"Peak_start"]<-focalList_short[index_temp[peak_region[j,1]],"BPstart"]
					MaxPeaks[k,"Peak_end"]<-focalList_short[index_temp[peak_region[j,2]],"BPend"]
					MaxPeaks[k,"Start_index"]<-focalList_short[index_temp[peak_region[j,1]],"index"]
					MaxPeaks[k,"End_index"]<-focalList_short[index_temp[peak_region[j,2]],"index"]
				}
				if (j>1){
					MaxPeaks<-rbind(MaxPeaks, rep(0,14))
					k<-k+1
					MaxPeaks[k,"index"]<-focalList_short[index_temp[1],"index"]
					MaxPeaks[k,"Chromo"]<-focalList_short[index_temp[1],"Chromo"]
					MaxPeaks[k,"BPstart"]<-focalList_short[index_temp[1],"BPstart"]
					MaxPeaks[k,"BPend"]<-focalList_short[index_temp[length(index_temp)],3]
					MaxPeaks[k,"MB"]<-(as.numeric(MaxPeaks[k,4])-as.numeric(MaxPeaks[k,3]))/1000000
					MaxPeaks[k,"Peak_start"]<-focalList_short[index_temp[peak_region[j,1]],2]
					MaxPeaks[k,"Peak_end"]<-focalList_short[index_temp[peak_region[j,2]],3]
					MaxPeaks[k,"Start_index"]<-focalList_short[index_temp[peak_region[j,1]],5]
					MaxPeaks[k,"End_index"]<-focalList_short[index_temp[peak_region[j,2]],5]
				}
			}
		}
	}	
}

## Delete previously added empty rows
MaxPeaks<-MaxPeaks[-((nrow(MaxPeaks)-1 ): nrow(MaxPeaks)),]


# Add CNV data to CGHset
if (CNVcalls==FALSE){
	fData(CGHset)$CNV <- .match_CNV2CGH(CGHset, CNVset)
}

for (i in 1:nrow(MaxPeaks)){
	MaxPeaks[i,"Freq_gain"]<-length(which(assayDataElement(CGHset, 'focal')[as.numeric(MaxPeaks[i,"Start_index"]),]>0))
	MaxPeaks[i,"Freq_loss"]<-length(which(assayDataElement(CGHset, 'focal')[as.numeric(MaxPeaks[i,"Start_index"]),]<0))
	MaxPeaks[i,"Total_gain"]<-length(which(calls(CGHset)[as.numeric(MaxPeaks[i,"Start_index"]),]>0))
	MaxPeaks[i,"Total_loss"]<-length(which(calls(CGHset)[as.numeric(MaxPeaks[i,"Start_index"]),]<0))
	# Add HDs and Amplification
	MaxPeaks[i,"freq_Amps"]<-length(which(calls(CGHset)[as.numeric(MaxPeaks[i,"Start_index"]),]==2))
	MaxPeaks[i,"Freq_HDs"]<-length(which(calls(CGHset)[as.numeric(MaxPeaks[i,"Start_index"]),]== (-2)))
	
	if (CNVcalls==TRUE){	
		MaxPeaks[i,"CNV_call"]<-length(which(assayDataElement(CGHset, 'matchedNormal')[as.numeric(MaxPeaks[i,"Start_index"]),which(assayDataElement(CGHset, 'focal')[as.numeric(MaxPeaks[i,"Start_index"]),]!=0)]!=0))}
	if (CNVcalls==FALSE){
		tmp.cnv<-sum(fData(CGHset)$CNV[MaxPeaks[i,"Start_index"]:MaxPeaks[i,"End_index"]])
		MaxPeaks[i,"CNV_call"]<-tmp.cnv/length(MaxPeaks[i,"Start_index"]:MaxPeaks[i,"End_index"])} 
}


cat("Matching CNV list to copynumber data ...", "\n")

# Match each probe to CNV locations from the dataframe:CNVset
if (CNVcalls==TRUE){
	CNV<-matrix(data=0, ncol=ncol(assayDataElement(CGHset, 'matchedNormal')), nrow=nrow(assayDataElement(CGHset, 'matchedNormal')))
	CNV[which(assayDataElement(CGHset, 'matchedNormal')!=0)]<- 1	
	CNV_sum<-rep(0,nrow(assayDataElement(CGHset, 'matchedNormal')))	
	CNV_sum[which(rowSums(CNV)>2)]<-1
	fData(CGHset)$CNV <- CNV_sum
}


#### Output (.txt & .RData)
cat("Generating output files...", "\n")

# Summary matrix file
write.table(MaxPeaks[,-1], file="Focals_list.txt", sep="\t", row.names=FALSE)
cat("Written text file to workDir...", "\n")


# Focal.call eSet (CGHset)
save(CGHset, file="focalCall.RData")
cat("Written 'focalCall.RData' to workDir...", "\n")


############### End analysis section ######	
	
# Tidy up nicely ;-)
	cat("FINISHED!\n")
	timeFinished <- round((proc.time() - timeStarted)[1] / 60)
	cat("Total time:", timeFinished, "minutes\n")	
	
return(CGHset)
}
# EOF
