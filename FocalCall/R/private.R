## .private
## Collection of functions used within focalCall
## No further documentation privided

## 
.getSegments<-function(Segs, Chromo){
	# Count segments
	Segment.count<-matrix(data=0, ncol=ncol(Segs), nrow=nrow(Segs))
	for (j in 1:ncol(Segs)){
		k=1
		test.data<-Segs[,j]
		for (i in 1:length(test.data)){
			ifelse(test.data[i-1]==test.data[i] & Chromo[i-1]==Chromo[i], k<-k, k<-k+1)
			Segment.count[i,j]<-k
		}
		#cat("Sample", j, "finished", "\n")
	}
	return(Segment.count)
}


##
.getSegments_increase<-function(Segs, Chromo){
	# Count segments
	Segment.count<-matrix(data=0, ncol=1, nrow=length(Segs))
	#for (j in 1:nrow(Segs)){
		k=1
		test.data<-Segs#[,j]
		for (i in 1:length(test.data)){
			ifelse(test.data[i-1]+1==test.data[i] & Chromo[i-1]==Chromo[i], k<-k, k<-k+1)
			Segment.count[i,1]<-k
		}
		#cat("Sample", j, "finished", "\n")
	#}
	return(Segment.count)
}

## Matching CNV regions to arrayCGHprobes or CN-seq bins
.match_CNV2CGH<-function(calls, CNVlist){
CNVlist[,1]<-paste(gsub("chr", "", CNVlist[,"Chromosome"]))
CNVlist[which(CNVlist[,"Chromosome"]=="X"),1]<-"23"
	CNV.probe<-rep(0,length(bpstart(calls)))
	for (i in 1:length(unique(chromosomes(calls)))){
		CNVlist.tmp<-CNVlist[which(CNVlist[,"Chromosome"]==i),]
		BPstart.tmp<-bpstart(calls)[chromosomes(calls)==i]
		BPend.tmp<-bpend(calls)[chromosomes(calls)==i]
		temp<-NULL
		for (j in 1:nrow(CNVlist.tmp)){
			tmp.temp<-which(BPstart.tmp<as.numeric(paste(CNVlist.tmp[j,"End"])) & BPend.tmp>as.numeric(paste(CNVlist.tmp[j,"Start"])))
			temp<-c(temp, tmp.temp)
		}
		CNV.probe[which(chromosomes(calls)==i)][temp]<-1
		#cat("Finished Chromosome", i, "\n")
	}
	return(CNV.probe)
}


## Identify peak region or multiple regions
.getPeaks<-function(input){
	nr<-nrow(input)
	max_height<-max(input[,"focal"])
	
# First part, will finish when 1 peak found
	for (i in 3:max_height){
		temp.input<-which(input[,"focal"]>i | input[,"focal"]==i)
		
		finished<-FALSE
		
		if(i==max_height & length(temp.input)-1==(max(temp.input)-min(temp.input))){
			finished<-TRUE
			peak_region<-rep(0,2)
			peak_region[1]<-min(temp.input)
			peak_region[2]<-max(temp.input)
		}
	}
# Second part when more than 1 peak found
	if (finished==FALSE){
		i=1
		temp.input<-which(input[,"focal"]>i | input[,"focal"]==i)		
		uni.regions<-rep(0,length(input[,"focal"]))
		z<-1
		uni.regions[1]<-1
		for (i in 2:length(input[,"focal"])){
			ifelse(input[i-1,"focal"]==input[i,"focal"], z<-z, z<-z+1)
			uni.regions[i]<-z
		}
		uni.reg<-unique(uni.regions)
		freq.call<-rep(0,length(input[,"focal"]))
		for (i in 1:max(uni.reg)){
			index.temp<-which(uni.regions==uni.reg[i])

			if (i==1){
				index.temp2<-which(uni.regions==uni.reg[i+1])
				if(input[index.temp[1],"focal"]<input[index.temp2[1],"focal"]) {freq.call[index.temp]<-1}
			}
			if (i>1 & i<max(uni.reg)){
				index.temp3<-which(uni.regions==uni.reg[i-1])
				index.temp2<-which(uni.regions==uni.reg[i+1])
				if(input[index.temp[1],"focal"]<input[index.temp2[1],"focal"] | input[index.temp[1],"focalSegment"]<input[index.temp3[1],"focalSegment"]) {freq.call[index.temp]<-1}
			}
			if (i==max(uni.reg)){
				index.temp3<-which(uni.regions==uni.reg[i-1])
				if(input[index.temp[1],"focal"]<input[index.temp3[1],"focal"]) {freq.call[index.temp]<-1}
			}
		}
		peaks<-unique(uni.regions[which(freq.call==0)])
		peak_region<-matrix(data=0, ncol=2, nrow=length(peaks))
		for (i in 1:length(peaks)){
			index.tmp<-which(peaks[i]==uni.regions)
			peak_region[i,1]<-min(index.tmp)
			peak_region[i,2]<-max(index.tmp)
		}
	}
	return(peak_region)
}
	

