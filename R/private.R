## .private
## Collection of functions used within focalCall
## No further documentation privided

## 
.getSegments<-function(Segs, Chromo){
	# Count segments
	Segment.count<-matrix(data=0, ncol=dim(Segs)[2], nrow=length(Segs[,1]))
	for (j in 1:dim(Segs)[2]){
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
	#for (j in 1:dim(Segs)[2]){
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
CNVlist[,1]<-paste(gsub("chr", "", CNVlist[,1]))
CNVlist[which(CNVlist[,1]=="X"),1]<-"23"
	CNV.probe<-rep(0,length(bpstart(calls)))
	for (i in 1:length(unique(chromosomes(calls)))){
		CNVlist.tmp<-CNVlist[which(CNVlist[,1]==i),]
		BPstart.tmp<-bpstart(calls)[chromosomes(calls)==i]
		BPend.tmp<-bpend(calls)[chromosomes(calls)==i]
		temp<-NULL
		for (j in 1:nrow(CNVlist.tmp)){
			tmp.temp<-which(BPstart.tmp<as.numeric(paste(CNVlist.tmp[j,3])) & BPend.tmp>as.numeric(paste(CNVlist.tmp[j,2])))
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
	max_height<-max(input[,4])
	
# First part, will finish when 1 peak found
	for (i in 3:max_height){
		temp.input<-which(input[,4]>i | input[,4]==i)
		
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
		temp.input<-which(input[,4]>i | input[,4]==i)		
		uni.regions<-rep(0,length(input[,4]))
		z<-1
		uni.regions[1]<-1
		for (i in 2:length(input[,4])){
			ifelse(input[i-1,4]==input[i,4], z<-z, z<-z+1)
			uni.regions[i]<-z
		}
		uni.reg<-unique(uni.regions)
		freq.call<-rep(0,length(input[,4]))
		for (i in 1:max(uni.reg)){
			index.temp<-which(uni.regions==uni.reg[i])

			if (i==1){
				index.temp2<-which(uni.regions==uni.reg[i+1])
				if(input[index.temp[1],4]<input[index.temp2[1],4]) {freq.call[index.temp]<-1}
			}
			if (i>1 & i<max(uni.reg)){
				index.temp3<-which(uni.regions==uni.reg[i-1])
				index.temp2<-which(uni.regions==uni.reg[i+1])
				if(input[index.temp[1],4]<input[index.temp2[1],4] | input[index.temp[1],6]<input[index.temp3[1],6]) {freq.call[index.temp]<-1}
			}
			if (i==max(uni.reg)){
				index.temp3<-which(uni.regions==uni.reg[i-1])
				if(input[index.temp[1],4]<input[index.temp3[1],4]) {freq.call[index.temp]<-1}
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
	

