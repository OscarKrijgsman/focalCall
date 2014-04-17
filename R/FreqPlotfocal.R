FreqPlotfocal<-function(calls, header="FrequencyPlot focal aberrations"){

	Chromo<-chromosomes(calls)
	focal.call<-assayDataElement(calls, 'focal')

## Add check for 'focal' within calls
	n_sam<-ncol(focal.call)
	n<-dim(focal.call)
	cat("n", "\t")
	Gains<-matrix(data=0, ncol=n[2], nrow=n[1])
	Losses<-matrix(data=0, ncol=n[2], nrow=n[1])

	Gains[which(focal.call>0)]<- 1
	Losses[which(focal.call<0)]<- 1

	Gain_sum<-rowSums(Gains)
	Loss_sum<-rowSums(Losses)

### Color focals and CNVs
	color.gain<-rep("grey", nrow(Gains))
	color.gain[which(fData(calls)$CNV==0 & Gain_sum!=0)]<-"red"
	color.loss<-rep("grey", nrow(Gains))
	color.loss[which(fData(calls)$CNV==0 & Loss_sum!=0)]<-"blue"

#### Plotting
	plot((Gain_sum/n_sam)*100, ylim=c(-100,100), type="h", col=color.gain, xlab="Chromosome", ylab="Frequency (Percentage)", main=paste(header),xaxt="n")

	points((-Loss_sum/n_sam)*100, type="h", col=color.loss)
	abline(h=0, cex=4)
	uni.chr<-unique(Chromo)
	temp<-rep(0,length(uni.chr))
	for (i in 1:length(uni.chr)){
		temp[i]<-max(which(uni.chr[i]==Chromo))
	}
	for (i in 1:length(temp)){
		abline(v=temp[i], col="black", lty="dashed")
	}

# ADD LABELS AND LINES
	nChroms<-length(unique(Chromo))
	begin<-c()
	for (d in 1:nChroms){
			chrom <- sum(Chromo == d)
     		begin <- append(begin,chrom)
     	   }
	uni.chr<-unique(Chromo)
	temp2<-rep(0,length(uni.chr))
	for (i in 1:length(uni.chr)){
        if(i==1){
                temp2[1]  <- (begin[1]*0.5)}
        else if(i>1){
        temp2[i]<- temp[i-1]+(begin[i]*0.5)}
	}
	for (i in 1:length(temp)){
        axis(1,at=temp2[i],labels=uni.chr[i],cex.axis=1)
	}
}




