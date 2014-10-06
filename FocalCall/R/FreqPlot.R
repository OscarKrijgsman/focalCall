FreqPlot<-function(calls, header="FrequencyPlot all aberrations"){
	Chromo<-chromosomes(calls)
	Calls<-calls(calls)

	n_sam<-ncol(Calls)
	Gains<-matrix(data=0, ncol=ncol(Calls), nrow=nrow(Calls))
	Losses<-matrix(data=0, ncol=ncol(Calls), nrow=nrow(Calls))

	Gains[which(Calls>0)]<- 1
	Losses[which(Calls<0)]<- 1

	Gain_sum<-rowSums(Gains)
	Loss_sum<-rowSums(Losses)

	#### Plotting
	plot((Gain_sum/n_sam)*100, ylim=c(-100,100), type="h", col="red", xlab="Chromosome", ylab="Frequency (Percentage)", main=paste(header),xaxt="n")
	points((-Loss_sum/n_sam)*100, type="h", col="blue")
	abline(h=0, cex=4)
	uni.chr<-unique(Chromo)
	temp<-rep(0,length(uni.chr))
	for (i in 1:length(uni.chr)){
		temp[i]<-max(which(uni.chr[i]==Chromo))}
	for (i in 1:length(temp)){
		abline(v=temp[i], col="black", lty="dashed")}

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




