#Loading script to workspace
args <- commandArgs(trailingOnly = TRUE)

source("spco-wt.r")
pop="Hui"    
sampling<-read.table(paste(pop,"_exc.txt",sep=""))
ind<-sampling$V1
##first chromosome
cat("============ Analysis on Chromosome 1 ============\n")
cat("loading data ...\n")
load(paste(pop,"_chr1.Rdata",sep=""))
#get number of individuals
numInds<-length(colnames(tbl))
#get chromosome length in Morgan
chromLen<-max(map[,2])/100
load(paste(pop,"_SPCO_chr1.Rdata",sep=""))
maxlvl<-7
wtstats_res<-spco.wtstats(spco_res, threshold=0.05, maxlevel=maxlvl)
centers<-NA
for(rep in 1:5000){
  wtss<-subset(wtstats_res$wt.summary,rownames(wtstats_res$wt.summary)!=ind[rep])
  threshold_val<-vector()
  for (i in 1:maxlvl){threshold_val[i]<-max(wtss[c(1:100),i],na.rm=TRUE)}
  for (i in 1:numInds-1){wtss[i,1:maxlvl]<-cutoff(wtss[i,1:maxlvl],thr=threshold_val[1:maxlvl])} #perform normalization
  for (i in 1:numInds-1){wtss[i,"center"] <- sum(wtss[i,1:maxlvl]*(1:maxlvl))/sum(wtss[i,1:maxlvl])}
  #Normalize subtract from each calculated "center" log2 of the length of chromosome in Morgans:
  wtss[,"center"] <- (wtss[,"center"]) - (log2(chromLen))
  centers<-c(centers,wtss[c(101:numInds-1),"center"])
}
centers<-centers[2:length(centers)]
cmat<-cbind(centers)

#recursive calculate for each chromosome
for ( chrom in 2:22){
  cat("============ Analysis on Chromosome",chrom," ============\n")
  cat("loading data ...\n")
  #load data from previous saved data
  load(paste(pop,"_chr",chrom,".Rdata",sep=""))
  
  #Start perform analysis
  
  #get chromosome length in Morgan
  chromLen<-max(map[,2])/100
  #Load SPCO results
  load(paste(pop,"SPCO_chr",chrom,".Rdata",sep=""))
  #Decide max level
  maxlvl<-0 #max level: 1-5, 7; 6-20, 6; 21-22, 5
  if(chrom<=5){
    maxlvl=7
  }
  else if(chrom>5 && chrom<=20){
    maxlvl=6
  }
  else{
    maxlvl=5
  }
  
  #Perform wavelet transformation
  wtstats_res<-spco.wtstats(spco_res, threshold=0.05, maxlevel=maxlvl)
  
  #Filter out noise and normalize:
  #Get threshold
  threshold_val<-vector()
  for (i in 1:maxlvl){threshold_val[i]<-max(wtstats_res$wt.summary[c(1:100),i],na.rm=TRUE)}
  for (i in 1:numInds){wtstats_res$wt.summary[i,1:maxlvl]<-cutoff(wtstats_res$wt.summary[i,1:maxlvl],thr=threshold_val[1:maxlvl])} #perform normalization
  for (i in 1:numInds){wtstats_res$wt.summary[i,"center"] <- sum(wtstats_res$wt.summary[i,1:maxlvl]*(1:maxlvl))/sum(wtstats_res$wt.summary[i,1:maxlvl])}
  #Normalize subtract from each calculated "center" log2 of the length of chromosome in Morgans:
  wtstats_res$wt.summary[,"center"] <- (wtstats_res$wt.summary[,"center"]) - (log2(chromLen))
  
  #Prepare output
  routfile=paste(pop,"_SPCO2_chr",chrom,".Rdata",sep="")#result output file
  #Save results
  save(pco_res,spco_res,wtstats_res,file=routfile)
  #Done
}
