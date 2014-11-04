#Loading script to workspace
args <- commandArgs(trailingOnly = TRUE)

source("spco-wt.r")
prefix="test" 
run_spco <- function(prefix, npop1, npop2, nadm){   
	##first chromosome
	#prepare data
	cat("============ Analysis on Chromosome 1 ============\n")
	cat("Start read and convert data ...\n")
	gfile <- paste(prefix, "Chr1.sgeno", sep="") #genotype file
    sfile <- paste(prefix, "Chr1.snp", sep="")  #SNPs position file
    ifile <- paste(prefix, ".ind", sep="")  #individual file
    #read individuals
    indList <- read.table(ifile)$V1
    numInds <- length(indList)
    #reading data and convert to tables
    tbl <- as.matrix(read.table(gfile)) #genotype table 
    snp <- read.table(sfile)
    colnames(tbl) <- as.character(indList)
    rownames(tbl) <- as.character(snp$V1)
    PhysPos <- snp$V4       #Physical position
    GenMap <- snp$V3*100    #Genetic map
    chromLen <- max(GenMap)/100   #in Morgan
    map <- matrix(c(PhysPos, GenMap), ncol=2)#matrix of map
    colnames(map) <- c("PhysPos", "GenMap") #Add column names
    rownames(map) <- PhysPos  #Add row names
    PhysPos<-structure(snp$V4,.Names=as.character(snp$V1)) # A List of Physical positions
    
    cat("Saving data ...\n")
    #Save data
    doutfile <- paste(prefix, "Chr1.Rdata", sep="")#data output file
    save(tbl, map, PhysPos, file=doutfile)
	
	cat("Start perform analysis ...\n")
    #Start perform analysis
    #Perform PCO, calculate PC axis on two parental populations
    #First 20 from parental 1, second 20 from parental 2, the rest from admixed
    pco_res<-pco(tbl,axis.who=colnames(tbl[,1:40])) 
    #Perform Step PCO, with first 20 and second 20 individuals as parental populations
    spco_res<-spco(tbl, pco=pco_res, genetic.map=map, Nbins=1024, window.size="3 sigma",pop1=colnames(tbl[,1:20]),pop2=colnames(tbl[,21:40]),max.window.size=Inf) 
    #Decide max level
    maxlvl<-5 #max level: 1-5, 7; 6-20, 6; 21-22, 5
    #Perform wavelet transformation
    wtstats_res<-spco.wtstats(spco_res, threshold=0.05, maxlevel=maxlvl)
    
    cat("Saving intermediate results ...\n")
    #Save wavelet transform summary
    wtoutfile<-paste(prefix, "_wtSummary.txt",sep="")
    wss<-wtstats_res$wt.summary
    write.table(wss,file=wtoutfile,sep="\t",quote=FALSE)
    #Save pco, spco, wavelet transform results
    routfile=paste(prefix, "_SPCO.Rdata",sep="")#result output file
    save(pco_res, spco_res, wtstats_res, file=routfile)
    #Done
	
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
