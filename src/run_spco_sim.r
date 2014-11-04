#Loading script to workspace
source("spco-wt.r")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0){
    prefix <- args[1]
    cat("Start read and convert data ... \n")
    #prepare input files
    gfile <- paste(prefix, ".sgeno", sep="") #genotype file
    sfile <- paste(prefix, ".snp", sep="")  #SNPs position file
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
    doutfile <- paste(prefix, ".Rdata", sep="")#data output file
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
    
    cat("Start bootstrap, randomly drop 1 individual each time ...\n")
    d <- read.table(wtoutfile)
    inds <- rownames(d)[41:length(d[,1])]
    sct<-NA
    for(rep in 1:5000){
        exc_ind<-sample(inds,1)
        sd<-subset(d,rownames(d)!=exc_ind)
        threshold_val<-vector()
        for (i in 1:maxlvl){threshold_val[i]<-max(sd[c(1:40),i],na.rm=TRUE)}
        for (i in 1:numInds){sd[i,1:maxlvl]<-cutoff(sd[i,1:maxlvl],thr=threshold_val[1:maxlvl])} #perform normalization
        for (i in 1:numInds){sd[i,"center"] <- sum(sd[i,1:maxlvl]*(1:maxlvl))/sum(sd[i,1:maxlvl])}
        #Normalize subtract from each calculated "center" log2 of the length of chromosome in Morgans:
        sd[,"center"] <- (sd[,"center"]) - (log2(chromLen))
        ct<-sd[,"center"]
        mct<-mean(ct,na.rm=TRUE)
        sct<-c(sct,mct)
    }
    sct<-sct[2:length(sct)]
    
    cat("Output bootstrap results ...\n")
    write.table(sct,file=paste(prefix,"_bs5000_center.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)
    cat("Plot bootstrap results ...\n")
    pdf(paste(prefix,"_bs5000_center.pdf", sep=""))
    plot(density(sct),xlab="center")
    dev.off()
    msct<-mean(sct,na.rm=TRUE)
    sdsct<-sd(sct,na.rm=TRUE)
    cat(prefix,"centers re-ampling\n")
    cat("mean:",msct,"\n")
    cat("sd:",sdsct,"\n")    
} else {
    cat("Usage: Rscript run_spco_sim.r prefix\n")
}


