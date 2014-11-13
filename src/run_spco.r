#############################################################################
#the function is used to perform StepPCO on 22 autosomal chromosome
#the input file should be named as $prefix_Chr$i.$suffix
#the file format support is eigenstrat, i.e. xx.geno, xx.snp and xx.ind
#before run the analysis, the genofile should convert into steppco format
#by using recode4spco.py, with extension sgeno
#Attention, the results was NOT normalized because I want to use the original
#ones to perform bootstrapping, if you don't care , just load the results and 
#perform normalizing as describe in the readme of spco-wt.r
#############################################################################
#Loading script to workspace
source("spco-wt.r") 
run_spco <- function(prefix, npop1, npop2){   
    for(chr in 1:22){
		#prepare data
		cat(paste("============ Analysis on Chromosome",chr,"============\n"))
		cat("Start read and convert data ...\n")
		gfile <- paste(prefix, "_Chr",chr,".sgeno", sep="") #genotype file
		sfile <- paste(prefix, "_Chr",chr,".snp", sep="")  #SNPs position file
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
		doutfile <- paste(prefix, "_Chr",chr,".Rdata", sep="")#data output file
		save(tbl, map, PhysPos, file=doutfile)
	
		cat("Start perform analysis ...\n")
		#Start perform analysis
		#Perform PCO, calculate PC axis on two parental populations
		#First npop1 from parental 1, second npop2 from parental 2, the rest from admixed
		pco_res<-pco(tbl,axis.who=colnames(tbl[,1:(npop1+npop2)])) 
		#Perform Step PCO, with first npop1 and second npop2 individuals as parental populations
		spco_res<-spco(tbl, pco=pco_res, genetic.map=map, Nbins=1024, window.size="3 sigma",pop1=colnames(tbl[,1:npop1]),pop2=colnames(tbl[,npop1+1:npop1+npop2]),max.window.size=Inf) 
		#Decide max level
		maxlvl <- 7	#max level: 1-5, 7; 6-20, 6; 21-22, 5
		if(chr <= 5){
			maxlvl <- 7
		}
		else if(chr > 5 && chr <= 20){
			maxlvl <- 6
		}
		else{
			maxlvl <- 5
		}
		
		#Perform wavelet transformation
		wtstats_res<-spco.wtstats(spco_res, threshold=0.05, maxlevel=maxlvl)
    
		cat("Saving intermediate results ...\n")
		#Save wavelet transform summary
		wtoutfile<-paste(prefix, "_Chr1_wtSummary.txt",sep="")
		wss<-wtstats_res$wt.summary
		write.table(wss,file=wtoutfile,sep="\t",quote=FALSE)
		#Save pco, spco, wavelet transform results
		routfile=paste(prefix, "_Chr1_SPCO.Rdata",sep="")#result output file
		save(pco_res, spco_res, wtstats_res, file=routfile)
		#Done
		cat(paste("Finished analyzing on Chromosome",chr))
	}
}
