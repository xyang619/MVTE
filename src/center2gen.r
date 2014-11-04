convert<-function(cmean, yval){
    x<-2:101
    y1<-cmean[1]
    x1<-x[1]
    xval<-0
    for(i in 2:100){
        y2<-cmean[i]
        x2<-x[i]
        if((yval-y1)*(yval-y2)<=0){
            xval<-x1+(x2-x1)*(yval-y1)/(y2-y1)
            break
        }else{
            x1<-x2
            y1<-y2
        }
    }
    return (xval)
} 

cvtall<-function(simfile, rsmpfile,pop="Hui"){
    data<-read.table(simfile)
    mat<-as.matrix(data)
    cmean<-apply(mat,2,mean)
    d<-read.table(rsmpfile)
    yvals<-d$V1
    gens<-convert(cmean, yvals[1])
    for(i in 2:length(yvals)){
        gen<-convert(cmean, yvals[i])
        gens<-c(gens,gen)
    }
    df<-data.frame(yvals, gens)
    outfile<-paste(pop,"_resamp_gen.txt",sep="")
    write.table(df,file=outfile,quote=F,row.names=F)
}


