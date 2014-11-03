source("mvt.r")
sim<-function(nsamp=50,gen=10,step=5,out="sim.pdf"){
  #nsamp	number of sample to be simulated
  #gen		base generation
  #step		by which step the generation increases
  #out		outputfile name
  sn<-round(nsamp/2)
  t<-sample(1:sn,1)
  for(i in 1:50){
    t<-c(t,sample(1:sn,1))
    if(sum(t)>=nsamp){
      t[length(t)]=nsamp-sum(t[1:length(t)-1])
      break
    }
  }
  cat(t,"\n")
  x<-rep(gen,t[1])
  for(i in 2:length(t)){
    x<-c(x,rep(gen+(i-1)*step,t[i]))
  }
  mx<-mean(sample(x,nsamp-1))
  #solveG(gens=mx,N=N)
  for(i in 2:10000){
    mx<-c(mx,mean(sample(x,nsamp-1)))
  }
  pdf(out)
  title=paste("mean generation N=",nsamp,"g=",gen,"k=",step)
  plot(density(mx),main=title)
  dev.off()
  return (mx)
}
#test case
#mx<-sim(N=100);plot(density(mx));solveG(gens=mx,N=100)
