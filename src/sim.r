source("mvt.r")
sim<-function(N=50,g=10,k=5,out="sim.pdf"){
  sn<-round(N/2)
  t<-sample(1:sn,1)
  for(i in 1:50){
    t<-c(t,sample(1:sn,1))
    if(sum(t)>=N){
      t[length(t)]=N-sum(t[1:length(t)-1])
      break
    }
  }
  cat(t,"\n")
  x<-rep(g,t[1])
  for(i in 2:length(t)){
    x<-c(x,rep(g+(i-1)*k,t[i]))
  }
  mx<-mean(sample(x,N-1))
  #solveG(gens=mx,N=N)
  for(i in 2:10000){
    mx<-c(mx,mean(sample(x,N-1)))
  }
  pdf(out)
  title=paste("mean generation N=",N,"g=",g,"k=",k)
  plot(density(mx),main=title)
  dev.off()
  return (mx)
}
#test case
#mx<-sim(N=100);plot(density(mx));solveG(gens=mx,N=100)
