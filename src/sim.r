source("mvt.r")
sim <- function(nsamp=50, gen=10, step=5, out="sim.pdf"){
  #nsamp	number of sample to be simulated
  #gen		base generation
  #step		by which step the generation increases
  #out		output file name
  sn <- round(nsamp/2)	#maximum number of sample in each wave
  t <- sample(1:sn, 1)		#number of individual in each wave
  #randomly generate a number between 1 and half of of the number of samples,
  #and add to the total number of sample, until the total number reach or greater than
  #the number of samples
  for(i in 1:50){
    t <- c(t, sample(1:sn, 1))
    if(sum(t) >= nsamp){
      t[length(t)] = nsamp - sum(t[1:length(t) - 1])
      break
    }
  }
  cat(t, "\n")
  #generate samples in each wave
  x <- rep(gen, t[1])
  for(i in 2:length(t)){
    x <- c(x, rep(gen + (i-1) * step, t[i]))
  }
  mx <- mean(sample(x, nsamp - 1))
  #solveG(gens=mx,N=nsamp)
  for(i in 2:10000){
    mx <- c(mx,mean(sample(x, nsamp-1)))
  }
  pdf(out)
  title=paste("mean generation N=", nsamp, "g=", gen, "k=", step)
  plot(density(mx), main = title)
  dev.off()
  return (mx)
}
#test case
#mx<-sim(nsamp=100);plot(density(mx));solveG(gens=mx,N=100)
