#A function to find the index of local maximums for a given vector x
localMaxima<-function(x){
  d1f <- diff(x)  #1 order differentiation 
  index <- 0
  for(i in 2:length(d1f)){
    #find critical point
    if((d1f[i-1] >= 0) && (d1f[i-1]*d1f[i] <= 0)){
      #check which one is larger
      if(x[i-1] > x[i]){
        index <- c(index, i-1)
      }else{
        index <- c(index, i)
      }
    } 
  }
  #remove the initial value
  index <- index[2:length(index)]
  return (index)
}

#A function to find the index of local minimums for a give vector x
localMinima<-function(x){
  d1f <- diff(x) #1 order differentiation 
  index <- 0
  for(i in 2:length(d1f)){
    #find critical point
    if((d1f[i-1] <= 0) && (d1f[i-1]*d1f[i] <= 0)){
      #check which one is closer to zero
      if(x[i-1] > x[i]){
        index <- c(index, i)
      }else{
        index <- c(index, i-1)
      }
    } 
  }
  #remove the initial value
  index <- index[2:length(index)]
  return (index)
}

solveG <- function(N, gens, breaks=2048,file="solve.pdf"){
  #N, the number of individuals sampled
  #gens, the mean generations bootstrapped, with 1 individual randomly dropped once
  #breaks, the number of breaks in estimate density function
  #file, the name of file output file, in which density of empirical and simulation 
  #with recovered parameters got from the density for the given mean generations
  f <- density(gens, n=breaks)
  maxs <- localMaxima(f$y)  #the index of peaks
  mins <- localMinima(f$y)  #the index of low peaks
  
  #calculate the probability for the first peak
  #simple integral over the left most point to the local minimum critical point
  prob <- 0
  index <- mins[1]
  for(i in 2:index){
    prob <- prob + f$y[i]*(f$x[i] - f$x[i-1])
  }
  #calculate the probabilities for the rest of peaks
  #integral over previous local minimum critical point to next local minimum critical point
  for(k in 2:length(mins)){
    indexp <- mins[k-1] #previous critical point
    indexc <- mins[k]   #current critical point
    tmp <- 0 #temp value for probability of k-th wave
    for(i in (indexp+1):indexc){
      tmp <- tmp + f$y[i]*(f$x[i] - f$x[i-1])
    }
    prob <- c(prob, tmp)
  }
  #the probability for last peak is the rest of value between 1 and the sum of previous p-values
  prob <- c(prob, 1-sum(prob))
  #the number of individuals in each wave, special attention should be paid to the approximation 
  numb <- round(N*prob) 
  
  #get the mean generations estimated for each peak
  index <- maxs[1]
  mgen <- f$x[index]
  for(i in 2:length(maxs)){
    index <- maxs[i]
    mgen <- c(mgen, f$x[index])
  }
  
  #write in matrix form Ag=b
  #g = (G1,G2,...,Gk)  #generations for k peaks
  #b = (N-1)*mg
  #construct the coefficients matrix A
  A <- rbind(numb)
  for(i in 2:length(maxs)){
    A <- rbind(A, numb)
  }
  #minus 1 from the digonal number
  for(i in 1:length(maxs)){
    A[i,i] <- A[i,i] - 1
  }
  
  b <- (N-1)*mgen
  #Solve for the generations
  gsolve <- solve(A, b)
  
  #Simulation with the recovered parameters, repeat 5000 times
  sim <- rep(gsolve[1], numb[1])
  for(i in 2:length(maxs)){
    sim <- c(sim, rep(gsolve[i], numb[i]))
  }
  #mean generations of simulations
  msim <- mean(sample(sim, N-1))
  for(i in 2:5000){
    msim <- c(msim, mean(sample(sim, N-1)))
  }
  sf <- density(msim)
  
  #plot empirical and simulated distributions
  pdf(file)
  xmin <- min(f$x, sf$x)
  xmax <- max(f$x, sf$x)
  ymax <- max(f$y, sf$y)
  plot(f$x, f$y, xlim = c(xmin, xmax), ylim = c(0, ymax), xlab = "Generation", ylab="Density", type="l")
  points(sf$x,sf$y,type="l",lty=2)
  legend("topright", c("Empirical", "Simulated"), lty = c(1, 2))
  dev.off()
  #prepare the return data frame
  svdf <- data.frame(prob, numb, mgen, gsolve)
  #return the generation
  return (svdf)
}
