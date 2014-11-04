################################################################################
# SPCO
# Copyright (C) 2011 Rostislav Matveev, Irina Pugach
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
################################################################################


################################################################################
# 
# AUX STUFF

msg <- function(...,sep="",new.line=FALSE){
  cat(paste(...,sep=sep))
  if(new.line) cat("\n")
}

msgn <- function(...,sep="",new.line=TRUE){
  msg(...,sep=sep,new.line=new.line)
}

tick <- function(point="."){
  cat(point)
}

pause <- function(){
  msg("PAUSE! press \"enter\" to continue")
  scan(nmax=1,quiet=TRUE,nlines=1)
}

make.array <- function(...,data=NA){
  DimNames   <- list(...)
  rank       <- length(DimNames)
  TotalSize  <- 1
  Dim        <- vector(mode="integer",length=rank)
  for(i in 1:rank){
    Dim[i]   <- length(DimNames[[i]])
    if(Dim[i] == 1 && is.double(DimNames[[i]])){
       Dim[i] <- DimNames[[i]]
       DimNames[[i]] <- 1:Dim[i]
    }      
    DimNames[[i]] <- as.character(DimNames[[i]])
    TotalSize <- TotalSize * Dim[i]
  }
  array(data=data, dim=Dim, dimnames = DimNames)
} 

how.many <- function(x){
  return(length(which(x)))
}

matches <- function(ptrn,str,ignore.case=TRUE){
  return(length(grep(ptrn,str,ignore.case=ignore.case))==1)
}

size <- function(object,newline=TRUE){
  S <- object.size(object)
  m <- floor(log2(S)/10)
  S <- floor(S*10/(1024^m))/10
  suffix=c("b","Kb","Mb","Gb","Tb")
  cat(paste(S,suffix[m+1],sep=""))
  if(newline) cat("\n")
}

# AUX STUFF
#
############################################################################

################################################################################
#
# TEMPLATES

# type=0 constant 0
# type=1 constant 1
# 0<type<1 mixed, step is exponetialy ditributed with the mean=1/rate.
# this would produce approximately "rate" steps
chromosome <- function(type=0,
                       rate=1024
                       ){
  if( type == 0 )
    return(c(-1,99))
  if( type == 1 )
    return(c(-1,0,1,99))
  if( abs(type-1/2) < 1/2 ){
    rate.white = rate/(1-type)
    rate.red   = rate/type
    chromo=double(1)
    i=0
    while(chromo[2*i+1] < 1){
      i <- i+1
      chromo[2*i]   <- chromo[2*i-1] + rexp(1,rate=rate.red)
      chromo[2*i+1] <- chromo[2*i]   + rexp(1,rate=rate.white)
    }
    return(c(-1,chromo[1:(2*i-1)],1,99))
  }
}

population <- function(N=0,
                       type=0
                       ){
  p <- list()
  class(p) <- "recsim.population"
  if( N < 1 ) return(p)
  chromo <- chromosome(type=type)
  for(i in 1:N) p[[i]] <- chromo
  return(p)
}

events <- function(){
  e <- list()
  class(e) <- "recsim.events"
  return(e)
}

event <- function(event="none",gen=0,...){
  e <- list(event=event,generation=gen,...)
  class(e) <- "recsim.event"
  return(e)
}

info <- function(){
  inf <- make.array(NULL,c("size"))
  class(inf) <- "recsim.info"
  return(inf)
}

simulation <- function(name="Population",
                       comment="Comment",
                       recmap=NULL
                      ){
  s <- list()
  class(s) <- "recsim"
  s$name <- name
  s$comment <- comment
  s$recmap <- recmap
  s$events <- events()
  s$info <- info()
  s$generations <- list()
  return(s)
}

recmap <- function(rate=1.5,
                   bkgd=0.2,
                   space=1/5000,
                   distance.shape=1.5,
                   intensity.shape=1.5,
                   width=1/100000,
                   width.shape=1){
  recmap <- list(rate=rate,
                  bkgd=bkgd,
                  space=space,
                  distance.shape=distance.shape,
                  intensity.shape=intensity.shape,
                  width=width,
                  width.shape=width.shape)
  class(recmap) <- "recsim.recmap"

  next.hs <- function()
    c<<-c+rgamma(1,shape=distance.shape,scale=space/distance.shape)

  c=0
  recmap$hotspots <- replicate(1.1/space,next.hs())
  recmap$hotspots <- recmap$hotspots[which(recmap$hotspots<1)]

  next.int <- function()
    c<<-c+rgamma(1,shape=intensity.shape,scale=space/intensity.shape)

  c=0
  recmap$intensity <- replicate(length(recmap$hotspots),next.int())
  recmap$intensity <-
    recmap$intensity/recmap$intensity[length(recmap$intensity)]
  
  recmap$sigma <- rgamma(length(recmap$hotspots),
                          shape=width.shape,
                          scale=width/width.shape)
  return(recmap)
}


# TEMPLATES
#
################################################################################

################################################################################
#
# OPERATIONS

# Add N(1) chromosomes of type(1) to pop
admix <- function(pop,N=1,type=1){
  chromo <- chromosome(type=type)
  l <- length(pop)
  for(i in (l+1):(l+N)) pop[[i]] <- chromo
  return(pop)
}

subsample <- function(...){
  UseMethod("subsample")
}

subsample.recsim.population <- function(pop,sample.size=100){
  l <- length(pop)
  sample.size <- min(sample.size,l)
  if(sample.size==0) return(NULL)
  pop <- pop[sample(l,sample.size)]
  class(pop) <- "recsim.population"
  return(pop)
}

# Find recombination points according to rec map
recpoints <- function(recmap){
  #num of rec events (at least 1)
  num.rec <- rpois(1,recmap$rate-1)+1
  points<-numeric(num.rec)
  for(i in 1:num.rec){
    #are we in a HS or BG
    if( runif(1) <= recmap$bkgd ){
      points[i] <- runif(1)
    }
    else{ 
      #choose in which HS are we located
      hs<-min(which(recmap$intensity >= runif(1)))
      points[i] <- rnorm(1,recmap$hotspot[hs],recmap$sigma[hs])
    }
  }
  return(sort(points))
}

add.event <- function(events,...){
  e <- list(...)
  class(e) <- "recsim.event"
  l <- length(events)
  events[[l+1]] <- e
  return(events)
}


offspring <- function(parents,recmap,r=1,max=8){
  N <- rbinom(n=1,size=max,p=2*r/max)
  if(N==0) return()
  deti <- list()
  for(j in 1:N){
    points <- recpoints(recmap)
    deti[[j]] <- c(-1)
    p0 <- -1
    i <- 1
    for(p1 in c(points,100)){
      a <- min(which( parents[[i]] >  p0 ))
      b <- max(which( parents[[i]] <= p1 ))
      if(a<=b) range=a:b
      else     range=c()
      if( (length(deti[[j]]) - a) %% 2 == 1){
        deti[[j]] <- c(deti[[j]],parents[[i]][range])
      }
      else{
        deti[[j]] <- c(deti[[j]],p0,parents[[i]][range])
      }
      i <- 3 - i
      p0 <- p1
    }
  }
  return(deti)
}

new.gen <- function(pop,
                    target.size=length(pop),
                    recmap,
                    max.children=8
                   ){
  N=length(pop)
  if(N<=1) return(population(N=0))
  N0 = N %/% 2
  pairs = sample(N,N0*2)
  r=target.size/N

  pop <- pop[sample(N,2*N0)]
  dim(pop) <- c(2,N0)
  pop.new <- apply(pop,2,offspring,recmap=recmap,max=max.children,r=r)
  pop.new <- unlist(pop.new,recursive=FALSE)
  class(pop.new) <- "recsim.population"
  return(pop.new)
}

# read command from stdin or file
read.action <- function(prompt="recsim> ",
                        file="",
                        what=list(action="c",
                                  value ="c",
                                  value2="c",
                                  value3="c",
                                  value4="c",
                                  value5="c",
                                  value6="c",
                                  value7="c",
                                  value8="c",
                                  value9="c"
                                  )
                          ){
  if(file=="") file="stdin"
  comment=TRUE
  while(comment){
    comment=FALSE
    if(file=="stdin") cat(prompt)
    scan(file=file,
         what=what,
         n=1,
         nlines=1,
         quiet=TRUE,
         fill=TRUE,
         comment.char = "",
         blank.lines.skip = FALSE) -> action
    if(!is.na(action$action))
      if(matches("^#",action$action)) comment=TRUE
  }
  for(i in 1:length(action)){
    if(length(action[[i]])==0) action[[i]]=NA
    if(is.na(action[[i]]) |
       action[[i]] == "" ) action[[i]]=NA
    else{
      suppressWarnings(as.numeric(action[[i]])) -> tmp
      if(! is.na(tmp)) action[[i]]<-tmp
    }
  }
  if(is.na(action$action)){
    action$action <- "n"
    action$value  <- 1
  }
  if(is.numeric(action$action)){
    action$value  <- floor(abs(action$action))
    action$action <- "n"
  }
  class(action)<-"recsim.action"
  return(action)
}

# OPERATIONS
#
################################################################################

################################################################################
#
# RECSIM


recsim <- function(name="Default Simulation",
                   comment="Size 100; type 0",
                   file="",
                   size=100,
                   type=0,
                   recmap=NA
                  ){
  
  ##############################
  # setup connection for reading commands
  if(file!="" & file!="stdin"){
    connection=file(description=file,open="r")
  }
  else connection=""
  #
  ##############################

  ##############################
  # INITIAL PARAMETERS
  
  help.message="
---------------------------------------------------------------------
c <Name> <Type> <Size> <Comment>
\t\t Create initial population for the simulation.

recmap <Rate> <Bkgd> <Space> <DistShape> <IntensShape> <Width> <WidthShape>
\t\t Create recombination map

n <N>\t\t Go forward N generations (default 1)
<enter>\t\t Go forward 1 generation
<N> \t\t Go forward N generations
g <N>\t\t Go to generation N.
a <N> <T>\t Admix N (default 10) chromosomes of type T (default 1)
r <R>\t\t Set growth rate to R (default 1)
bn <N>\t\t Bottle neck, discard all but N chromosomes
q\t\t End simulation.

s <N>\t\t Set sample size to display, if 0 or no value - don't display
svd <N>\t\t Size of the sample to save ea generation (default 0)
sv  <N>\t\t Sample size to save for current generation (default 100)
sn <N> <filename>
\t\t Plot N (default 10) samples to the file
\t\t <filename> (pdf,png or jpg)(default is from pop name)

d\t\t Redisplay all the data.
h\t\t Show population history up to this timepoint.
v <N>\t\t Set verbosity level to N (default 2)
?\t\t Help
---------------------------------------------------------------------
"
  if(file=="" | file=="stdin")
    v <- 2
  else
    v <- 0
  samples      <- 0
  max.children <- 8
  r            <- 1
  save.default <- 0
  
  #
  ##############################
  
  
  ##############################
  # create new default population
  gen <- 1
  s<-simulation(name=name,comment=comment)
  cur <- population(N=size,
                    type=type
                   )
 
  s$events <- add.event(events=s$events,
                        event="create",
                        generation=gen,
                        name=name,
                        comment=comment,
                        type=type,
                        size=size
                       )
                       
  if(is.na(recmap))
    s$recmap <- recmap(rate=2.7,
                       bkgd=1,
                       space=0.002,
                       distance.shape=1.5,
                       intensity.shape=1.5,
                       width=1/100000,
                       width.shape=1
                      )
  else
    s$recmap=recmap
  
  
  s$events <- add.event(events=s$events,
                        event="recmap",
                        generation=gen,
                        recmap=recmap
                       )
    
  s$generations[[as.character(gen)]] <- subsample(cur,save.default)
  s$info <- rbind(s$info,c(size))
  target.size=size
  #
  ##############################

  ##############################
  # START
  finish <- FALSE
  while(!finish){
    if(v>=2){
      cat("\n")
      msgn("Simulation:\t",s$name)
      msgn("\t\t", s$comment)
      msgn("Generation:\t",gen)
      msgn("Size:\t\t",length(cur))
      msgn("Growth rate:\t",r)
      msgn("Show:\t\t", samples," samples")
      msgn("Save:\t\t", save.default," samples")
      msgn("Verbosity:\t",v)
      msgn("Target.size:\t",round(target.size*100)/100)
    }
    read.action(file=connection)->action
    if(v>=1) print(action)

    # Go to next generation
    if(action$action=="g" | action$action=="goto" ){
      action$action <- "n"
      action$value <- action$value - gen
      if(action$value<1){
        if(v>=0) msgn("Cannot go back in time.")
        next()
      }
    }
    if(action$action=="n" | action$action=="next" ){
      if(is.na(action$value) | !is.numeric(action$value))
        action$value=1
    
      for(k in 1:action$value){
        target.size=r*target.size
        if(v>=1) msg("Generation:",gen+1)
        cur <- new.gen(cur,
                       target.size=target.size,
                       recmap=s$recmap
                      )
        if(v>=1) msgn("; ",length(cur),"/",round(target.size))
        gen <- gen + 1
        s$info <- rbind(s$info,c(length(cur)))
        s$generations[[as.character(gen)]] <- subsample(cur,save.default)
        if(length(cur)==0){
          if(v>=1) cat("Population is extinct.\n")
          finish=TRUE
          s$events <- add.event(events=s$events,
                                event="extinct",
                                generation=gen
                               )
          break()
        }
      }
      plot(cur,samples=samples,main=paste(s$name,"; generation ",gen,sep=""))
      next()
    }

    # Create
    if(action$action=="c" | action$action=="create"){
      gen <- 1
      s <- simulation(name=action$value,comment=action$value4)
      cur <- population(N=action$value3,
                        type=action$value2
                       )
      s$events <- add.event(events=s$events,
                            event="create",
                            generation=gen,
                            name=action$value,
                            comment=action$value4,
                            type=action$value2,
                            size=action$value3
                           )
      #print(s)
      s$generations[[as.character(gen)]] <- subsample(cur,save.default)
      s$info <- rbind(s$info,c(size))
      target.size=action$value3
      next()
    }
      
    # Make recombination map
    if(action$action=="recmap"){
      s$recmap <- recmap(rate=action$value,
                         bkgd=action$value2,
                         space=action$value3,
                         distance.shape=action$value4,
                         intensity.shape=action$value5,
                         width=action$value6,
                         width.shape=action$value7
                        )
      s$events <- add.event(events=s$events,
                            event="recmap",
                            generation=gen,
                            recmap=s$recmap
                           )
      next()
    }


    # Bottleneck
    if(action$action=="bn" | action$action=="bottleneck" ){
      if(is.na(action$value) | !is.numeric(action$value)) action$value <- 10
      s$events <- add.event(events=s$events,
                            event="bottleneck",
                            generation=gen,
                            size=action$value
                           )
      target.size=action$value
      next()
    }
      
    # Redisplay
    if(action$action=="d" | action$action=="display"){
      if(samples > 0) plot(cur,samples=samples,main=paste(s$name,"; generation ",gen,sep=""))
      next()
    }

    # Save
    if(action$action=="sv" | action$action=="save"){
      if(is.na(action$value)) action$value=100
      s$generations[[as.character(gen)]] <- subsample(cur,action$value)
      next()
    }

    # Set default sample size to save
    if(action$action=="svd" | action$action=="savedefault" ){
      if(is.na(action$value)) save.default <- 0
      else save.default <- action$value
      next()
    }

    # Snapshot
    if(action$action=="sn" | action$action=="snapshot"){
      if(is.na(action$value)) action$value=10
      if(is.na(action$value2)) action$value2=paste(sub(" ","_",s$name),"_Gen",gen,".png",sep="")
      filetype=sub(".*[.](...)","\\1",action$value2)
      msgn("file=",action$value2)
      msgn("filetype=",filetype)
      if(filetype=="png") png(filename=action$value2,width=1024,height=768)
      else if(filetype=="pdf") pdf(file=action$value2,width=1024,height=768)
      else if(filetype=="jpg") jpeg(filename=action$value2,width=1024,height=768)
      else next()
      plot(cur,samples=action$value,main=paste(s$name,"; generation ",gen,sep=""))
      dev.off()
      next()
    }
    
    # Show history
    if(action$action=="h" | action$action=="history"){
      msgn("==============================================================================")
      print(s$events)
      msgn("==============================================================================")
      next()
    }
    
    # End simulation
    if(action$action=="q" | action$action=="quit"){
      finish=TRUE
      s$events <- add.event(events=s$events,
                            event="finish",
                            generation=gen
                           )
      next()
    }

    # Change verbosity level
    if(action$action=="v" | action$action=="verbosity"){
      if(is.na(action$value) | !is.numeric(action$value)) v <- 2
      else  v <- action$value
      next()
    }
    
    # Set growth rate
    if(action$action=="r" | action$action=="growthrate"){
      if(is.na(action$value) | !is.numeric(action$value)) r <- 1
      else r=action$value
      s$events <- add.event(events=s$events,
                            event="setrate",
                            generation=gen,
                            rate=r)
      next()
    }

    # Admix
    if(action$action=="a" | action$action=="admix"){
      if(is.na(action$value) | !is.numeric(action$value)) action$value=10
      if(is.na(action$value2)| !is.numeric(action$value2)) action$value2=1
      cur <- admix(cur,N=action$value,type=action$value2)
      s$events <- add.event(events=s$events,
                            event="admix",
                            generation=gen,
                            type=action$value2,
                            size=action$value)
      target.size=target.size+action$value
      s$generations[[as.character(gen)]] <- subsample(cur,save.default)
      plot(cur,samples=samples,main=paste(s$name,"; generation ",gen,sep=""))
      next()
    }

    # Set sample size
    if(action$action=="s"){
      if(is.na(action$value) | !is.numeric(action$value)) action$value=0
      samples=action$value
      plot(cur,samples=samples,main=paste(s$name,"; generation ",gen,sep=""))
      next()
    }

    # show help
    if(action$action=="?" | TRUE){
      cat(help.message)
      next()
    }
  }
  if(file != "") close(connection)
  msg("Simulation:",s$name," ")
  size(s)
  return(s)
}

# RECSIM
#
################################################################################

################################################################################
#
# PRINT

print.recsim.individual <- function(ind){
  msg("Individual: ")
  if( length(ind) == 2 ){
    msgn("Pure type: ",0)
    return(invisible())
  }
  if( length(ind) == 4 &
      ind[2] == 0 &
      ind[3] == 1){
    msgn("Pure type: ",1)
    return(invisible())
  }
  n <- length(ind)
  msg("Type: ")
  for(i in 2*(1:((n - 2)/2))){
    msg(ind$genome[i],"--",ind$genome[i+1],", ")
  }
  msgn("")
}

print.recsim.population <- function(p){
  msgn("Population: Sample size: ",length(p))
}

print.ranges <- function(x){
  y <- ""
  x <- sort(x)
  last <- length(x)
  xb <- x[1]
  if(length(x)>1){
    for(i in 2:length(x)){
      if(x[i]-x[i-1] > 1){
        if(xb==x[i-1])
          y <- paste(y,xb,"; ",sep="")
        else
          y <- paste(y,xb,"--",x[i-1],"; ",sep="")
        xb <- x[i]
      }
    }
      if(xb==x[last])
        y <- paste(y,xb,sep="")
      else
        y <- paste(y,xb,"--",x[last],sep="")
      xb <- x[i]
  }
  else y <- paste(xb,sep="")
  return(y)
}
   
print.recsim <- function(s){
  msgn("Simulation:\n\t\"",s$name,"\"; \"",s$comment,"\"")
  msg("Size:\n\t")
  size(s)
  print(s$recmap)
  msg("Generations:\n\t")
  msgn(dim(s$info)[1])
  msg("Sampled generations:\n\t")
  msgn(print.ranges(as.integer(names(s$generations))))
  msgn()
  print(s$events)
}



print.recsim.action <- function(action){
  msg("action: ")
  for(i in 1:length(action)){
    if(is.na(action[[i]])) break()
    if(is.character(action[[i]]) & matches("[[:space:]]",action[[i]]))
      msg("\"",action[[i]],"\" ")
    else
      msg(action[[i]]," ")
  }
  msgn("")
}

print.recsim.recmap <- function(recmap){
  msgn("Recombination map:")
  msgn("\t","Parameters:\tGiven:\tEmpirical:")
  msgn("\t","rate\t\t",recmap$rate,"\t",recmap$rate)
  msgn("\t","Bkgd\t\t",recmap$bkgd,"\t",recmap$bkgd)
  l <- length(recmap$hotspots)
  sp <- recmap$hotspots - c(0,recmap$hotspots[1:(l-1)])
  sp.mean <- mean(sp)
  msgn("\t","space\t\t",recmap$space,"\t",sp.mean)
  msgn("\t","distance.shape\t",recmap$distance.shape,"\t",sp.mean^2/var(sp))
  int <- recmap$intensity - c(0,recmap$intensity[1:(l-1)])
  msgn("\t","intensity.shape\t",recmap$intensity.shape,"\t",mean(int)^2/var(int))
  w.mean <- mean(recmap$sigma)
  msgn("\t","width\t\t",recmap$width,"\t",w.mean)
  msgn("\t","width.shape\t",recmap$width.shape,"\t",w.mean^2/var(recmap$sigma))
  msgn("\t","no.hotspots\t",floor(1/recmap$space),"\t",length(recmap$hotspots))
}

print.recsim.event <- function(e,j="",tab="\t\t"){
  if(j=="") prefix=tab
  else prefix=paste(tab,j,". ",sep="")
  if(length(e$event)==0){
    print(unclass(e))
    return(invisible(0))
  }
  if(e$event == "create"){
    msgn(prefix,
         "Create population \"",
         e$name,
         "\" of size ",
         e$size,
         " and type ",
         e$type
        )
    msgn(tab,"\t\"",e$comment,"\"")
    return(invisible())
  }

  if(e$event == "admix"){
    if(e$size %% 10 == 1) suffix=""
    else suffix="s"
    msgn(prefix,
         "Admix ",
         e$size,
         " individual",
         suffix,
         " of type ",
         e$type
        )
    return(invisible())
  }
    
  if(e$event == "setrate"){
    msgn(prefix,
         "Set growth rate to ",
         e$rate
        )
    return(invisible())
  }

  if(e$event == "recmap"){
    msgn(prefix,
         "Create recombination map")
    return(invisible())
  }

  if(e$event == "bottleneck"){
    msgn(prefix,
         "Bottleneck: ",
         e$size,
         " individuals"
        )
    return(invisible())
  }
  
  if(e$event == "extinct"){
    msgn(prefix,
         "Population is extinct.")
    return(invisible())
  }

  if(e$event == "finish"){
    msgn(prefix,
         "End simulation")
    return(invisible())
  }
  
  msg(prefix)
  for(n in names(e)) msg(n,": ",e[[n]],"; ")
  msgn("")
}

print.recsim.events <- function(events){
  #  msgn("Generations:\t",length(events))
  msgn("EVENTS:")
  current.generation <- -99
  j=1
  l <- length(events)
  if(l==0) return(invisible())
  for(i in 1:length(events)){
    if(length(events[[i]]) == 0){
      warning("Event ",i," is empty.")
      next()
    }
    if(events[[i]]$generation != current.generation){
      if(events[[i]]$generation < current.generation)
        warning("Events are not in order!")
      current.generation <- events[[i]]$generation
      msgn("\t","Generation: ",current.generation)
      j <- 1
    }
    print(events[[i]],j=j,tab="\t\t")
    j <- j+1
  }
}

# PRINT
#
################################################################################

################################################################################
#
# PLOT SIM

plot.chromosome <- function(chromo,new=TRUE){
  fill.color=rgb(255,0,0,alpha=100,maxColorValue=255)
  line.color=rgb(255,0,0,alpha=200,maxColorValue=255)
  color0="green"
  color1="blue"
  if(new) plot.new()
  plot.window(c(0,1),c(-0.05,1.05))
  lines(c(0,1),c(-0.02,-0.02),lwd=3,col=color0)
  lines(c(0,1),c(1.02,1.02),  lwd=3,col=color1)
  l<-length(chromo)
  oldx=0
  if(l>2){
    for( i in 2*(1:((l-2)/2)) ){
      polygon(c(chromo[i],chromo[i],chromo[i+1],chromo[i+1]),
              c(0,1,1,0),
              col=fill.color,
              border=NA
              )
      lines(x=c(oldx,chromo[i],chromo[i],chromo[i+1],chromo[i+1]),
            y=c(0,0,1,1,0),
            col=line.color,
            lwd=1
            )
      oldx <- chromo[i+1]
    }
  }
  lines(c(oldx,1),c(0,0),lwd=1,col=line.color)
  box()
}

plot.recsim.population <- function(p,
                                   samples=10,
                                   main=""
                                  ){
  if(samples<=0) return(invisible())
  samples=min(samples,length(p))
  s <- sample(length(p),samples)
  par(mfrow=c(samples + 1,1),mar=c(0,0,0,0))
  plot.new()
  plot.window(c(0,1),c(0,1))
  text(x=0.5,y=0.5,adj=c(0.5,0.5),col="blue",cex=3,labels=main)
  for(c in s) plot.chromosome(p[[c]])
}
  

# PLOT SIM
#
################################################################################

################################################################################
#
# STATS

#find lists of lengths of white and red intervals for an individual genome
intervals <- function(chromo,check=FALSE){
  l <- length(chromo)
  chromo[1] <- 0
  chromo[l] <- 1
  gaps <- chromo[2:l] - chromo[1:(l-1)]
  h <- gaps[2*(1:(l/2))-1]
  h <- h[h!=0]
  if(l>2) f <- gaps[2*(1:(l/2-1))]
  else f <- numeric(0)
  if(check) if(sum(h)+sum(f)!=1) warning("Intervals: inconsistency")
  return(list(h=h,f=f))
}

recsim.stats <- function(s,
                         blockwidththreshold=0,
                         admixed.only=TRUE,
                         check=FALSE,
                         verbose=1
                        ){
  result <- list()
  result$name <- s$name
  result$comment <- s$comment
  result$recmap <- s$recmap
  result$info <- s$info
  result$events <- s$events
  result$sampled.generations <- as.numeric(names(s$generations))
  total.gen <- dim(result$info)[1]
  result$ind.stats <- make.array(numeric(0),
                                c("generation",
                                  "chromosome",
                                  "breakpoints",
                                  "hosttotal",
                                  "hostblockwidth",
                                  "foreigntotal",
                                  "foreignblockwidth"
                                 )
                               )
  result$gen.stats <- make.array(result$sampled.generations,
                                 c("generation",
                                   "sampled",
                                   "admixed",
                                   "breakpoints",
                                   "hosttotal",
                                   "hostblockwidth",
                                   "foreigntotal",
                                   "foreignblockwidth"
                                  )
                                ) 

  for(g in result$sampled.generations){
    gg <- as.character(g)
    sampled <- length(s$generations[[gg]])
    ind.stats <- make.array(1:sampled,colnames(result$ind.stats))
    for(c in 1:sampled){
      chromo <- s$generations[[gg]][[c]]
      l <- length(chromo)
      if(check){
        if(l%%2 != 0) warning("Length of chromo is not even")
        if(is.unsorted(chromo)) warning("Chromo is unsorted")
      }
      inter <- intervals(chromo)
      inter$f <- inter$f[inter$f >= blockwidththreshold]
      inter$h <- inter$h[inter$h >= blockwidththreshold]
      ind.stats[c,"generation"] <- g
      ind.stats[c,"chromosome"] <- c
      ind.stats[c,"breakpoints"] <- how.many( chromo!=0 & chromo!=1 ) - 2
      ind.stats[c,"hosttotal"] <- sum(inter$h)
      ind.stats[c,"hostblockwidth"] <- if(length(inter$h)==0) 0 else mean( inter$h )
      ind.stats[c,"foreigntotal"] <- sum(inter$f)
      ind.stats[c,"foreignblockwidth"] <- if(length(inter$f)==0) 0 else mean( inter$f )
    }
    result$ind.stats <- rbind(result$ind.stats,ind.stats)
    result$gen.stats[gg,"generation"] <- g 
    result$gen.stats[gg,"sampled"]    <- sampled
    result$gen.stats[gg,"admixed"]    <- how.many(ind.stats[,"hosttotal"] < 1)/sampled
    if(admixed.only) range <- (ind.stats[,"breakpoints"] != 0)
    else             range <- 1:sampled
    result$gen.stats[gg,"breakpoints"] <- mean(ind.stats[range,"breakpoints"])
    result$gen.stats[gg,"hosttotal"] <- mean(ind.stats[range,"hosttotal"])
    result$gen.stats[gg,"hostblockwidth"] <- mean(ind.stats[range,"hostblockwidth"])
    result$gen.stats[gg,"foreigntotal"] <- mean(ind.stats[range,"foreigntotal"])
    result$gen.stats[gg,"foreignblockwidth"] <- mean(ind.stats[range,"foreignblockwidth"])
  }
  class(result) <- "recsim.stats"
  return(result)
}

plot.recsim.stats <- function(st,
                             main=paste("\"",st$name,"\" \"",st$comment,"\"",sep=""),
                             sub="",
                             log=TRUE,
                             what="?"
                            ){
  if(what=="all") what <- c("breakpoints", "hostblockwidth", "foreigntotal", "foreignblockwidth")

  l <- length(what)
  if(l > 1){
    plot.new()
    par(mfrow=c((l+1)%/%2,2),mar=c(2,3.5,1.9,0.5)) #bot,left,top,right
    for(w in what) plot(st,main=paste(main,what,sep="; "),sub=sub,log=log,what=w)
    return(invisible())
  }
  
  if(log) log="x"
  else    log=""

  
  if(what == "sampled"){
    plot(x=st$gen.stat[,"generation"],
         y=st$gen.stat[,"sampled"],
         col="red",
         type="l",
         log=log,
         xlab="",
         ylab=""
        )
    title(main=main,sub=sub)
    title(xlab="generations",line=1.9)
    title(ylab="sampled",line=1.9)
    return(invisible())
  }
  if(what == "admixed"){
    plot(x=st$gen.stat[,"generation"],
         y=st$gen.stat[,"admixed"]*100,
         col="red",
         type="l",
         log=log,
         xlab="",
         ylab=""
        )
    title(main=main,sub=sub)
    title(xlab="generations",line=1.9)
    title(ylab="admixed",line=1.9)
    return(invisible())
  }
  if(what == "breakpoints"){
    plot(x=st$ind.stat[,"generation"],
         y=st$ind.stat[,"breakpoints"],
         log=log,
         xlab="",
         ylab="")
    title(main=main,sub=sub)
    title(xlab="generations",line=1.9)
    title(ylab="breakpoints",line=1.9)
    lines(x=st$gen.stat[,"generation"],
          y=st$gen.stat[,"breakpoints"],
          col="red")
    return(invisible())
  }
  if(what == "hosttotal"){
    plot(x=st$ind.stat[,"generation"],
         y=st$ind.stat[,"hosttotal"],
         log=log,
         xlab="",
         ylab="")
    title(main=main,sub=sub)
    title(xlab="generations",line=1.9)
    title(ylab="Total amount of foreign genome",line=1.9)
    lines(x=st$gen.stat[,"generation"],
          y=st$gen.stat[,"hosttotal"],
          col="red")
    return(invisible())
  }
  if(what == "hostblockwidth"){
    plot(x=st$ind.stat[,"generation"],
         y=st$ind.stat[,"hostblockwidth"],
         log=log,
         xlab="",
         ylab="")
    title(main=main,sub=sub)
    title(xlab="generations",line=1.9)
    title(ylab="Ave width of host blocks",line=1.9)
    lines(x=st$gen.stat[,"generation"],
          y=st$gen.stat[,"hostblockwidth"],
          col="red")
    return(invisible())
  }
  if(what == "foreigntotal"){
    plot(x=st$ind.stat[,"generation"],
         y=st$ind.stat[,"foreigntotal"],
         log=log,
         xlab="",
         ylab="")
    title(main=main,sub=sub)
    title(xlab="generations",line=1.9)
    title(ylab="Total amount of foreign genome",line=1.9)
    lines(x=st$gen.stat[,"generation"],
          y=st$gen.stat[,"foreigntotal"],
          col="red")
    return(invisible())
  }
  if(what == "foreignblockwidth"){
    plot(x=st$ind.stat[,"generation"],
         y=st$ind.stat[,"foreignblockwidth"],
         log=log,
         xlab="",
         ylab="")
    title(main=main,sub=sub)
    title(xlab="generations",line=1.9)
    title(ylab="Ave width of foreign blocks",line=1.9)
    lines(x=st$gen.stat[,"generation"],
          y=st$gen.stat[,"foreignblockwidth"],
          col="red")
    return(invisible())
  }
  msgn("Say what to plot!\nplot(st,...,what=what,...)
         what = \"sampled\", \"admixed\", \"breakpoints\", \"hosttotal\",
         \"hostblockwidth\", \"foreigntotal\", \"foreignblockwidth\" or \"all\""
      )
  
}

# STATS
#
################################################################################

################################################################################
#
# WT

wt.chromo <- function(chromo,levels=17){
  # levels should be at least 2
  levels <- max(levels,2)

  # fill empty wt
  wt <- list()
  class(wt) <- "wavelet"
  power2=1
  for(i in 1:levels){
    wt[[i]] <- rep(0,power2)
    power2 <- power2*2
  }
  #strip boundaries (-1 und 99)
  l=length(chromo)
  if( l==2 ){
    wt$average <- 0
    return(wt)
  }
  chromo <- chromo[2:(l-1)]
  l <- l-2
  evenrange <- 2*(1:(l/2))
  wt$average <- sum(chromo[evenrange]-chromo[evenrange-1])

  #the last 1 should be removed
  if(chromo[l]==1) chromo <- chromo[1:(l-1)]
  
  # cycle through bpts and wt levels.
  # each bpt contributes to only one coeff in each level
  sign <- +1
  for(bpt in chromo){
    power2 <- 1 # =(1/2)^(level-1)
    for(level in 1:levels){
      #which coeff
      k=floor(bpt/power2)+1
      wt[[level]][k] <-
        wt[[level]][k] +
        sign * ( abs(bpt/power2 - k + 1/2) - 1/2 )
      power2 <- power2/2
    }
    sign <- -sign
  }
  
  return(wt)
}

wt.cutoff <- function(wt,threshold){
  level=length(wt)-1
  for (i in 1:level){
    wt[[i]][abs(wt[[i]])<=threshold] <- 0
  }
  return(wt)
}


wt.summary <- function(wt,p=1,maxlevel=length(wt)-1){
  l <- length(wt)-1
  s <- make.array(c(1:maxlevel,"average","center"))
  for(i in 1:maxlevel) s[i] <- ( sum(abs(wt[[i]])^p)^(1/p) ) / (2^(i-1))
  s["average"] <- wt$average
  s["center"] <- sum(s[1:maxlevel]*(1:maxlevel))/sum(s[1:maxlevel])
  class(s) <- "wavelet.summary"
  return(s)
}

plot.wavelet.summary <- function(s,main=""){
  l <- length(s)-1
  barplot(s[1:l],main=main)
  title(sub=paste("center=",round(s["center"]*10000)/10000,sep=""),line=0)
  lines(c(s["center"],s["center"]),c(0,1),col="red",lwd=3)
}


recsim.wtstats <- function(s,
                           verbose=1,
                           threshold=0,
                           maxlevel
                          ){
  result <- list()
  class(result) <- "recsim.wtstats"
  result$name <- s$name
  result$comment <- s$comment
  result$recmap <- s$recmap
  result$info <- s$info
  result$events <- s$events
  result$sampled.generations <- as.numeric(names(s$generations))
  total.gen <- dim(result$info)[1]
  result$ind.wtstats <- make.array(c(),
                                     c("generation","No",1:maxlevel,"average","center"))
  result$gen.wtstats <- make.array(result$sampled.generations,
                                   c("generation",
                                     1:maxlevel,
                                     "average",
                                     "center")
                                  )
  for(g in result$sampled.generations){
    if( verbose>=1 ) msg("Generation: ",g)
    gg <- as.character(g)
    sampled <- length(s$generations[[gg]])
    ind.wtstats <- make.array(1:sampled,colnames(result$ind.wtstats))
    for(c in 1:sampled){
      if(verbose >=2) tick()
      chromo <- s$generations[[gg]][[c]]
      wt <- wt.chromo(chromo,levels=maxlevel)
      wt <- wt.cutoff(wt,threshold=threshold)
      wt.summary <- wt.summary(wt,maxlevel=maxlevel)
      ind.wtstats[c,] <- c(g,c,wt.summary)
    }
    result$ind.wtstats <- rbind(result$ind.wtstats,ind.wtstats)
    gen.wtstats <- apply(ind.wtstats[,c(as.character(1:maxlevel),"average","center")],
                         2,
                         mean,
                         na.rm=TRUE
                        )
    result$gen.wtstats[gg,] <- c(g,gen.wtstats)
    if( verbose >=1 ) msgn()
  }
  return(result)
}

plot.recsim.wtstats <- function(wts,
                         main=paste("\"",wts$name,"\" \"",wts$comment,"\"",sep=""),
                         sub="",
                         log=TRUE,
                         what="?"
                        ){
  if(what=="all") what <- c("centers","averages")
  l <- length(what)
  if(l > 1){
    plot.new()
    par(mfrow=c((l+1)%/%2,2),mar=c(2,3.5,1.9,0.5)) #bot,left,top,right
    for(w in what) plot(wts,main=paste(main,w,sep="; "),sub=sub,log=log,what=w)
    return(invisible())
  }

  if(log) log="x"
  else    log=""

  if(what=="centers"){
    plot(x=wts$ind.wtstats[,"generation"],
         y=wts$ind.wtstats[,"center"],
         log=log,
         xlab="",
         ylab=""
        )
        title(main=main,sub=sub)
    title(xlab="generations",line=1.9)
    title(ylab="wt centers",line=1.9)
    lines(x=wts$gen.wtstats[,"generation"],
          y=wts$gen.wtstats[,"center"],
          col="red"
         )
    lines(x=wts$gen.wtstats[,"generation"],
          y=wts$gen.wtstats[,"a.center"],
          col="blue"
         )
  }
}
# WT
#
################################################################################
