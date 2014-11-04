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
# Auxilary stuff
#

msg <- function(...,sep="",new.line=FALSE){
  cat(paste(...,sep=sep))
  if(new.line) cat("\n")
}

msgn <- function(...,sep="",new.line=TRUE){
  msg(...,sep=sep,new.line=new.line)
}

msgc <- function(...,sep="",width=80,fill=" ",new.line=FALSE){
  space=""
  for( i in 1:width) space=paste(space,fill,sep="")
  msg("\r",space,"\r")
  msg(...,sep=sep,new.line=new.line)
}

tick <- function(point="."){
  cat(point)
}

pause <- function(){
  msg("PAUSE! press \"enter\" to continue:")
  scan(file="",nmax=1,quiet=TRUE,nlines=1)
}

######

make.array <- function(...,data=NA){
  DimNames   <- list(...)
  rank       <- length(DimNames)
  TotalSize  <- 1
  Dim        <- vector(mode="integer",length=rank)
  for(i in 1:rank){
    Dim[i]   <- length(DimNames[[i]])
    DimNames[[i]] <- as.character(DimNames[[i]])
    TotalSize <- TotalSize * Dim[i]
  }
  array(data=data, dim=Dim, dimnames = DimNames)
} 

######

how.many <- function(x){
  return(length(which(x)))
}

matches <- function(ptrn,str,ignore.case=TRUE){
  return(1:length(str) %in% grep(ptrn,str,ignore.case=ignore.case))
}


size <- function(object,newline=TRUE){
  S <- object.size(object)
  m <- floor(log2(S)/10)
  S <- floor(S*10/(1024^m))/10
  suffix=c("b","Kb","Mb","Gb","Tb")
  cat(paste(S,suffix[m+1],sep=""))
  if(newline) cat("\n")
}

add.alpha <- function(col,alpha=1)
  return(rgb(t(col2rgb(col)),alpha=alpha*255,maxColorValue=255))


# AUX STUFF
#
################################################################################

################################################################################
#
# AUX CALCULATIONS




Norm <- function(v,p=2){
  return(( if(p==Inf) max(abs(v)) else (sum(abs(v)^p))^(1/p) )* length(v)^(1-1/p))
}


# Stats used in windowsize decisions
# T - table, A - axis
# returns mean, var, table cols averaged with weights in axis
  
stats <- function(T,A){
  stats <- list()
  people <- colnames(T)
  values <- numeric(length(people))
  names(values) <- people
  for(p in 1:length(people)){
      values[p] <- sum(T[,p]*A,na.rm=TRUE)/Norm(A[which(!is.na(T[,p]))],p=1)
    }
  
  stats$mean  <- mean(values,na.rm=TRUE) 
  stats$sigma <- var (values,na.rm=TRUE)
  v <- values - stats$mean
  vp <- v[v>0]
  vm <- v[v<0]
  lp <- length(vp)
  lm <- length(vm)
  stats$values <- values
  return(stats)
}

# Nice graphical representation of stats()
plot.stats <- function(s1,s2,main="",sub=""){
  solid1="red"
  solid2="blue"
  alpha=0.1
  transp1=add.alpha(solid1,alpha)
  transp2=add.alpha(solid2,alpha)
  plot.new()
  plot.window(c(-1,1),c(-1,1))
  m1 <- s1$mean
  sig1 <- sqrt(s1$sigma)
  m2 <- s2$mean
  sig2 <- sqrt(s2$sigma)
  v1 <- s1$values
  v2 <- s2$values
  points(rep( 0.02,length(v1)),v1,pch=15,col=solid1)
  points(rep(-0.02,length(v2)),v2,pch=20,col=solid2)
  lines(c(-0.1,0.6),c(m1,m1),col=solid1)
  lines(c(0.5,0.5),c(m1-sig1,m1+sig1),col=solid1,lwd=5)
  rect(xleft=-0.1,xright=0.7,ybottom=m1-sig1,ytop=m1+sig1,border=NA,col=transp1)
  lines(c(0.5,0.5),c(m1-2*sig1,m1+2*sig1),col="red",lwd=3)
  rect(xleft=-0.1,xright=0.7,ybottom=m1-2*sig1,ytop=m1+2*sig1,border=NA,col=transp1)
  lines(c(0.5,0.5),c(m1-3*sig1,m1+3*sig1),col="red",lwd=1)
  rect(xleft=-0.1,xright=0.7,ybottom=m1-3*sig1,ytop=m1+3*sig1,border=NA,col=transp1)
  lines(c(-0.6,0.1),c(m2,m2),col="blue")
  lines(c(-0.5,-0.5),c(m2-sig2,m2+sig2),col="blue",lwd=5)
  rect(xleft=-0.6,xright=0.1,ybottom=m2-sig2,ytop=m2+sig2,border=NA,col=transp2)
  lines(c(-0.5,-0.5),c(m2-2*sig2,m2+2*sig2),col="blue",lwd=3)
  rect(xleft=-0.6,xright=0.1,ybottom=m2-2*sig2,ytop=m2+2*sig2,border=NA,col=transp2)
  lines(c(-0.5,-0.5),c(m2-3*sig2,m2+3*sig2),col="blue",lwd=1)
  rect(xleft=-0.6,xright=0.1,ybottom=m2-3*sig2,ytop=m2+3*sig2,border=NA,col=transp2)
  title(main=main,sub=sub)
}
  

######
# x,X,y,Y are min/MAX of the ellipse wrt x-, y-coords
# returns a list of vertices of approximating polygon with 4n vertices
ellipse <- function(x,X,y,Y,n=16){
  cx <- (x+X)/2
  cy <- (y+Y)/2
  wx <- abs(x-X)/2
  wy <- abs(y-Y)/2
  angles <- (pi/2)*((1:n)/n)
  xx <- wx*cos(angles)
  yy <- wy*sin(angles)
  ell <- list(x=c(cx + xx, cx - xx[(n-1):1], cx - xx[2:n], cx + xx[(n-1):1]),
              y=c(cy + yy, cy + yy[(n-1):1], cy - yy[2:n], cy - yy[(n-1):1])
             )
  class(ell) <- "ellipse"
  return(ell)
}

plot.ellipse <- function(ell,col="red",alpha=0.5){
  col <- add.alpha(col,alpha)
  polygon(ell$x,ell$y,col=col,border=col)
}
  
# AUX CALCULATIONS
#
################################################################################


################################################################################
#
# PCO

pco <- function(table, # SNPS x PEOPLE matrix
                who=colnames(table),
                axis.who=who,
                normalize=FALSE,
                no.axes=2,
                nas="0",   # what to do with NAs when calculating axis in SNP coords; another possibility is "ignore"
                estimates=TRUE,  #to find estimated coords for all possible NA values
                verbose=2
               ){
  # Prepare variables
  people <- colnames(table)
  names(people) <- people
  axis.who <- people[axis.who]
  who <- people[who]
  all <- unique(c(who,axis.who))
  if(length(all)<length(people)) table <- table[,all]
  snps <- rownames(table)
  no.snps   <- length(snps)
  result <- list()
  class(result) <- "pco"

  # check whether we do NA estimates
  nnas <- how.many(is.na(table))
  if( estimates & (nnas == 0) & (verbose >= 2) )
    msgn("No NA's. Will not do estimates.")
  estimates <- estimates & (nnas > 0)
  if(estimates) etable <- table[,all]
  
  # prepare return object
  result$estimates <- estimates
  result$normalize <- normalize
  result$call <- sys.call()
  result$date <- date()
  result$table <- match.call()$table
  result$people <- who
  result$axis.people <- axis.who
  result$no.snps <- no.snps
  result$pap  <- make.array(people,1:no.axes)
  result$pa   <- make.array(snps,1:no.axes)
  result$li   <- make.array(people,1:no.axes)
  if(estimates) result$li.upper <- result$li.lower <- make.array(people,1:no.axes)
  result$info <- "
    call      --  function call.
    table     --  name of R object containing table
                  from which PCO was calculated.
    date      --  date
    people    --  colnames of columns from which projections on PC
                  where calculated.
    axis.people -- these colnames where taken into account for finding PC
                   axis.
    no.snps   -- number of SNPs in the table
    estimates --  whether li.upper and li.lower were calculated
    normalize -- whether table rows were normalized to have variance 1
    pap       --  columns are principal axes in \"people coordinates\"
    pa        --  columns are principal axes in \"snp coordinates\"
    li        --  columns are principal coordinates
    li.upper
    li.lower  
    G         --  matrix of quadratic form in \"people coordinates\"
    ev        --  vector of eigenvalues
"

  # START
  if(verbose >= 2){
    msgn("Number of people:\t", length(who))
    msgn("Number of people\n for finding axes:\t", length(axis.who))
    msgn("Number of SNPs:\t\t",   no.snps)
  }

  #Put origin at the center of mass
  if(verbose>=2) msg("Centering...")
  no.checks <- 5
  center.threshold <- 1/(10*length(axis.who))
  checks <- sample(no.snps,no.checks)
  a=0
  for(i in checks) a <- a + abs(mean(table[i,axis.who],na.rm=TRUE))/no.checks

  if( a <= center.threshold ){
    if(verbose>=2) msg("Table is already centered...")
  }
  else{
    table <- t(apply(table,1,function(x) return(x-mean(x[axis.who],na.rm=TRUE))))
  }
  if(verbose>=2) msgn("Done")
  
  # Normalize
  if(normalize){
    if(verbose>=2) msg("Normalizing...")
    table <- t(apply(table,1,function(x) return(x/sqrt(var(x[axis.who],na.rm=TRUE)))))
    if(verbose>=2) msgn("Done")
  }
  
  #Find matrix of scalar products
  if(verbose>=2) msg("Calculate scalar products...")
  result$G <- var(table[,axis.who],na.rm=TRUE)
  if(verbose>=2) msgn("Done")
  
  # Find eigenvectors in "People coords"
  if(verbose>=2) msg("Find spectral data...")
  spec       <- eigen(result$G,symmetric=TRUE)
  result$ev  <- spec$values[1:no.axes]
  result$pap <- spec$vectors[,1:no.axes]
  if(verbose>=2) msgn("done")
  
  # Find eigenvectors in "snp coords"
  # replace na's by 0
  if(verbose>=2) msg("Eigenvectors in SNPs coordinates.")
  tbl.a <- table[,axis.who]
  if(as.character(nas) == "0") tbl.a[which(is.na(tbl.a))] <- 0
  for(i in 1:no.axes){
    result$pa[,i] <- apply(tbl.a,
                           1,
                           function(x) return(sum(x*result$pap[,i],na.rm=TRUE))
                           )
    result$pa[,i] <- result$pa[,i]/sqrt(sum(result$pa[,i]^2))
    if(verbose>=1) tick()
  }
  rm(tbl.a)
  if(verbose>=2) msgn("Done")

  # Find principal components
  if(verbose>=2) msg("Find principal components.")
  for(i in 1:no.axes){
    result$li[,i] <- apply(table[,who],
                           2,
                           function(x) return(sum(x*result$pa[,i],na.rm=TRUE))
                          )/sqrt(no.snps)
    if(verbose>=1) tick()
  }
  if(verbose>=2) msgn("Done")

  
  # Upper/lower estimates
  if(estimates){
    if(verbose>=2) msg("Find upper/lower estimates for missing values..")
    no.who <- length(who)
    for(i in 1:no.axes){
      table.u <- t(apply(cbind(table[,who],result$pa[,i]),
                         1,
                         function(x){
                           if(x[no.who+1]>0) fcn <- max
                           else              fcn <- min
                           x[is.na(x)] <- fcn(x[1:no.who],na.rm=TRUE)
                           return(x[1:no.who])
                         }
                        ))
      if(verbose>=1) tick()
      table.l <- t(apply(cbind(table[,who],result$pa[,i]),
                         1,
                         function(x){
                           if(x[no.who+1]>0) fcn <- min
                           else              fcn <- max
                           x[is.na(x)] <- fcn(x[1:no.who],na.rm=TRUE)
                           return(x[1:no.who])
                         }
                        ))
      if(verbose>=1) tick()
      result$li.upper[,i] <- apply(table.u,
                                   2,
                                   function(x) return(sum(x*result$pa[,i]))
                                  )/sqrt(no.snps)
      if(verbose>=1) tick()
      result$li.lower[,i] <- apply(table.l,
                                   2,
                                   function(x) return(sum(x*result$pa[,i]))
                                  )/sqrt(no.snps)

      if(verbose>=1) tick()
    }
    if(verbose>=2) msgn("Done")
  }

  if(verbose>=2){
    msgn("Done")
    msgn(result$info)
    msg("The resulting object is ")
    size(result)
  }
 
  return(result)
}


#
#
################################################################################

################################################################################
#
# Plot pco

plot.pco <- function(pco,
                     who=rownames(pco$li),
                     names="full",
                     dots=TRUE,
                     spots=TRUE,
                     level=1,
                     alpha=0.2,
                     main="", #paste(pco$table,pco$date),
                     sub="",
                     colors="colorfile.txt"  #filename, or table of colors and pch
                    ){
  people <- rownames(pco$li)
  names(people) <- people
  people <- people[who]
  estimates <- spots & pco$estimates
  xmax <- max(pco$li[people,level],
              pco$li.upper[people,level],
              pco$li.lower[people,level])
  xmin <- min(pco$li[people,level],
              pco$li.upper[people,level],
              pco$li.lower[people,level])
  ymax <- max(pco$li[people,level+1],
              pco$li.upper[people,level+1],
              pco$li.lower[people,level+1])
  ymin <- min(pco$li[people,level+1],
              pco$li.upper[people,level+1],
              pco$li.lower[people,level+1])


  if(length(colors)==1){
    as.vector(read.table(colors)[people,"color"])->my.colors
  }
  else{
    colors[people,"color"]->my.colors
  }

  names(my.colors) <- people
  text=TRUE

  if(class(names)=="function") shortpeople <- names(people)
  else{
    if(names=="empty" | names=="no" | names=="none") text=FALSE
    if(names=="full")  shortpeople <- people
    if(names=="short") shortpeople <- sub("([a-Z]*)_(.*)","\\2",people)
  }
  
  plot.new()
  plot.window(c(xmin,xmax),c(ymin,ymax))
  box()
  title(main=main,
        sub=sub,
        xlab=paste("PC",level,sep=""),
        ylab=paste("PC",level+1,sep="")
        )
  axis(1)
  axis(2)
  axis(3)
  axis(4)
  if(text){
    text(pco$li[people,level],
         pco$li[people,level+1],
         labels=shortpeople,             
         adj=c(0.5,0.5),                 
         cex=0.8,                        
         col=my.colors                   
	)
  }
  if(dots){
    points(pco$li[people,level],pco$li[people,level+1],  
           pch=19,
           col=my.colors,
          )
  }
  if(estimates){
    for(p in people){
      plot(ellipse(pco$li.lower[p,level],
                   pco$li.upper[p,level],
                   pco$li.lower[p,level+1],
                   pco$li.upper[p,level+1]
                  ),
           col=my.colors[p],
           alpha=alpha
          )
    }
  }
}


print.pco <- function(pco){
  msgn("Principal Component Analysis")
  msgn("----------------------------")
  msgn("Table:\t\t",pco$table)
  msgn("Created:\t",pco$date)
  msgn("Number of axes:\t",length(pco$ev))
  msgn("People:\t\t",length(pco$people))
  msgn("People used\n   for axis:\t",length(pco$axis.people))
  msgn("Number of SNPs:\t", pco$no.snps)
  msgn("Estimates:\t",pco$estimates)
  msgn("Normalize:\t",pco$normalize)
  msgn("----------------------------")
}
#
#
################################################################################





################################################################################
#
# GENETIC MAP

# genetic.map - Nx2 table, phys.pos - vector

get.contiguous.regions <- function(v){
  regions <- list()
  i <- 1
  V <- cbind(1:length(v),v)
  while( (s <- suppressWarnings(min(which(is.na(V[,2]))))) != Inf ){
    V <- V[s:dim(V)[1],,drop=FALSE]
    e <- suppressWarnings(min(which(!is.na(V[,2]))))
    if(e==Inf){
      regions[[i]] <- c(V[1,1]-1,Inf)
      break()
    }
    s <- V[1,1]-1
    V <- V[e:dim(V)[1],,drop=FALSE]
    e <- V[1,1]
    regions[[i]] <- c(s,e)
    i <- i+1
  }
  return(regions)
}

interpolate.genetic.map <- function(phys.pos,
                                    gen.map,
                                    verbose=2,
                                    method=linear){
  if(is.unsorted(phys.pos)) {
    if(verbose>=2) msgn("Physical positions are unsorted. Sorting...")
    phys.pos <- sort(phys.pos)
  }
  if(is.unsorted(gen.map[,1])){
    if(verbose>=2) msgn("Physical positions in genetic map are unsorted. Sorting...")
    gen.map <- gen.map[order(gen.map[,1]),]
  }
  if(is.unsorted(gen.map[,2])) stop("Genetic map is not monotone. Exiting...")

  new.gen.map <- make.array(phys.pos,c("PhysPos","GenMap"),data=NA)
  new.gen.map[,1] <- phys.pos

  which(phys.pos %in% gen.map[,1])->definedpp
  which( gen.map[,1] %in% phys.pos)->definedgm
  new.gen.map[definedpp,2] <- gen.map[definedgm,2]
  regions <- get.contiguous.regions(new.gen.map[,2])
  if(verbose>=2) msgn("Found ",length(phys.pos)-length(definedpp),
                      " undefined values in ",
                      length(regions)," contiguous regions"
                     )

  #interpolate
  NAS <- which(!is.na(new.gen.map[,2]))
  if( length(NAS) < 2 ) stop("Too few values defined. Cannot interpolate. Exiting...")
  S <- min(NAS)
  E <- max(NAS)
  general.slope <- (new.gen.map[E,2]-new.gen.map[S,2])/(new.gen.map[E,1]-new.gen.map[S,1])
  for(i in 1:length(regions)){
    if(verbose>=1) tick()
    s <- regions[[i]][1]
    e <- regions[[i]][2]
    if(s==0){
      reg <- (s+1):(e-1)
      slope <- new.gen.map[e,2]/new.gen.map[e,1]
      new.gen.map[reg,2] <- (new.gen.map[reg,1])*slope
      next()
    }
    if(e==Inf){
      reg <- (s+1):length(phys.pos)
      e <- length(phys.pos)
      new.gen.map[reg,2] <- (new.gen.map[reg,1]-new.gen.map[s,1])*general.slope + new.gen.map[s,2]
      break()
    }
    reg <- (s+1):(e-1)
    slope=(new.gen.map[e,2]-new.gen.map[s,2])/(new.gen.map[e,1]-new.gen.map[s,1])
    new.gen.map[reg,2] <- (new.gen.map[reg,1]-new.gen.map[s,1])*slope + new.gen.map[s,2]
  }
  return(new.gen.map)
}

plot.genetic.map <- function(gm,main="",sub="")
  plot(x=gm[,1],
       y=gm[,2],
       type="l",
       main=main,sub=sub,
       xlab="Chromosome (bp)",
       ylab="Chromosome (cM)"
       )
plot.genetic.map.density <- function(gm,main="",sub=""){
  a=gm[,2]
  y=(a[2:length(a)]-a[1:(length(a)-1)])/(gm[2:length(a),1]-gm[1:(length(a)-1),1])
  x=(gm[2:length(a)]+a[1:(length(a)-1)])/2
  plot(x=x,
       y=y,
       type="l",
       main=main,sub=sub,
       xlab="Chromosome (bp)",
       ylab="Density (cM/bp)"
       )
}

# GENETIC MAP
#
################################################################################



################################################################################
#
# SPCO
# 
# table        -- table of SNP calls, rows - snps, cols - people
# pco          -- pco ran on the table, if not supplied, then calculated automatically
# phys.pos     -- vector of physical positions of SNPs along a chromosome,
#                 if not given then snps are considered to be uniformly distributed
# who          -- the set of individuals for whom spco is to be calculated, 
#                 if not given, then everybody
# window.size  -- window size, posible values:
#                 numeric -- fixed window size in bp 
#                 "Number bp" -- fixed window size in bp
#                 "Number snp" -- window size chosen to have these many snps
#                 "Number sigma" -- window size chosen to expect a distribution with
#                                 a given variance
#                 "Number wt"  -- window chosen using wavelet transform
# max.window.size -- maximal allowable window size
# step -- step
# remove.fixed.snps -- TRUE/FALSE
# density.threshold -- minimal no of SNPs not to reject a given window
# quality.threshold -- minimal acceptable quality of calls within a window
#                      for an individual
# na.estimates      -- whether estimates for the missing calls should be made
#
#



spco <- function(table,                    # table of calls
                 pco=NA,                   # pco on the table
                 phys.pos=NA,              # vector of phys. positions of SNPs
                 genetic.map=NA,
                 who=colnames(table),      # list of people for whom calc. spco
                 pop1=NA,
                 pop2=NA,
                 window.size = 1000000,    # size of the sliding window in basepairs
                 max.window.size = Inf,
                 step = NA,            # step in basepairs, or cM
                 Nbins = 1024,
                 remove.fixed.snps=TRUE,   # 
                 density.threshold = 1,    # min no. of SNPs in a window
                 quality.threshold = 10,    # min % of no NAs for each person
                 na.estimates=TRUE,        # shall we find upper and lower estimates?
                 normalize=TRUE,          #
                 verbose=2
                 ){

  if(verbose>=2)
    msgn("\n================================================================================")

  # Is genetic map supplied?
  if(!is.na(genetic.map[1])){
    if(verbose>=2) msgn("Use genetic map for physical positions of SNPs.")
    use.gm <- TRUE
    phys.pos <- genetic.map[,2]
    pp.suffix="cM"
  }
  else{
    use.gm <- FALSE
    # Is physical position vector supplied?
    # If not, set it to 1:N
    if(is.na(phys.pos)[1]){
      if(verbose>=2)
        msgn("No physical position vector supplied.")
      phys.pos <- make.array(rownames(table),Data=0)
      phys.pos <- 1:length(phys.pos)
    }
    pp.suffix="bp"
  }
  
  # Calculate principal axis, if not given
  if(is.na(pco[1])){
    if(verbose>=2) msg("Calculating principal axis...")
    pco=pco(table,nas=0,estimates=FALSE,verbose=verbose)
    if(verbose>=2) msgn("Done")
  }
  else if(verbose>=2) msgn("Using supplied principal axis")
  axis=pco$pa[,1]

  # Check whether data is OK
  if(verbose>=2) msg("Check consistency...")
  if( dim(table)[1] != length(axis)     ) stop("Data do no match 1. Exiting...")
  if( length(axis)  != length(phys.pos) ) stop("Data do no match 2. Exiting...")
  if(verbose>=2) msgn("Done")

  # remove fixed snps
  if(remove.fixed.snps){
    if(verbose>=2) msg("Removing fixed snps...")
    not.fixed <- which(axis!=0)
    no.snps.fixed <- length(axis)-length(not.fixed)
    table <- table[not.fixed,,drop=FALSE]
    axis  <- axis[not.fixed,drop=FALSE]
    phys.pos <- phys.pos[not.fixed,drop=FALSE]
    rm(not.fixed)
    if(verbose>=2) msgn("Done")
  }
  
  # Check whether we do NA estimates
  if( na.estimates & (how.many(is.na(table)) == 0) ){
    if(verbose>=2) msgn("No NA's. Will not do upper/lower estimates.")
    na.estimates=FALSE
  }

  
  # Find upper/lower estimates    
  if( na.estimates ){
    if(verbose>=2) msg("Find upper and lower estimates for NA's...")
    pluses    <- which(axis >= 0)
    minuses   <- which(axis <  0)
    nasplus   <- nasminus  <- is.na(table)
    table.u   <- table.l <- table
    nasplus[minuses,] <- FALSE
    nasminus[pluses,] <- FALSE
    table.u[nasplus]  <- 1
    table.u[nasminus] <- -1
    table.l[nasplus]  <- -1
    table.l[nasminus] <- 1
    rm(nasplus,nasminus,pluses,minuses)
    if(verbose>=2) msgn("Done")
  }


  # Prepare aux and output data
  if(verbose>=2) msg("Preparing data...")
  people        <- colnames(table)
  names(people) <- people
  people        <- people[who]
  no.people     <- length(people)
  no.snps       <- length(phys.pos)
  chromolength  <- max(phys.pos)
  if(!is.na(Nbins)) step <- chromolength/(Nbins-1/2)
  no.bins       <- floor(chromolength/step)+1
  coords        <- make.array(people,1:no.bins,data=NA)
  if(na.estimates){
    coords.upper  <- make.array(people,1:no.bins,data=NA)
    coords.lower  <- make.array(people,1:no.bins,data=NA)
  }
  no.md.bins    <- 0
  no.mq.bins    <- 0
  quality       <- make.array(c(people,"Average"),1:no.bins,data=NA)
  density       <- make.array(1:no.bins,data=NA)
  windows       <- make.array(1:no.bins,data=NA)
  windows.snp   <- make.array(1:no.bins,data=NA)
  if(verbose>=2) msgn("Done")


  
  # Figure out what window size method we are using
  if(verbose>=2) msgn("Figure out window size")
  window.size <- window.size[1]
  if(is.numeric(window.size)){
    if(verbose>=2) msgn("Window method: fixed")
    window.method="fixed"
  }else{
    if(is.character(window.size)){
      if( matches("^[[:blank:]]*[0-9]+[[:blank:]]*(bp?)|(cm)[[:blank:]]*$",window.size) ){
        if(verbose>=2) msgn("Window method: fixed")
        window.method="fixed"
        window.size <- as.numeric(sub("(bp?)|(cm)","",window.size,ignore.case=TRUE))
      }else
      if( matches("^[[:blank:]]*[0-9]+[[:blank:]]*snps?[[:blank:]]*$",window.size) ){
        if(verbose>=2) msgn("Window method: SNPs")
        window.method="snps"
        window.size <- as.integer(sub("snps?","",window.size,ignore.case=TRUE))
      }else
      if( matches("^[[:blank:]]*[0-9]+[.]?[0-9]*[[:blank:]]*sig?m?a?[[:blank:]]*$",window.size) ){
        if(verbose>=2) msgn("Window method: sigma")
        window.method="sigma"
        if( is.na(pop1[1]) | is.na(pop2[1]))
          stop("I need pop1 and pop2 for window.method=\"sigma\"")
        window.size=as.numeric(sub("sig?m?a?","",window.size,ignore.case=TRUE))
      }else
      # WT not supported yet!!!!!!
      if( matches("^[[:blank:]]*[0-9]+[.]?[0-9]*[[:blank:]]*wt[[:blank:]]*$",window.size) ){
        if(verbose>=2) msgn("Window method: wt")
        window.method="wt"
        window.size=as.numeric(sub("[wW][tT]","",window.size,ignore.case=TRUE))
      }else
      stop("Window method", window.size, "is not supported")
    }else
    stop("Window method", window.size, "is not supported")
  }

  if(verbose>=2) msgn("Done")

  # POP1 POP2 business
  if(verbose>=2) msg("Parental populations...")

  if(!is.na(pop1[1]) & !is.na(pop2[1])){
    collect.stats=TRUE
    stats <- make.array(c("mean.l","mean.u","sigma.l","sigma.u"),
                        1:no.bins,data=NA
                       )
    Pop1 <- Pop2 <- character(0)
    for(i in 1:length(pop1))
      Pop1 <- union(Pop1,grep(pop1[i],people,ignore.case=TRUE,value=TRUE))
    for(i in 1:length(pop2))
      Pop2 <- union(Pop2,grep(pop2[i],people,ignore.case=TRUE,value=TRUE))
    if(verbose>=3){
      msg("Population1: ")
      for(p in Pop1) msg(p,", ")
      msgn()
      msg("Population2: ")
      for(p in Pop2) msg(p,", ")
      msgn()
    }
    if(mean(pco$li[Pop1,1])<mean(pco$li[Pop2,1])){
      Pop.l=Pop1
      Pop.u=Pop2
    }
    else{
      Pop.l=Pop2
      Pop.u=Pop1
    }
  }
  msgn("Done")
  
  # define get window procs
  get.window.fixed <- function(bin){
    start   <- (bin-1)*step - window.size/2
    end     <- (bin-1)*step + window.size/2
    window  <- which( (phys.pos >= start) & (phys.pos <= end) )
    return(window)
  }
  
  get.window.snps <- function(bin){
    center  <- (bin-1)*step
    c <- max(which(phys.pos<= center),0)
    s <- max(c-window.size/2+1,1)
    e <- min(c+window.size/2,no.snps)
    window <- s:e
    if(verbose >= 4) msgn("center=",center," c=",c," s=",s," e=",e)
    return(window)
   }

  get.window.sigma <- function(bin){
    center <- (bin-1)*step+1
    c <- max(which(phys.pos<= center),1)
    #
    max.size.reached <- FALSE
    size.bad  <- density.threshold
    size.good <- Inf
    while((abs(size.good-size.bad))>1){
      size <- (size.bad+min(3*size.bad,size.good))/2
      if(verbose>=3) msgn("size.bad=",size.bad,", size.good=",size.good,", size=",size)
      s <- max(c-size/2+1,1)
      e <- min(c+size/2,no.snps)
      if(phys.pos[e]-phys.pos[s]>max.window.size){
        max.size.reached <- TRUE
        if(verbose>=3) msgn("Max window size is reached!")
        start   <- (bin-1)*step - max.window.size/2
        end     <- (bin-1)*step + max.window.size/2
        window  <- which( (phys.pos >= start) & (phys.pos <= end) )
      }
      else{
        window <- s:e
      }
      actual.size <- length(window)
      # check new window
      if(actual.size==0) return(integer(0))
      A <- axis[window,drop=FALSE]
      T <- table[window,people,drop=FALSE]
      quality <- min(100*how.many(!is.na(T[,Pop.l,drop=FALSE]))/(length(Pop.l)*actual.size),
                     100*how.many(!is.na(T[,Pop.u,drop=FALSE]))/(length(Pop.u)*actual.size)
                    )
      if(verbose>=3) msgn("Q=",quality)
      if(quality<quality.threshold){ # window is bad for quality reasons
        if(max.size.reached) return(integer(0))
        size.bad <- size
        next()
      }
      stats.l <- stats(T[,Pop.l,drop=FALSE],A)
      stats.u <- stats(T[,Pop.u,drop=FALSE],A)
      spread <- abs(stats.u$mean-stats.l$mean)/(sqrt(stats.u$sigma)+sqrt(stats.l$sigma))
      window.ok <- spread >= window.size
      if(is.na(window.ok)) window.ok <- FALSE
      if(verbose>=3){
        msgn("m.l=",prettyNum(stats.l$mean,digits=2),
             ", m.u=",prettyNum(stats.u$mean,digits=2),
             ", sig.l=",prettyNum(sqrt(stats.l$sigma),digits=2),
             ", sig.u=",prettyNum(sqrt(stats.u$sigma),digits=2),
             ", spread=",prettyNum(spread,digits=2),
             ", OK=",window.ok)
        
        if(verbose>=4){
          plot.stats(stats.l,
                     stats.u,
                     main=paste("m.l=",prettyNum(stats.l$mean,digits=2),
                       ", m.u=",prettyNum(stats.u$mean,digits=2),
                       ", sig.l=",prettyNum(sqrt(stats.l$sigma),digits=2),
                       ", sig.u=",prettyNum(sqrt(stats.u$sigma),digits=2),
                       ", spread=",prettyNum(spread,digits=2),
                       ", OK=",window.ok),
                     sub=paste("BIN=",bin,
                       ", size.bad=",round(size.bad),
                       ", size.good=",round(size.good),
                       ", size=",round(size))
                     )
          pause()
        }
      }
      if(window.ok) size.good <- size
      else{
        size.bad  <- size
        if(max.size.reached) break()
      }
    }
    return(window)
  }

  
  if(window.method=="fixed"){
    get.window <- get.window.fixed
  } else
  if(window.method=="snps"){
    get.window <- get.window.snps
  } else
  if(window.method=="sigma"){
    get.window <- get.window.sigma
  } else
  if(window.method=="wt"){
    get.window <- get.window.wt
  }
  
  #
  if(verbose>=2){
    msgn("\n================================================================================")
    msgn("Apparent length of the chromosome: ",
         format(chromolength,scientific=FALSE,big.mark=","),
         pp.suffix) 
    msgn("Step:\t\t\t",format(step,scientific=FALSE,big.mark=","),pp.suffix)
    msgn("Window method:\t\t",window.method)
    msgn("Window size:\t\t",format(window.size,scientific=FALSE,big.mark=","))
    msgn("Number of bins:\t\t",format(no.bins,scientific=FALSE,big.mark=","))
    msg("Number of people:\t",format(no.people,scientific=FALSE,big.mark=","))
    if(window.method=="sigma") msgn("/",length(Pop.l),"/",length(Pop.u))
    else msgn(" ")
    msgn("Number of SNPs:\t\t",format(no.snps,scientific=FALSE,big.mark=","))
    msgn("Removed fixed SNPs:\t",remove.fixed.snps)
    msgn("Do NoCall estimates:\t",na.estimates)
    msgn("Collect stats:\t\t",collect.stats)
    msgn("================================================================================")
  }
    
  # Start calculations
  if(verbose>=2){
    msgn("Starting calculation...")
    startingtime=as.integer(Sys.time())
    tpb <- 2
    eta <- 0
  }
  for(bin in 1:no.bins){
    #find window
    window <- get.window(bin)
    density[bin]  <- length(window)
    if(density[bin]>0){
      windows[bin]  <- phys.pos[max(window)]-phys.pos[min(window)]
      windows.snp[bin] <- length(window)
    }
    
    if(verbose>=2) msgc(" bin=",bin,"/",no.bins,";\twin=",density[bin],"snps/",windows[bin],pp.suffix)
    if(verbose>=3) msgn()
    if(verbose>=4) {
      msgn("window beg: ",phys.pos[max(window)], ";  window end: ", phys.pos[min(window)])
    }
    # If there are few SNPs in the window, go to the next bin
    if(density[bin] < density.threshold){
      no.md.bins         <- no.md.bins + 1
      next()
    }
    
    #find projected axis and table
    A <- axis[window,drop=FALSE]
    T <- table[window,,drop=FALSE]
    if(na.estimates){
      T.u <- table.u[window,,drop=FALSE]
      T.l <- table.l[window,,drop=FALSE]
    }

    #Check quality of call within a window
    quality[people,bin] <- 100*(rep(1,density[bin]) %*% !is.na(T))/density[bin]
    quality["Average",bin] <- mean(quality[people,bin],na.rm=TRUE)
    good.people <- people[quality[people,bin] >= quality.threshold]
    if(length(good.people) == 0){
      if(verbose>=2)
        msgc(" bin=",bin,"/",no.bins,";\twin=",density[bin],"snps/",windows[bin],"bp")
      no.mq.bins           <- no.mq.bins + 1
      next()
    }

    # Calculate values for the given bin
    for(gp in good.people){
      coords[gp,bin] <- sum(T[,gp]*A,na.rm=TRUE)/Norm(A[which(!is.na(T[,gp]))],p=1)
      if(na.estimates){
        N <- Norm(A,p=1)
        coords.upper[gp,bin] <- sum(T.u[,gp]*A,na.rm=TRUE)/N
        coords.lower[gp,bin] <- sum(T.l[,gp]*A,na.rm=TRUE)/N
      }
    }

    # Center and normalize
    UPPER=mean(coords[Pop.u,bin],na.rm=TRUE)
    LOWER=mean(coords[Pop.l,bin],na.rm=TRUE)
    coords[,bin]=(2*coords[,bin]-(UPPER+LOWER))/(UPPER-LOWER)
    
    #Fill in the statistics
    if(collect.stats){
      stats["mean.l",bin]         <- mean(coords[Pop.l,bin],na.rm=TRUE)
      stats["mean.u",bin]         <- mean(coords[Pop.u,bin],na.rm=TRUE)
      stats["sigma.l",bin]         <- var(coords[Pop.l,bin],na.rm=TRUE)
      stats["sigma.u",bin]         <- var(coords[Pop.u,bin],na.rm=TRUE)
    }

    # Centering and normalization for estimated statistics. 
    if(na.estimates){
      coords.upper[,bin] <- coords.upper[,bin] - mean(coords.upper[,bin],na.rm=TRUE)
      coords.lower[,bin] <- coords.lower[,bin] - mean(coords.lower[,bin],na.rm=TRUE)
      if(normalize){
        uppers <- which(coords.upper[,bin]+coords.lower[,bin]>0)
        M <- abs(mean(coords.upper[uppers,bin]+coords.lower[uppers,bin]))/2
        coords.upper[uppers,bin] <- coords.upper[uppers,bin]/M
        coords.lower[uppers,bin] <- coords.lower[uppers,bin]/M
        lowers <- which(coords.upper[,bin]+coords.lower[,bin]<0)
        M <- abs(mean(coords.upper[lowers,bin]+coords.lower[lowers,bin]))/2
        coords.upper[lowers,bin] <- coords.upper[lowers,bin]/M
        coords.lower[lowers,bin] <- coords.lower[lowers,bin]/M
      }
    }
    if(verbose>=3) pause()
  }
  if(verbose>=2) msgn("\nDone")
  
  # Put everything together and return
  result <- list()
  class(result) <- "spco"
  result$call        <- sys.call()
  result$date        <- date()
  result$tbl         <- match.call()$table
  result$window.size <- window.size
  result$window.method <- window.method
  result$use.gm <- use.gm
  result$step        <- step
  result$people      <- people
  if(collect.stats){
    result$pop.l       <- Pop.l
    result$pop.u       <- Pop.u
    result$stats       <- stats 
  }
  result$na.estimates<- na.estimates
  result$remove.fixed.snps<- remove.fixed.snps
  result$normalized <- normalize
  result$collect.stats <- collect.stats
  result$chromosome.length <- chromolength
  result$density.threshold <- density.threshold
  result$quality.threshold <- quality.threshold
  result$no.bins     <- no.bins
  result$no.snps     <- no.snps
  result$no.snps.fixed <- no.snps.fixed
  result$no.missed.bins <- make.array(c("Density","Quality"))
  result$no.missed.bins["Density"] <- no.md.bins
  result$no.missed.bins["Quality"] <- no.mq.bins
  result$density     <- density
  result$quality     <- quality
  result$windows     <- windows
  result$windows.snp <- windows.snp
  result$coords      <- coords
  if(na.estimates){
    result$coords.upper <- coords.upper
    result$coords.lower <- coords.lower
  }
  result$info         <- "
  call               -- call to the function that created this object.
  date               -- time and date when spco was created.
  tbl                -- name of the R object containing the table of SNP calls.
  window.method      -- Method to determine the size of the sliding window.
  window.size        -- Size of the sliding window
  step               -- Step
  people             -- Individuals in the analysis.
  pop1,pop2          -- parental populations.
  na.estimates       -- indicates whether upper/lower estimates were made for missing data,
                        those are saved in coords.upper, coords.lower entries.
                        if FALSE, then only coords entry will be present
  remove.fixed.snps  -- indicates whether fixed SNPs were excluded from calculation.
  normalize          -- indicates whether normalization of coordinates was performed.
  snps               -- SNPs in the analysis (not stored)
  chromosome.length  -- Length of the chromosome
  density.threshold  -- Minimum allowed number of SNPs in a window
  quality.threshold  -- Minimum allowed percentage of not NA's for each individual
  no.bins            -- Number of bins
  no.snps            -- number of SNPs
  no.snps.fixed      -- number of SNPs that are fixed.
  no.missed.bins     -- Number of bins skipped for quality or density reasons
  density            -- Density of SNPs for each window
  quality            -- Quality (percentage of valid calls for each window and
                        individual. Also average across individuals for each window.
  windows            -- window size for each bin.
  coords             -- coordinates 
  coords.upper       -- Upper readings for step.pco coordinates
  coords.lower       -- Lower readings for step.pco coordinates
  stats              -- statistics for each bin.
  \n"

  if(verbose>=2){
    msg("The size of the resulting object is ")
    size(result,newline=TRUE)
  }
  return(result)
}

# SPCO
#
################################################################################

################################################################################
#
# STATS


wt <- function(x){
  n <- floor(log2(length(x)))
  if(length(x)!=2^n){
    warning("length of x is not 2^n")
    return(invisible())
  }

  w <- list()
  i = n
  while (length(x)>1){
    l <- length(x)
    even <- 2*(1:(l/2)) #to generate all even indeces
    odd  <- even - 1    #to generate all odd indeces
    w[[i]] <- (x[odd]-x[even])/2 # wavelet coefficients
    x <- (x[odd]+x[even])/2 # new vector x, half the length of previous one
    i=i-1
  }
  w$average <- x # remaining x, (which correponds to average)
  class(w)<-"wavelet"
  return(w)
}


#inverse wavelet transform, recover wave from wavelet coefficients
iwt <- function(w,level=length(w)-1){ #level specifies a filtering cutoff
  x<-w$average
  for (L in 1:level){
    odd <- 2*(1:2^(L-1))-1
    even <- 2*(1:2^(L-1))
    newx <- numeric(2^L)
    newx[odd] <- x+w[[L]]#reconstructing x from wavelet coefficients and old x
    newx[even] <- x-w[[L]]
    x <- newx
  }
  return(x)
}

wt.cutoff <- function(wt,threshold){
  level=length(wt)-1
  for (i in 1:level){
    wt[[i]][abs(wt[[i]])<=threshold] <- 0
  }
  return(wt)
}

cutoff <- function(x,thr=0.05){
  x[which(abs(x)<=thr)] <- 0
  return(x)
}

# as an argument takes a sparse vector of diploid chromosome.
wt.diploid <- function(diploid,levels=17){
  # levels should be at least 2
  levels <- max(levels,2)

  # fill empty wt
  wt <- list()
  class(wt) <- "wavelet"
  power2=1
  for(i in 1:levels){
    wt[[i]] <- numeric(power2)
    power2 <- power2*2
  }
  wt$average <- 0

  l=dim(diploid)[2]

  # if no breakpoints, exit
  if( l==0 ) return(wt)

  # find average
  for( i in 1:l) wt$average <- wt$average + (1-diploid[1,i])*diploid[2,i]

  #the last 1 should be removed
  if(diploid[1,l]==1) l <- l-1

  if( l==0 ) return(wt)
  
  # cycle through bpts and wt levels.
  # each bpt contributes to only one coeff in each level
  for(bpt in 1:l){
    power2 <- 1 # =(1/2)^(level-1)
    for(level in 1:levels){
      #which coeff
      k=floor(diploid[1,bpt]/power2)+1
      wt[[level]][k] <-
        wt[[level]][k] +
        diploid[2,bpt]*( abs(diploid[1,bpt]/power2 - k + 1/2) - 1/2 )
      power2 <- power2/2
    }
  }
  class(wt) <- "wavelet"
  return(wt)
}

wt.summary <- function(wt,p=1,maxlevel=length(wt)-1){
  l <- length(wt)-1
  s <- make.array(c(1:maxlevel,"average","center"))
  for(i in 1:maxlevel){
    s[i] <- ( sum(abs(wt[[i]])^p)^(1/p) ) / (2^(i-1))
  }
  s["average"] <- wt$average
  s["center"] <- sum(s[1:maxlevel]*(1:maxlevel))/sum(s[1:maxlevel])
  class(s) <- "wavelet.summary"
  return(s)
}

plot.wavelet.summary <- function(s,main=""){
  l <- length(s)-1
  barplot(s[1:l],main=main,ylim=c(0,1))
  title(sub=paste("center=",round(s["center"]*10000)/10000,sep=""),line=0)
  lines(c(s["center"],s["center"]),c(0,1),col="red",lwd=3)
}


make.sparse <- function(v){
  
  jumps <- c(v,0) - c(0,v)
  up1    <- which(jumps==2)
  uphalf <- which(jumps==1)
  down1    <- which(jumps== -2)
  downhalf <- which(jumps== -1)
  result <- rbind(c(up1,
                    uphalf,
                    down1,
                    downhalf),
                  c(rep(2,length(up1)),
                    rep(1,length(uphalf)),
                    rep(-2,length(down1)),
                    rep(-1,length(downhalf)))
                 )
  result[1,] <- (result[1,]-1)/length(v)
  result <- result[,order(result[1,])]
  return(result)
}


spco.stats <- function(spco,
                       who=spco$people,
                       lambda=2,
                       breakpointsthreshold=0.5,
                       wtlevels=17,
                       wt.summary.threshold=0.1
                       ){
  result <- list()
  class(result) <- "spco.stats"
  ml <- spco$stats["mean.l",]
  mu <- spco$stats["mean.u",]
  sl <- sqrt(spco$stats["sigma.l",])
  su <- sqrt(spco$stats["sigma.u",])
  uu <- mu-(lambda-breakpointsthreshold)*su
  ul <- mu-(lambda+breakpointsthreshold)*su
  ll <- ml+(lambda-breakpointsthreshold)*sl
  lu <- ml+(lambda+breakpointsthreshold)*sl
  dcoords <- make.array(who,1:spco$no.bins)
  tmp <- make.array(who,1:spco$no.bins,1:8,data=NA)
  level <- function(p,i)
    return((spco$coords[p,i]>=ll[i]) +
           (spco$coords[p,i]>=lu[i]) +
           (spco$coords[p,i]>=ul[i]) +
           (spco$coords[p,i]>=uu[i])
          )

  msgn("Detect admixed regions")
  for(p in who){
    old.state <- 0
    for( i in 1:spco$no.bins ){
       # NA keep oldstate go to next bin
      if(is.na(spco$coords[p,i])){
        dcoords[p,1] <- NA
        next()
      }
      L <- level(p,i)
      
      if(old.state==0){
        if     ( L==4 ) { old.state <- dcoords[p,i] <-  1 }
        else if( L==0 ) { old.state <- dcoords[p,i] <- -1 }
        else            {              dcoords[p,i] <-  0 }
      }

      if(old.state==-1){
        if     ( L==4 ) { old.state <- dcoords[p,i] <-  1 }
        else if( L >1 ) { old.state <- dcoords[p,i] <-  0 }
        else            {              dcoords[p,i] <- -1 }
      }

      if(old.state==1) {
        if     ( L==0 ) { old.state <- dcoords[p,i] <- -1 }
        else if( L< 3 ) { old.state <- dcoords[p,i] <-  0 }
        else            {              dcoords[p,i] <-  1 }
      }
      tmp[p,i,] <- c(old.state, L, ll[i], lu[i], ul[i], uu[i], spco$coords[p,i], dcoords[p,i])
    }
    msgc(p)
  }
  msgn("\nDone")
  result$dcoords <- dcoords
  result$coords.sparse <- list()
  result$breakpoints <- make.array(who)
  result$nb <- result$fg <- result$foreignblockwidth <- make.array(who)
  result$wt <- make.array(who,c(1:10,"average","center"))
  result$wt.raw <- make.array(who,c(1:10,"average","center"))
  result$wt.sp <- make.array(who,c(1:10,"average","center"))
  msgn("Gather stats")
  for(p in who){
    msg(p," ")    
    result$coords.sparse[[p]] <- make.sparse(result$dcoords[p,])
    # breakpoints
    result$breakpoints[p] <- how.many(result$coords.sparse[[p]][1,] != 0 &
                                      result$coords.sparse[[p]][1,] !=1
                                     )/2
    result$fg[p] <- mean(result$dcoords[p,] + 1,na.rm=TRUE)/2
    result$nb[p] <- (
                      how.many(result$coords.sparse[[p]][2,]==1) +
                      2*how.many(result$coords.sparse[[p]][2,]==2)
                    )
    msg(result$nb[p])
    #corrections at the beginning and end
    if(result$coords.sparse[[p]][1,1]==0 & result$coords.sparse[[p]][2,1]>0){
      result$nb[p] <- result$nb[p] + 1
      msg("+1")
    }
    l <- length(result$coords.sparse[[p]][1,])
    if(result$coords.sparse[[p]][1,l]==1 & result$coords.sparse[[p]][2,l]>0){
      result$nb[p] <- result$nb[p] - 1
      msg("-1")
    }
    msg("=",result$nb[p])
    result$foreignblockwidth[p] <- result$fg[p]/result$nb[p]
    wtc <- wt.diploid(result$coords.sparse[[p]],levels=10)
    result$wt.sp[p,] <- wt.summary(wtc)
    wtc <- wt(spco$coords[p,])
    result$wt.raw[p,] <- wt.summary(wtc)
    wtc <- wt(result$dcoords[p,])
    result$wt[p,] <- wt.summary(wtc)
    msgn()
  }
  msgn("\nDone")
  result$tmp <- tmp
  return(result)
}

spco.wtstats <- function(spco,
                         who=spco$people,
                         threshold=0,
                         p=1,
                         maxlevel=log2(length(spco$coords[1,]))
                        ){
  result <- list()
  class(result) <- "spco.wtstats"
  result$spco <- match.call()$spco
  result$threshold <- threshold
  result$p <- p
  result$wt <- list()
  result$wt.summary <- make.array(who,c(1:maxlevel,"average","center")) 
  for(ind in who){
    result$wt[[ind]] <- wt(spco$coords[ind,])
 
    result$wt.summary[ind,] <- wt.summary(wt.cutoff(result$wt[[ind]],threshold=threshold),maxlevel=maxlevel)
  }
  return(result)
}



# STATS
#
################################################################################

################################################################################
# PLOT SPCO
#
# AUX functions                          

slog <- function(x,coeff=1){
  return(sign(x)*log(exp(coeff)*abs(x)+1)/coeff)
}

#########################
print.spco <- function(spco){
  msgn("SPCO:")
  msgn("Created:\t\t",spco$date)
  msgn("table:\t\t\t\"",spco$tbl,"\"")
  if(spco$window.method=="fixed") sfx=" bp"
  if(spco$window.method=="snps")  sfx=" snps"
  if(spco$window.method=="sigma")  sfx=" sigma"
  if(spco$window.method=="wavelets")  sfx=" WT"
  msgn("Window method:\t\t",spco$window.method)
  msgn("Window:\t\t\t",format(spco$window.size,scientific=FALSE,big.mark=","),sfx)
  msgn("Ave. win. size:\t\t",
       format(round(mean(spco$windows,na.rm=TRUE)),scientific=FALSE,big.mark=","),
       " bp")
  msgn("Step:\t\t\t",format(spco$step,scientific=FALSE,big.mark=",")," bp")
  msgn("No. people:\t\t",format(length(spco$people),scientific=FALSE,big.mark=","))
  msgn("No. snps:\t\t",
       format(spco$no.snps,scientific=FALSE,big.mark=","),
       " + ",
       format(spco$no.snps.fixed,scientific=FALSE,big.mark=","),
       "(fixed)"
       )
  msgn("Chromo. length:\t\t",
       format(spco$chromosome.length,scientific=FALSE,big.mark=","),
       " bp")
  msgn("Density threshold:\t",
       format(spco$density.threshold,scientific=FALSE,big.mark=","),
       " bp/window")
  msgn("Ave. density:\t\t",
       format(floor(mean(spco$density[!is.na(spco$quality["Average",])])),
              scientific=FALSE,big.mark=","
              ),
       " snps/window"
       )
  msgn("Quality thereshold:\t",round(spco$quality.threshold*100),"%")
  msgn("Ave. quality:\t\t",
       format(floor(mean(spco$quality["Average",],na.rm=TRUE)),
              scientific=FALSE,
              big.mark=","
              ),
       " %")
  msgn("No. bins:\t\t", format(spco$no.bins,scientific=FALSE,big.mark=","))
  msgn("No. of missed bins:\t",
       format(sum(spco$no.missed.bins),scientific=FALSE,big.mark=","),
       " = ",
       format(spco$no.missed.bins["Density"],scientific=FALSE,big.mark=","),
       "(d) + ",
       format(spco$no.missed.bins["Quality"],scientific=FALSE,big.mark=","),
       "(q)"
       )
  msgn("NA estimates:\t\t",spco$na.estimates)
  msgn("Removed fixed snps:\t",spco$remove.fixed.snps)
 # find plotting ranges
  range <- 1:spco$no.bins
  ranges <- list()
  i=0
  while(length(range)>0){
    # remove NAs at the beginning
    b <- min(which(!is.na(spco$quality["Average",range])))
    if(b==Inf) break
    range <- range[b:length(range)]
    i=i+1
    # find the end of the plotting range
    e <- min(which(is.na(c(spco$quality["Average",range],NA))))-1
    ranges[[i]] <- range[1:e]
    if( e==length(range) ) break
    range<-range[(e+1):length(range)]
  }
  msgn("No. contig. ranges:\t",format(i,,scientific=FALSE,big.mark=","))
  reg  <- ""
  for(i in 1:length(ranges)){
    reg <-paste(reg,"",min(ranges[[i]]),"--",max(ranges[[i]]),";  ",sep="")
  }
  msgn("Ranges:\t\t\t",reg)
}

plot.spco.info <- function(spco){
  gap1<-0.045
  gap2<-0.052
  plot.new()
  y <- 1
  EE <- environment()
  plottext <- function(...){
    text(0,y,paste(...,sep=""),pos=4)
    assign("y",y-gap2,envir=EE)
  }
  
  title(main=paste("SPCO:",spco$date,spco$tbl),col="green")

  call <- format(spco$call)
  y    <- 1
  text(0,y,"CALL:  ",pos=4)
  y <- y+gap1
  for(i in 1:length(call)){
    y<-y-gap1
    text(0.1,y,call[i],pos=4,col="blue")
  }
  y<-y-gap2
  
  if(spco$window.method=="fixed") sfx=" bp"
  if(spco$window.method=="snps")  sfx=" snps"
  if(spco$window.method=="sigma")  sfx=" sigma"
  plottext("Window method:\t\t",spco$window.method)
  plottext("Window:\t\t\t",format(spco$window.size,scientific=FALSE,big.mark=","),sfx)
  plottext("Ave. win. size:\t\t",
       format(round(mean(spco$windows,na.rm=TRUE)),scientific=FALSE,big.mark=","),
       " bp")
  plottext("Step:\t\t\t",format(spco$step,scientific=FALSE,big.mark=",")," bp")
  plottext("No. people:\t\t",format(length(spco$people),scientific=FALSE,big.mark=","))
  plottext("No. snps:\t\t",
       format(spco$no.snps,scientific=FALSE,big.mark=","),
       " + ",
       format(spco$no.snps.fixed,scientific=FALSE,big.mark=","),
       "(fixed)"
       )
  plottext("Chromo. length:\t\t",
       format(spco$chromosome.length,scientific=FALSE,big.mark=","),
       " bp")
  plottext("Density threshold:\t",
       format(spco$density.threshold,scientific=FALSE,big.mark=","),
       " bp/window")
  plottext("Ave. density:\t\t",
       format(floor(mean(spco$density[!is.na(spco$quality["Average",])])),
              scientific=FALSE,big.mark=","
              ),
       " snps/winfow"
       )
  plottext("Quality thereshold:\t",round(spco$quality.threshold*100),"%")
  plottext("Ave. quality:\t\t",
       format(floor(mean(spco$quality["Average",],na.rm=TRUE)),
              scientific=FALSE,
              big.mark=","
              ),
       " %")
  plottext("No. bins:\t\t", format(spco$no.bins,scientific=FALSE,big.mark=","))
  plottext("No. of missed bins:\t",
       format(sum(spco$no.missed.bins),scientific=FALSE,big.mark=","),
       " = ",
       format(spco$no.missed.bins["Density"],scientific=FALSE,big.mark=","),
       "(d) + ",
       format(spco$no.missed.bins["Quality"],scientific=FALSE,big.mark=","),
       "(q)"
       )
  plottext("NA estimates:\t\t",spco$na.estimates)
  plottext("Removed fixed SNPs:\t",spco$remove.fixed.snps)
 # find plotting ranges
  range <- 1:spco$no.bins
  ranges <- list()
  i=0
  while(length(range)>0){
    # remove NAs at the beginning
    b <- min(which(!is.na(spco$quality["Average",range])))
    if(b==Inf) break
    range <- range[b:length(range)]
    i=i+1
    # find the end of the plotting range
    e <- min(which(is.na(c(spco$quality["Average",range],NA))))-1
    ranges[[i]] <- range[1:e]
    if( e==length(range) ) break
    range<-range[(e+1):length(range)]
  }
  plottext("No. contig. ranges:\t",format(i,,scientific=FALSE,big.mark=","))
  reg  <- ""
  for(i in 1:length(ranges)){
    reg <-paste(reg,"",min(ranges[[i]]),"--",max(ranges[[i]]),";  ",sep="")
  }
  plottext("Ranges:\t\t\t",reg)
}


plot.spco.coords <-function(spco,                     # spco as returned by step.pco
                            who=spco$people,          # vector of individuals
                            range=1:spco$no.bins,     # range within chromosome
                            slog=FALSE,               # Log scale?
                            v.range.total=FALSE,      # how to choose vert. scale
                            type="c",
                               # what to plot, c=center,a=average,u=upper,l=lower
                            colors="colorfile.txt",   # file to get colors from
                            alpha=0.3,
                            b.alpha=alpha,
                            new=TRUE                  # Initiate new plot?
                            ){
  # Initiate new plot
  if(new)       plot.new()

  # WHO???
  people <- spco$people
  names(people) <- people
  people <- people[who]

  # Get the colors
  if(length(colors)==1){
    as.vector(read.table(colors)[people,"color"])->my.colors}
  else{
    colors[people,"color"]->my.colors}
  names(my.colors) <- people
  
  # Shall we slog?
  slog.tmp <- function(x,...) return(x)
  if(class(slog)=="function"){
    slog.tmp=slog
  } else
  if(slog)      slog.tmp <- function(x) return(slog(x,coeff=slog))

  # What to plot?
  if(!spco$na.estimates) type="c"
  if(type=="ul" | type=="lu"){
    type="ul"
    cu <-  slog.tmp(spco$coords.upper[people,])
    cl <-  slog.tmp(spco$coords.lower[people,])
  }
  if(type=="c") co <-  slog.tmp(spco$coords[people,,drop=FALSE])
  if(type=="a") co <- slog.tmp((spco$coords.upper[people,,drop=FALSE] +
                                spco$coords.lower[people,,drop=FALSE])/2)
  if(type=="u") co <-  slog.tmp(spco$coords.upper[people,,drop=FALSE])
  if(type=="l") co <-  slog.tmp(spco$coords.lower[people,,drop=FALSE])

  # Find window dimensions
  if(v.range.total) range1 <- 1:spco$no.bins
  else range1 <- range
  w=min(range)
  W=max(range)
  if(type=="ul"){
    M=max(cu[,range1],na.rm=TRUE)
    m=min(cl[,range1],na.rm=TRUE)
  }
  else{
    M=max(co[,range1],na.rm=TRUE)
    m=min(co[,range1],na.rm=TRUE)
  }

  # Plot
  plot.window(c(w,W),c(m,M))
  title(ylab="PC1",line=2)
  axis(2,tick=TRUE,labels=TRUE)
  
  if(type=="ul"){
    for(p in people){
      # Find plotting ranges
      range1 <- range
      ranges <- list()
      i=0
      while(length(range1)>0){
          # remove NAs at the beginning
        b <- min(which(!is.na(cu[p,range1]) & !is.na(cl[p,range1])))
        if(b==Inf) break
        range1 <- range1[b:length(range1)]
        i=i+1
          # find the end of the plotting range
        e <- min(which(is.na(c(cu[1,range1],NA))))-1
        ranges[[i]] <- range1[1:e]
        if( e==length(range1) ) break
        range1<-range1[(e+1):length(range1)]  
      }

      # Plot 
      if(length(ranges) == 0) next
      col=t(col2rgb(my.colors[p]))
      b.col=rgb(col,alpha=b.alpha*255,maxColorValue=255)
      f.col=rgb(col,alpha=alpha*255,maxColorValue=255)
      for( i in 1:length(ranges)){
        r=ranges[[i]]
        rr=r[length(r):1]
        polygon(c(r,rr,r[1]),c(cu[p,r],cl[p,rr],cu[p,r[1]]),
                col=f.col,
                border=b.col
                )
      }
    }
  }
  else{
    for(p in people){
      lines(range,co[p,range],col=my.colors[p])
    }
  }
}

plot.spco.stats <- function(spco,                     # spco as returnd by step.pco
                            who=integer(0),           # vector of individuals
                            range=1:spco$no.bins,     # range within chromosome
                            slog=FALSE,               # Log scale?
                            v.range.total=FALSE,      # how to choose vert. scale
                            colors="colorfile.txt",   # file to get colors from
                            alpha=0.2,
                            new=TRUE,                 # Initiate new plot?
                            main=paste("SPCO:",spco$tbl,spco$date),
                            sub=paste("w=",
                                format(spco$window.size,scientific=FALSE,big.mark=","),
                                if(spco$window.method=="fixed") "bp, " else
                                if(spco$window.method=="sigma") "sigma " else "snps, ",
                                "  s=",
                                format(spco$step,scientific=FALSE,big.mark=","),
                                "bp",
                                sep="")
                            ){
  # Initiate new plot
  if(new)       plot.new()

  # WHO???
  people <- spco$people
  names(people) <- people
  people <- people[who]

  # Get the colors
  if(length(colors)==1){
    as.matrix(read.table(colors))[,"color"]->my.colors
  }
  else{
    colors[,"color"]->my.colors
  }
  
  m.l <- spco$stats["mean.l",range]
  m.u <- spco$stats["mean.u",range]
  s.l <- sqrt(spco$stats["sigma.l",range])
  s.u <- sqrt(spco$stats["sigma.u",range])
  
  # Find window dimensions, setup window
  w=min(range)
  W=max(range)
  M=max(m.u+(3*s.u),spco$coords[people,],na.rm=TRUE)
  m=min(m.l-(3*s.l),spco$coords[people,],na.rm=TRUE)
  plot.window(c(w,W),c(m,M))
  title(ylab="PC1",line=2)
  title(main=main,sub=sub)
  axis(2,tick=TRUE,labels=TRUE)
  draw.axis=TRUE
  if(draw.axis){
    no.lab=7
    lstep=length(range)/no.lab
    lbl=range[1+(0:no.lab)*lstep]
    lbl=(lbl-1)*spco$step
    axis(side=1,
         at=range[1+(0:no.lab)*lstep],
         labels=lbl)
  }

  # Find plotting ranges
  range1 <- range
  ranges <- list()
  i=0
  while(length(range1)>0){
    # remove NAs at the beginning
    b <- min(which(!is.na(m.l[range1])))
    if(b==Inf) break()
    range1 <- range1[b:length(range1)]
    i=i+1
    # find the end of the plotting range
    e <- min(which(is.na(c(m.l[range1],NA))))-1
    ranges[[i]] <- range1[1:e]
    if( e==length(range1) ) break()
    range1<-range1[(e+1):length(range1)]  
  }

  # Plot 
  if(length(ranges) >= 0){
    solid.l=my.colors[spco$pop.l[1]]#"forestgreen"
    solid.u=my.colors[spco$pop.u[1]]#"blue"
    #print(my.colors)
    #msgn("solid.l=",solid.l)
    transp.l=add.alpha(solid.l,alpha)
    transp.u=add.alpha(solid.u,alpha)
    for( i in 1:length(ranges)){
      r=ranges[[i]]
      rr=r[length(r):1]
      xs <- c(r,rr,r[1])
      for(k in 1:3){
        ys <- c((m.l + k*s.l)[r],
                (m.l - k*s.l)[rr],
                (m.l + k*s.l)[r[1]]
               )
        polygon(x=xs,y=ys, col=transp.l,border=NA)
        ys <- c((m.u + k*s.u)[r],
                (m.u - k*s.u)[rr],
                (m.u + k*s.u)[r[1]]
               )
        polygon(x=xs,y=ys, col=transp.u,border=NA)
      }
    }
  }
  lines(range,m.l,col=solid.l,lwd=2)
  lines(range,m.u,col=solid.u,lwd=2)
  for(p in people){
      lines(range,spco$coords[p,range],col=my.colors[p])
    }
}



plot.spco.quality <- function(spco,                      # spco as returnd by step.pco
                              who=spco$people,           # vector of individuals
                              range=1:spco$no.bins,      # range within chromosome
                              v.range.total=FALSE,       # how to choose vert. scale
                              colors="colorfile.txt",     # file to get colors from
                              new=TRUE                   # Initiate new plot?
                              ){
  # Initiate new plot
  if(new)       plot.new()

  # WHO???
  people <- spco$people
  names(people) <- people
  people <- people[who]

  # Get the colors
  if(length(colors)==1){
    as.vector(read.table(colors)[c(people,"Average"),"color"])->my.colors}
  else{
    colors[c(people,"Average"),"color"]->my.colors}
  names(my.colors) <- c(people,"Average")
  
  q <- spco$quality[c(people,"Average"),]

  if(v.range.total) range1 <- 1:spco$no.bins
  else range1 <- range
  M=max(q[,range1],na.rm=TRUE)
  m=min(q[,range1],na.rm=TRUE)
  w=min(range)
  W=max(range)
  plot.window(c(w,W),c(m,M))
  title(ylab="Quality\n%",
        line=2
        )
  axis(2,col=my.colors["Average"])
  for(p in people){
    lines(range,q[p,range],col=my.colors[p])
  }
  lines(range,q["Average",range],col=my.colors["Average"],lwd=3)
}


plot.spco.density <- function(spco,                  # spco as returned by step.pco
                              range=1:spco$no.bins,  # range within chromosome
                              v.range.total=TRUE,    # how to choose vert. scale
                              col="green",           # color for the density graph
                              t.col="red",           # color for the threshold line
                              a.col="yellow",        # average
                              palette=                # palette
                              rainbow(32,start=4/6,end=0),
                              type="g",
                              new=TRUE               # Initiate new plot?
                              ){
  # Initiate new plot
  if(new)       plot.new()

  den <- spco$density

  if(v.range.total) range1 <- 1:spco$no.bins
  else range1 <- range
  M=max(den[range1],na.rm=TRUE)
  m=min(den[range1],na.rm=TRUE)
  w=min(range)
  W=max(range)

  if(type=="g"){
    plot.window(c(w,W),c(m,M))
    title(ylab="Density\nsnp/window",
          line=2
          )
    axis(2,col=col)
    lines(range,den[range],col=col)
    notmissed<-intersect(which(!is.na(spco$quality["Average",])),range)
    ad<-make.array(1:spco$no.bins,data=mean(den[notmissed]))
    ad[is.na(spco$quality["Average",])]<-NA
    lines(range,ad[range],col=a.col,lwd=2)
    lines(c(w,W),c(spco$density.threshold,spco$density.threshold),col=t.col,lwd=2)
  }
  else{
    colorf<-function(value){
      value<-min(M,max(value,m))
      return(palette[floor((value-m)/(M-m+1) * length(palette))+1])
    }
    plot.window(c(w,W),c(0,1))
    for(i in range){
      rect(i-1,0,i,1, col=colorf(den[i]),border=NA)
    }
  }
}

plot.spco.windowsize <- function(spco,                  # spco as returned by step.pco
                                range=1:spco$no.bins,  # range within chromosome
                                v.range.total=TRUE,    # how to choose vert. scale
                                col="blue",           # color for the graph
                                a.col="yellow",        # average
                                f.col="violet",
                                new=TRUE               # Initiate new plot?
                              ){
  # Initiate new plot
  if(new)       plot.new()

  win <- spco$windows

  if(v.range.total) range1 <- 1:spco$no.bins
  else range1 <- range
  M=max(win[range1],na.rm=TRUE)
  m=min(win[range1],na.rm=TRUE)
  w=min(range)
  W=max(range)

  plot.window(c(w,W),c(m,M))
  title(ylab="Window\nsize",
        line=2
        )
  axis(2,col=col)
  lines(range,win[range],col=col)
  notmissed <- intersect(which(!is.na(spco$quality["Average",])),range)
  aw <- make.array(1:spco$no.bins,data=mean(win[notmissed]))
  aw[is.na(spco$quality["Average",])] <- NA 
  lines(range,aw[range],col=a.col,lwd=2)
  if(spco$window.method=="fixed"){
    lines(c(w,W),c(spco$window.size,spco$window.size),col=f.col,lwd=2)
  }
}

plot.spco <- function(spco,                     # spco as returnd by step.pco
                      who=spco$people,          # vector of individuals
                      what=NA,
                      range=1:spco$no.bins,     # range within chromosome
                      slog=FALSE,               # Log scale?
                      colors="colorfile.txt",   # file to get colors from
                      proportions=c(7,1,1,1),   # vertical proportion of the parts
                      coords="coords" %in% what,# plot coords?
                      quality="quality" %in% what,             # plot quality graph?
                      density="density" %in% what,             # plot density graph?
                      windows="windows" %in% what,
                      density.type="g",         # type of density graph ("graph"/"rainbow")
                      palette=rainbow(32,start=4/6,end=0),
                      coords.type="ul",          # type of coords graph, cloud=ul/average/upper/lower
                      draw.axis=TRUE,
                      labels="chromo",
                      main=paste("SPCO:",spco$tbl,spco$date),              # main title
                      sub=paste("w=",
                                format(spco$window.size,scientific=FALSE,big.mark=","),
                                if(spco$window.method=="fixed") {
                                  if(!spco$use.gm) "bp, "
                                  else "cM, "
                                  }
                                else 
                                if(spco$window.method=="sigma") "sigma " else "snps, ",
                                "  s=",
                                format(spco$step,scientific=FALSE,big.mark=","),
                                if(!spco$use.gm) "bp"
                                else "cM" ,
                                sep="")                 # subtitle
                      ){
  leftmar<-0.1
  rightmar<-0.95
  botmar <-0.2
  topmar <-0.1
  if(!coords)  proportions[1] <- 0
  if(!density) proportions[2] <- 0
  if(!quality) proportions[3] <- 0
  if(!windows) proportions[4] <- 0
  qt <- (1-botmar-topmar)/sum(proportions)
  qc <- qt * proportions[1]
  qd <- qt * proportions[2]
  qq <- qt * proportions[3]
  qw <- qt * proportions[4]

  plot.new()
  plot.window(c(0,1),c(0,1))
  title(main=main,sub=sub)
  
  if(draw.axis){
    no.lab=7
    par(plt=c(leftmar,rightmar,botmar-0.0001,botmar))
    w=min(range)
    W=max(range)
    plot.window(c(w,W),c(0,0.1))
    lstep=length(range)/no.lab
    lbl=range[1+(0:no.lab)*lstep]
    if(labels=="chromo")
      lbl=prettyNum((lbl-1)*spco$step,digits=4)
    axis(side=1,
         at=range[1+(0:no.lab)*lstep],
         labels=lbl)
  }

  if(windows){
    par(plt=c(leftmar,rightmar,botmar,botmar+qw))
    box()
    plot.spco.windowsize(spco,
                         range=range,
                         v.range.total=FALSE,
                         new=FALSE)
  }

  if(quality){
    par(plt=c(leftmar,rightmar,botmar+qw,botmar+qw+qq))
    box()
    plot.spco.quality(spco,
                      who=who,
                      range=range,
                      v.range.total=FALSE,
                      colors=colors,
                      new=FALSE)
  }

  if(density){
    par(plt=c(leftmar,rightmar,botmar+qw+qq,botmar+qw+qq+qd))
    box()
    plot.spco.density(spco,
                      range=range,
                      v.range.total=FALSE,
                      col="green",
                      t.col="red",
                      a.col="yellow",
                      type=density.type,
                      palette=palette,
                      new=FALSE)
  }

  if(coords){
    par(plt=c(leftmar,rightmar,botmar+qw+qd+qq,botmar+qw+qd+qq+qc))
    box()
    plot.spco.coords(spco,
                     who=who,
                     range=range,
                     v.range.total=FALSE,
                     slog=slog,
                     new=FALSE,
                     colors=colors,
                     type=coords.type)
  }
}

#
#
################################################################################
