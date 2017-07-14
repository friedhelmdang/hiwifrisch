
## Functions


"sm.start.timer" <- function(depth=1)
  {
    st.start.timer(depth)
  }

"sm.stop.timer" <- function(depth=1,info.level=0)
{
  st.stop.timer(depth,info.level)
}

mab.df <- function(d,f){paste (d, "/",f,sep="")}

mab.simulate <- function (
   simulation.name  ,
   simulation.run   ,
   repetitions      = 1000,      
   recurrent.parent = "H",
   donor.parent     = "H",
   linkage.map      ,
   target.loci      = NULL,
   flanking.loci    = NULL,
   recipient.loci   = NULL,
   gen.type         = NULL,
   population.size  = NULL,
   sel.strategy     = NULL,
   no.selected      = NULL,
   no.preselected   = NULL,
   reg.chr          = NULL,                            
   reg.begin        = NULL,                            
   reg.end          = NULL,                            
   missing.allele   = "9",
   result.file      = NULL,
   success.factor   = 1   ,
   recode.infiles   = FALSE                       
){

  if ( is.null(no.selected) ) no.selected <- rep(1,length(population.size)) 
  l1 <- length(gen.type)
  l2 <- length(population.size)
  l3 <- length(sel.strategy)
  l4 <- length(no.selected)
  l5 <- length(no.preselected)
  if ( ( (l1 != l2) || (l2 != l3) ) || (l3 != l4) ) {
    info(-2,"Definition of generations incorrect")
    return (invisible(NULL))
  }
  if ( (l5 > 0)  && (l5 != l4) ) {
    info(-2,"Definition of generations incorrect")
    return (invisible(NULL))
  }
  l1 <- length(reg.chr  )
  l2 <- length(reg.begin)
  l3 <- length(reg.end  )
  if (  (l1 != l2) || (l2 != l3)  ) {
    info(-2,"Definition of target regions  incorrect")
    return (invisible(NULL))
  }
  # Alias names for selection strategies
  sel.strategy <- as.character(sel.strategy)
  if ( length (sel.strategy) > 0 )
    for (i in 1:length(sel.strategy))
    {
      if ( sel.strategy[i]=="0")  sel.strategy[i] <- "n" else
      if ( sel.strategy[i]=="1")  sel.strategy[i] <- "t" else
      if ( sel.strategy[i]=="2")  sel.strategy[i] <- "tb" else
      if ( sel.strategy[i]=="3")  sel.strategy[i] <- "tfb" 
    }
#
  sm.start.timer()
# Constants
  t.name <- "t"
  f.name <- "fm"
  r.name <- "r"
  m.name <- "m"
# Load map file and prepare for plabsim
  map <- read.table(
                   file=mab.df(st.input.dir,linkage.map),
                   stringsAsFactors=F
                   )
  names(map) <- c("chrom","pos","name")
  map$pos <- map$pos/100
  map$class <- rep(m.name ,nrow(map))
# Sort map for chromosomes
  info(1,"Processing mapfile")
  repeat{
    change <- FALSE
    for (i in 2:nrow(map)){
      if (map[i-1,1] > map[i,1]) {
        tm <- map[i,];
        map[i,] <- map[i-1,];
        map[i-1,] <- tm
        change <- TRUE
      }
    }
    if (!change) break; 
  }
# Sort map for positions
  repeat{
    change <- FALSE
    for (i in 2:nrow(map)){
      if ( (map[i-1,2] > map[i,2]) && (map[i-1,1] == map[i,1])) {
        tm <- map[i,];
        map[i,] <- map[i-1,];
        map[i-1,] <- tm
        change <- TRUE
      }
    }
    if (!change) break; 
  }
 # Reduce map and marker data to common markers
  if (  (recurrent.parent != "H") && (donor.parent != "H") )
    { 
      p1.raw <- read.table(
                           file=mab.df(st.input.dir,donor.parent),
                           stringsAsFactors=F,
                           skip=1
                           )
      p1.name <- scan(
                      file=mab.df(st.input.dir,donor.parent),
                      what="character",
                      n=1,
                      quiet=T
                      )
      if (2==length(p1.raw)) p1.raw <- cbind(p1.raw ,p1.raw[,2] )
      names(p1.raw) <- c("name","allele1","allele2")
      p2.raw <- read.table(
                           file=mab.df(st.input.dir,recurrent.parent),
                           stringsAsFactors=F,
                           skip=1
                           )
      p2.name <- scan(
                      file=mab.df(st.input.dir,recurrent.parent),
                      what="character",
                      n=1,
                      quiet=T
                      )
      if (2==length(p2.raw)) p2.raw <- cbind(p2.raw ,p2.raw[,2] )
      names(p2.raw) <- c("name","allele3","allele4")
      a1 <- merge(map,p1.raw,sort=F)
      a2 <- merge(map,p2.raw,sort=F)
      index <- ( (!(a1$allele1==missing.allele)) &&
                        (!(a1$allele2==missing.allele)) )
      a1 <- a1[index,]
      index <- ( (!(a2$allele3==missing.allele)) &&
                        (!(a2$allele4==missing.allele)) )
      a2 <- a2[index,]
      c.m <- data.frame(name=merge(a1,a2,sort=F)$name)
      map    <- merge(map,c.m,sort=F)
      map <- data.frame(chrom=map$chrom,
                        pos=map$pos    ,
                        name=as.character(map$name)  ,
                        class=as.character(map$class),
                        stringsAsFactors=F)
      p1.raw <- merge(p1.raw,c.m,sort=F)        
      p2.raw <- merge(p2.raw,c.m,sort=F)        
    }  
  # Add the target loci and flanking markers to the map
  if (length(target.loci) !=0 )
    {
      found <- FALSE
      for (i in 1:nrow(map))
        {
          for (j in 1:length(target.loci))
            {
              if( map[i,3] == target.loci[j] )
                {
                  map[i,4] <- t.name
                  found <- TRUE
                }
            }
        }
      if (found==FALSE)
        {
          info(-2,"Map position or marker data for target locus missing")
          return (invisible(NULL))
        }
    }
  # Add the flanking markers to the map
  if (length(flanking.loci) !=0 )
    {
      found <- FALSE
      for (i in 1:nrow(map))
        {
          for (k in 1:length(flanking.loci))
            {
              if( map[i,3] == flanking.loci[k] )
                {
                  map[i,4] <- f.name
                  found <- TRUE
                }
            }
        }
      if (found==FALSE)
        {
          info(-2,"Map position or marker data for flanking locus missing")
          return (invisible(NULL))
        }
    }
    # Add the recipient loci to the map
   if (length(recipient.loci) !=0 )
    {
      found <- FALSE
      for (i in 1:nrow(map))
        {
          for (k in 1:length(recipient.loci))
            {
              if( map[i,3] == recipient.loci[k] )
                {
                  map[i,4] <- r.name
                  found <- TRUE
                }
            }
        }
      if (found==FALSE)
        {
          info(-2,"Map position or marker data for recipient locus missing")
          return (invisible(NULL))
        }
    }
  # 
   maptmp <- mab.df(st.data.dir,
                  as.character(sprintf("%05i",round(runif(1,1,10000)))))
   write.table (map,file=maptmp, quote=F,col.names = F,row.names = F)
   no.chroms <- max(map$chrom)
   chrom.length <- 0.001+tapply(map$pos,map$chrom,max)
   # initialize simulations
   info(0,paste(simulation.name,"-",simulation.run,sep=""))
   genome.parameter.set(           
     no.chrom  =  no.chroms,          
     no.hom    = 2,          
     chrom.len = chrom.length
    ) 
   load.linkage.map(maptmp)
   unlink(maptmp)
   define.effects(t.name, paste("0,",t.name," 1 uniform 1"))
   define.effects(f.name, paste("0,",f.name," 8 uniform 1"))
   define.effects(r.name, paste("0,",r.name," 8 uniform 1"))
   define.effects(m.name, paste("0,",m.name," 8 uniform 1"))
 # Load parents
   init.population("P1",homozygote(1)); do.name <- "H1"
   init.population("P2",homozygote(8)); rp.name <- "H8"
   if ( (recurrent.parent != "H") && (donor.parent == "H") ) {
     info(-2,"Marker data for donor parent missing")
     return (invisible(NULL))
   }
   if ( (recurrent.parent == "H") && (donor.parent != "H") ) {
    info(-2,"Marker data for recuttent parent missing")
    return (invisible(NULL))
   }
#
 if (  (recurrent.parent != "H") && (donor.parent != "H") )
  {
        x <- return.population ("P1 P2")
        remove.population("P1 P2")
        x$name <- splitdt(row.names(x))[,1]
        x$row.label <- row.names(x)
        y <- merge (x,p1.raw)
        y <- merge (y,p2.raw)
        for (j in (-1+2*(1:(nrow(y)/2))) ){
          rp.alleles <- as.numeric ((y$allele1[j] == y$allele3[j]) ||
                                    (y$allele1[j] == y$allele4[j]) ) +
                       as.numeric ( (y$allele2[j] == y$allele3[j]) ||
                                    (y$allele2[j] == y$allele4[j]) ) 
          if ( rp.alleles == 0)  y[j:(j+1),2] <- c (1,0) else
          if ( rp.alleles == 1)  y[j:(j+1),2] <- c (1,1) else
          if ( rp.alleles == 2)  y[j:(j+1),2] <- c (0,1)
        }
#
        if (length(target.loci) !=0 )
          {
            for (i in (-1+2*(1:(nrow(y)/2))) ) 
              {
                for (j in 1:length(target.loci))
                  {
                    if( y[i,1] == target.loci[j] )
                      {
                         if ( ( y[i,2]   == y[i,3] )  &&
                              ( y[i+1,2] == y[i+1,3] )   ) 
                           {
                             info(-2,
                                  paste("Monomorphic target locus '",
                                        target.loci[j], "'", sep="")
                                  )
                             return(invisible(NULL))
                           }
                          if ( ( ( y[i,2]  == y[i+1,3] ) &&
                               (  y[i+1,2] == y[i  ,3] )    ) &&
                               (y$allele1[i] != y$allele2[i])    )
                           {
                             y[i:(i+1),2] <- c (1,1)
                             info(0,
                                  paste("Heterozygous target locus '",
                                        target.loci[j], "'", sep="")
                                  )
                           }
                      }
                  }
              }
          }
        z <- data.frame(P1.1=y[,2],P2.1=y[,3])
        rownames(z) <- y$row.label
        # replace the x matrix !
        if (recode.infiles)
          {
            zz <- z
            colnames (zz) <- c(p1.name,p2.name)
            outfname <- mab.df(st.output.dir,
                               paste ("recoded", donor.parent,
                                      recurrent.parent, sep="-"))
            write.table ( zz,
                         file = outfname,
                         row.names = TRUE,
                         col.names = TRUE,
                         quote     = FALSE
                         )
          }
        generate.population(
                            z,
                            backcross=TRUE
                            )
  }
# Simulation
   p.s <- population.size
   g.t <- gen.type
   s.s <- sel.strategy
   n.s <- no.selected 
   n.p <- no.preselected 
   n.gen <- length (g.t)
   if (n.gen > 0 )
     {
       # population names
       pop <- paste("GEN",1:n.gen,sep="")
       s1 <- paste(pop,"s1",sep=""); s2 <- paste(pop,"s2",sep="")
       s3 <- paste(pop,"s3",sep="")
       sl <- paste(pop,"se",sep=""); st <- paste(pop,"st",sep="")
       # breeding scheme
       remove.population(st)
       ot <- round(repetitions/10)
       ht <- matrix(0,nrow=n.gen,ncol=repetitions)
       sm <- matrix(0,nrow=n.gen,ncol=repetitions)
       for (j in 1:repetitions){
         reset.mdp()
         if ( (j/ot) == (round(j/ot)) ) info.cat(0,paste(j,"of",repetitions))
         for (i in 1:n.gen ){
           # crossing
           if ( (g.t[i] == "f1")||(i ==1) ){
             cross(pop[i],"P1","P2",p.s[i])
           } else
           if (g.t[i] == "bc"){
             cross(pop[i],sl[i-1],"P2",p.s[i])
           } else
           if (g.t[i] == "s"){
             cross(pop[i],sl[i-1],sl[i-1],p.s[i],self=1)     
           } else
           if (g.t[i] == "dh"){
             dh(pop[i],sl[i-1],p.s[i-1])     
           } 
           # selection
           if (i==1) nrp <- "P1" else nrp <-  sl[i-1]
           if (s.s[i]=="t") { 
             select.all.best(s1[i],pop[i],t.name )
             population.sample(sl[i],s1[i],n.s[i],replace=TRUE)
           } else
           if (s.s[i]=="tb") {
             select.all.best(s1[i],pop[i],t.name )
             # Count genome wide HT
             ht[i,j] <- get.population.size(s1[i])$NoInds
             # Count genome wide HT
             select.n.best  (sl[i],s1[i], m.name,n.s[i])    
           } else
           if (s.s[i]=="tf") {
             select.all.best(s1[i],pop[i],t.name )
             # Count flanking SM
             reset.mdp()
             evaluate.mdp(s1[i], nrp,"P2", f.name)
             sm[i,j] <- get.mdp()
             reset.mdp()
             # Count flanking SM
             select.n.best  (sl[i],s1[i], f.name,n.s[i])    
           } else
           if (s.s[i]=="tfb") {
             select.all.best(s1[i],pop[i],t.name )
             # Count flanking SM
             reset.mdp()
             evaluate.mdp(s1[i], nrp,"P2", f.name)
             sm[i,j] <- get.mdp()
             reset.mdp()
             select.all.best ( s2[i], s1[i], f.name )
             # Minimum number preselected
             if ( ! is.null(n.p) ) {
               nnbest <- get.population.size(s2[i])$NoInds
               if ( nnbest < n.p[i] ) {
                 select.n.best ( "tmp-fillup", s1[i], f.name, (n.p[i]-nnbest) )
                 population.concat (s2[i],"tmp-fillup")
               }
             }
             # Count genome wide HT
             ht[i,j] <- get.population.size(s2[i])$NoInds
             select.n.best ( sl[i], s2[i], m.name, n.s[i] )    
           } else
           if (s.s[i]=="tr") {
             select.all.best(s1[i],pop[i],t.name )
             # Count genome wide HT
             ht[i,j] <- get.population.size(s1[i])$NoInds
             # Count genome wide HT
             select.n.best  (sl[i],s1[i], r.name,n.s[i])    
           } else
           if (s.s[i]=="trb") {
             select.all.best(s1[i],pop[i],t.name )
             # Count genome wide HT
             ht[i,j] <- get.population.size(s1[i])$NoInds
             # Count genome wide HT
             select.all.best(s2[i],s1[i], r.name )
             select.n.best  (sl[i],s2[i], m.name,n.s[i])    
           } else
           if (s.s[i]=="trfb") {
             select.all.best(s1[i],pop[i],t.name )
             # Count genome wide HT
             ht[i,j] <- get.population.size(s1[i])$NoInds
             # Count genome wide HT
             select.all.best(s2[i],s1[i], r.name )
             select.all.best(s3[i],s2[i], f.name )
             # Minimum number preselected
             if ( ! is.null(n.p) ) {
               nnbest <- get.population.size(s3[i])$NoInds
               if ( nnbest < n.p[i] ) {
                 select.n.best ( "tmp-fillup", s2[i], f.name, (n.p[i]-nnbest) )
                 population.concat (s3[i],"tmp-fillup")
               }
             }
             select.n.best  (sl[i],s3[i], m.name,n.s[i])    
           } else
           if (s.s[i]=="tfrb") {
             select.all.best(s1[i],pop[i],t.name )
             # Count flanking SM
             reset.mdp()
             evaluate.mdp(s1[i], nrp,"P2", f.name)
             sm[i,j] <- get.mdp()
             reset.mdp()
             select.all.best(s2[i],s1[i], f.name )
             # Minimum number preselected
             if ( ! is.null(n.p) ) {
               nnbest <- get.population.size(s2[i])$NoInds
               if ( nnbest < n.p[i] ) {
                 select.n.best ( "tmp-fillup", s1[i], f.name, (n.p[i]-nnbest) )
                 population.concat (s2[i],"tmp-fillup")
               }
             }
             # Count genome wide HT
             ht[i,j] <- get.population.size(s2[i])$NoInds
             select.all.best(s3[i],s2[i], r.name )
             select.n.best  (sl[i],s3[i], m.name,n.s[i])    
           } else
           if (s.s[i]=="ts") {
             select.all.best(s1[i],pop[i],t.name )
             # Count genome wide HT
             ht[i,j] <- get.population.size(s1[i])$NoInds
             # Count genome wide HT
             select.n.best.segments  (sl[i],s1[i],n.s[i])    
           } else
           if (s.s[i]=="tfs") {
             select.all.best(s1[i],pop[i],t.name )
             # Count flanking SM
             reset.mdp()
             evaluate.mdp(s1[i], nrp,"P2", f.name)
             sm[i,j] <- get.mdp()
             reset.mdp()
             # Count flanking SM
             select.all.best(s2[i],s1[i], f.name )
             if ( ! is.null(n.p) ) {
               nnbest <- get.population.size(s2[i])$NoInds
               if ( nnbest < n.p[i] ) {
                 select.n.best ( "tmp-fillup", s1[i], f.name, (n.p[i]-nnbest) )
                 population.concat (s2[i],"tmp-fillup")
               }
             }
             # Count genome wide HT
             ht[i,j] <- get.population.size(s2[i])$NoInds
             # Count genome wide HT
             select.n.best.segments  (sl[i],s2[i],n.s[i])    
           } else
           if (s.s[i]=="n") {
             population.sample(sl[i],pop[i],n.s[i])
           } else
           {
             info(-2,
              paste("Selection strategy '",
                    s.s[i], "' not implemented ", sep="")
              )
             return(invisible(NULL))
           }
           # success factor
           if ( success.factor < 1) {
             succ <- rbinom( 1 , n.s[i], success.factor )
             if ( succ < 1 ) succ <- 1
             population.sample ( "temp" , sl[i] , succ )
             swap.population.name ( "temp" , sl[i] )
           }
           # success factor
           append.population (st[i],sl[i])
         } # i in 1:n.gen
       } # j in 1:repetitions
     }
   # Analyse simulated data
   # Save output data
   outp <- NULL
   # simulation parameters
   outp$simulation.name  <-    simulation.name 
   outp$simulation.run   <-    simulation.run  
   outp$repetitions      <-    repetitions     
   outp$recurrent.parent <-    recurrent.parent
   outp$donor.parent     <-    donor.parent    
   outp$linkage.map      <-    linkage.map     
   outp$target.loci      <-    target.loci     
   outp$flanking.loci    <-    flanking.loci   
   outp$recipient.loci   <-    recipient.loci     
   outp$gen.type         <-    gen.type        
   outp$population.size  <-    population.size 
   outp$sel.strategy     <-    sel.strategy    
   outp$no.selected      <-    no.selected
   outp$no.preselected   <-    no.preselected
   outp$reg.chr          <-    reg.chr     
   outp$reg.begin        <-    reg.begin   
   outp$reg.end          <-    reg.end     
   outp$success.factor   <-    success.factor

   # Populations and Labels
  if (n.gen > 0 ){
    pare <- c("P1")
    a.pop <- c(pare,st)
    gen.label <- c(pare,paste("G",0:(n.gen-1),"-",g.t,sep=""))
  } else {
    pare <- c("P1")
    a.pop <- c(pare)
    gen.label <- c(pare)
  }
  
   # Target frequency
   if (!is.null(target.loci))
     {
       res <- summarize.gvalue (a.pop,"t")
       rownames(res) <-gen.label
       outp$TargetAlleles <- res
     }
  
   # Genomwide RPG
   res <-  genome.contribution(
                               a.pop,
                               8
                               )
   rownames(res) <-gen.label
   outp$RPG.gw <- res

  # Regionwide RPG
  if (!is.null(reg.chr))
    for(i in 1:length(reg.chr)) {
      if (reg.end[i] ==0)
        reg.end[i] <- genome.parameter.get()$chrom.len[reg.chr[i]]
      reg.begin <- reg.begin*100

      res <- genome.contribution(
                       a.pop,
                       8    ,
                       chromosome = reg.chr[i]   ,
                       begin      = reg.begin[i] ,
                       end        = reg.end[i]
                       )

      cmd <- paste ("outp$RPG.reg",i,"<- res",sep="")
      eval(parse(text=cmd))  
    }

   # Number of donor segments
   res <-  genome.segments(
                           a.pop,
                           1
                          )
   rownames(res) <-gen.label
   outp$DonorSegments.gw <- res

   # Regionwide number of donor segments
  if (!is.null(reg.chr))
    for(i in 1:length(reg.chr)) {
      if (reg.end[i] ==0)
        reg.end[i] <- genome.parameter.get()$chrom.len[reg.chr[i]]
      reg.begin <- reg.begin*100

      res <- genome.segments(
                       a.pop,
                       1    ,
                       chromosome = reg.chr[i]   ,
                       begin      = reg.begin[i] ,
                       end        = reg.end[i]
                       )

      cmd <- paste ("outp$DonorSegments.reg",i,"<- res",sep="")
      eval(parse(text=cmd))  
    }

  
   # Linkage drag
   if (!is.null(target.loci ))
    for(i in 1:length(target.loci)) {

      mm <- get.map()
      ml <- mm[mm$name==target.loci[i],]
      chrom <- ml$chrom
      pos   <- ml$pos
      res <- linkage.drag(
                       a.pop,
                       1,
                       chrom,
                       pos
                       )
      rownames(res) <-gen.label
      cmd <- paste ("outp$LinkageDrag",i,"<- res",sep="")
      eval(parse(text=cmd))  
    }

   # SM Marker data points
   if (n.gen > 0 )
    {
      mean.sm <- matrix(0,nrow=1+n.gen,ncol=1)
      for (i in 1:n.gen) { mean.sm[1+i] <- mean(sm[i,]) }
      rownames(mean.sm) <-c("Total",gen.label[2:length(gen.label)])
      colnames(mean.sm) <-"Mean"
      mean.sm[1] <- sum(mean.sm)
      outp$SM <- mean.sm
    }
  
   # HT Marker data points
   if (n.gen > 0 )
    {
      mean.ht <- matrix(0,nrow=1+n.gen,ncol=1)
      for (i in 1:n.gen) { mean.ht[1+i] <- mean(ht[i,]) }
      rownames(mean.ht) <-c("Total",gen.label[2:length(gen.label)])
      colnames(mean.ht) <-"Mean"
      mean.ht[1] <- sum(mean.ht)
      outp$HT <- mean.ht
    }
  
   # save
   if ( is.null(result.file) ) result.file <-
            (paste(simulation.name,"-",simulation.run,".rda",sep=""))
   save(outp,file=mab.df(st.data.dir,result.file))
  sm.stop.timer()
  return(invisible(outp))

}



mab.load.data <- function(
   simulation.name  ,
   simulation.run   ,
   result.file  = NULL                      
){
   if ( is.null(result.file) ) result.file <-
            (paste(simulation.name,"-",simulation.run,".rda",sep=""))
   load(file=mab.df(st.data.dir,result.file))
   outp
}  



mab.tabulate<- function(
               simulation.names.runs,
               par.summary       = T,
               ext.summary       = F,
               RPG.gw.mean       = F,
               TargetAlleles     = F,
               RPG.gw.Q10        = F,
               RPG.reg.Q10       = F,
               RPG.reg.mean      = F,
               DonorSegments.gw  = F,
               DonorSegments.reg = F,
               LinkageDrag       = F,
               MarkerAnalyses    = F,
               digits            = 1,
               percent           = TRUE
                            )
{
 if (percent) pr <- 100 else pr <- 1
 # Read data
  snr <- matrix(simulation.names.runs,ncol=2,byrow=T)
  n.sc <- nrow(snr)
  outp <- vector("list",n.sc)
  sc.names <- vector("character",n.sc)
  max.gen <- 0
  for (i in 1:n.sc){
    outp[[i]] <-  mab.load.data (
                     simulation.name  = snr[i,1],
                     simulation.run   = snr[i,2] )
    sc.names[i] <- paste(snr[i,1],"-",snr[i,2],sep="")
    ngen <- 1+length(outp[[i]]$gen.type)
    if  (ngen > max.gen) max.gen <- ngen
  }

  # Check whether the target regions are identical
 {
   if (n.sc > 1)
     for (i in 2:n.sc){
       if (
           !(identical( outp[[i]]$target.loci,outp[[i-1]]$target.loci))
           )
         {
           info(-2,"Target loci are not identical in all scenarios")
           return(invisible(NULL))
         }
     }
 }
  # Check whether RPG regions are identical
  if (RPG.reg.mean || RPG.reg.Q10){
      if (n.sc > 1)
        for (i in 2:n.sc){
          if (
                !((identical( outp[[i]]$reg.chr,outp[[i-1]]$reg.chr)
                &&
                identical( outp[[i]]$reg.begin,outp[[i-1]]$reg.begin) )
                &&
                identical( outp[[i]]$reg.begin,outp[[i-1]]$reg.begin) )
              )
            {
              info(-2,"Chromosome regions are not identical in all scenarios")
              return(invisible(NULL))
            }
        }
    }
  if (max.gen>1){
    ge.names <- c("P1",paste("G",0:(max.gen-2),sep=""))
  } else {
    ge.names <- c("P1")
  }
 ret <- NULL
  # Parameter summary
  if ((par.summary) & (!ext.summary))
    {
    n.lines <- 6
    summ <- matrix(nrow=n.lines, ncol=n.sc)
    colnames(summ) <- sc.names
    rn <- vector("character",n.lines)
    for (i in 1:n.sc){
     summ[1,i] <- outp[[i]]$repetitions
     rn[1] <- "Reps"
     summ[2,i] <- outp[[i]]$donor.parent
     rn[2] <- "DP" 
     summ[3,i] <- paste(outp[[i]]$gen.type ,collapse="/")
     rn[3] <- "Pop type" 
     summ[4,i] <- paste(outp[[i]]$population.size,collapse="/")
     rn[4] <- "Pop size"
     summ[5,i] <- paste(outp[[i]]$sel.strategy,collapse="/")
     rn[5] <- "Sel str"
     summ[6,i] <- paste(outp[[i]]$no.selected,collapse="/")
     rn[6] <- "n Sel"
    }
    rownames(summ) <- rn
    ret$par.summary <- summ
  }
  # Extended summary
  if (ext.summary){
    if ( !is.null(outp[[1]]$reg.chr) ) n.lines <- 15 else n.lines <- 12
    summ <- matrix(nrow=n.lines, ncol=n.sc)
    colnames(summ) <- sc.names
    rn <- vector("character",n.lines)
    for (i in 1:n.sc){
     summ[1,i] <- outp[[i]]$repetitions
     rn[1] <- "Reps"
     summ[2,i] <- outp[[i]]$recurrent.parent
     rn[2] <- "RP" 
     summ[3,i] <- outp[[i]]$donor.parent
     rn[3] <- "DP" 
     summ[4,i] <- outp[[i]]$linkage.map
     rn[4] <- "Map" 
     summ[5,i] <- paste(outp[[i]]$target.loci,collapse="/")
     rn[5] <- "Target" 
     summ[6,i] <- paste(outp[[i]]$flanking.loci,collapse="/")
     rn[6] <- "Flm" 
     summ[7,i] <- paste(outp[[i]]$gen.type ,collapse="/")
     rn[7] <- "Pop type" 
     summ[8,i] <- paste(outp[[i]]$population.size,collapse="/")
     rn[8] <- "Pop size"
     summ[9,i] <- paste(outp[[i]]$sel.strategy,collapse="/")
     rn[9] <- "Sel str"
     summ[10,i] <- paste(outp[[i]]$no.selected,collapse="/")
     rn[10] <- "n Sel"
     summ[11,i] <- paste(outp[[i]]$no.preselected,collapse="/")
     rn[11] <- "n PreSel"
     summ[12,i] <- paste(outp[[i]]$success.factor,collapse="/")
     rn[12] <- "Succ f"
     if ( !is.null(outp[[1]]$reg.chr) ){
       summ[13,i] <- paste(outp[[i]]$reg.chr,collapse="/")
       rn[13] <- "Reg chr"
       summ[14,i] <- paste(outp[[i]]$reg.begin,collapse="/")
       rn[14] <- "Reg beg"
       summ[15,i] <- paste(outp[[i]]$reg.end,collapse="/")
       rn[15] <- "Reg end"
     }
   }
    rownames(summ) <- rn
    ret$ext.summary <- summ
  }
  # Tabulate TargetAlleles
  if (TargetAlleles){
    TA <- matrix(nrow=max.gen, ncol=n.sc,
                           dimnames = list(ge.names,sc.names))
    for (i in 1:n.sc){
      ddd <- data.frame(outp[[i]]$TargetAlleles)$Mean
      TA[1:length(ddd),i] <- ddd
      
    }
    ret$TargetAlleles.mean <- round(TA,digits=digits)
  }
  # Tabulate RPG.gw.mean
  if (RPG.gw.mean){
    RPG.mean <- matrix(nrow=max.gen, ncol=n.sc,
                           dimnames = list(ge.names,sc.names))
    for (i in 1:n.sc){
      ddd <-  data.frame(outp[[i]]$RPG.gw)$Mean
      RPG.mean[1:length(ddd),i] <- ddd
    }
      ret$RPG.gw.mean <- round(pr*RPG.mean,digits=digits)
  }
  # Tabulate RPG.gw.Q10
  if (RPG.gw.Q10){
    RPG.Q10 <- matrix(nrow=max.gen, ncol=n.sc,
                           dimnames = list(ge.names,sc.names))
    for (i in 1:n.sc){
      ddd <- data.frame(outp[[i]]$RPG.gw)$Q10
      RPG.Q10[1:length(ddd),i] <- ddd
    }
    ret$RPG.gw.Q10 <- round(pr*RPG.Q10,digits=digits)
  }
  # Tabulate RPG.reg.mean
  if (RPG.reg.mean){
     n.reg <- length(outp[[i]]$reg.chr)
     for (rr in 1:n.reg){
       RPG.mean <- matrix(nrow=max.gen, ncol=n.sc,
                             dimnames = list(ge.names,sc.names))
       for (i in 1:n.sc)
         {
           cmd <- paste("ddd <- data.frame(outp[[",i,"]]$RPG.reg",
                        rr,
                        ")$Mean",
                        sep="")
           eval(parse(text=cmd))
           cmd <- paste("RPG.mean[1:length(ddd),",
                        i,
                        "] <- ddd",
                        sep="")
           eval(parse(text=cmd))
         }
       cmd <- paste("ret$RPG.reg",
                    rr,
                    ".mean <- round(",pr,"*RPG.mean,digits=",digits,")",
                    sep="")
       eval(parse(text=cmd))
     }
  }
  # Tabulate RPG.reg.Q10
  if (RPG.reg.Q10){
     n.reg <- length(outp[[i]]$reg.chr)
     if (n.reg > 0 ) for (rr in 1:n.reg){
       RPG.Q10 <- matrix(nrow=max.gen, ncol=n.sc,
                             dimnames = list(ge.names,sc.names))
       for (i in 1:n.sc)
         {
           cmd <- paste("ddd <- data.frame(outp[[",i,"]]$RPG.reg",
                        rr,
                        ")$Q10"
                        ,sep="")
           eval(parse(text=cmd))
           cmd <- paste("RPG.Q10[length(ddd),",
                        i,
                        "] <- ddd"
                        ,sep="")
           eval(parse(text=cmd))
         }
       cmd <- paste( "ret$RPG.reg",
                      rr,
                      ".Q10 <- round(",pr,"*RPG.Q10,digits=",digits,")",
                    ,sep="")
       eval(parse(text=cmd))
     }
  }
  # Tabulate DonorSegments.gw.mean
  if (DonorSegments.gw){
    DS <- matrix(nrow=max.gen, ncol=n.sc,
                             dimnames = list(ge.names,sc.names))
    for (i in 1:n.sc){
      ddd <- data.frame(outp[[i]]$DonorSegments.gw)$Mean
      DS[1:length(ddd),i] <- ddd
    }
    ret$DonorSegments.gw.mean <- round(DS,digits=digits)
  }

  # Tabulate DonorSegments.reg.mean
  if (DonorSegments.reg){
     n.reg <- length(outp[[i]]$reg.chr)
     if (n.reg > 0 ) for (rr in 1:n.reg){
       DonorSegments <- matrix(nrow=max.gen, ncol=n.sc,
                             dimnames = list(ge.names,sc.names))
       for (i in 1:n.sc)
         {
           cmd <- paste("ddd <- data.frame(outp[[",i,"]]$DonorSegments.reg",
                        rr,
                        ")$Mean",
                        sep="")
           eval(parse(text=cmd))
           cmd <- paste("DonorSegments[1:length(ddd),",
                        i,
                        "] <- ddd",
                        sep="")
           eval(parse(text=cmd))
         }
       cmd <- paste("ret$DonorSegments.reg",
                    rr,
                    ".mean <- round(DonorSegments,digits=",digits,")",
                    sep="")
       eval(parse(text=cmd))
     }
  }

  # Tabulate linkage drag
   if (LinkageDrag){
     n.tar <- length(outp[[i]]$target.loci)
     for (rr in 1:n.tar){
       LD <- matrix(nrow=max.gen, ncol=n.sc,
                             dimnames = list(ge.names,sc.names))
       for (i in 1:n.sc)
         {
           cmd <- paste("ddd <- data.frame(outp[[",i,"]]$LinkageDrag",
                        rr,
                        ")$Mean"
                        ,sep="")
           eval(parse(text=cmd))
           cmd <- paste("LD[1:length(ddd),",
                        i,
                        "] <- ddd"
                        ,sep="")
           eval(parse(text=cmd))

         }
       cmd <- paste( "ret$LinkageDrag",
                      rr,
                      ".mean <- round(100*LD,digits=",digits,")"
                    ,sep="")
       eval(parse(text=cmd))
     }
  }
  # Tabulate MarkerAnalyses
  ma.names <- c("Total",ge.names[2:length(ge.names)])
  if (MarkerAnalyses){
    SM <- matrix(nrow=max.gen, ncol=n.sc,
                           dimnames = list(ma.names,sc.names))
    for (i in 1:n.sc){
      ddd <- data.frame(outp[[i]]$SM)$Mean
      SM[1:length(ddd),i] <- ddd
      
    }
    ret$SM.mean <- round(SM,digits=digits)
    HT <- matrix(nrow=max.gen, ncol=n.sc,
                           dimnames = list(ma.names,sc.names))
    for (i in 1:n.sc){
      ddd <- data.frame(outp[[i]]$HT)$Mean
      HT[1:length(ddd),i] <- ddd
      
    }
    ret$HT.mean <- round(HT,digits=digits)
  }

   ret
}

mab.save.inds <- function(pop,pname){
   x <- return.population(pop)
   for (j in 1:ncol(x) ){
     mdta <- NULL  
     i <- 1
     while(i < nrow(x)) {
            n1 <- splitdt(row.names(x)[i])[1]
            n2 <- splitdt(row.names(x)[i+1])[1]
            l1 <- splitdt(row.names(x)[i])[2]
            l2 <- splitdt(row.names(x)[i+1])[2]
            if( n1 == n2){
              a1 <- x[i,j]
              a2 <- x[i+1,j]
              if ( (a1 == 1)  && (a2 == 0) ) mgt <- c(l1,l1)  else
              if ( (a1 == 1)  && (a2 == 1) ) mgt <- c(l1,l2)  else
              if ( (a1 == 0)  && (a2 == 1) ) mgt <- c(l2,l2)
              mdta <- rbind(mdta,c(n1,mgt))
              i <- i +2
            } else {
              mgt <- c(l1,l1)  
              mdta <- rbind(mdta,c(n1,mgt))
              i <- i+1 
            }
          }
     if (i==nrow(x)){
       n1 <- splitdt(row.names(x)[i])[1]
       l1 <- splitdt(row.names(x)[i])[2]
       mgt <- c(l1,l1)  
       mdta <- rbind(mdta,c(n1,mgt))
       i <- i+1 
     }
   iname = paste( pname,  sprintf("%03d",j), sep="")
   fname = paste(iname, ".mda", sep="")
   write(iname,file=mab.df(st.output.dir,fname),append=FALSE)
   write.table(
               mdta,
               file=mab.df(st.output.dir,fname),
               quote=F,
               row.names=F,
               col.names=F,
               append=TRUE
     )
   }
}


mab.compare <- function (
   a.simulation.name  ,
   a.simulation.run   ,
   a.repetitions      = 1000,      
   a.recurrent.parent = "H",
   a.donor.parent     = "H",
   a.linkage.map      ,
   a.target.loci      = NULL,
   a.flanking.loci    = NULL,
   a.recipient.loci   = NULL,
   a.gen.type         = NULL,
   a.population.size  = NULL,
   a.sel.strategy     = NULL,
   a.no.selected      = NULL,
   a.no.preselected   = NULL,
   a.reg.chr          = NULL,                            
   a.reg.begin        = NULL,                            
   a.reg.end          = NULL,                            
   a.missing.allele   = "9" ,
   a.success.factor   = 1   ,                      
#
   b.simulation.name =NULL, c.simulation.name =NULL, d.simulation.name =NULL,
   b.simulation.run  =NULL, c.simulation.run  =NULL, d.simulation.run  =NULL,
   b.repetitions     =NULL, c.repetitions     =NULL, d.repetitions     =NULL,      
   b.recurrent.parent=NULL, c.recurrent.parent=NULL, d.recurrent.parent=NULL,
   b.donor.parent    =NULL, c.donor.parent    =NULL, d.donor.parent    =NULL,
   b.linkage.map     =NULL, c.linkage.map     =NULL, d.linkage.map     =NULL,
   b.target.loci     =NULL, c.target.loci     =NULL, d.target.loci     =NULL,
   b.flanking.loci   =NULL, c.flanking.loci   =NULL, d.flanking.loci   =NULL,
   b.recipient.loci  =NULL, c.recipient.loci  =NULL, d.recipient.loci  =NULL,
   b.gen.type        =NULL, c.gen.type        =NULL, d.gen.type        =NULL,
   b.population.size =NULL, c.population.size =NULL, d.population.size =NULL,
   b.sel.strategy    =NULL, c.sel.strategy    =NULL, d.sel.strategy    =NULL,
   b.no.selected     =NULL, c.no.selected     =NULL, d.no.selected     =NULL,
   b.no.preselected  =NULL, c.no.preselected  =NULL, d.no.preselected  =NULL,
   b.reg.chr         =NULL, c.reg.chr         =NULL, d.reg.chr         =NULL,
   b.reg.begin       =NULL, c.reg.begin       =NULL, d.reg.begin       =NULL,
   b.reg.end         =NULL, c.reg.end         =NULL, d.reg.end         =NULL,
   b.missing.allele  =NULL, c.missing.allele  =NULL, d.missing.allele  =NULL,
   b.success.factor  =NULL, c.success.factor  =NULL, d.success.factor  =NULL,
#
   e.simulation.name =NULL, f.simulation.name =NULL, g.simulation.name =NULL,
   e.simulation.run  =NULL, f.simulation.run  =NULL, g.simulation.run  =NULL,
   e.repetitions     =NULL, f.repetitions     =NULL, g.repetitions     =NULL,      
   e.recurrent.parent=NULL, f.recurrent.parent=NULL, g.recurrent.parent=NULL,
   e.donor.parent    =NULL, f.donor.parent    =NULL, g.donor.parent    =NULL,
   e.linkage.map     =NULL, f.linkage.map     =NULL, g.linkage.map     =NULL,
   e.target.loci     =NULL, f.target.loci     =NULL, g.target.loci     =NULL,
   e.flanking.loci   =NULL, f.flanking.loci   =NULL, g.flanking.loci   =NULL,
   e.recipient.loci  =NULL, f.recipient.loci  =NULL, g.recipient.loci  =NULL,
   e.gen.type        =NULL, f.gen.type        =NULL, g.gen.type        =NULL,
   e.population.size =NULL, f.population.size =NULL, g.population.size =NULL,
   e.sel.strategy    =NULL, f.sel.strategy    =NULL, g.sel.strategy    =NULL,
   e.no.selected     =NULL, f.no.selected     =NULL, g.no.selected     =NULL,
   e.no.preselected  =NULL, f.no.preselected  =NULL, g.no.preselected  =NULL,
   e.reg.chr         =NULL, f.reg.chr         =NULL, g.reg.chr         =NULL,
   e.reg.begin       =NULL, f.reg.begin       =NULL, g.reg.begin       =NULL,
   e.reg.end         =NULL, f.reg.end         =NULL, g.reg.end         =NULL,
   e.missing.allele  =NULL, f.missing.allele  =NULL, g.missing.allele  =NULL,
   e.success.factor  =NULL, f.success.factor  =NULL, g.success.factor  =NULL,
#
   h.simulation.name =NULL, i.simulation.name =NULL,
   h.simulation.run  =NULL, i.simulation.run  =NULL,
   h.repetitions     =NULL, i.repetitions     =NULL,     
   h.recurrent.parent=NULL, i.recurrent.parent=NULL,
   h.donor.parent    =NULL, i.donor.parent    =NULL,
   h.linkage.map     =NULL, i.linkage.map     =NULL,
   h.target.loci     =NULL, i.target.loci     =NULL,
   h.flanking.loci   =NULL, i.flanking.loci   =NULL,
   h.recipient.loci  =NULL, i.recipient.loci  =NULL,
   h.gen.type        =NULL, i.gen.type        =NULL,
   h.population.size =NULL, i.population.size =NULL,
   h.sel.strategy    =NULL, i.sel.strategy    =NULL,
   h.no.selected     =NULL, i.no.selected     =NULL,
   h.no.preselected  =NULL, i.no.preselected  =NULL, 
   h.reg.chr         =NULL, i.reg.chr         =NULL,
   h.reg.begin       =NULL, i.reg.begin       =NULL,
   h.reg.end         =NULL, i.reg.end         =NULL,
   h.missing.allele  =NULL, i.missing.allele  =NULL,
   h.success.factor  =NULL, i.success.factor  =NULL
)
{
  n.scn <- 0

  prefix <- c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
  suffix <- c("simulation.name", 
              "simulation.run",  
              "repetitions",     
              "recurrent.parent",
              "donor.parent",    
              "linkage.map",     
              "target.loci",     
              "flanking.loci",
              "recipient.loci",
              "gen.type",        
              "population.size", 
              "sel.strategy",    
              "no.selected", 
              "no.preselected", 
              "reg.chr",     
              "reg.begin",  
              "reg.end",    
              "missing.allele",
              "success.factor" )

  an <- c( "sn",
           "sr",
           "re",
           "rp",
           "dp",
           "li",
           "tl",
           "fl",
           "rl",
           "gt",
           "ps",
           "ss",
           "ns",
           "np",
           "gh",
           "gb",
           "ge",
           "am",
           "sf"
          )

 for (i in 1:length(an)){
   cmd <- paste(an[i], " <- vector('list',9)",sep="")
   eval(parse(text=cmd))
 }


  for (i in 1:9)
    {
      if (
          eval(parse(text=paste(
           " !is.null(",
           prefix[i],"simulation.name",
          ")" ,
           sep="" )))
          )
        {
  
          n.scn   <- n.scn+1

          sn[[i]] <-
            eval(parse(text=paste(
                         prefix[i],"simulation.name",
                         sep="" )))
          
          for (k in 2:length(an)){
            cmd <- paste( "if( !is.null( a.",suffix[k],") ) ",
                         an[k],"[[",i,"]]","<- a.", suffix[k],
                         sep="" )
            eval(parse(text=cmd))
            cmd <- paste( "if( !is.null(",prefix[i],suffix[k],") ) ", 
                         an[k],"[[",i,"]]", "<-", prefix[i],suffix[k],
                         sep="" )
            eval(parse(text=cmd))
          }
          
        }
    }
 
  
  for (i in 1:n.scn){
    mab.simulate(
                   simulation.name   =  sn[[i]] , 
                   simulation.run    =  sr[[i]] ,
                   repetitions       =  re[[i]] ,
                   recurrent.parent  =  rp[[i]] ,
                   donor.parent      =  dp[[i]] ,
                   linkage.map       =  li[[i]] ,
                   target.loci       =  tl[[i]] ,
                   flanking.loci     =  fl[[i]] ,
                   recipient.loci    =  rl[[i]] ,
                   gen.type          =  gt[[i]] ,
                   population.size   =  ps[[i]] ,
                   sel.strategy      =  ss[[i]] ,
                   no.selected       =  ns[[i]] ,
                   no.preselected    =  np[[i]] ,
                   reg.chr           =  gh[[i]] ,                            
                   reg.begin         =  gb[[i]] ,                            
                   reg.end           =  ge[[i]] ,                            
                   missing.allele    =  am[[i]] ,
                   success.factor    =  sf[[i]] 
                   )
  }
}

