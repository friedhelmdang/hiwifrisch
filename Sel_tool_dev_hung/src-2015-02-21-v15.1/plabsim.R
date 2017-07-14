plabsim.R.version <- "15.1"
plabsim.R.date    <- "2015/01/13"

plabsim <- list()
plabsim$deprecated.info <- FALSE
plabsim$allele.missing.indicator <- 3333

"info" <-  function(lev,msg){
    c <- .C("r_info",
            as.integer(lev),
            as.character(msg)
            );
  }

"info.cat" <-  function(lev,msg){

  level = 0
  c  <-  .C(  "get_info_level",
               level = as.integer(level)   )

  if ( lev >= c$level) cat ( sprintf( "M: %s\r", msg ))
  
  }


"write.version.2" <-
  function(){
    c <- .C("print_version_2",
            msg = as.character(plabsim.R.date)
            );
  }

"plabsim.init" <- function() {
  rng.init()
}


"population.append" <- function(
  NameP1,
  NameP2
) {
  c  <- .C(
    "append_population",
    as.character(NameP1),
    as.character(NameP2)
  )
}

"population.concat" <- function(
  NameP1,
  NameP2
) {
  c  <- .C(
    "concat_population",
    as.character(NameP1),
    as.character(NameP2)
  )
}

"population.copy" <- function(
  NameP1,
  NameP2,
  start=1,
  n=-1
) {
  c  <- .C(
    "copy_population",
    as.character(NameP1),
    as.character(NameP2),
    format(start,scientific=FALSE),
    format(n,scientific=FALSE)
  )
}

"population.divide" <- function(
  NameP1,
  NameP2,
  NoI=-1
) {
  c  <- .C(
    "divide_population",
    as.character(NameP1),
    as.character(NameP2),
    format(NoI,scientific=FALSE)
  )
}

"population.exist" <- function(PopName) {
   invisible(as.logical(sum(as.character(population.list()$PopName)==PopName)))
}

"population.info.get" <- function(
  PopName,
  ind
) {
  if ((ind > 0) && (ind <= population.size.get(PopName)$NoInds)) {
    info <- paste(rep("a",36),collapse="")
    c  <- .C(
      "get_population_info",
      as.character(PopName),
      format(ind,scientific=FALSE),
      info = as.character(info)
    )
    c$info
  }
}

"population.name.swap" <- function(
  NameP1,
  NameP2
) { 
  c  <- .C(
    "swap_population_name",
    as.character(NameP1),
    as.character(NameP2)
  )
}

"population.optimize" <-
  function(PopNames){
    c  <- .C("optimize_population",
             as.character(PopNames))
  }

"population.remove" <-  function(
  PopNames
) {
  if (length(PopNames) > 1) {
    PopNames <- paste(PopNames, collapse=" ")
  }
  c  <- .C(
    "remove_population",
    as.character(PopNames)
  )
}

"population.remove.all" <- function(){
  y <- as.vector((population.list()$PopName))
  allPops <- paste(y,collapse = " ")
  population.remove(allPops)
  info(1,"All populations removed")
}

"population.rename" <- function(
  OldName,
  NewName
) {
  c <- .C(
    "rename_population",
    as.character(OldName),
    as.character(NewName)
  )
}

"population.resize" <- function(
  PopName,
  newSize
) {
  c <- .C(
    "resize_population",
    as.character(PopName),
    format(newSize,scientific=FALSE)
  )
}

"population.sample" <-
  function(NameP1, NameP2, size=-1, replace=FALSE) {
    if (replace) {
      rep <- 1
    } else {
      rep <- 0
    }
    c  <- .C("sample_population",
             as.character(NameP1),
             as.character(NameP2),
             format(size,scientific=FALSE),
             as.integer(rep))
  }

"population.size.get" <- function(
  PopName
) {
  if (length(PopName) > 1) {
    PopName <- paste(PopName, collapse=" ")
  }
  digits <- 50;
  sNoInds <- paste(rep("a",digits),collapse="");
  sNoIndsAlloc <- paste(rep("a",digits),collapse="");
  c  <- .C(
    "get_population_size",
    as.character(PopName),
    NoInds=as.character(sNoInds),
    NoIndsAlloc=as.character(sNoIndsAlloc)
  )
  if( sNoInds==c$NoInds ) return ( invisible(NULL) )
  data.frame(
    NoInds=as.integer(c$NoInds),
    NoIndsAlloc=as.integer(c$NoIndsAlloc)
  )
  }

"population.transfer" <- function(
  NameP1,
  NameP2,
  population.NameP2.delete=TRUE
) {
  if (population.NameP2.delete) {
    population.append(NameP1, NameP2)
  } else {
    population.concat(NameP1, NameP2)
  }
}

"be.quiet" <-
  function(){
    set.info.level(0);
  }

"cross" <- function(
  NamePg,
  NameP1,
  NameP2,
  NoPg = -1,
  self = 0,
  mode = 1,
  crossing.scheme = 0
) {
  significant <- -1
  c  <- .C(
    "cross",
    as.character(NamePg),
    as.character(NameP1),
    as.character(NameP2),
    as.double(self),
    format(NoPg,scientific=FALSE),
    as.double(significant),
    as.integer(mode),
    as.integer(crossing.scheme)
  )
}

"define.effects" <- function (
  name,
  description
) {
  tmp.file <- tempfile()
  generate.effect.file(tmp.file, description)
  load.effmap(name, tmp.file)
  unlink(tmp.file)
}

"define.map" <- function (
  description
) {
  tmp.file <- tempfile()
  generate.map.file(tmp.file,description);
  linkage.map.load(tmp.file);
  unlink(tmp.file)
}

"dh" <- function (
  NamePg,
  NameP1,
  NoPg=-1
) {
  significant <- -1;
  c  <- .C(
    "dh",
    as.character(NamePg),
    as.character(NameP1),
    format(NoPg,scientific=FALSE),
    as.double(significant)
  )
}

"evaluate.allele" <- function(
  PopName,
  Loci,
  missing=0
) {
  if (1==missing) {
    evallo <- evaluate.locus(PopName,Loci)
    summe <- evallo[,1] + evallo[,2]
    locivec <- unlist(strsplit(Loci," "))
    Loci <- paste(locivec[which(summe>0)],collapse=" ")
  }
  if (0 < nchar(Loci)) {
    bg <- plabsim$allele.missing.indicator
    result <- .Call("allele_evaluate",
                    as.character(PopName),
                    as.character(unlist(strsplit(Loci," "))),
                    as.integer(missing),
                    as.integer(bg));
    if ( is.null(result) ) return (invisible( NULL) )
    result <- as.data.frame(result)
    result <- result[,-(dim(result)[2])]
  } else {
    return(invisible( NULL))
  }
  result
}

"evaluate.allele.freq" <- function( PopName,
                                    Loci,
                                    missing=0 )
{
  if ( 1 <length(Loci) ){
    info(-2,"Only available for single loci")
    return(invisible(NULL))
  }
  e.a <- evaluate.allele(PopName,Loci,missing)
  if (is.null(e.a)) return(invisible(NULL))
  e.a$freq <-e.a$count/sum(e.a$count)
  return(e.a)
}


"evaluate.genome" <- function (
  PopName,
  allele,
  chromosome=NULL,
  begin=NULL,
  end=NULL,
  verbose=FALSE
) {
  empty <- rep(0,as.integer(population.size.get(PopName)$NoInds));
  if (is.null(chromosome)) {
    c  <- .C("evaluate_genome",
             as.character(PopName),
             as.integer(allele),
             ratio  = as.double(empty),
             length = as.double(empty),
             total  = as.double(empty),
             blocks = as.integer(empty))
  } else {
    if (is.null(begin)) {
      begin <- 0
    }
    if (is.null(end)) {
      end <- genome.parameter.get()$chrom.len[chromosome]
    }
    c  <- .C("evaluate_genome_region",
             as.character(PopName),
             as.integer(allele),
             as.integer(chromosome-1),
             as.double(begin),
             as.double(end),
             ratio  = as.double(empty),
             length = as.double(empty),
             total  = as.double(empty),
             blocks = as.integer(empty))
  }
  ratio  <- data.frame(ratio=c$ratio)
  length <- data.frame(length=c$length)
  total  <- data.frame(total=c$total)
  blocks <- data.frame(noBlocks=c$blocks)
  if (verbose) {
    return.value <- cbind(ratio, length)
    return.value <- cbind(return.value, total)
    return.value <- cbind(return.value, blocks)
  } else {
    return.value <- ratio
  }   
  return.value
}

"evaluate.genotype" <- function(PopName, Loci, alleles=0, mode=2, bg=0) {
  result <- .Call("genotype_evaluate",
                  as.character(PopName),
                  as.character(unlist(strsplit(Loci," "))),
                  as.integer(mode));
  if ( is.null(result) ) return( invisible(NULL) )
  if (bg == 1) {
    col.name <- colnames(result)
    res <- matrix(result[,1:(dim(result)[2]-2)],ncol=dim(result)[2]-2)
    result <- matrix(result[(0==rowSums(plabsim$allele.missing.indicator == res)),],ncol=dim(result)[2])
    result[,dim(result)[2]] <- result[,dim(result)[2]-1]/sum(result[,dim(result)[2]-1])
    colnames(result) <- col.name
  } else {
    result[,dim(result)[2]] <- result[,dim(result)[2]-1]/result[,dim(result)[2]]
  }
  if (alleles[1] != 0) {
    ng <- nrow(result)
    na <- length(alleles)
    passt <- 0; i <- 0
    while ( (i < ng) && (passt==0) ){
      passt <- 1 ; i <- i+1
      for (j in 1:na){
        if (result[i,j] != alleles[j]) {passt <- 0}
      }
      if (passt==1) {select <- i} 
    }
    if (passt==1) {result <- result[select,]}
  }  
  as.data.frame(result)
}

"evaluate.haplotype" <-
  function(PopName,Loci,missing=0,bg=plabsim$allele.missing.indicator){
    genotype.population(PopName);
    y <- .Call("evaluate_haplotype",
             as.character(PopName),
             as.character(Loci),
             as.integer(missing),
             as.integer(bg));
    y <- y[-1,]
    if (!is.null(dim(y))) {
      y <- as.data.frame(y);
    } else {
      y <- as.data.frame(t(y));
    }
    y
  }

"evaluate.hap.init" <-
  function(PopName, doubleHetero=0, missing=0){
# doubleHetero 0: double heterozygotes are not removed
# doubleHetero 1: double heterozygotes are removed
# missing 0: missing alleles are removed
# missing 1: missing alleles are not removed
    bg <- plabsim$allele.missing.indicator  
    genotype.population(PopName);
    c <- .C("haplotype_init",
            as.character(PopName),
            as.integer(doubleHetero),
            as.integer(missing),
            as.integer(bg));
  }

"evaluate.hap.calc" <-
  function(locusA,locusB){
    c <- .C("haplotype_next",
            as.integer(locusA),
            as.integer(locusB),
            noTI = as.integer(0),
            noI  = as.integer(0),
            nodH = as.integer(0),
            dLocA= as.integer(0),
            dLocB= as.integer(0));
    c(c$noTI,c$noI,c$nodH,c$dLocA,c$dLocB);
  }

"evaluate.hap.calc.d" <-
  function(locusA,locusB){
    c <- .C("haplotype_calc_d",
            dd = as.double(0),
            dp = as.double(0),
            rr = as.double(0))
            
    invisible(c(c$dd, c$dp, c$rr))
  }

"evaluate.hap.return.genotype" <-
  function(){
    y <- .Call("haplotype_return_genotype");
    as.data.frame(y);  
  }

"evaluate.hap.return.table" <-
  function(){
    y <- .Call("haplotype_return_table");
    as.data.frame(y);  
  }

"evaluate.hap.free" <-
  function(){
    c <- .C("haplotype_free");
  }

"evaluate.locus" <- function(
  PopName,
  eLoci
) {
  bg <- plabsim$allele.missing.indicator  
  genotype.population(PopName);
  noLoci <- length(unlist(strsplit(eLoci," ")));
  empty <- rep(0,3 * noLoci)
  c  <- .C(
    "evaluate_locus",
    as.character(PopName),  
    as.character(eLoci),
    matrix=as.integer(empty),
    as.integer(noLoci),
    as.integer(bg)
  )
  dim(c$matrix) <- c(noLoci, 3)
  result <- data.frame(c$matrix);
  names(result) <- unlist(c("homozygote","heterozygote","missing"));
  rownames(result) <- unlist(strsplit(eLoci," "));
  result;
}

"evaluate.all.loci" <- function(
  PopName,
  loci = paste(as.vector(unlist(linkage.map.get()[,3])), collapse=" ")
) {
  evaluate.locus(
    PopName,
    loci
  )    
}

"evaluate.mdp" <- function(
  NamePg,
  NameP1,
  NameP2,
  effectfile=NULL
) {
  genotype.population(NameP1)
  genotype.population(NameP2)
  if (is.null(effectfile)) {
    abc  <- list.effects()
    for (i in 1:nrow(abc)) { 
      if (0 != abc[i,2]) {
        evaluate.mdp(NamePg, NameP1, NameP2, abc[i,1])
      }
    }
  } else {
    c <- .C(
      "evaluate_mdp",
      as.character(NamePg),
      as.character(NameP1),
      as.character(NameP2),
      as.character(effectfile)
    )
  }
}

"evaluate.population" <- function(
  PopNames,
  effName=NULL
) {
  c <- .C(
    "evaluate_population",
    as.character(PopNames),
    as.character(effName)
  )
}

"generate.effect.file" <- function(
  fName,
  description
) { 
  c <- .C(
    "generate_effect_file",
    as.character(fName),
    as.character(description)
  )
}

"generate.map.file" <- function(
  fName,
  description
) {
  if (length(description)>1) {
    description <- paste(description,collapse=", ")
  }
  c <- .C(
    "generate_map_file",
    as.character(fName),
    as.character(description)
  )
}

"generate.population" <- function(
  dta,
  backcross=FALSE
) {
  bg <- plabsim$allele.missing.indicator
  p <- data.params(dta)
  no.pop.u <- p$no.pop
  pop.u    <- 1:p$no.pop
  no.ind.u <- p$no.ind
  ind.u    <- p$ind.list
  no.mar.u <- p$no.mar
  mar.u    <- 1:p$no.mar
  DTA  <- as.matrix(dta)
  FREQ <- matrix(16,nrow=sum(p$no.all[mar.u]),ncol=p$no.pop)
  pop.list <- paste(as.vector(p$pop.list),collapse = " ")

  #old <- options(warn = -1)
  #on.exit(options(old))
  b <- as.data.frame(as.numeric(p$rowhead[,2]));
  colnames(b) <- c("allele")
  if ( (any(is.na(b$allele))) || (any(b$allele==0)) ) {
    info(-2,"Error in allele encoding: only objects of type integer allowed")
  } else {
    c  <- .C(
             "generate_population",
             as.integer(p$no.pop),
             as.integer(p$no.ind),
             as.integer(p$no.mar),
             as.integer(p$no.all),
             as.integer(no.pop.u),
             as.integer(pop.u),
             as.integer(no.ind.u),
             as.integer(ind.u),
             as.integer(no.mar.u),
             as.integer(mar.u),
             DTA =as.double(DTA),
             FREQ=as.double(FREQ),
    as.character(paste(as.vector(p$pop.list),collapse = " ")),
    as.character(paste(as.vector(p$mar.list),collapse = " ")),
             as.integer(b$allele),
             as.integer(bg),
             NAOK=TRUE,
             as.integer(backcross)
             )
  }
}

"genome.contribution" <-
  function (
            pops,
            allele,
            chromosome=NULL,
            begin=NULL,
            end=NULL
            ) 
{
  p <- unlist(strsplit(pops, " "))
  n <- length(p)
  result <- matrix(ncol = 7, nrow = n)
  for (i in 1:n) {
    g <- evaluate.genome(
            p[i],
            allele,
            chromosome = chromosome,
            begin      = begin,
            end        = end
          )
    result[i, 1] <- nrow(g)
    result[i, 2] <- mean(g$ratio)
    result[i, 3] <- sdev(g$ratio)
    result[i, 4] <- min(g$ratio)
    result[i, 5] <- quantile(g$ratio, 0.1)
    result[i, 6] <- quantile(g$ratio, 0.5)
    result[i, 7] <- max(g$ratio)
  }
  colnames(result) <- c("Obs", "Mean", "SDev", "Min", "Q10", 
                        "Med", "Max")
  rownames(result) <- p
  result
}


"genotype.population" <-
  function(PopNames){
    c  <- .C("genotype_population",
             as.character(PopNames))
  }

"genome.parameter.get" <- function() {
  no.chrom <- 0;
  no.hom   <- 0;
  c  <- .C("get_genome_par",
           no.chrom=as.integer(no.chrom),
           no.hom  =as.integer(no.hom))
  d  <- .C("get_genome_par_b",
           chrom.len = as.double(rep(0,c$no.chrom)))
  list (
    no.chrom  = c$no.chrom,
    no.hom    = c$no.hom,
    chrom.len = d$chrom.len 
  )
}

"linkage.map.get" <- function() {
  digits <- 50;
  no.mappoints <- 0;
  maxname      <- 0;
  classlgth <- namelgth <- paste(rep("a",digits),collapse="");
  c  <- .C("get_map_a",
           no.mappoints = as.integer(no.mappoints),
           maxname      = as.integer(maxname),
           classlgth    = as.character(classlgth),
           namelgth     = as.character(namelgth))
  chrom <- rep(0,c$no.mappoints)
  pos   <- rep(0,c$no.mappoints)
  class <- paste(rep(" ", c$classlgth), collapse="")
  name  <- paste(rep(" ", c$namelgth), collapse="")
  d  <- .C("get_map_b",
           chrom = as.double(chrom),
           pos   = as.double(pos),
           name  = as.character(name),
           class = as.character(class))
  name <- unlist(strsplit(d$name," "))
  class <- unlist(strsplit(d$class," "))
  data.frame(
    chrom = d$chrom,
    pos   = d$pos,
    name  = name,
    class = class
  )
}

"get.mdp" <- function() {
  digits <- 50;
  cmdp <- paste(rep("a",digits),collapse="");
  c  <- .C(
    "get_mdp",
    mdp = as.character(cmdp)
  )
  as.integer(c$mdp)
}

"set.mdp" <- function(
  MDP=0
) {
  c <- .C(
    "set_mdp",
    format(MDP,scientific=FALSE)
  )
}

"get.population" <- function(
  PopName
) {
  count <- 0;
  c  <- .C(
    "get_allele_number",
    as.character(PopName),
    count=as.integer(count)
  )
  empty <- rep(0,c$count)
  c  <- .C(
    "get_population",
    as.character(PopName),
    ind   = as.integer(empty),
    chrom = as.integer(empty),
    hom   = as.integer(empty),
    pos   = as.double (empty),
    all   = as.integer(empty)
  )
  data.frame(
    ind   = c$ind,
    chrom = c$chrom,
    hom   = c$hom,
    pos   = c$pos,
    all   = c$all
  )
}

"get.population.gvalue" <-function(
  name,
  EffName=NULL
) {
  empty <- rep(0,as.integer(population.size.get(name)$NoInds))
  c  <- .C(
    "get_population_gvalue",
    as.character(name),
    gvalue = as.double(empty),
    as.character(EffName)
  )
  data.frame(gvalue=c$gvalue)
}

"get.score" <-
  function(pop,effectfile=NULL){
    genotype.population(pop);
    evaluate.population(pop);
    g <- get.population.gvalue(pop,effectfile)
    g
  }


"homozygote" <-  function (
  allele,
  NoInd = 1
) {
  c  <- .C(
    "get_genome_par",
    no.chrom = as.integer(0),
    no.hom   = as.integer(0)
  )  
  pop  <-  NULL
  for (ind in 1:NoInd) {
    for (chrom in 1:c$no.chrom) {
      for (hom in 1:c$no.hom) {
        pop <- rbind(
          pop,
          c(ind, chrom, hom, 0, allele)
        )
      } #hom
    } #chrom
  } #ind
  pop <- as.data.frame(pop)
  names(pop) <- c("ind","chrom","hom","pos","all")
  pop
}

"plant" <-  function (
  genotype,
  NoInd=1
) {
  c  <- .C(
    "get_genome_par",
    no.chrom = as.integer(0),
    no.hom   = as.integer(0)
  )  
  pop  <-  NULL
  for (ind in 1:NoInd) {
    for (chrom in 1:c$no.chrom) {
      for (hom in 1:c$no.hom) {
        pop <- rbind(
          pop,
          c(ind, chrom, hom, 0, genotype[hom])
        )
      } #hom
    } #chrom
  } #ind
  pop <- as.data.frame(pop)
  names(pop) <- c("ind","chrom","hom","pos","all")
  pop
}

"init.population" <- function(
  name,
  data,
  delete=TRUE
) {
  if (5 != dim(data)[2]) stop("Wrong data format.")
  if (delete == TRUE) {
    population.remove(name)
  }
  c <- .C(
    "init_population",
    as.character(name),
    format(nrow(data),scientific=FALSE),
    as.integer(data$ind),
    as.integer(data$chrom),
    as.integer(data$hom),
    as.double(data$pos),
    as.integer(data$all)
  )
}

"list.effects" <-
  function(){
    lEffName <- 0;
    nrEffects   <-  0;
    b <- .C("list_eff_parameter",
            lEffName=as.integer(lEffName),
            nrEffects=as.integer(nrEffects));
    empty  <- rep(0,b$nrEffects);
    name <- paste(rep("a",b$lEffName),collapse="")
    c  <- .C("list_eff",
             name=as.character(name),
             weight=as.double(empty))
    x <- data.frame(strsplit(c$name," "),c$weight);
    names(x)  <- c("effect","weight");
    x
  }

"population.list" <-
  function(){
    lPopName <- 0
    NoPops <- 0
    b <- .C("list_populations_parameter",
            lPopName=as.integer(lPopName),
            NoPops=as.integer(NoPops))
    empty <- rep(0,b$NoPops)
    name <- paste(rep("a",b$lPopName),collapse="")
    c  <- .C("list_populations",
             name=as.character(name),
             inds=as.integer(empty) )
    PopName <- strsplit(c$name," ")
    x <- data.frame(PopName,c$inds)
    names (x) <- c("PopName","count")
    x
  }

"load.effmap" <- function(
  name,
  file = NA
) {
  if (is.na(file)) {
    file <- name
  }
  allPops <- paste(as.vector(population.list()$PopName),collapse = " ")
  remove.evaluate.population(allPops)
  c  <- .C(
    "load_effmap",
    as.character(file),
    as.character(name)
  )
}

"load.internal.effmap" <-
  function(name, spec, weight){
    c <- population.list();    
    y <- as.vector((c$PopName));
    allPops <- paste(y,collapse = " ");
    remove.evaluate.population(allPops)
    b  <- .C("load_effmap_mod",
             as.character(name),
             as.character(spec),
             as.double(weight));
  }

"linkage.map.load" <-
  function(file, disperse=FALSE, file.disperse=NA, disperse.factor=100){
    allPops <- paste(as.vector((population.list()$PopName)),collapse = " ");
    remove.evaluate.population(allPops);  
    remove.effmaps();
    remove.genotype.population(allPops);
    remove.map();
    if (disperse) {
      map     <- read.table(file=file)
      left    <- 1
      right   <- 1;
      counter <- 0;
      if (dim(map)[1]>1) { # more than one marker locus
        for (line in 2:dim(map)[1]) {
          if ((map[line-1,2] != map[line,2])||((map[line-1,1] != map[line,1]))) {
            if (left<right) {
              counter <- counter + right - left
              if ((map[line-1,1] == map[line,1])) { #Marker Loci befinden sich nicht am Ende des Chromosoms
                delta <- (map[line,2]-map[left,2])/(disperse.factor*(right-left+1))
                map[left:right,2] <- map[left:right,2] + (0:(right-left))*delta
              } else {
                if (map[left,2] == get.genome.par()$chrom.len[map[line-1,1]]) { #vorher | marker loci sitzen am Ende des Chromosomes | hinter ihnen ist kein Platz
                  if ((left>1)&&(map[left,1]==map[left-1,1])) { # es existiert noch ein Marker auf dem Chromosom vor Marker left
                    delta <- (map[left,2]-map[left-1,2])/(disperse.factor*(right-left+1))
                    map[left:right,2] <- map[left:right,2] - ((right-left):0)*delta
                  } else { # Marker ist erster Marker auf dem Chromosom
                    if (0 < map[left,2]) { #zwischen Marker left und chromosomanfang ist noch "platz"
                      delta <- (map[left,2])/(disperse.factor*(right-left+1))
                      map[left:right,2] <- map[left:right,2] - ((right-left):0)*delta
                    } else {
                      info(-1,paste("Dispersion not possible for marker loci:",paste(as.vector(map[left:right,3]), collapse =" "),"!"))
                    }
                  }
                } else {# marker left right sind am ende des Chromosoms, doch dahinter ist noch platz
                  delta <- (get.genome.par()$chrom.len[map[line-1,1]]-map[left,2])/(disperse.factor*(right-left+1))
                  map[left:right,2] <- map[left:right,2] + (0:(right-left))*delta
                }
              }
            }
            left <- line;
          }
          right <- line;
        } #line
      }
      if (left<right) {
        counter <- counter + right - left
        if (map[left,2] == get.genome.par()$chrom.len[map[line-1,1]]) {
          if ((left>1)&&(map[left,1]==map[left-1,1])) {
            delta <- (map[left,2]-map[left-1,2])/(disperse.factor*(right-left+1))
            map[left:right,2] <- map[left:right,2] - ((right-left):0)*delta
          } else {
            if (0 < map[left,2]) {
              delta <- (map[left,2])/(disperse.factor*(right-left+1))
              map[left:right,2] <- map[left:right,2] - ((right-left):0)*delta
            } else {
              info(-1,paste("Dispersion not possible for marker loci:",paste(as.vector(map[left:right,3]), collapse =" "),"!"))
            }
          }
        } else {
          delta <- (get.genome.par()$chrom.len[map[line-1,1]]-map[left,2])/(disperse.factor*(right-left+1))
          map[left:right,2] <- map[left:right,2] + (0:(right-left))*delta
        }
      }
      info(0,paste(counter,"loci positions modified"))
      file.disperse <- ifelse(is.na(file.disperse),paste(file,".disperse",sep=""),file.disperse)
      info(0,paste("Dispersed linkage map saved to file:",file.disperse))
      write.table(map,file=file.disperse,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
      file <- file.disperse
    }
    c  <- .C("define_map",
             as.character(file))
  }

"linkage.map.save" <- function(file, dta.map=linkage.map.get()) {
  write.table(
    dta.map, 
    file, 
    row.names = FALSE,
    col.names = FALSE,
    quote     = FALSE
  )  
}

"population.individual.remove" <- function (
  PopName,
  IndPos
) {
  ## testen ob Population xxpop vorhanden? und loeschen der genotype info!!
  PopSize <- population.size.get(PopName)$NoInds
  if ((IndPos>=1)&&(IndPos<=PopSize)) {
    if ((1 < PopSize)) {
      if (1==IndPos) {
        population.copy("xxPOP",PopName,2)
        population.swap.name("xxPOP",PopName)
        population.remove("xxPOP")
      } else {
        if (PopSize==IndPos) {
          population.copy("xxPOP",PopName,1,PopSize-1)
          population.swap.name("xxPOP",PopName)
          population.remove("xxPOP")
        } else {
          population.copy("xxPOP",PopName,1,IndPos-1)
          population.copy("xyPOP",PopName,IndPos+1)
          population.swap.name("xxPOP",PopName)
          population.remove("xxPOP")
          concat.population(PopName,"xyPOP")
        }
      }
    } else {
      remove.population(PopName)
    } 
  }
}

"effmap.remove" <-
  function(Name){ 
    c  <- .C("effmap_remove",
             as.character(Name))
  }

"effmap.remove.all" <-
  function(){ 
    c  <- .C("effmap_remove_all");
  }

"remove.effmaps" <-
  function(){ 
    effmap.remove.all()
  }

"remove.evaluate.population" <-
  function(PopNames){ 
    c  <- .C("remove_evaluate_population",
             as.character(PopNames))
  }

"remove.genotype.population" <-
  function(PopNames){
    c  <- .C("remove_genotype_population",
             as.character(PopNames))
  }

"remove.map" <-
  function(){
    c  <- .C("remove_map")
  }

"reset.all" <-
  function(){
    population.remove.all();
    remove.effmaps();
    remove.map();
    reset.mdp();
  }

"reset.mdp" <-
  function(){
    c  <- .C("reset_mdp");
  }

"resources" <-
  function(expr) {
    expr <- substitute(expr)
    stime <- proc.time()
    w <- eval(expr)
    etime <- proc.time()
    on.exit()
    time <- etime - stime
    time[3] <- max(time[3], time[1] + time[2])
    info(0, c(CPU = time[1] + time[2],
            Elapsed = time[3],
            "% CPU" = round((100 * (time[1] + time[2]))/time[3], 1),
           Child = time[4] + time[5]))
    invisible(w)
  }


"return.population" <- function(
  PopNames,
  missing.rm=TRUE
) {
  if (1<length(PopNames)) {
    PopNames <- paste(PopNames, collapse=" ")
  }
  ## hier noch abfrage ob loci verloren gehen!!!
    
  if (missing.rm) {
    bg <- plabsim$allele.missing.indicator
  } else {
    bg <- 3334
  }
  genotype.population(PopNames);
  y <- .Call(
    "return_population",
    as.character(PopNames),
    as.integer(bg)
  )
  if ( is.null(y) ) return ( invisible(NULL) )
  as.data.frame(y)
}

"select.all.best" <- function(
  newPop,
  oldPop,
  effect = NULL,
  x = 1,
  decreasing = TRUE
) {
  popsize <- population.size.get(oldPop)$NoInds
  if ( is.null(popsize) || (0 == popsize)) {
    info(-2,paste("Population",oldPop,"is of size 0"))
    return (invisible(NULL))
  } else {
    genotype.population(oldPop)
    evaluate.population(oldPop, effect)
    population.sort(oldPop, decreasing)
    c <- select.all.best.intern(oldPop,x)
    population.divide(newPop,oldPop,c$n)
    invisible(c)
  }
}

"select.all.best.intern" <- function(
  NamePop,
  n=1
) {
  digits   <- 50
  minscore <- 0
  maxscore <- 0
  indexbest <- paste(rep("a",digits),collapse="")
  c <- .C("select_all_n_best",
          as.character(NamePop),
          indexbest=as.character(indexbest),
          as.integer(n),
          minscore=as.double(minscore),
          maxscore=as.double(maxscore))
  data.frame(n        = as.integer(c$indexbest),
             minscore = as.double(c$minscore),
             maxscore = as.double(c$maxscore))
}

"select.n.best" <- function(
  newPop,
  oldPop,
  effect=NULL,
  n=1,
  decreasing = TRUE
) {
  popsize <- population.size.get(oldPop)$NoInds

  if (0 == popsize) {
    info(-2,paste("Population",oldPop,"is of size 0"))
    return (invisible(NULL))
  } else {
    if (n > popsize) {
      info(1,
           paste( "Population size: ", popsize, ". Selection of ",
                 n, " individuals not possible.",sep=""
                 )
           )
    }
  }
  genotype.population(oldPop)
  evaluate.population(oldPop,effect)
  population.sort(oldPop, decreasing)
  population.divide(newPop,oldPop,n)
}

"set.co.freq" <-
  function(cofreq=1){ 
    c  <- .C("set_co_freq",
             as.double(cofreq))
  }

"set.NoLociInit" <-
  function(NoLociInit){ 
    c  <- .C("set_nolociinit",
             as.integer(NoLociInit))
  }


"effect.weight.set.all" <- function(
  weight
) {
  c <- .C(
    "effect_weight_set_all",
    as.double(weight)
  )
}


"set.eff.weight" <- function(
  fname,
  weight
) {
  c <- .C(
    "set_eff_weight",
    as.character(fname),
    as.double(weight)
  )
}

"genome.parameter.set" <- function(
  no.chrom,
  no.hom=2,
  chrom.len
) {
  if (length(chrom.len) != no.chrom) {
    info(-2, "More chromosomes than chromosome lenghts specified")
  } else {
    population.remove.all();
    remove.map();
    c  <- .C(
             "set_genome_par",
             as.integer(no.chrom),
             as.integer(no.hom),
             as.double(chrom.len)
             )
  }
}

"set.info.level" <- function(
  level
) {
  c  <-  .C(
    "set_info_level",
    as.integer(level)
  )
}

"set.population.gvalue" <-
  function(name,gvalue){
    c  <- .C(
      "set_population_gvalue",
      as.character(name),
      as.double(gvalue)
    )
  }

"population.info.set" <-
  function(PopName,ind,info){
    if ((ind > 0) && (ind <= population.size.get(PopName)$NoInds)) {
      c  <- .C("set_population_info",
               as.character(PopName),
               format(ind,scientific=FALSE),
               as.character(info))
    }
  }

"single.cross" <-
  function(PgName,PopName,classSize=1){
    population.remove(PgName);
    b  <- population.size.get(PopName);
    gtypes <- b$NoInds / classSize;
    for(i in 1:(as.integer(gtypes)-1)){
      copy.population("single-cross-b",PopName,i*classSize+1);
      b  <- population.size.get("single-cross-b");
      for(j in 1:classSize){
        copy.population("single-cross-a",PopName,(i-1)*classSize+j,1);
        cross("single-cross-tmp","single-cross-a","single-cross-b",b$NoInds,crossing.scheme=3)
        concat.population(PgName,"single-cross-tmp");
      }
      population.remove("single-cross-b");
    }
    population.remove("single-cross-a");
  }

"population.sort" <-  function(
  PopName,
  decreasing = TRUE
) {
  c  <- .C("sort_population",
           as.character(paste(PopName,sep=" ")),
           as.integer(decreasing)
           )
}

"ssd.mating" <- function(
  NamePg,
  NameP,
  NoPg=1,
  maxcycles=1
) {
  significant <- -1;
  c  <- .C(
    "ssd_mating",
    as.character(NamePg),
    as.character(NameP),
    format(NoPg,scientific=FALSE),
    as.integer(maxcycles),
    as.double(significant)
  )
}

"summarize.gvalue" <-
  function(pops,effectfile=NULL){
    p <- unlist(strsplit(pops," "))
    n <- length(p) 
    result <- matrix(ncol=7,nrow=n)
    for (i in 1:n) {
      genotype.population(p[i]);
      evaluate.population(p[i]);
      g <- as.matrix(get.population.gvalue(p[i],effectfile))
      result[i,1] <- nrow(g) 
      result[i,2] <- mean(g)
      result[i,3] <- sdev(g)
      result[i,4] <- min (g) 
      result[i,5] <- quantile(g,0.1)
      result[i,6] <- quantile(g,0.5)
      result[i,7] <- max(g)
    }
    colnames(result) <- c("Obs", "Mean", "SDev","Min","Q10","Med","Max")
    rownames(result) <- p
    result
  }

"talk.to.me" <-
  function(){
    set.info.level(1);
  }

"data.params" <-
  function (dta){
    colhead <- splitdt(colnames(dta))
    rowhead <- splitdt(row.names(dta))
    # Number of markers
    no.mar  <- length(levels(as.factor(rowhead[,1])))
    # Number of alleles per marker and list of markers
    no.all      <- matrix(0,ncol=no.mar)
    no.all[1]   <- 1
    mar.list    <- matrix("",nrow=no.mar)
    mar.list[1] <- rowhead[1,1] 
    no.rows <- nrow(rowhead)
    if (no.rows>1){
      k <- 1
      for (i in 2:no.rows) {
        if ( rowhead[i,1] == rowhead[i-1,1] ) {
          no.all[k] <- no.all[k] + 1
        }
        else {
          k            <- k+1
          no.all[k]    <- 1
          mar.list[k]  <- rowhead[i,1]
        }
      }
    }
    # Number of populations
    no.pop  <- length(levels(as.factor(colhead[,1])))
    # Number of individuals and list of populations
    no.ind      <- matrix(0,ncol=no.pop)
    no.ind[1]   <- 1
    pop.list    <- matrix("",nrow=no.pop)
    pop.list[1] <- colhead[1,1]
    no.cols <- nrow(colhead)
    if (no.cols>1){
      k <- 1
      for (i in 2:no.cols) {
        if ( colhead[i,1] == colhead[i-1,1] ) {
          no.ind[k] <- no.ind[k] + 1
        }
        else {
          k           <- k+1
          no.ind[k]   <- 1
          pop.list[k] <- colhead[i,1] 
        }
      }
    }
    # List of individuals
    ind.list <- 1:(no.ind[1])
    if (no.pop>1){
      for (i in 2:no.pop ) ind.list <- c(ind.list,1:(no.ind[i]))
    }
    # Return values
    list(colhead=colhead,rowhead=rowhead,no.mar=no.mar,
         no.all=no.all,no.pop=no.pop,no.ind=no.ind,ind.list=ind.list,
         mar.list=mar.list, pop.list=pop.list)
  }


"splitdt" <-
  function(m.a) {
    .Call("splitdt", as.character(m.a) )
  }

"plabsim.version" <- function() {
  write.version.2()
}


"population.matrix.load" <- function(
  file,
  backcross=FALSE,
  keep.IndID=TRUE
) {
  dta <- read.table(file,header=TRUE)
  generate.population(dta,backcross=backcross)

  if (keep.IndID) {
    dta.param <- gd.data.parameters(dta)
    line <- 1
    for (p in 1:dta.param$no.pop) {
      for (i in 1:dta.param$no.ind[1,p]) {
        set.population.info(
          dta.param$pop.list[p],
          i,
          dta.param$colhead[line,2]
        )
        line <- line + 1
      } #i
    } #p
  }
}

"population.matrix.save" <- function (file,
                                      PopNames,
                                      missing.rm=TRUE,
                                      keep.IndID=TRUE) {
  dta <- return.population(PopNames, missing.rm)

  if (keep.IndID) {
    dta.param <- gd.data.parameters(dta)
    col.names <- NULL
    for (p in 1:dta.param$no.pop) {
      for (i in 1:dta.param$no.ind[1, p]) {
        indid <- as.integer(population.info.get(dta.param$pop.list[p], i))
        if (is.na(indid)) {
          info(-2,"IndID not numeric")
        } else {
        col.names <- c(col.names, 
                       paste(dta.param$pop.list[p],
                             indid,
                             sep=".")
                       )
      }
      } #i     
    } #p
    colnames(dta) <- col.names
  }
  
  write.table(dta,
              file = file,
              quote = FALSE)
}

"population.plabsim.load" <- function(file, PopName=NA) {
  dta.description <- read.table(
    file,
    comment.char = "",
    nrows        = 1
  )
  PopName <- ifelse(
    is.na(PopName),
    as.character(dta.description[1,2]),
    PopName
  )
  dta.info <- read.table(
    file,
    comment.char = "",
    skip         = 1,
    nrows        = dta.description[1,3]
  )
  dta.info <- as.character(dta.info[,3])
  col.classes <- c("integer",      # 1
                   "integer",      # 2
                   "integer",      # 3
                   "numeric",      # 4
                   "integer")      #12

  dta <- read.table(
    file,
    colClasses=col.classes,
    skip = 1+dta.description[1,3]
  )
  colnames(dta) <- c("ind","chrom","hom","pos","all")
  init.population(PopName,dta)
  for (ind in 1:dta.description[1,3]) {
    set.population.info(PopName,ind,dta.info[ind])
  } #ind
}


"population.plabsim.save" <- function(file, PopName) {
  dta <- get.population(PopName)
  write(
    file = file, 
    paste("#",PopName,population.size.get(PopName)$NoInds,sep=" ")
  )
  for (i in 1:population.size.get(PopName)$NoInds) {
    write(
      file   = file, 
      paste("# ",i," \"",population.info.get(PopName,i),"\"",sep=""),
      append = TRUE
    )
  } #i
  write.table(
    file     = file,
    dta,
    quote    = FALSE,
    col.name = FALSE,
    row.name = FALSE,
    append   = TRUE
  )
}

# Code from plabsim.R
# End

# Random numbers
# Start

"rng.choose" <-
  function(rngName="") {
    c <- .C("rng_choose",
            as.character(rngName));
  }

"rng.info" <-
  function() {
    c <- .C("rng_info")
  }

"rng.list" <-
  function() {
    c <- .C("rng_list")
  }

"rng.init" <-
  function() {
    c <- .C("rng_init")
  }

# Random numbers
# End

# Deprecated functions from Plabsim
# Start

"deprecated" <-  function(NameOld,NameNew) {}

"append.population" <-
  function(NameP1, NameP2){
    deprecated("append.population()","population.append()")
    population.append(NameP1, NameP2)
  }

"concat.population" <-
  function(NameP1, NameP2){
    deprecated("concat.population()","population.concat()")
    population.concat(NameP1, NameP2)
  }

"copy.population" <-
  function(NameP1, NameP2, start=1, n=-1) {
    deprecated("copy.population()","population.copy()")
    population.copy(NameP1, NameP2, start, n)
  }

"divide.population" <-
  function(NameP1, NameP2, NoI=-1) {
    deprecated("divide.population()","population.divide()")
    population.divide(NameP1, NameP2, NoI)
  }

"evaluate.allele.ver1" <-
  function(PopName, eLoci, missing=0){
    deprecated("evaluate.allele.ver1()","evaluate.allele()")
    bg     <- 3333
    nrGt   <- 0
    nrLoci <- 0
    llName <- 0
    separator <- "$GPS";
    tail <- "count";
    c  <- .C("evaluate_allele",
             as.character(PopName),  
             eLoci   = as.character(eLoci),
             nrGt    = as.integer(nrGt),
             nrLoci  = as.integer(nrLoci),
             llName = as.integer(llName),
             as.character(separator),
             as.character(tail),
             as.integer(missing),
             as.integer(bg))
    if (c$nrGt<1) return (NA)
    empty <- rep(0,c$nrGt * ((c$nrLoci)+1))
    lName <- paste(rep("a",c$llName),collapse="")
    mode <- 0;
    d  <- .C("print_genotype",
             as.integer(c$nrLoci),
             as.integer(c$nrGt),
             matrix=as.integer(empty),
             as.integer(mode),
             as.character(c$eLoci),
             lName=as.character(lName),
             as.character(separator),
             as.character(tail))
    dim(d$matrix) <- c(c$nrGt, (c$nrLoci)+1)
    lName <- strsplit(d$lName,"\\$GPS");
    result <- data.frame(d$matrix);
    names(result) <- unlist(lName[1]);
    result;
  }


"evaluate.genotype.ver1" <-
  function(PopName, eLoci, alleles=0, mode=2,bg=0) {
    deprecated("evaluate.genotype.ver1()","evaluate.genotype()")
    nrGt   <- 0;
    nrLoci <- 0;
    llName <- 0;
    separator <- "$GPS";
    tail <- "count";
    c  <- .C("evaluate_genotype",
             as.character(PopName),  
             eLoci   = as.character(eLoci),
             as.integer(mode),
             nrGt    = as.integer(nrGt),
             nrLoci  = as.integer(nrLoci),
             llName  = as.integer(llName),
             as.character(separator),
             as.character(tail))
    empty <- rep(0,c$nrGt * ((c$nrLoci)+1))
    lName <- paste(rep("a",c$llName),collapse="")
    mode <- 1;
    d  <- .C("print_genotype",
             as.integer(c$nrLoci),
             as.integer(c$nrGt),
             matrix = as.integer(empty),
             as.integer(mode),
             as.character(c$eLoci),
             lName  = as.character(lName),
             as.character(separator),
             as.character(tail))
    dim(d$matrix) <- c(c$nrGt, (c$nrLoci)+1)
    lName <- strsplit(d$lName,"\\$GPS");
    result <- data.frame(d$matrix);
    n <- length(result)
    summe <- sum(result[,n])
    freq  <- result[,n] / summe
    result <- data.frame(result,freq)
    names(result) <- c(unlist(lName[1]),"Frequency");
    if (alleles[1] != 0) {
      ng <- nrow(result)
      na <- length(alleles)
      passt <- 0; i <- 0
      while ( (i < ng) && (passt==0) ){
        passt <- 1 ; i <- i+1
        for (j in 1:na){
          if (result[i,j] != alleles[j]) {passt <- 0}
        }
        if (passt==1) {select <- i} 
      }
      if (passt==1) {result <- result[select,]}
    }

    if (bg == 1) {
      for (row in dim(result)[1]:1) {
        weg <- 0;
        for (col in 1:(dim(result)[2]-2)) {
          if (3333==result[row,col]) weg <- 1;
        } # col
        if (weg==1) result <- result[-row,];
      } # row
    }
    result
  }

"get.genome.par" <-
  function() {
    deprecated("get.genome.par()","genome.parameter.get()")
    genome.parameter.get()
  }

"get.map" <-
  function(){
    deprecated("get.map()","linkage.map.get()")
    linkage.map.get()
  }

"get.population.info" <-
  function(PopName,ind){
    deprecated("get.population.info()","population.info.get()")
    population.info.get(PopName,ind)
  }

"get.population.size" <-
  function(PopName){
    deprecated("get.population.size()","population.size.get()")
    population.size.get(PopName)
  }

"load.linkage.map" <-
  function(file, disperse=FALSE, file.disperse=NA, disperse.factor=100){
    deprecated("load.linkage.map()","linkage.map.load()")
    linkage.map.load(file, disperse, file.disperse, disperse.factor)
  }

"list.populations" <-
  function(){
    deprecated("list.populations()","population.list()")
    population.list()    
  }

"optimize.population" <-
  function(PopNames){
    deprecated("optimize.population()","population.optimize()")
    population.optimize(PopNames)
    
  }

"population.swap.name" <-
  function(NameP1,NameP2){
    deprecated("population.swap.name()","population.name.swap()")
    population.name.swap(NameP1,NameP2)
  }

"remove.all.populations" <-
  function(){
    deprecated("remove.all.populations()","population.remove.all()")   
    population.remove.all()
  }

"remove.population" <-
  function(PopNames) {
    deprecated("remove.population()","population.remove()")
    population.remove(PopNames)
  }

"rename.population" <-
  function(OldName, NewName){
    deprecated("rename.population()","population.rename()")
    population.rename(OldName, NewName)
  }

"resize.population" <-
  function(PopName, newSize){
    deprecated("rename.population()","population.rename()")
    population.resize(PopName, newSize)
  }

"sample.population" <-
  function(NameP1, NameP2, NoI=-1, rep=1){
    deprecated("sample.population()","population.sample()")
    population.sample(NameP1, NameP2, NoI, rep)
  }

"save.linkage.map" <- function(file, dta.map=get.map()) {
    deprecated("save.linkage.map()","linkage.map.save()")
    linkage.map.save(file, dta.map)
  }


"set.genome.par" <-
  function(no.chrom,no.hom,chrom.len) {
    deprecated("set.genome.par()","genome.parameter.set()")
    genome.parameter.set(no.chrom,no.hom,chrom.len)
  }

"set.population.info" <-
  function(PopName,ind,info){
    deprecated("set.population.info()","population.info.set()")
    population.info.set(PopName,ind,info)
  }

"sort.population" <-
  function( PopNames,
            decreasing = TRUE ) {
    deprecated("sort.population()","population.sort()")
    population.sort(PopNames,decreasing)
  }

"swap.population.name" <-
  function(NameP1,NameP2) {
    deprecated("swap.population.name()","population.swap.name()")
    population.swap.name(NameP1,NameP2)
  }

# Code from plabsim.deprecated.R
# End

## Begin Code 2010

genome.segments <-
function (
            pops,
            allele,
            chromosome=NULL,
            begin=NULL,
            end=NULL
            ) 
{
 p <- unlist (strsplit(pops, " "))
  n <- length(p)
  result <- matrix(ncol = 7, nrow = n)
  for (i in 1:n) {
    g <- evaluate.genome(
            p[i],
            as.integer(allele),
            chromosome = chromosome,
            begin      = begin,
            end        = end,
            verbose    = T
          )
    result[i, 1] <- nrow(g)
    result[i, 2] <- mean(g$noBlocks)
    result[i, 3] <- sdev(g$noBlocks)
    result[i, 4] <- min(g$noBlocks)
    result[i, 5] <- quantile(g$noBlocks, 0.1)
    result[i, 6] <- quantile(g$noBlocks, 0.5)
    result[i, 7] <- max(g$noBlocks)
  }
  colnames(result) <- c("Obs", "Mean", "SDev", "Min", "Q10", 
                        "Med", "Max")
  rownames(result) <- p
  result
}


"select.n.best.segments" <- function(
  newPop,
  oldPop,
  n=1,
  allele = 1 ,                                  
  decreasing = FALSE
) {
  popsize <- population.size.get(oldPop)$NoInds
  if ((n > popsize) || (0 == popsize)) {
    info(-2,
      paste(
        "Population",oldPop,
        "has only",popsize,
        "individuals. Selection of\n",n,
        "individuals is not possible!"
      )
    )
    return (invisible(NULL))
  }
  g <- evaluate.genome(
            oldPop,
            as.integer(allele),
            chromosome = NULL,
            begin      = NULL,
            end        = NULL,
            verbose    = T
          )$noBlocks
  set.population.gvalue(oldPop,g)
  population.sort(oldPop, decreasing)
  population.divide(newPop,oldPop,n)
}

linkage.drag <-
function (
            pops,
            allele,
            chrom,
            pos
            ) 
{
  p <- unlist (strsplit(pops, " "))
  n <- length(p)
  result <- matrix(ncol = 7, nrow = n)
  for (i in 1:n) {
    g <- evaluate.ld(
            p[i],
            as.integer(allele),
            chrom = chrom,
            pos   = pos
          )
    result[i, 1] <- nrow(g)
    result[i, 2] <- mean(g$length)
    result[i, 3] <- sdev(g$length)
    result[i, 4] <- min(g$length)
    result[i, 5] <- quantile(g$length, 0.1)
    result[i, 6] <- quantile(g$length, 0.5)
    result[i, 7] <- max(g$length)
  }
  colnames(result) <- c("Obs", "Mean", "SDev", "Min", "Q10", 
                        "Med", "Max")
  rownames(result) <- p
  result
}

"evaluate.ld" <- function (
  PopName,
  allele,
  chrom,
  pos
)
{

  empty <- rep(0,as.integer(population.size.get(PopName)$NoInds));

  c  <- .C("evaluate_ld",
             as.character(PopName),
             as.integer(allele),
             as.integer(chrom),
             as.double(pos),
             ld  = as.double(empty)
           )
  
  ld  <- data.frame(length=c$ld)
  return(ld)
}


## End Code 2010


###############################################################################
# Genetic Distances BEGIN
###############################################################################

"gd.allow.zero.frequencies" <-
  function(allow=0){
    c  <- .C("set_allow_zeros",
             as.integer(allow)
             )
  }

"gd.status.zero.frequencies" <-
  function(){
    allow=99;
    c  <- .C("get_allow_zeros",
             allow=as.integer(allow)
             )
    c$allow
  }

"gd.allele.frequencies" <-
  function(dta){
    p <- gd.data.parameters(dta)
    no.pop.u <- p$no.pop
    pop.u    <- 1:p$no.pop
    no.ind.u <- p$no.ind
    ind.u    <- p$ind.list
    no.mar.u <- p$no.mar
    mar.u    <- 1:p$no.mar
    DTA  <- as.matrix(dta)
    DTA[is.na(DTA)] <- -1
    FREQ <- matrix(16,nrow=sum(p$no.all[mar.u]),ncol=p$no.pop)
    c  <- .C("allele_freq",
             as.integer(p$no.pop),
             as.integer(p$no.ind),
             as.integer(p$no.mar),
             as.integer(p$no.all),
             as.integer(no.ind.u),
             as.integer(ind.u),
             DTA =as.double(DTA),
             FREQ=as.double(FREQ)
            )
    c$FREQ[c$FREQ==-1] <- NA
    FREQ <- data.frame(matrix(c$FREQ,ncol=p$no.pop))
    names(FREQ) <- c(as.character(p$pop.list))
    row.names(FREQ) <- row.names(dta)
    FREQ
  }

"gd.correct.missing" <-
  function(dta){
    p <- gd.data.parameters(dta)
    no.pop.u <- p$no.pop
    pop.u    <- 1:p$no.pop
    no.ind.u <- p$no.ind
    ind.u    <- p$ind.list
    no.mar.u <- p$no.mar
    mar.u    <- 1:p$no.mar
    DTA  <- as.matrix(dta)
    DTA[is.na(DTA)] <- -1
    FREQ <- matrix(16,nrow=sum(p$no.all[mar.u]),ncol=p$no.pop)
    c  <- .C("correct_missing",
             as.integer(p$no.pop),
             as.integer(p$no.ind),
             as.integer(p$no.mar),
             as.integer(p$no.all),
             as.integer(no.ind.u),
             as.integer(ind.u),
             DTA =as.double(DTA),
             FREQ=as.double(FREQ)
            )
    c$DTA[c$DTA==-1] <- NA
    c$DTA[c$DTA==0.5] <- 1
    result <- data.frame(matrix(c$DTA,ncol=length(p$ind.list)))
    names(result) <- names (dta)
    rownames(result) <- rownames (dta)
    result
 }

"gd.list.missing" <-
  function (dta){
    r <- data.frame(ind=" ",marker=" ")
    y <- gd.splitdot(rownames(dta))[,1]
    for (i in 1:length(dta)){
      m <- levels(as.factor(as.character(y[is.na(dta[,i])])))
      n <-  names(dta)[i]
      r <- rbind(r,data.frame(ind=rep(n,length(m)),marker=m))
      remove(m,n)
    }
    r <- r[-1,]
    r
  }

"gd.list.irregular" <-
  function (dta1){
    dta <- gd.correct.missing(dta1)
    r <- data.frame(ind=" ",marker=" ")
    y <- gd.splitdot(rownames(dta))[,1]
    for (i in 1:length(dta)){
      index <-  (dta[,i] != 0) & (dta[,i] != 1) & (dta[,i] != 0.5)
      m <- levels(as.factor(as.character(y[index])))
      n <-  names(dta)[i]
      r <- rbind(r,data.frame(ind=rep(n,length(m)),marker=m))
      remove(m,n)
    }
    r <- r[-1,]
    r
  }

"gd.distance.similarity" <-
  function(dta,measure,par,boot=0,
           no.pop.u=0,pop.u=0,no.ind.u=0,ind.u=0,no.mar.u=0,mar.u=0){
    # Distance measures
    if (measure=="rd")  m <- 1
    if (measure=="mrd") m <- 2
    if (measure=="euc") m <- 3
    # Similarity coefficients
    if (measure=="dic") m <- 11
    if (measure=="jac") m <- 12
    if (measure=="sma") m <- 13
    if (pmatch("sdev",par,nomatch=0)) s <- 1 else s <- 0
    if (pmatch("sdev/jack",par,nomatch=0)) s <- 1 else s <- 0
    if (pmatch("sdev/boot",par,nomatch=0)) s <- boot
    p <- gd.data.parameters(dta)
    if(no.pop.u==0) no.pop.u <- p$no.pop
    if(pop.u   ==0) pop.u    <- 1:p$no.pop
    if(no.ind.u==0) no.ind.u <- p$no.ind
    if(ind.u   ==0) ind.u    <- p$ind.list
    if(no.mar.u==0) no.mar.u <- p$no.mar
    if(mar.u   ==0) mar.u    <- 1:p$no.mar
    DTA  <- as.matrix(dta)
    DTA[is.na(DTA)] <- -1
    FREQ <- matrix(16,nrow=sum(p$no.all[mar.u]),ncol=p$no.pop)
    u    <- (((p$no.pop^2)-p$no.pop)/2 )
    DIST <- as.double(1:u)
    if (s) SDEV  <- DIST
    else SDEV  <- 1
    c  <- .C("gen_dist",
             as.integer(m),
             as.integer(p$no.pop),
             as.integer(p$no.ind),
             as.integer(p$no.mar),
             as.integer(p$no.all),
             as.integer(no.pop.u),
             as.integer(pop.u),
             as.integer(no.ind.u),
             as.integer(ind.u),
             as.integer(no.mar.u),
             as.integer(mar.u),
             DTA =as.double(DTA),
             FREQ=as.double(FREQ),
             DIST=as.double(DIST),
             as.integer(s),  
             SDEV=as.double(SDEV),
             as.character(par)
            )
    r1 <- r2 <- array("",c(u,1))
    i3 <- 0;
    for (i1 in 1:(p$no.pop-1)){ 
      for (i2 in (i1+1):p$no.pop){
        i3 <- i3 + 1;
        r1[i3] <- p$pop.list[i1]; 
        r2[i3] <- p$pop.list[i2];
      }
    }
    result <- data.frame(OTU1=r1,OTU2=r2,Measure=c$DIST)
    if (s)  result <- data.frame(result,SDev=c$SDEV);
    result
  }

"gd.similarity.coefficient" <-
  function(dta,measure="jac",par="dist",boot=0,
           no.pop.u=0,pop.u=0,no.ind.u=0,ind.u=0,no.mar.u=0,mar.u=0,
           resample.primers=FALSE){

    if ( (measure != "dic") && (measure != "jac") && (measure != "sma") ){
      info(-2,"Distance measure not defined")
      return (invisible(NULL))
    }
    p <- gd.data.parameters(dta)
    if ( sum(p$no.all)==p$no.mar){
      if ( gd.status.zero.frequencies()==0 ){
        info (-1,"Dominant marker data?")
        info (-1,"Consider 'gd.allow.zero.frequencies(1)'")
      }
    }
    if ( sum(p$no.all)>p$no.mar){
      if ( gd.status.zero.frequencies()==1 ){
        info (-1,"Co-dominant marker data?")
        info (-1,"Consider 'gd.allow.zero.frequencies(0)'")
      }
    }
    if ( sum(p$no.ind)>p$no.pop){
      info (-1,"Each population should consist of only one individual")
      info (-1,"when using similarity coefficients.")
    }
    y <- dta
    if (resample.primers){
      # Respample over primers
      x <- rownames(dta)
      m <- regexpr("\\.",x)
      good <- m > 0
      lens <- attr(m,"match.length")
      x <- paste(substring(x,1,m[good]-1),substring(x,m[good]+lens),sep="")
      m <- regexpr("\\+",x)
      good <- m > 0
      lens <- attr(m,"match.length")
      x <- paste(substring(x,1,m[good]-1),".",substring(x,m[good]+lens),sep="")
      rownames(y) <- x
    }
    
    gd.distance.similarity(y,measure,par,boot,
           no.pop.u,pop.u,no.ind.u,ind.u,no.mar.u,mar.u)
  }

"gd.genetic.distance" <-
  function(dta,measure="euc",par="dist",boot=0,
           no.pop.u=0,pop.u=0,no.ind.u=0,ind.u=0,no.mar.u=0,mar.u=0){

    if ( (measure != "euc") && (measure != "rd") && (measure != "mrd") ){
      info(-2,"Distance measure not defined")
      return (invisible(NULL))
    }
    
    p <- gd.data.parameters(dta)

    if ( sum(p$no.all)==p$no.mar){
      if ( gd.status.zero.frequencies()==0 ){
        info (-1,"Dominant marker data?")
        info (-1,"Consider 'gd.allow.zero.frequencies(1)'")
      }
      if (measure!="euc") {
        info (-1,"Consider the Euclidean distance for dominant data")
      }
      if ( sum(p$no.ind)>p$no.pop){
        info (-1,"analysis on population level for dominant data")
      }
    }

    if ( sum(p$no.all)>p$no.mar){
      if ( gd.status.zero.frequencies()==1 ){
        info (-1,"Co-dominant marker data?")
        info (-1,"Consider 'gd.allow.zero.frequencies(0)'")
      }
    }
    gd.distance.similarity(dta,measure,par,boot,
           no.pop.u,pop.u,no.ind.u,ind.u,no.mar.u,mar.u)
  }

"gd.data.parameters" <-
  function (dta){
    colhead <- gd.splitdot(colnames(dta))
    rowhead <- gd.splitdot(row.names(dta))
    # Number of markers
    no.mar  <- length(levels(as.factor(rowhead[,1])))
    # Number of alleles per marker and list of markers
    no.all      <- matrix(0,ncol=no.mar)
    no.all[1]   <- 1
    mar.list    <- matrix("",nrow=no.mar)
    mar.list[1] <- rowhead[1,1] 
    no.rows <- nrow(rowhead)
    if (no.rows>1){
      k <- 1
      for (i in 2:no.rows) {
        if ( rowhead[i,1] == rowhead[i-1,1] ) {
          no.all[k] <- no.all[k] + 1
        }
        else {
          k            <- k+1
          no.all[k]    <- 1
          mar.list[k]  <- rowhead[i,1]
        }
      }
    }
    # Number of populations
    no.pop  <- length(levels(as.factor(colhead[,1])))
    # Number of individuals and list of populations
    no.ind      <- matrix(0,ncol=no.pop)
    no.ind[1]   <- 1
    pop.list    <- matrix("",nrow=no.pop)
    pop.list[1] <- colhead[1,1]
    no.cols <- nrow(colhead)
    if (no.cols>1){
      k <- 1
      for (i in 2:no.cols) {
        if ( colhead[i,1] == colhead[i-1,1] ) {
          no.ind[k] <- no.ind[k] + 1
        }
        else {
          k           <- k+1
          no.ind[k]   <- 1
          pop.list[k] <- colhead[i,1] 
        }
      }
    }
    # List of individuals
    ind.list <- 1:(no.ind[1])
    if (no.pop>1){
      for (i in 2:no.pop ) ind.list <- c(ind.list,1:(no.ind[i]))
    }
    # Return values
    list(colhead=colhead,rowhead=rowhead,no.mar=no.mar,
         no.all=no.all,no.pop=no.pop,no.ind=no.ind,ind.list=ind.list,
         mar.list=mar.list, pop.list=pop.list)
  }

"gd.splitdot" <-
  function(m.a) {
    .Call("splitdot", as.character(m.a) )
  }

"gd.mk.matr" <- function (dta,distance){
  p <- gd.data.parameters(dta)
  aa <- distance[,3]
  bb <- p$no.pop
  cc <- p$pop.list
  dd <- structure(aa, Size = bb, Labels = cc,
                  Diag = FALSE, Upper = FALSE, method = "euclidean", class = "dist")
  as.matrix(dd)         
}

"gd.pcoa" <-  function(dta,k=3){
  distance <- gd.genetic.distance(dta,measure="mrd")
  D <- gd.mk.matr(dta,distance)
  Y <-  cmdscale( D, k=k ) 
  Y
}

###############################################################################
# Genetic Distances End
###############################################################################

sdev <- function (x, na.rm = FALSE) 
{
    if (is.matrix(x)) {
        apply(x, 2, sd, na.rm = na.rm)
    }
    else if (is.vector(x)) 
        sqrt(var(x, na.rm = na.rm))
    else if (is.data.frame(x)) {
        sapply(x, sd, na.rm = na.rm)
    }
    else sqrt(var(as.vector(x), na.rm = na.rm))
}

