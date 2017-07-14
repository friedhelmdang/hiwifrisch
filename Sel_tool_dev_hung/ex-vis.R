###########################################################################
# LD based haplotypes with SelectionTools
###########################################################################

{
   st.input.dir  <- "input"
   st.output.dir <- "output"
   st.script.dir <- "src/"
   source(paste(st.script.dir,"/SelectionTools.R",sep=""))
   version <- .C("print_version") 
   library("GenABEL")
   library("LDheatmap")
 }

# Load maize data set of Crossa et al. 2010

{ 
  st.read.marker.data ("ex-crossa.pop" ,format="m")
  st.read.map         ("ex-crossa.map", format="mcp", skip=1)
  #
  cl.s <- read.table(file="input/ex-crossa.clusters")  #  Clustering results
  st.copy.marker.data ( "q3", "default" )              # Only cluster #3
  #
  st.restrict.marker.data   ( ind.list=row.names(cl.s)[cl.s==3] , data.set="q3") 
  st.restrict.marker.data   ( MaMis.MAX = 0.1, data.set="q3" ) 
  st.restrict.marker.data   ( InMis.MAX = 0.1, data.set="q3" ) 
  st.restrict.marker.data   ( ExHet.MIN = 0.1, data.set="q3" ) 
}

# LD heatmaps for chromosome 3 in cluster 3

CHR <- 8

{
  ld <- st.calc.ld  ( ld.measure="r2", data.set="q3" )
  l  <- st.LDplot.r ( ld, chrom=CHR )
  m  <- st.LDplot.m ( chrom=CHR, data.set="q3")
  dev.new()
  LDheatmap( l, m,
             distances = "genetical", LDmeasure = "r",
             title=paste("LD calculated with SelectionTools, Chrom.",CHR),
             col=heat.colors(20) )
}

{
  q3 <- st.conv.genable("q3")
  idx.marker <- (gtdata(q3)@chromosome == CHR)
  c <- gtdata(q3) [,idx.marker]
  l <- r2fast(c)
  m <- c@map/10000
  dev.new() 
  LDheatmap( l, m,
            distances = "genetical", LDmeasure = "r",
            title=paste("LD calculated with GenABEL, Chrom.",CHR),
            col=heat.colors(20) )
}






{
  st.copy.marker.data ( "t","q3" ) 
  h <- st.def.hblocks ( hap=3, data.set="t" ) 
  st.recode.hil  (data.set="t")        
  dev.new()
  st.plot.ggt ( data.set="t" )
}

nrow(h)

{
  st.copy.marker.data ( "t","q3" ) 
  h <- st.def.hblocks ( hap.unit=0, data.set="t" ) 
  st.recode.hil  (data.set="t")        
  dev.new()
  st.plot.ggt (data.set="t" )
}

{

  st.copy.marker.data ( "q","p" ) 

  h <- st.def.hblocks ( ld.adjacent=0,
                        ls.multilocus=0,
                        data.set="q" )
  
  st.recode.hil  (data.set="q")        
  dev.new()
  st.plot.ggt ( i.list=cl.x, data.set="q" )
}



{ # Load sugar beet data set
  st.read.marker.data ( "geno_train.csv", format = "t" )
  st.read.map("ZR_INT_1101_MABC-corrected.map",m.stretch=1)
  st.restrict.marker.data   ( NoAll.MAX = 2 )  
  st.restrict.marker.data   ( MaMis.MAX = 0.05 )  # Max missing at a marker
  st.restrict.marker.data   ( InMis.MAX = 0.05 )  # Max missing per individual
  st.restrict.marker.data   ( ExHet.MIN = 0.05 )  # Minimum gene diversity
  st.copy.marker.data ( "q3", "default" )        
  #
  x <- st.marker.data.statistics(data.set="q3")$marker.list$Name
  keep <- x[! x %in% c("s39g24xxs1","s3mp0166xx")]
  st.restrict.marker.data ( mar.list=keep, data.set="q3") 
}


#######################################################################
# Compare with Genable
#######################################################################





nsnps(q3); nids(q3)


{ 
  idx.marker <- (gtdata(q3)@chromosome == 3)
  c <- gtdata(q3) [,idx.marker]
  r <- r2fast(c)
  m <- c@map/1000
  dev.new(); LDheatmap( r, m,  col=heat.colors(20))
  r <- dprfast(c)
  dev.new();LDheatmap( r, m,  col=heat.colors(20))
}

r[1:10,1:10]



{ 
  idx.marker <- (gtdata(p)@chromosome == 3)
  c <- gtdata(p) [,idx.marker]
  r <- r2fast(c)
  m <- c@map/1000
}

LDheatmap( r, m,  col=heat.colors(20))


  r <- r2fast(p)
  m <- p@map/1000
  LDheatmap( r, m,  col=heat.colors(20))

r[1:5,1:5]
m[1:5]

dprfast(p)
