###########################################################################
# LD based haplotypes with SelectionTools
###########################################################################

{
   st.input.dir  <- "input"
   st.output.dir <- "output"
   #
   st.script.dir <- "src/"
   source(paste(st.script.dir,"/SelectionTools.R",sep=""))
   version <- .C("print_version")
   # library("SelectionTools")
   #
   library("LDheatmap")            
   LDplot <- function ( ld.list, chrom, data.set,
                       ld.measure="r2", title="", no.colors=20)
     {
       l  <- st.LDplot.ld  ( ld.list, chrom )
       m  <- st.LDplot.map ( chrom, data.set)
       if ("Dp"==ld.measure) LDmeasure <- "Dprime" else  LDmeasure <- "r"
       LDheatmap ( l, m, distances = "genetical", LDmeasure = LDmeasure,
                 title=paste("Chrom.",chrom,title), col=heat.colors(no.colors) )
     }
 }


{ # Load maize data set of Crossa et al. 2010
  st.read.marker.data ("ex-crossa.pop" ,format="m")
  st.read.map         ("ex-crossa.map", format="mcp", skip=1)
  st.restrict.marker.data ( MaMis.MAX = 0.01 ) 
  st.restrict.marker.data ( InMis.MAX = 0.01 ) 
  st.restrict.marker.data ( ExHet.MIN = 0.1 ) 
  #
  cl.s <- read.table(file="input/ex-hapLD.clusters")  #  Clustering results
  st.copy.marker.data ( "q2", "default" )              # Only cluster #3
  st.copy.marker.data ( "q3", "default" )              # Only cluster #3
  #
  st.restrict.marker.data ( ind.list=row.names(cl.s)[cl.s==2], data.set="q2") 
  st.restrict.marker.data ( ind.list=row.names(cl.s)[cl.s==3], data.set="q3") 
}

  ld1 <- st.calc.ld  ( ld.measure="r2")

  ld2 <- st.calc.ld  ( ld.measure="r2-estimate-phases")

head(data.frame(ld1,ld2[,6]))

  options(width=120)

  st.copy.marker.data ( "t2", "q2" )       
  st.copy.marker.data ( "t3", "q3" )       

  ld <- st.calc.ld  ( ld.measure="r2",
                      data.set="t3")

  h <- st.def.hblocks ( ld.threshold = 0.6,
                        tolerance = 3,
                        ld.criterion = "flanking", 
                        data.set="t3" )

  h2 <- st.set.hblocks (h,
                        hap.symbol="c", 
                        data.set="t"  )

  pdf(file="ex-hapLD-fig01.pdf",width=4, height=4,pointsize=10)

{
  st.recode.hil  (data.set="t3" )
  ld <- st.calc.ld  ( ld.measure="r2",data.set="t3" )
  x11();LDplot ( ld, 8, "t3")
}

{
  st.recode.hil  (data.set="t3" )
  ld <- st.calc.ld  ( ld.measure="r2",data.set="t3" )
  x11();LDplot ( ld, 8, "t3")
}

  dev.off()
