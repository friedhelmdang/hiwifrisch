###########################################################################
# Visualizing DNA-marker data
###########################################################################

{
   st.input.dir  <- "input"
   st.output.dir <- "output"
   st.script.dir <- "src/"
   source(paste(st.script.dir,"/SelectionTools.R",sep=""))
   version <- .C("print_version") 
}

library("SelectionTools")

{ st.read.marker.data ("ex-crossa.pop" ,format="m")
  st.restrict.marker.data   ( NoAll.MAX = 2   ) 
  st.restrict.marker.data   ( InMis.MAX = 0.1 ) 
  st.restrict.marker.data   ( ExHet.MIN = 0.1 ) 
  st.read.map         ("ex-crossa.map", format="mcp", skip=1)  }

st.restrict.marker.data   ( MaMis.MAX = 0.001 )

a <- st.calc.ld.2(ld.measure="r2")
b <- st.calc.ld(ld.measure="r2-estimate-phases")

head(a)
head(b)

a <- st.calc.ld.2(ld.measure="Dp")
b <- st.calc.ld(ld.measure="Dp-estimate-phases")

head(a)
head(b)
  
###########################################################################
# Genetic distances
###########################################################################

st.set.info.level(1)
options(width=256)

st.read.marker.data ("ex-crossa.pop" ,format="m")

dist.mat <- st.genetic.distances.fct ( measure="rd", split=3, format="m" )  
as.matrix(dist.mat)[1:10,1:10]


dist.mat.c <- st.genetic.distances ( measure="rd",  format="m" )  
as.matrix(dist.mat.c)[1:10,1:10]

