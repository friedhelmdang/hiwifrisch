###########################################################################
# Speedtests
###########################################################################


# Load portable version 
  st.script.dir <- "~/1-daten/p111-selection-tools/SelectionTools-dev/src/"
  source(paste(st.script.dir,"/SelectionTools.R",sep=""))
  version <- .C("print_version")
  st.set.info.level(1)
#

library(SelectionTools)
st.set.info.level(1)


###########################################################################
# Matrix calculations with buildin R functions / BLAS of R
###########################################################################

st.set.num.threads(4)
system.time({ x <- replicate(5e3, rnorm(5e3)); tcrossprod(x) })

# x8890 @2.9 MKL2016
#_4T_ : 2.7s 

# x5670 2x6c@2.9 MKL2016
#_4T_ : 6.2s 
#_8T_ : 4.85s 

# i7 6920 4c@2.9 
# Rblas: 38s
# Openblas: 2.67s

# i7 4910 4c@2.9 
# Openblas: 3.05s

st.set.num.threads(4)
system.time({ x <- replicate(1e4, rnorm(1e4)); tcrossprod(x) })

# i7 4910 4c@2.9 
# Openblas: 15.59s

# i5 6440 4c@2.6
# Openblas: 15.28

###########################################################################
# Matrix calculations
###########################################################################

setwd("x:/1-daten/p111-selection-tools/SelectionTools-data/fk-mais")
setwd("e:/1-daten/p111-selection-tools/SelectionTools-data/fk-mais")

setwd("~/1-daten/p111-selection-tools/SelectionTools-data/fk-mais")

st.input.dir  <- "./"
st.data.dir   <- "./"
st.output.dir <- "./output"
st.set.info.level(1)

# { # Small data set
#   st.read.marker.data ("small-1.mpo",format="m",auxfiles=F)
#   st.read.performance.data ( "large-1-gdc.p" )
# }

{
  st.read.marker.data ("large-1.g",format="m",auxfiles=F)
  st.restrict.marker.data   ( NoAll.MAX = 2   )  # Max missing at a marker
  st.restrict.marker.data   ( MaMis.MAX = 0.1 )  # Max missing at a marker
  st.restrict.marker.data   ( InMis.MAX  =0.1 )  # Max missing per individual
  st.restrict.marker.data   ( ExHet.MIN = 0.1 )  # Minimum gene diversity
  st.read.performance.data ( "large-1-gdc.p" )
}

##########################################################################

st.set.num.threads(24)
st.set.num.threads(12)
st.set.num.threads(8)
st.set.num.threads(4)

st.get.num.threads()


st.set.matrix.ops("st") 
st.set.matrix.ops("mt") 
st.set.matrix.ops("gsl") 
st.set.matrix.ops("blaslapack") 

{
  gs.build.Z    ( auxfiles=F )   
  gs.lambda.hsq ( hsq=0.9, auxfiles=F )    
  gs.mme.coeff  ( auxfiles=F )
  gs.mme.rhs    ( auxfiles=F ) 
  gs.mme.solve  ( auxfiles=F ) 
  gs.plot.model.fit()
}

# x8890, @2.2ghz, lin mkl SelectionTools 16.1
# _4_ threads
# gs.mme.coeff: 11.3 
# gs.mme.solve: 57.57s
# _8_ threads
# gs.mme.coeff: 7.53 
# gs.mme.solve: 32.29

# x5670 2x6c@2.9 (s302) lin MKL2016
# _4_ threads
# gs.mme.coeff: 26.72s
# gs.mme.solve: 3.27
# _8_ threads
# gs.mme.coeff: 18.6s
# gs.mme.solve: 1.69m, 1.4m

# x7542, 4x6c@2,6ghz, 256G, lin mkl
# gs.mme.coeff: 24.4s, 19.63s
# gs.mme.solve: 1.39m, 1.4m

# i7-4190,   4c@2.9ghz,     32G, win openblasp2.19
# gs.mme.coeff: 9.95s
# gs.mme.solve: 1.25m

###########################################################################
# OMP parallel performance GD
###########################################################################

{
  st.read.marker.data ("large-1.g",format="m",auxfiles=F)
  st.read.performance.data ( "large-1-gdc.p" )
}

st.set.num.threads(48)
st.set.num.threads(12)
st.set.num.threads(8)
st.set.num.threads(4)

st.get.num.threads()

st.genetic.distances (auxfiles=FALSE)  


# x8890, @2.2ghz, (s302) 
# _4_  threads: 2.53m
# _8_  threads: 1.21m
# _12_ threads: 1.08m

# x5670 2x6c@2.9 (s302) 
# _4_ threads:  5.97m
# _8_ threads:  4.48m 
# _12_ threads: 3.98m



###########################################################################
# OMP parallel performance LD
###########################################################################


setwd("~/1-daten/p111-selection-tools/SelectionTools-data/2013-aph-cr")
st.input.dir  <- "input"
st.output.dir <- "output"

{ # Read in marker and phenotypic data 
  st.read.marker.data ( "zr.mpo", format = "m")
  st.read.map("zr.map",format="mcp")
}

  st.restrict.marker.data   ( NoAll.MAX = 2   )  # Max missing at a marker
  st.restrict.marker.data   ( ExHet.MIN = 0.1 )  # Minimum gene diversity


##########################################################################

st.set.num.threads(24)
st.set.num.threads(12)
st.set.num.threads(8)
st.set.num.threads(4)
st.set.num.threads(1)

st.get.num.threads()

system.time()

ld <- st.calc.ld(ld.measure = "r2", auxfiles = FALSE)
ld <- st.calc.ld(ld.measure = "r2-estimate-phases", auxfiles = FALSE)

